from __future__ import annotations

import csv
import json
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Sequence, Tuple

_LOGGER = logging.getLogger("magpie.rrna")

FASTA_EXTS = (".fa", ".fna", ".fasta")


@dataclass(frozen=True)
class FastaRecord:
    header: str
    seq: str

    @property
    def id(self) -> str:
        return self.header.split()[0]


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _strip_fasta_suffix(name: str) -> str:
    s = name.strip()
    low = s.lower()
    for ext in FASTA_EXTS:
        if low.endswith(ext):
            return s[: -len(ext)]
    return s


def _genome_id_from_16s_filename(path: Path) -> str:
    base = _strip_fasta_suffix(path.name)
    if base.endswith("_16S"):
        return base[:-4]
    return base


def _iter_16s_fastas(dir_: Path) -> Iterable[Path]:
    if not dir_.exists():
        return
    for p in sorted(dir_.iterdir(), key=lambda x: x.name):
        if not p.is_file():
            continue
        low = p.name.lower()
        if any(low.endswith(f"_16s{ext}") for ext in FASTA_EXTS):
            yield p


def _read_fasta(path: Path) -> List[FastaRecord]:
    records: List[FastaRecord] = []
    header: str | None = None
    seq_parts: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header=header, seq="".join(seq_parts)))
                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA without header in {path}")
                seq_parts.append(line)

    if header is not None:
        records.append(FastaRecord(header=header, seq="".join(seq_parts)))

    return records


def _write_fasta(records: Sequence[FastaRecord], path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        for rec in records:
            out.write(f">{rec.header}\n")
            out.write(f"{rec.seq}\n")


def _copy_file(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def _pick_longest(records: Sequence[FastaRecord]) -> FastaRecord:
    if not records:
        raise ValueError("Cannot pick longest from empty record list.")
    return sorted(records, key=lambda r: (len(r.seq), r.header), reverse=True)[0]


def _run_vsearch_cluster_fast(
    *,
    in_fasta: Path,
    centroids_out: Path,
    cluster_id: float,
    threads: int,
    vsearch_bin: str,
    uc_out: Path | None = None,
) -> None:
    cmd = [
        vsearch_bin,
        "--cluster_fast",
        str(in_fasta),
        "--id",
        str(cluster_id),
        "--centroids",
        str(centroids_out),
        "--threads",
        str(threads),
    ]
    if uc_out is not None:
        cmd.extend(["--uc", str(uc_out)])

    _LOGGER.debug("Running VSEARCH: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"VSEARCH failed for {in_fasta}\n"
            f"Command: {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )


def _domain_paths(out: Path, domain: str) -> Dict[str, Path]:
    dom = out / domain
    return {
        "root": dom,
        "single": dom / "single",
        "multiple": dom / "multiple",
        "clustered": dom / "clustered",
        "copies_tsv": dom / f"{domain}_16S_copies.tsv",
        "genes_fasta": dom / f"{domain}_16S_genes.fasta",
        "genes_map": dom / f"{domain}_16S_genes.map.tsv",
        "centroids": dom / f"{domain}_16S_centroids.fasta",
        "clusters_uc": dom / f"{domain}_16S_clusters.uc",
    }


def _initialise_domain_dirs(paths: Dict[str, Path]) -> None:
    for key in ("root", "single", "multiple", "clustered"):
        _ensure_dir(paths[key])


def _count_and_split_domain(
    *,
    barrnap_domain_dir: Path,
    paths: Dict[str, Path],
) -> Tuple[List[List[str]], Dict[str, int]]:
    rows: List[List[str]] = []
    counts = {
        "total_fastas": 0,
        "with_16S": 0,
        "single": 0,
        "multiple": 0,
        "none": 0,
    }

    for fasta_fp in _iter_16s_fastas(barrnap_domain_dir):
        genome_id = _genome_id_from_16s_filename(fasta_fp)
        records = _read_fasta(fasta_fp)
        n = len(records)

        counts["total_fastas"] += 1
        category = "none"
        out_name = f"{genome_id}.fna"

        if n == 1:
            category = "single"
            counts["with_16S"] += 1
            counts["single"] += 1
            _write_fasta(records, paths["single"] / out_name)
        elif n > 1:
            category = "multiple"
            counts["with_16S"] += 1
            counts["multiple"] += 1
            _write_fasta(records, paths["multiple"] / out_name)
        else:
            counts["none"] += 1

        rows.append([genome_id, str(n), category, str(fasta_fp)])

    with paths["copies_tsv"].open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "copies", "category", "source_fasta"])
        w.writerows(rows)

    return rows, counts


def _cluster_multi_copy_domain(
    *,
    paths: Dict[str, Path],
    vsearch_bin: str,
    multi_cluster_id: float,
    threads: int,
) -> int:
    clustered_files = 0

    multi_files = sorted(
        p for p in paths["multiple"].iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    )

    for fp in multi_files:
        out_fp = paths["clustered"] / fp.name
        _run_vsearch_cluster_fast(
            in_fasta=fp,
            centroids_out=out_fp,
            cluster_id=multi_cluster_id,
            threads=threads,
            vsearch_bin=vsearch_bin,
        )
        clustered_files += 1

    return clustered_files


def _build_genome_representatives(
    *,
    domain: str,
    paths: Dict[str, Path],
) -> Dict[str, int]:
    picked: Dict[str, FastaRecord] = {}
    source: Dict[str, str] = {}

    clustered_files = sorted(
        p for p in paths["clustered"].iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    )
    single_files = sorted(
        p for p in paths["single"].iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTS
    )

    # First, take clustered multi-copy genomes (one representative per genome)
    for fp in clustered_files:
        genome_id = _strip_fasta_suffix(fp.name)
        recs = _read_fasta(fp)
        if not recs:
            continue
        rep = _pick_longest(recs)
        picked[genome_id] = FastaRecord(header=genome_id, seq=rep.seq)
        source[genome_id] = "clustered"

    # Then, add single-copy genomes not already present
    for fp in single_files:
        genome_id = _strip_fasta_suffix(fp.name)
        if genome_id in picked:
            continue
        recs = _read_fasta(fp)
        if not recs:
            continue
        rep = _pick_longest(recs)
        picked[genome_id] = FastaRecord(header=genome_id, seq=rep.seq)
        source[genome_id] = "single"

    ordered_ids = sorted(picked)
    _write_fasta([picked[g] for g in ordered_ids], paths["genes_fasta"])

    with paths["genes_map"].open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "source"])
        for gid in ordered_ids:
            w.writerow([gid, source[gid]])

    lengths = [len(picked[g].seq) for g in ordered_ids]
    return {
        "representatives": len(ordered_ids),
        "min_len": min(lengths) if lengths else 0,
        "max_len": max(lengths) if lengths else 0,
        "mean_len": (sum(lengths) / len(lengths)) if lengths else 0.0,
    }


def _final_cluster_domain(
    *,
    paths: Dict[str, Path],
    vsearch_bin: str,
    final_cluster_id: float,
    threads: int,
) -> None:
    recs = _read_fasta(paths["genes_fasta"])
    if not recs:
        paths["centroids"].write_text("", encoding="utf-8")
        paths["clusters_uc"].write_text("", encoding="utf-8")
        _LOGGER.warning("No representative 16S sequences found for %s", paths["root"].name)
        return

    _run_vsearch_cluster_fast(
        in_fasta=paths["genes_fasta"],
        centroids_out=paths["centroids"],
        uc_out=paths["clusters_uc"],
        cluster_id=final_cluster_id,
        threads=threads,
        vsearch_bin=vsearch_bin,
    )


def rrna_step(
    *,
    barrnap_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    vsearch_bin: str,
    multi_cluster_id: float = 0.90,
    final_cluster_id: float = 1.0,
) -> None:
    """
    Count, cluster, and select representative 16S genes from Barrnap outputs.

    Expected inputs:
      barrnap_dir/bacteria/*.gff
      barrnap_dir/bacteria/*_16S.fna
      barrnap_dir/archaea/*.gff
      barrnap_dir/archaea/*_16S.fna

    Outputs:
      out/bacteria/
        single/
        multiple/
        clustered/
        bacteria_16S_copies.tsv
        bacteria_16S_genes.fasta
        bacteria_16S_genes.map.tsv
        bacteria_16S_centroids.fasta
        bacteria_16S_clusters.uc
      out/archaea/
        ...
      out/summary.tsv
      out/report.json
    """
    bac_in = barrnap_dir / "bacteria"
    arc_in = barrnap_dir / "archaea"

    if not bac_in.exists() and not arc_in.exists():
        raise FileNotFoundError(
            f"Could not find Barrnap domain directories under {barrnap_dir}. "
            f"Expected at least one of: {bac_in} or {arc_in}"
        )

    if out.exists() and any(out.iterdir()) and not force:
        raise FileExistsError(f"rrna outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    domains = {
        "bacteria": bac_in,
        "archaea": arc_in,
    }

    summary_rows: List[List[str]] = []
    report: Dict[str, object] = {
        "inputs": {
            "barrnap_dir": str(barrnap_dir),
            "vsearch_bin": vsearch_bin,
            "cpus": cpus,
            "multi_cluster_id": multi_cluster_id,
            "final_cluster_id": final_cluster_id,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    for domain, src_dir in domains.items():
        paths = _domain_paths(out, domain)
        _initialise_domain_dirs(paths)

        if not src_dir.exists():
            _LOGGER.warning("Skipping missing Barrnap domain dir: %s", src_dir)
            report["domains"][domain] = {
                "present": False,
                "source_dir": str(src_dir),
            }
            continue

        copy_rows, copy_counts = _count_and_split_domain(
            barrnap_domain_dir=src_dir,
            paths=paths,
        )

        n_clustered = _cluster_multi_copy_domain(
            paths=paths,
            vsearch_bin=vsearch_bin,
            multi_cluster_id=multi_cluster_id,
            threads=cpus,
        )

        rep_stats = _build_genome_representatives(
            domain=domain,
            paths=paths,
        )

        _final_cluster_domain(
            paths=paths,
            vsearch_bin=vsearch_bin,
            final_cluster_id=final_cluster_id,
            threads=cpus,
        )

        summary_rows.append([
            domain,
            str(copy_counts["total_fastas"]),
            str(copy_counts["with_16S"]),
            str(copy_counts["single"]),
            str(copy_counts["multiple"]),
            str(copy_counts["none"]),
            str(n_clustered),
            str(rep_stats["representatives"]),
            str(rep_stats["min_len"]),
            str(rep_stats["max_len"]),
            f"{rep_stats['mean_len']:.3f}",
        ])

        report["domains"][domain] = {
            "present": True,
            "source_dir": str(src_dir),
            "outputs": {
                "root": str(paths["root"]),
                "single_dir": str(paths["single"]),
                "multiple_dir": str(paths["multiple"]),
                "clustered_dir": str(paths["clustered"]),
                "copies_tsv": str(paths["copies_tsv"]),
                "genes_fasta": str(paths["genes_fasta"]),
                "genes_map": str(paths["genes_map"]),
                "centroids_fasta": str(paths["centroids"]),
                "clusters_uc": str(paths["clusters_uc"]),
            },
            "counts": {
                **copy_counts,
                "clustered_multi_files": n_clustered,
                **rep_stats,
            },
        }

        _LOGGER.info(
            "%s: total=%d with_16S=%d single=%d multiple=%d none=%d representatives=%d",
            domain,
            copy_counts["total_fastas"],
            copy_counts["with_16S"],
            copy_counts["single"],
            copy_counts["multiple"],
            copy_counts["none"],
            rep_stats["representatives"],
        )

    with (out / "summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "domain",
            "total_fastas",
            "with_16S",
            "single",
            "multiple",
            "none",
            "clustered_multi_files",
            "representatives",
            "min_len",
            "max_len",
            "mean_len",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")