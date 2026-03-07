from __future__ import annotations

import csv
import gzip
import json
import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

_LOGGER = logging.getLogger("magpie.barrnap")

FASTA_EXTS = (".fa", ".fna", ".fasta")
GZ_EXT = ".gz"


@dataclass(frozen=True)
class RrnaHit:
    seqid: str
    start: int
    end: int
    strand: str
    product: str
    attributes: str


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _strip_fasta_suffix(name: str) -> str:
    s = name.strip()
    if s.lower().endswith(GZ_EXT):
        s = s[:-len(GZ_EXT)]
    for ext in FASTA_EXTS:
        if s.lower().endswith(ext):
            s = s[:-len(ext)]
            break
    return s


def _iter_genome_files(dir_: Path) -> Iterable[Path]:
    if not dir_.exists():
        return
    for p in dir_.iterdir():
        if not p.is_file():
            continue
        low = p.name.lower()
        if any(low.endswith(ext) for ext in FASTA_EXTS) or any(low.endswith(ext + GZ_EXT) for ext in FASTA_EXTS):
            yield p


def _read_text_fasta(path: Path) -> str:
    if path.name.lower().endswith(GZ_EXT):
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            return fh.read()
    return path.read_text(encoding="utf-8", errors="replace")


def _read_fasta_dict(path: Path) -> Dict[str, str]:
    """
    Read FASTA into dict keyed by first token of header.
    """
    text = _read_text_fasta(path)
    seqs: Dict[str, List[str]] = {}
    current = None

    for line in text.splitlines():
        if not line:
            continue
        if line.startswith(">"):
            current = line[1:].strip().split()[0]
            if current in seqs:
                raise ValueError(f"Duplicate FASTA header '{current}' in {path}")
            seqs[current] = []
        else:
            if current is None:
                raise ValueError(f"Malformed FASTA without header in {path}")
            seqs[current].append(line.strip())

    return {k: "".join(v).upper() for k, v in seqs.items()}


def _revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def _parse_gff_16s(gff_path: Path) -> List[RrnaHit]:
    hits: List[RrnaHit] = []

    with gff_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) != 9:
                continue

            seqid, source, feature, start, end, score, strand, phase, attrs = parts

            # Barrnap usually labels 16S in attributes, e.g. Name=16S_rRNA;product=16S ribosomal RNA
            attrs_low = attrs.lower()
            feature_low = feature.lower()

            is_rrna = feature_low == "rrna" or "rrna" in attrs_low
            is_16s = "16s" in attrs_low or "16s_rRNA".lower() in attrs_low

            if not (is_rrna and is_16s):
                continue

            product = ""
            for field in attrs.split(";"):
                if field.startswith("product="):
                    product = field.split("=", 1)[1]
                    break
                if field.startswith("Name="):
                    product = field.split("=", 1)[1]

            hits.append(
                RrnaHit(
                    seqid=seqid,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    product=product,
                    attributes=attrs,
                )
            )

    return hits


def _extract_hits_to_fasta(
    *,
    genome_path: Path,
    genome_id: str,
    hits: List[RrnaHit],
    out_fasta: Path,
) -> int:
    seqs = _read_fasta_dict(genome_path)

    with out_fasta.open("w", encoding="utf-8") as out:
        for i, hit in enumerate(hits, start=1):
            if hit.seqid not in seqs:
                raise ValueError(
                    f"GFF seqid '{hit.seqid}' not found in FASTA for genome '{genome_id}' ({genome_path})"
                )

            full = seqs[hit.seqid]
            subseq = full[hit.start - 1:hit.end]
            if hit.strand == "-":
                subseq = _revcomp(subseq)

            header = f">{genome_id}|{hit.seqid}:{hit.start}-{hit.end}({hit.strand})|16S_{i}"
            out.write(header + "\n")
            out.write(subseq + "\n")

    return len(hits)


def _run_barrnap_one(
    *,
    genome_path: Path,
    kingdom: str,
    threads: int,
    barrnap_bin: str,
    reject: float,
    gff_out: Path,
) -> None:
    """
    Run barrnap and write GFF to gff_out.
    Supports plain or gzipped FASTA input.
    """
    if genome_path.name.lower().endswith(GZ_EXT):
        with gzip.open(genome_path, "rt", encoding="utf-8", errors="replace") as fin, \
             gff_out.open("w", encoding="utf-8") as fout:
            proc = subprocess.run(
                [
                    barrnap_bin,
                    "--kingdom", kingdom,
                    "--threads", str(threads),
                    "--reject", str(reject),
                ],
                stdin=fin,
                stdout=fout,
                stderr=subprocess.PIPE,
                text=True,
            )
    else:
        with genome_path.open("r", encoding="utf-8", errors="replace") as fin, \
             gff_out.open("w", encoding="utf-8") as fout:
            proc = subprocess.run(
                [
                    barrnap_bin,
                    "--kingdom", kingdom,
                    "--threads", str(threads),
                    "--reject", str(reject),
                ],
                stdin=fin,
                stdout=fout,
                stderr=subprocess.PIPE,
                text=True,
            )

    if proc.returncode != 0:
        raise RuntimeError(
            f"Barrnap failed for {genome_path}\n"
            f"stderr:\n{proc.stderr}"
        )


def _copy_into_bucket(src: Path, dst_dir: Path) -> None:
    _ensure_dir(dst_dir)
    shutil.copy2(src, dst_dir / src.name)


def barrnap_step(
    *,
    qc_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    barrnap_bin: str,
    reject: float = 0.8,
) -> None:
    """
    Predict 16S genes from QC-filtered MAG bins.

    Expected inputs:
      qc_dir/bins/bacteria
      qc_dir/bins/archaea
    """
    bins_root = qc_dir / "bins"
    bac_in = bins_root / "bacteria"
    arc_in = bins_root / "archaea"

    if not bac_in.exists() and not arc_in.exists():
        raise FileNotFoundError(
            f"Could not find QC bins under {bins_root}. "
            f"Expected at least one of: {bac_in} or {arc_in}"
        )

    out_bac = out / "bacteria"
    out_arc = out / "archaea"
    out_bac_single = out / "bacteria_16S_single"
    out_arc_single = out / "archaea_16S_single"
    out_bac_multi = out / "bacteria_16S_multiple"
    out_arc_multi = out / "archaea_16S_multiple"

    for d in (out, out_bac, out_arc, out_bac_single, out_arc_single, out_bac_multi, out_arc_multi):
        _ensure_dir(d)

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"Barrnap outputs already exist under {out}. Use --force to overwrite.")

    summary_rows: List[List[str]] = []

    def process_domain(domain: str, kingdom: str, src_dir: Path, dst_dir: Path, single_dir: Path, multi_dir: Path) -> None:
        if not src_dir.exists():
            _LOGGER.warning("Skipping missing domain dir: %s", src_dir)
            return

        for genome_fp in sorted(_iter_genome_files(src_dir), key=lambda p: p.name):
            genome_id = _strip_fasta_suffix(genome_fp.name)

            gff_out = dst_dir / f"{genome_id}.gff"
            fa_out = dst_dir / f"{genome_id}_16S.fna"

            _run_barrnap_one(
                genome_path=genome_fp,
                kingdom=kingdom,
                threads=cpus,
                barrnap_bin=barrnap_bin,
                reject=reject,
                gff_out=gff_out,
            )

            hits = _parse_gff_16s(gff_out)

            if hits:
                n_hits = _extract_hits_to_fasta(
                    genome_path=genome_fp,
                    genome_id=genome_id,
                    hits=hits,
                    out_fasta=fa_out,
                )
            else:
                fa_out.write_text("", encoding="utf-8")
                n_hits = 0

            category = "none"
            if n_hits == 1:
                category = "single"
                _copy_into_bucket(fa_out, single_dir)
            elif n_hits > 1:
                category = "multiple"
                _copy_into_bucket(fa_out, multi_dir)

            summary_rows.append([
                genome_id,
                domain,
                str(genome_fp),
                str(gff_out),
                str(fa_out),
                str(n_hits),
                category,
            ])

    process_domain("Bacteria", "bac", bac_in, out_bac, out_bac_single, out_bac_multi)
    process_domain("Archaea", "arc", arc_in, out_arc, out_arc_single, out_arc_multi)

    with (out / "summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "genome_id",
            "domain",
            "genome_path",
            "gff_path",
            "fasta_16S_path",
            "n_16S_hits",
            "category",
        ])
        w.writerows(summary_rows)

    bac_ids = sorted(row[0] for row in summary_rows if row[1] == "Bacteria" and row[5] != "0")
    arc_ids = sorted(row[0] for row in summary_rows if row[1] == "Archaea" and row[5] != "0")

    (out / "bacteria.txt").write_text("\n".join(bac_ids) + ("\n" if bac_ids else ""), encoding="utf-8")
    (out / "archaea.txt").write_text("\n".join(arc_ids) + ("\n" if arc_ids else ""), encoding="utf-8")

    report = {
        "inputs": {
            "qc_dir": str(qc_dir),
            "bacteria_bins": str(bac_in),
            "archaea_bins": str(arc_in),
            "barrnap_bin": barrnap_bin,
            "reject": reject,
            "cpus": cpus,
        },
        "counts": {
            "total_genomes": len(summary_rows),
            "bacteria_genomes": sum(1 for r in summary_rows if r[1] == "Bacteria"),
            "archaea_genomes": sum(1 for r in summary_rows if r[1] == "Archaea"),
            "with_16S": sum(1 for r in summary_rows if r[5] != "0"),
            "single_16S": sum(1 for r in summary_rows if r[6] == "single"),
            "multiple_16S": sum(1 for r in summary_rows if r[6] == "multiple"),
            "no_16S": sum(1 for r in summary_rows if r[6] == "none"),
        },
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "bacteria_dir": str(out_bac),
            "archaea_dir": str(out_arc),
            "bacteria_16S_single": str(out_bac_single),
            "archaea_16S_single": str(out_arc_single),
            "bacteria_16S_multiple": str(out_bac_multi),
            "archaea_16S_multiple": str(out_arc_multi),
        },
    }
    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info(
        "Barrnap complete. Total=%d, with_16S=%d, single=%d, multiple=%d, none=%d",
        len(summary_rows),
        sum(1 for r in summary_rows if r[5] != "0"),
        sum(1 for r in summary_rows if r[6] == "single"),
        sum(1 for r in summary_rows if r[6] == "multiple"),
        sum(1 for r in summary_rows if r[6] == "none"),
    )