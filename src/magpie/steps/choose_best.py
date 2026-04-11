from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

_LOGGER = logging.getLogger("magpie.choose_best")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _has_fasta_records(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                return True
    return False


def _read_id_map(path: Path) -> Dict[str, str]:
    mp: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            return mp

        cols = {c.lower(): c for c in r.fieldnames}
        of = cols.get("original_filename") or cols.get("original") or r.fieldnames[0]
        ni = cols.get("new_id") or cols.get("id") or r.fieldnames[1]

        for row in r:
            old = row.get(of, "").strip()
            new = row.get(ni, "").strip()
            if old and new:
                mp[old] = new
    return mp


def _strip_genome_suffix(genome_id: str) -> str:
    x = genome_id.strip()
    if x.endswith("_genomic"):
        x = x[:-8]
    return x


def _read_metadata(
    *,
    metadata_tsv: Path,
    domain: str,
    id_map_tsv: Path | None,
) -> Dict[str, Dict[str, str]]:
    """
    Read QC metadata keyed by harmonised genome ID.
    """
    id_map = _read_id_map(id_map_tsv) if (id_map_tsv is not None and id_map_tsv.exists()) else {}

    rows: Dict[str, Dict[str, str]] = {}

    with metadata_tsv.open("r", encoding="utf-8", newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        if not r.fieldnames:
            return rows

        first_col = r.fieldnames[0]

        for row in r:
            raw_id = row.get(first_col, "").strip()
            if not raw_id:
                continue

            if "domain" in row and row["domain"].strip():
                if row["domain"].strip().lower() != domain.lower():
                    continue

            mapped = id_map.get(raw_id, raw_id)
            mapped = _strip_genome_suffix(mapped)

            rows[mapped] = row

    return rows


def _parse_clusters(path: Path, genes_in_alignment: set[str]) -> Dict[str, List[str]]:
    """
    Parse vsearch .uc file into centroid -> members list.
    """
    clusters: Dict[str, List[str]] = defaultdict(list)

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            if not raw.strip() or raw.startswith("#"):
                continue

            line = raw.rstrip("\n")
            parts = line.split("\t")
            rec_type = parts[0] if parts else ""

            if rec_type in {"S", "C"}:
                q = parts[8] if len(parts) > 8 else None
                if q and q in genes_in_alignment:
                    clusters.setdefault(q, [])
                continue

            if rec_type == "H":
                q = parts[8] if len(parts) > 8 else None
                t = parts[9] if len(parts) > 9 else None
                if t in genes_in_alignment:
                    clusters.setdefault(t, [])
                    if q is not None:
                        clusters[t].append(q)
                continue

            if len(parts) >= 2 and parts[-1] != "*":
                centroid, member = parts[-2], parts[-1]
                if centroid in genes_in_alignment:
                    clusters.setdefault(centroid, [])
                    clusters[centroid].append(member)

    for gid in genes_in_alignment:
        clusters.setdefault(gid, [])

    return clusters


def _to_float_or_none(value: str | None) -> float | None:
    if value is None:
        return None
    s = value.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def _choose_best(
    metadata: Dict[str, Dict[str, str]],
    candidates: List[str],
) -> str:
    """
    Rank by:
      1. completeness descending
      2. contamination ascending
      3. genome_id ascending (stable fallback)
    """
    unique_candidates: List[str] = []
    for c in candidates:
        if c not in unique_candidates:
            unique_candidates.append(c)

    avail = [g for g in unique_candidates if g in metadata]
    if not avail:
        return sorted(unique_candidates)[0]

    def sort_key(genome_id: str) -> Tuple[float, float, str]:
        row = metadata[genome_id]
        comp = _to_float_or_none(row.get("checkm_completeness"))
        cont = _to_float_or_none(row.get("checkm_contamination"))

        comp_key = comp if comp is not None else float("-inf")
        cont_key = cont if cont is not None else float("inf")

        return (-comp_key, cont_key, genome_id)

    return sorted(avail, key=sort_key)[0]


def _choose_best_one_domain(
    *,
    domain: str,
    clusters_uc: Path,
    aligned_fasta: Path,
    metadata_tsv: Path,
    id_map_tsv: Path | None,
    outdir: Path,
) -> Dict[str, object]:
    _ensure_dir(outdir)

    out_fa = outdir / f"{domain}_16S_centroids_ssu_align_best.fna"
    proc_path = outdir / f"{domain}_16S_clusters_processed.tsv"
    md_out = outdir / f"{domain}_metadata_clusters_ssu_align_centroids.tsv"

    if not aligned_fasta.exists():
        raise FileNotFoundError(f"Aligned FASTA not found for {domain}: {aligned_fasta}")

    if not clusters_uc.exists():
        raise FileNotFoundError(f"Cluster file not found for {domain}: {clusters_uc}")

    if not _has_fasta_records(aligned_fasta):
        _LOGGER.warning("Skipping %s choose-best: aligned FASTA contains no sequences", domain)

        out_fa.write_text("", encoding="utf-8")
        with proc_path.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["centroid", "best", "all_genomes"])
        with md_out.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow([])

        return {
            "domain": domain,
            "status": "skipped_empty_alignment",
            "aligned_fasta": str(aligned_fasta),
            "clusters_uc": str(clusters_uc),
            "best_fasta": str(out_fa),
            "processed_clusters_tsv": str(proc_path),
            "metadata_tsv": str(md_out),
            "n_alignment_records": 0,
            "n_output_records": 0,
        }

    metadata = _read_metadata(
        metadata_tsv=metadata_tsv,
        domain=domain,
        id_map_tsv=id_map_tsv,
    )

    records = list(SeqIO.parse(str(aligned_fasta), "fasta"))
    for r in records:
        r.id = _strip_genome_suffix(r.id)
        r.name = r.id
        r.description = ""

    genes_16s = [r.id for r in records]
    genes_set = set(genes_16s)

    cluster_map = _parse_clusters(clusters_uc, genes_set)

    best_map: Dict[str, str] = {}
    processed_rows: List[Tuple[str, str, str]] = []

    for centroid, members in cluster_map.items():
        candidates = [centroid] + [m for m in members if m != centroid]

        ordered: List[str] = []
        for c in candidates:
            if c not in ordered:
                ordered.append(c)

        best = _choose_best(metadata, ordered)
        best_map[centroid] = best
        processed_rows.append((centroid, best, ",".join(ordered)))

    new_records: List[SeqRecord] = []
    out_ids: List[str] = []
    for r in records:
        new_id = best_map.get(r.id, r.id)
        r.id = new_id
        r.name = new_id
        r.description = ""
        new_records.append(r)
        out_ids.append(new_id)

    SeqIO.write(new_records, str(out_fa), "fasta")

    with proc_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["centroid", "best", "all_genomes"])
        for row in processed_rows:
            w.writerow(row)

    chosen_ids = sorted(set(out_ids))
    if metadata:
        fieldnames = list(next(iter(metadata.values())).keys())
        genome_col = fieldnames[0] if fieldnames else "genome_id"

        with md_out.open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for gid in chosen_ids:
                if gid in metadata:
                    row = dict(metadata[gid])
                    row[genome_col] = gid
                    w.writerow(row)
    else:
        with md_out.open("w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["genome_id"])

    return {
        "domain": domain,
        "status": "processed",
        "aligned_fasta": str(aligned_fasta),
        "clusters_uc": str(clusters_uc),
        "best_fasta": str(out_fa),
        "processed_clusters_tsv": str(proc_path),
        "metadata_tsv": str(md_out),
        "n_alignment_records": len(records),
        "n_output_records": len(new_records),
        "n_unique_output_genomes": len(set(out_ids)),
    }


def choose_best_step(
    *,
    prep_dir: Path,
    qc_dir: Path,
    rrna_dir: Path,
    align_dir: Path,
    out: Path,
    force: bool,
) -> None:
    """
    Choose the best genome per 16S cluster for archaea and bacteria.

    Inputs:
      prep_dir/id_map.tsv
      qc_dir/checkm_filtered.tsv
      rrna_dir/{domain}/{domain}_16S_clusters.uc
      align_dir/{domain}_16S_centroids_ssu_align.fna
    """
    id_map_tsv = prep_dir / "id_map.tsv"
    metadata_tsv = qc_dir / "checkm_filtered.tsv"

    if not metadata_tsv.exists():
        raise FileNotFoundError(f"Missing QC metadata: {metadata_tsv}")

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"Choose-best outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    domains = {
        "archaea": {
            "clusters_uc": rrna_dir / "archaea" / "archaea_16S_clusters.uc",
            "aligned_fasta": align_dir / "archaea_16S_centroids_ssu_align.fna",
        },
        "bacteria": {
            "clusters_uc": rrna_dir / "bacteria" / "bacteria_16S_clusters.uc",
            "aligned_fasta": align_dir / "bacteria_16S_centroids_ssu_align.fna",
        },
    }

    report = {
        "inputs": {
            "prep_dir": str(prep_dir),
            "qc_dir": str(qc_dir),
            "rrna_dir": str(rrna_dir),
            "align_dir": str(align_dir),
            "id_map_tsv": str(id_map_tsv),
            "metadata_tsv": str(metadata_tsv),
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, paths in domains.items():
        domain_out = out / domain

        stats = _choose_best_one_domain(
            domain=domain,
            clusters_uc=paths["clusters_uc"],
            aligned_fasta=paths["aligned_fasta"],
            metadata_tsv=metadata_tsv,
            id_map_tsv=id_map_tsv if id_map_tsv.exists() else None,
            outdir=domain_out,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats.get("n_alignment_records", 0)),
            str(stats.get("n_output_records", 0)),
            str(stats.get("n_unique_output_genomes", 0)),
            str(stats["best_fasta"]),
            str(stats["processed_clusters_tsv"]),
            str(stats["metadata_tsv"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_alignment_records",
            "n_output_records",
            "n_unique_output_genomes",
            "best_fasta",
            "processed_clusters_tsv",
            "metadata_tsv",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("Choose-best complete.")