from __future__ import annotations

import csv
import gzip
import json
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

_LOGGER = logging.getLogger("magpie.ec_table")

EC_FULL_RE = re.compile(r"^\d+\.\d+\.\d+\.\d+$")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _open_text_maybe_gzip(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _read_fasta_ids(path: Path) -> List[str]:
    ids: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                ids.append(raw[1:].strip().split()[0])

    if not ids:
        raise ValueError(f"No FASTA IDs found in: {path}")

    return ids


def _extract_genome_id_from_annotation_filename(path: Path) -> str:
    """
    Examples:
      MAG0001.emapper.annotations     -> MAG0001
      MAG0001.emapper.annotations.gz  -> MAG0001
    """
    name = path.name

    if name.endswith(".gz"):
        name = name[:-3]

    if ".emapper" in name:
        return name.split(".emapper", 1)[0]

    return name.split(".", 1)[0]


def _iter_annotation_files(ann_dir: Path) -> Iterable[Path]:
    patterns = [
        "*.emapper.annotations",
        "*.emapper.annotations.gz",
    ]

    seen: set[Path] = set()
    for pattern in patterns:
        for p in sorted(ann_dir.glob(pattern), key=lambda x: x.name):
            if p.is_file() and p not in seen:
                seen.add(p)
                yield p


def _split_multi_value_field(field: str) -> List[str]:
    field = field.strip()

    if not field or field == "-" or field.upper() == "NA":
        return []

    parts = re.split(r"[,\|]", field)
    return [p.strip() for p in parts if p.strip()]


def _normalise_ec(ec: str) -> str | None:
    ec = ec.strip()

    if not ec or ec == "-" or ec.upper() == "NA":
        return None

    if ec.lower().startswith("ec:"):
        ec = ec[3:].strip()

    if EC_FULL_RE.match(ec):
        return ec

    return None


def _find_ec_column(cols: List[str]) -> int | None:
    """
    Try common eggNOG EC column names.
    """
    candidates = ["EC", "EC_number", "ECs"]

    for c in candidates:
        if c in cols:
            return cols.index(c)

    return None


def _parse_one_eggnog_annotation(path: Path) -> Counter:
    """
    Parse one eggNOG-mapper .emapper.annotations file and count EC occurrences.
    """
    counts: Counter = Counter()
    ec_idx = None
    header_found = False

    with _open_text_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                text = line[1:]
                cols = text.split("\t")

                maybe_idx = _find_ec_column(cols)
                if maybe_idx is not None:
                    ec_idx = maybe_idx
                    header_found = True

                continue

            if not header_found or ec_idx is None:
                continue

            cols = line.split("\t")
            if len(cols) <= ec_idx:
                continue

            for raw_ec in _split_multi_value_field(cols[ec_idx]):
                ec = _normalise_ec(raw_ec)
                if ec is not None:
                    counts[ec] += 1

    if not header_found:
        raise ValueError(f"No eggNOG header containing an EC column found in: {path}")

    return counts


def _parse_eggnog_dir(ann_dir: Path) -> Dict[str, Counter]:
    """
    Return:
      genome_id -> Counter(EC -> copy_number)
    """
    ann_files = list(_iter_annotation_files(ann_dir))

    if not ann_files:
        raise FileNotFoundError(f"No *.emapper.annotations files found in: {ann_dir}")

    out: Dict[str, Counter] = {}

    for ann in ann_files:
        gid = _extract_genome_id_from_annotation_filename(ann)

        if gid in out:
            raise ValueError(
                f"Duplicate eggNOG annotation prefix '{gid}' in {ann_dir}. "
                f"MAGPIE requires one annotation file per genome ID."
            )

        _LOGGER.info("Parsing eggNOG EC annotations: %s -> %s", ann, gid)
        out[gid] = _parse_one_eggnog_annotation(ann)

    return out


def _write_ec_table(
    *,
    genome_ids: List[str],
    all_counts: Dict[str, Counter],
    out_path: Path,
    allow_missing_annotations: bool,
) -> Tuple[int, int, int]:
    """
    Write PICRUSt2-style EC table:
      assembly    1.1.1.1    1.1.1.2 ...

    Returns:
      n_genomes_written, n_ecs, n_missing_annotations
    """
    missing = [gid for gid in genome_ids if gid not in all_counts]

    if missing and not allow_missing_annotations:
        preview = ", ".join(missing[:10])
        raise ValueError(
            f"{len(missing)} genomes in reference FASTA had no matching eggNOG annotation file. "
            f"First missing IDs: {preview}. "
            f"Use --allow-missing-eggnog if you intentionally want zero-filled rows."
        )

    ecs: Set[str] = set()
    for gid in genome_ids:
        if gid in all_counts:
            ecs.update(all_counts[gid].keys())

    sorted_ecs = sorted(ecs)

    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        header = ["assembly"] + sorted_ecs
        out.write("\t".join(header) + "\n")

        for gid in genome_ids:
            counter = all_counts.get(gid, Counter())
            row = [gid] + [str(counter.get(ec, 0)) for ec in sorted_ecs]
            out.write("\t".join(row) + "\n")

    return len(genome_ids), len(sorted_ecs), len(missing)


def ec_table_step(
    *,
    package_ref_dir: Path,
    eggnog_existing_dir: Path,
    out: Path,
    force: bool,
    allow_missing_eggnog: bool,
) -> None:
    """
    Build PICRUSt2-style EC trait tables from existing eggNOG annotations.

    Inputs:
      package_ref_dir/bac_ref/bac_ref.fna
      package_ref_dir/arc_ref/arc_ref.fna
      eggnog_existing_dir/*.emapper.annotations

    Outputs:
      package_ref_dir/bac_ref/ec.txt.gz
      package_ref_dir/arc_ref/ec.txt.gz
      out/summary.tsv
      out/report.json
    """
    bac_fasta = package_ref_dir / "bac_ref" / "bac_ref.fna"
    arc_fasta = package_ref_dir / "arc_ref" / "arc_ref.fna"

    bac_out = package_ref_dir / "bac_ref" / "ec.txt.gz"
    arc_out = package_ref_dir / "arc_ref" / "ec.txt.gz"

    for p in [bac_fasta, arc_fasta]:
        if not p.exists():
            raise FileNotFoundError(f"Required reference FASTA missing: {p}")

    if not eggnog_existing_dir.exists():
        raise FileNotFoundError(f"eggNOG annotation directory not found: {eggnog_existing_dir}")

    key_outputs = [out / "summary.tsv", out / "report.json", bac_out, arc_out]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(
            f"EC table outputs already exist. Use --force to overwrite.\n"
            f"Checked: {', '.join(str(p) for p in key_outputs)}"
        )

    _ensure_dir(out)

    bac_ids = _read_fasta_ids(bac_fasta)
    arc_ids = _read_fasta_ids(arc_fasta)

    overlap = set(bac_ids) & set(arc_ids)
    if overlap:
        preview = ", ".join(sorted(overlap)[:10])
        raise ValueError(
            f"Genome IDs overlap between bacterial and archaeal references "
            f"({len(overlap)} IDs). First overlaps: {preview}"
        )

    all_counts = _parse_eggnog_dir(eggnog_existing_dir)

    bac_n, bac_ecs, bac_missing = _write_ec_table(
        genome_ids=bac_ids,
        all_counts=all_counts,
        out_path=bac_out,
        allow_missing_annotations=allow_missing_eggnog,
    )

    arc_n, arc_ecs, arc_missing = _write_ec_table(
        genome_ids=arc_ids,
        all_counts=all_counts,
        out_path=arc_out,
        allow_missing_annotations=allow_missing_eggnog,
    )

    rows = [
        ["bacteria", str(bac_n), str(bac_ecs), str(bac_missing), str(bac_fasta), str(bac_out)],
        ["archaea", str(arc_n), str(arc_ecs), str(arc_missing), str(arc_fasta), str(arc_out)],
    ]

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "domain",
            "n_genomes",
            "n_ecs",
            "n_missing_annotations",
            "reference_fasta",
            "ec_table",
        ])
        w.writerows(rows)

    report = {
        "inputs": {
            "package_ref_dir": str(package_ref_dir),
            "eggnog_existing_dir": str(eggnog_existing_dir),
            "allow_missing_eggnog": allow_missing_eggnog,
        },
        "outputs": {
            "bacteria_ec": str(bac_out),
            "archaea_ec": str(arc_out),
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
        "counts": {
            "annotation_files_parsed": len(all_counts),
            "bacteria_genomes": bac_n,
            "archaea_genomes": arc_n,
            "bacteria_ecs": bac_ecs,
            "archaea_ecs": arc_ecs,
            "bacteria_missing_annotations": bac_missing,
            "archaea_missing_annotations": arc_missing,
        },
    }

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info(
        "EC table step complete. Bacteria: %d genomes, %d ECs. Archaea: %d genomes, %d ECs.",
        bac_n,
        bac_ecs,
        arc_n,
        arc_ecs,
    )