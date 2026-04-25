from __future__ import annotations

import csv
import gzip
import json
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

_LOGGER = logging.getLogger("magpie.ko_table")


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


def _split_kegg_ko_field(kfield: str) -> List[str]:
    if not kfield:
        return []

    kfield = kfield.strip()
    if not kfield or kfield == "-" or kfield.upper() == "NA":
        return []

    if "," in kfield:
        raw = kfield.split(",")
    elif "|" in kfield:
        raw = kfield.split("|")
    else:
        raw = [kfield]

    out: List[str] = []

    for ko in raw:
        ko = ko.strip()
        if not ko or ko == "-" or ko.upper() == "NA":
            continue

        if ko.startswith("ko:"):
            ko = ko[3:]

        if ko.startswith("K") and len(ko) == 6 and ko[1:].isdigit():
            out.append(ko)

    return out


def _parse_one_eggnog_annotation(path: Path) -> Counter:
    """
    Parse one eggNOG-mapper .emapper.annotations file and count KO occurrences.
    """
    counts: Counter = Counter()
    ko_idx = None
    header_found = False

    with _open_text_maybe_gzip(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            if not line:
                continue

            if line.startswith("#"):
                text = line[1:]
                cols = text.split("\t")

                if "KEGG_ko" in cols:
                    ko_idx = cols.index("KEGG_ko")
                    header_found = True

                continue

            if not header_found or ko_idx is None:
                continue

            cols = line.split("\t")
            if len(cols) <= ko_idx:
                continue

            for ko in _split_kegg_ko_field(cols[ko_idx]):
                counts[ko] += 1

    if not header_found:
        raise ValueError(f"No eggNOG header containing 'KEGG_ko' found in: {path}")

    return counts


def _parse_eggnog_dir(ann_dir: Path) -> Dict[str, Counter]:
    """
    Return:
      genome_id -> Counter(KO -> copy_number)
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

        _LOGGER.info("Parsing eggNOG annotation: %s -> %s", ann, gid)
        out[gid] = _parse_one_eggnog_annotation(ann)

    return out


def _write_ko_table(
    *,
    genome_ids: List[str],
    all_counts: Dict[str, Counter],
    out_path: Path,
    allow_missing_annotations: bool,
) -> Tuple[int, int, int]:
    """
    Write PICRUSt2-style KO table:
      assembly    ko:K00001    ko:K00002 ...

    Returns:
      n_genomes_written, n_kos, n_missing_annotations
    """
    missing = [gid for gid in genome_ids if gid not in all_counts]

    if missing and not allow_missing_annotations:
        preview = ", ".join(missing[:10])
        raise ValueError(
            f"{len(missing)} genomes in reference FASTA had no matching eggNOG annotation file. "
            f"First missing IDs: {preview}. "
            f"Use --allow-missing-eggnog if you intentionally want zero-filled rows."
        )

    kos: Set[str] = set()
    for gid in genome_ids:
        if gid in all_counts:
            kos.update(all_counts[gid].keys())

    sorted_kos = sorted(kos)

    with gzip.open(out_path, "wt", encoding="utf-8") as out:
        header = ["assembly"] + [f"ko:{ko}" for ko in sorted_kos]
        out.write("\t".join(header) + "\n")

        for gid in genome_ids:
            counter = all_counts.get(gid, Counter())
            row = [gid] + [str(counter.get(ko, 0)) for ko in sorted_kos]
            out.write("\t".join(row) + "\n")

    return len(genome_ids), len(sorted_kos), len(missing)


def ko_table_step(
    *,
    package_ref_dir: Path,
    eggnog_existing_dir: Path,
    out: Path,
    force: bool,
    allow_missing_eggnog: bool,
) -> None:
    """
    Build PICRUSt2-style KO trait tables from existing eggNOG annotations.

    Inputs:
      package_ref_dir/bac_ref/bac_ref.fna
      package_ref_dir/arc_ref/arc_ref.fna
      eggnog_existing_dir/*.emapper.annotations

    Outputs:
      package_ref_dir/bac_ref/ko.txt.gz
      package_ref_dir/arc_ref/ko.txt.gz
      out/summary.tsv
      out/report.json
    """
    bac_fasta = package_ref_dir / "bac_ref" / "bac_ref.fna"
    arc_fasta = package_ref_dir / "arc_ref" / "arc_ref.fna"

    bac_out = package_ref_dir / "bac_ref" / "ko.txt.gz"
    arc_out = package_ref_dir / "arc_ref" / "ko.txt.gz"

    for p in [bac_fasta, arc_fasta]:
        if not p.exists():
            raise FileNotFoundError(f"Required reference FASTA missing: {p}")

    if not eggnog_existing_dir.exists():
        raise FileNotFoundError(f"eggNOG annotation directory not found: {eggnog_existing_dir}")

    key_outputs = [out / "summary.tsv", out / "report.json", bac_out, arc_out]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(
            f"KO table outputs already exist. Use --force to overwrite.\n"
            f"Checked: {', '.join(str(p) for p in key_outputs)}"
        )

    _ensure_dir(out)

    _LOGGER.info("Loading final bacterial reference IDs from: %s", bac_fasta)
    bac_ids = _read_fasta_ids(bac_fasta)

    _LOGGER.info("Loading final archaeal reference IDs from: %s", arc_fasta)
    arc_ids = _read_fasta_ids(arc_fasta)

    overlap = set(bac_ids) & set(arc_ids)
    if overlap:
        preview = ", ".join(sorted(overlap)[:10])
        raise ValueError(
            f"Genome IDs overlap between bacterial and archaeal references "
            f"({len(overlap)} IDs). First overlaps: {preview}"
        )

    _LOGGER.info("Parsing eggNOG annotations from: %s", eggnog_existing_dir)
    all_counts = _parse_eggnog_dir(eggnog_existing_dir)

    _LOGGER.info("Writing bacterial KO table: %s", bac_out)
    bac_n, bac_kos, bac_missing = _write_ko_table(
        genome_ids=bac_ids,
        all_counts=all_counts,
        out_path=bac_out,
        allow_missing_annotations=allow_missing_eggnog,
    )

    _LOGGER.info("Writing archaeal KO table: %s", arc_out)
    arc_n, arc_kos, arc_missing = _write_ko_table(
        genome_ids=arc_ids,
        all_counts=all_counts,
        out_path=arc_out,
        allow_missing_annotations=allow_missing_eggnog,
    )

    rows = [
        ["bacteria", str(bac_n), str(bac_kos), str(bac_missing), str(bac_fasta), str(bac_out)],
        ["archaea", str(arc_n), str(arc_kos), str(arc_missing), str(arc_fasta), str(arc_out)],
    ]

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "domain",
            "n_genomes",
            "n_kos",
            "n_missing_annotations",
            "reference_fasta",
            "ko_table",
        ])
        w.writerows(rows)

    report = {
        "inputs": {
            "package_ref_dir": str(package_ref_dir),
            "eggnog_existing_dir": str(eggnog_existing_dir),
            "allow_missing_eggnog": allow_missing_eggnog,
        },
        "outputs": {
            "bacteria_ko": str(bac_out),
            "archaea_ko": str(arc_out),
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
        "counts": {
            "annotation_files_parsed": len(all_counts),
            "bacteria_genomes": bac_n,
            "archaea_genomes": arc_n,
            "bacteria_kos": bac_kos,
            "archaea_kos": arc_kos,
            "bacteria_missing_annotations": bac_missing,
            "archaea_missing_annotations": arc_missing,
        },
    }

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info(
        "KO table step complete. Bacteria: %d genomes, %d KOs. Archaea: %d genomes, %d KOs.",
        bac_n,
        bac_kos,
        arc_n,
        arc_kos,
    )