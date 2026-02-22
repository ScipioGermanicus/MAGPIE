from __future__ import annotations

import csv
import json
import shutil
from pathlib import Path
from typing import Final

SUMMARY_FILES: Final[tuple[str, str]] = (
    "gtdbtk.bac120.summary.tsv",
    "gtdbtk.ar53.summary.tsv",
)

# File suffixes we may encounter for genomes
SUFFIXES: Final[tuple[str, ...]] = (
    "_genomic.fna.gz",
    "_genomic.fna",
    ".fa.gz",
    ".fna.gz",
    ".fasta.gz",
    ".fa",
    ".fna",
    ".fasta",
)


def _strip_suffix(name: str) -> str:
    """
    Strip known FASTA/GTDB filename suffixes, but do NOT treat dot-separated genome IDs
    (e.g. BC13.bin.1.603) as having a file extension.
    """
    s = str(name).strip()

    # Remove known suffixes (for filenames)
    for suf in SUFFIXES:
        if s.endswith(suf):
            return s[: -len(suf)]

    # If none of the known FASTA suffixes matched, assume this is already an ID.
    # Important: do NOT do Path(s).stem here, because it would truncate IDs with dots.
    return s


def _infer_domain(classification: str | None) -> str:
    if isinstance(classification, str):
        if classification.startswith("d__Bacteria"):
            return "Bacteria"
        if classification.startswith("d__Archaea"):
            return "Archaea"
    return "Unknown"


def _read_id_map(id_map_fp: Path) -> dict[str, str]:
    """
    Reads the prep id_map.tsv. We tolerate different column header casings.
    Expected in MAGPIE prep: original_filename<tab>new_id
    Returns mapping from *original stem* -> new_id (MAG0001, ...)

    Note: This is mainly useful if GTDB-Tk used original IDs. If you run GTDB-Tk
    on the prepared MAG0001_genomic.fna.gz files, user_genome will already be MAG0001-ish.
    """
    if not id_map_fp.exists():
        raise FileNotFoundError(f"Expected id_map.tsv at: {id_map_fp}")

    with id_map_fp.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"id_map.tsv has no header: {id_map_fp}")

        cols = [c.strip().lower() for c in reader.fieldnames]
        # map lowercase -> original
        colmap = {c.strip().lower(): c for c in reader.fieldnames}

        old_col = None
        new_col = None
        for c in ("original_filename", "original", "raw_id", "old_id", "user_genome"):
            if c in cols:
                old_col = colmap[c]
                break
        for c in ("new_id", "formatted_id", "mag_id"):
            if c in cols:
                new_col = colmap[c]
                break
        if old_col is None or new_col is None:
            raise ValueError(
                f"id_map.tsv must have columns like original_filename/new_id. Found: {reader.fieldnames}"
            )

        mapping: dict[str, str] = {}
        for row in reader:
            old = _strip_suffix(str(row[old_col]).strip())
            new = str(row[new_col]).strip()
            if old and new:
                mapping[old] = new
        return mapping


def _read_gtdb_summaries(classify_dir: Path) -> dict[str, str]:
    """
    Read bac120 + ar53 summary TSVs from classify_dir.
    Returns mapping user_genome -> classification (deduplicated).
    """
    result: dict[str, str] = {}
    for fn in SUMMARY_FILES:
        fp = classify_dir / fn
        if not fp.exists():
            continue
        with fp.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError(f"Empty/invalid TSV: {fp}")
            required = {"user_genome", "classification"}
            if not required.issubset(set(reader.fieldnames)):
                raise ValueError(f"{fp} missing required columns: user_genome/classification")

            for row in reader:
                ug = str(row["user_genome"]).strip()
                cl = str(row["classification"]).strip()
                if ug and ug not in result:
                    result[ug] = cl
    if not result:
        raise FileNotFoundError(f"No GTDB-Tk summaries read from: {classify_dir}")
    return result


def _find_prepped_genome(prep_mags_dir: Path, genome_id: str) -> Path | None:
    """
    Find a genome file in prep_mags_dir that corresponds to genome_id.
    We try genome_id + known suffix patterns.
    """
    candidates = []
    # Most likely in MAGPIE rename-default:
    candidates.append(prep_mags_dir / f"{genome_id}_genomic.fna.gz")
    candidates.append(prep_mags_dir / f"{genome_id}_genomic.fna")

    # Fallback generic suffixes
    for suf in SUFFIXES:
        candidates.append(prep_mags_dir / f"{genome_id}{suf}")

    for p in candidates:
        if p.exists():
            return p
    return None


def taxonomy_step(
    *,
    prep_dir: Path,
    classify_dir: Path,
    out: Path,
    force: bool,
    move_files: bool = False,
) -> None:
    """
    Infer domain taxonomy from GTDB-Tk summaries and split prepared genomes by domain.

    Inputs:
      - prep_dir: MAGPIE prep output directory (expects prep_dir/mags and prep_dir/id_map.tsv)
      - classify_dir: GTDB-Tk classify directory containing summary files
      - out: output directory for taxonomy artefacts and split genomes

    Outputs (in out/):
      - domain_map.tsv (id <tab> domain)
      - bacteria.txt / archaea.txt (IDs)
      - genomes_to_search_barrnap/bacteria/  (genome files)
      - genomes_to_search_barrnap/archaea/   (genome files)
      - taxonomy_report.json
    """
    out.mkdir(parents=True, exist_ok=True)
    report_fp = out / "taxonomy_report.json"
    if report_fp.exists() and not force:
        raise FileExistsError(f"{report_fp} exists. Use --force to overwrite.")

    prep_mags_dir = prep_dir / "mags"
    if not prep_mags_dir.exists():
        raise FileNotFoundError(f"Expected prepared MAGs at: {prep_mags_dir}")

    id_map_fp = prep_dir / "id_map.tsv"
    rename_map = _read_id_map(id_map_fp)

    # Read GTDB summaries
    gtdb = _read_gtdb_summaries(classify_dir)

    # Infer domain per user_genome, map to MAGPIE IDs when needed
    domain_by_id: dict[str, str] = {}

    unknown = 0
    for user_genome, classification in gtdb.items():
        dom = _infer_domain(classification)
        if dom == "Unknown":
            unknown += 1
            continue

        ug_norm = _strip_suffix(user_genome)

        # If GTDB-Tk ran on prepared files, ug_norm is already MAG0001 (or MAG0001_genomic -> MAG0001).
        # If it ran on originals, ug_norm may be original stem; map it to MAG id if possible.
        mag_id = rename_map.get(ug_norm, ug_norm)

        domain_by_id[mag_id] = dom

    # Prepare output dirs
    genome_dir = out  # keep outputs grouped under this step
    split_root = genome_dir / "genomes_to_search_barrnap"
    bac_dir = split_root / "bacteria"
    arc_dir = split_root / "archaea"
    bac_dir.mkdir(parents=True, exist_ok=True)
    arc_dir.mkdir(parents=True, exist_ok=True)

    # Write domain_map.tsv
    domain_map_fp = genome_dir / "domain_map.tsv"
    if domain_map_fp.exists() and not force:
        raise FileExistsError(f"{domain_map_fp} exists. Use --force to overwrite.")

    with domain_map_fp.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["genome_id", "domain"])
        for gid in sorted(domain_by_id):
            w.writerow([gid, domain_by_id[gid]])

    # Copy/move genomes into bacteria/archaea
    moved_or_copied = 0
    missing: list[str] = []

    for gid, dom in domain_by_id.items():
        src = _find_prepped_genome(prep_mags_dir, gid)
        if src is None:
            missing.append(gid)
            continue
        dst_dir = bac_dir if dom == "Bacteria" else arc_dir
        dst = dst_dir / src.name

        if dst.exists() and not force:
            raise FileExistsError(f"{dst} exists. Use --force to overwrite.")

        if move_files:
            shutil.move(str(src), str(dst))
        else:
            shutil.copy2(str(src), str(dst))
        moved_or_copied += 1

    if missing:
        raise FileNotFoundError(
            f"{len(missing)} genomes listed by GTDB-Tk were not found in prepared MAGs: {prep_mags_dir}\n"
            f"Examples: {missing[:10]}"
        )

    # Regenerate lists from actual split dirs (guarantees consistency)
    def write_id_list(dir_path: Path, out_fp: Path) -> None:
        ids = []
        for fn in sorted(p.name for p in dir_path.iterdir() if p.is_file() and not p.name.startswith(".")):
            ids.append(_strip_suffix(fn))
        out_fp.write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")

    bac_list_fp = genome_dir / "bacteria.txt"
    arc_list_fp = genome_dir / "archaea.txt"
    write_id_list(bac_dir, bac_list_fp)
    write_id_list(arc_dir, arc_list_fp)

    report = {
        "prep_dir": str(prep_dir),
        "prep_mags_dir": str(prep_mags_dir),
        "classify_dir": str(classify_dir),
        "out_dir": str(out),
        "move_files": move_files,
        "n_gtdb_genomes": len(gtdb),
        "n_domain_assigned": len(domain_by_id),
        "n_unknown": unknown,
        "outputs": {
            "domain_map_tsv": str(domain_map_fp),
            "bacteria_txt": str(bac_list_fp),
            "archaea_txt": str(arc_list_fp),
            "bacteria_dir": str(bac_dir),
            "archaea_dir": str(arc_dir),
        },
        "n_split_genomes": moved_or_copied,
    }
    report_fp.write_text(json.dumps(report, indent=2), encoding="utf-8")