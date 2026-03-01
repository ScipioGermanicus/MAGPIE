# src/magpie/steps/qc.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import json
import re
import shutil
from typing import Dict, Iterable, List, Optional, Tuple

from magpie.util.shell import run as shell_run


SEP_RE = re.compile(r"^\s*-{5,}\s*$")
SPLIT_RE = re.compile(r"\s{2,}")  # CheckM QA tables align with >=2 spaces


@dataclass(frozen=True)
class CheckMRow:
    genome_id: str
    marker_lineage: str
    completeness: float
    contamination: float


def _read_domain_map(domain_map_tsv: Path) -> Dict[str, str]:
    """
    domain_map.tsv written by taxonomy step.
    Expected columns: genome_id \t domain
    (If yours uses different headers, update here.)
    """
    domain_by_id: Dict[str, str] = {}
    with domain_map_tsv.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # try common keys
        genome_key = "genome_id" if "genome_id" in reader.fieldnames else "Genome"
        domain_key = "domain" if "domain" in reader.fieldnames else "Domain"
        if genome_key not in reader.fieldnames or domain_key not in reader.fieldnames:
            raise ValueError(f"domain_map.tsv missing expected columns. Found: {reader.fieldnames}")
        for row in reader:
            gid = (row.get(genome_key) or "").strip()
            dom = (row.get(domain_key) or "").strip()
            if gid:
                domain_by_id[gid] = dom
    return domain_by_id


def parse_checkm_qa_table(tsv_like: Path) -> List[CheckMRow]:
    """
    Parse CheckM QA (-o 2) wide table (space-aligned).
    Works with the style you pasted (Bin Id, Marker lineage, ..., Completeness, Contamination, ...).
    """
    header_cols: Optional[List[str]] = None
    idx_comp = idx_cont = None

    rows: List[CheckMRow] = []

    with tsv_like.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            s = line.strip()
            if not s:
                continue
            if SEP_RE.match(s):
                continue

            # Header line
            if s.startswith("Bin Id"):
                header_cols = SPLIT_RE.split(s)
                # locate indices defensively
                try:
                    idx_comp = header_cols.index("Completeness")
                    idx_cont = header_cols.index("Contamination")
                except ValueError as e:
                    raise ValueError(f"Could not locate Completeness/Contamination in header: {header_cols}") from e
                continue

            # Data lines should start with a genome id token
            if header_cols is None or idx_comp is None or idx_cont is None:
                # Not yet seen header; skip until we do
                continue

            parts = SPLIT_RE.split(s)
            # parts should have at least up to contamination index
            if len(parts) <= max(idx_comp, idx_cont):
                continue

            genome_id = parts[0].strip()
            marker_lineage = parts[1].strip() if len(parts) > 1 else ""

            try:
                completeness = float(parts[idx_comp])
                contamination = float(parts[idx_cont])
            except ValueError:
                # If the row is malformed, skip it rather than crashing
                continue

            rows.append(CheckMRow(
                genome_id=genome_id,
                marker_lineage=marker_lineage,
                completeness=completeness,
                contamination=contamination
            ))

    # de-dup by genome_id (keep first)
    seen = set()
    dedup: List[CheckMRow] = []
    for r in rows:
        if r.genome_id in seen:
            continue
        seen.add(r.genome_id)
        dedup.append(r)

    if not dedup:
        raise ValueError(f"No rows parsed from CheckM QA table: {tsv_like}")

    return dedup


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def place_bins(
    kept_ids: List[str],
    src_dir: Path,
    dst_dir: Path,
    mode: str = "symlink",  # symlink | copy | hardlink
) -> Tuple[int, List[str]]:
    """
    Place genomes by EXACT filename stem match inside src_dir.
    No fuzzy matching. If a genome isn't found, it is reported missing.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)

    # Build an index of exact stems in src_dir (non-recursive, because taxonomy already split)
    index: Dict[str, Path] = {}
    for p in src_dir.iterdir():
        if p.is_file():
            index[p.stem] = p

    placed = 0
    missing: List[str] = []

    for gid in kept_ids:
        src = index.get(gid)
        if src is None:
            missing.append(gid)
            continue

        dst = dst_dir / src.name
        if dst.exists():
            continue

        try:
            if mode == "symlink":
                dst.symlink_to(src)
            elif mode == "hardlink":
                dst.hardlink_to(src)
            elif mode == "copy":
                shutil.copy2(src, dst)
            else:
                raise ValueError(f"Unknown placement mode: {mode}")
        except OSError:
            # On HPC / Windows, symlinks may fail; fall back to copy for link modes
            if mode in ("symlink", "hardlink"):
                shutil.copy2(src, dst)
            else:
                raise

        placed += 1

    return placed, missing


def qc_step(
    *,
    checkm_qa: Path,
    domain_map: Path,
    bacteria_dir: Path,
    archaea_dir: Path,
    out: Path,
    completeness_min: float = 90.0,
    contamination_max: float = 10.0,
    place_mode: str = "symlink",
    force: bool = False,
) -> None:
    out.mkdir(parents=True, exist_ok=True)

    # Basic overwrite policy
    out_all = out / "checkm_clean_all.tsv"
    out_filt = out / "checkm_filtered.tsv"
    if (out_all.exists() or out_filt.exists()) and not force:
        raise FileExistsError(f"QC outputs exist under {out}. Use --force to overwrite.")

    domain_by_id = _read_domain_map(domain_map)

    rows = parse_checkm_qa_table(checkm_qa)

    # Strict join: require every CheckM genome_id to exist in domain_map (fail early on ambiguity)
    merged = []
    missing_in_domain = []
    for r in rows:
        dom = domain_by_id.get(r.genome_id)
        if dom is None:
            missing_in_domain.append(r.genome_id)
            continue
        merged.append((r, dom))

    if missing_in_domain:
        # This is consistent with your "fail early on ambiguity"
        missing_preview = ", ".join(missing_in_domain[:10])
        raise ValueError(
            f"{len(missing_in_domain)} CheckM IDs not present in domain_map.tsv "
            f"(first 10: {missing_preview}). This usually means an ID mismatch between "
            f"prepared genomes and CheckM input."
        )

    # Write all parsed (clean)
    write_tsv(
        out_all,
        header=["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"],
        rows=[
            [r.genome_id, r.marker_lineage, f"{r.completeness:.2f}", f"{r.contamination:.2f}", dom]
            for (r, dom) in merged
        ],
    )

    # Filter
    kept = [
        (r, dom) for (r, dom) in merged
        if r.completeness >= completeness_min and r.contamination <= contamination_max
    ]

    write_tsv(
        out_filt,
        header=["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"],
        rows=[
            [r.genome_id, r.marker_lineage, f"{r.completeness:.2f}", f"{r.contamination:.2f}", dom]
            for (r, dom) in kept
        ],
    )

    # ID lists
    bac_ids = sorted({r.genome_id for (r, dom) in kept if dom.lower().startswith("bact")})
    arc_ids = sorted({r.genome_id for (r, dom) in kept if dom.lower().startswith("arch")})

    (out / "bacteria.txt").write_text("\n".join(bac_ids) + ("\n" if bac_ids else ""), encoding="utf-8")
    (out / "archaea.txt").write_text("\n".join(arc_ids) + ("\n" if arc_ids else ""), encoding="utf-8")

    # Optional bins placement (exact IDs only)
    bins_root = out / "bins"
    placed_bac, miss_bac = place_bins(bac_ids, bacteria_dir, bins_root / "bacteria", mode=place_mode)
    placed_arc, miss_arc = place_bins(arc_ids, archaea_dir, bins_root / "archaea", mode=place_mode)

    # Write a small report/manifest
    report = {
        "inputs": {
            "checkm_qa": str(checkm_qa),
            "domain_map": str(domain_map),
            "bacteria_dir": str(bacteria_dir),
            "archaea_dir": str(archaea_dir),
        },
        "thresholds": {
            "completeness_min": completeness_min,
            "contamination_max": contamination_max,
        },
        "kept": {"bacteria": len(bac_ids), "archaea": len(arc_ids), "total": len(kept)},
        "placed": {"bacteria": placed_bac, "archaea": placed_arc},
        "missing_files": {"bacteria": miss_bac[:20], "archaea": miss_arc[:20]},
        "placement_mode": place_mode,
    }
    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")