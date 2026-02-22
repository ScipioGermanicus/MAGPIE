from __future__ import annotations

import csv
import logging
import os
import shutil
from pathlib import Path
from typing import Literal

from .checkm import checkm_step

_LOGGER = logging.getLogger("magpie.qc")

PlaceMode = Literal["copy", "symlink", "hardlink"]


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _read_domain_map(domain_map_tsv: Path) -> dict[str, str]:
    """
    Expect a TSV with at least: genome_id, domain
    (If your taxonomy uses user_genome instead of genome_id, we accept that.)
    """
    mapping: dict[str, str] = {}
    with open(domain_map_tsv, "r", encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames is None:
            raise ValueError(f"domain_map.tsv has no header: {domain_map_tsv}")

        gid_key = "genome_id" if "genome_id" in r.fieldnames else ("user_genome" if "user_genome" in r.fieldnames else None)
        dom_key = "domain" if "domain" in r.fieldnames else None
        if gid_key is None or dom_key is None:
            raise ValueError(f"domain_map.tsv must contain genome_id (or user_genome) and domain. Found: {r.fieldnames}")

        for row in r:
            gid = (row.get(gid_key) or "").strip()
            dom = (row.get(dom_key) or "").strip()
            if gid and dom:
                mapping[gid] = dom

    if not mapping:
        raise ValueError(f"No entries read from domain_map.tsv: {domain_map_tsv}")
    return mapping


def _read_checkm_min(checkm_min_tsv: Path) -> dict[str, tuple[float, float]]:
    out: dict[str, tuple[float, float]] = {}
    with open(checkm_min_tsv, "r", encoding="utf-8", newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames is None:
            raise ValueError(f"CheckM TSV has no header: {checkm_min_tsv}")
        need = {"genome_id", "checkm_completeness", "checkm_contamination"}
        if not need.issubset(set(r.fieldnames)):
            raise ValueError(f"CheckM TSV must contain {sorted(need)}. Found: {r.fieldnames}")

        for row in r:
            gid = (row.get("genome_id") or "").strip()
            if not gid:
                continue
            try:
                cpl = float((row.get("checkm_completeness") or "").strip())
                cnt = float((row.get("checkm_contamination") or "").strip())
            except ValueError:
                continue
            out[gid] = (cpl, cnt)

    if not out:
        raise ValueError(f"No parsable rows in CheckM TSV: {checkm_min_tsv}")
    return out


def _place_file(src: Path, dst: Path, mode: PlaceMode) -> None:
    _ensure_dir(dst.parent)
    if dst.exists():
        return
    try:
        if mode == "copy":
            shutil.copy2(src, dst)
        elif mode == "symlink":
            os.symlink(src, dst)
        elif mode == "hardlink":
            os.link(src, dst)
        else:
            raise ValueError(f"Unknown place mode: {mode}")
    except OSError as e:
        _LOGGER.warning("Link failed (%s). Falling back to copy: %s -> %s", e, src, dst)
        shutil.copy2(src, dst)


def _index_split_dir(d: Path) -> dict[str, Path]:
    """
    Map genome_id -> file path, assuming filename stem == genome_id.
    Handles dotted IDs correctly (Path(stem) preserves dots).
    """
    idx: dict[str, Path] = {}
    for p in d.iterdir():
        if not p.is_file():
            continue
        name = p.name
        if name.lower().endswith(".gz"):
            name = name[:-3]
        gid = Path(name).stem
        idx[gid] = p
    return idx


def qc_step(
    *,
    tax_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    checkm_results: Path | None,
    completeness_min: float,
    contamination_max: float,
    place_mode: str,
    checkm_bin: str = "checkm",
) -> None:
    tax_dir = tax_dir.resolve()
    out = out.resolve()

    if place_mode not in ("copy", "symlink", "hardlink"):
        raise ValueError("--place-mode must be one of: copy, symlink, hardlink")

    # Your actual taxonomy outputs
    domain_map_fp = tax_dir / "domain_map.tsv"
    split_root = tax_dir / "genomes_to_search_barrnap"
    bac_dir = split_root / "bacteria"
    arc_dir = split_root / "archaea"

    if not domain_map_fp.exists():
        raise FileNotFoundError(f"Expected taxonomy output not found: {domain_map_fp}")
    if not bac_dir.is_dir():
        raise FileNotFoundError(f"Expected taxonomy split dir not found: {bac_dir}")
    if not arc_dir.is_dir():
        raise FileNotFoundError(f"Expected taxonomy split dir not found: {arc_dir}")

    _ensure_dir(out)

    # Run CheckM unless user provided merged results
    if checkm_results is None:
        merged_min = checkm_step(
            bacteria_dir=bac_dir,
            archaea_dir=arc_dir,
            out=out / "checkm",
            cpus=cpus,
            force=force,
            checkm_bin=checkm_bin,
        )
    else:
        merged_min = checkm_results

    domain_map = _read_domain_map(domain_map_fp)
    checkm = _read_checkm_min(merged_min)

    # Strict ID policy: CheckM IDs must exist in domain_map
    missing_in_domain_map = sorted(set(checkm.keys()) - set(domain_map.keys()))
    if missing_in_domain_map:
        raise ValueError(
            "CheckM produced genome IDs that are missing from domain_map.tsv (ID mismatch). "
            f"First few: {missing_in_domain_map[:10]}"
        )

    # Build rows + filter
    rows_all: list[dict[str, str]] = []
    rows_keep: list[dict[str, str]] = []

    for gid, (cpl, cnt) in checkm.items():
        dom = domain_map[gid]
        rec = {
            "genome_id": gid,
            "domain": dom,
            "checkm_completeness": f"{cpl:.3f}",
            "checkm_contamination": f"{cnt:.3f}",
        }
        rows_all.append(rec)
        if cpl >= completeness_min and cnt <= contamination_max:
            rows_keep.append(rec)

    _LOGGER.info(
        "QC filter kept %d/%d genomes (completeness>=%.1f, contamination<=%.1f).",
        len(rows_keep),
        len(rows_all),
        completeness_min,
        contamination_max,
    )

    # Write tables
    out_all = out / "checkm_clean_all.tsv"
    out_keep = out / "checkm_filtered.tsv"
    for fp, rows in ((out_all, rows_all), (out_keep, rows_keep)):
        with open(fp, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(
                fh,
                delimiter="\t",
                fieldnames=["genome_id", "domain", "checkm_completeness", "checkm_contamination"],
            )
            w.writeheader()
            w.writerows(rows)

    # Write ID lists
    bac_ids = sorted({r["genome_id"] for r in rows_keep if r["domain"].lower().startswith("bact")})
    arc_ids = sorted({r["genome_id"] for r in rows_keep if r["domain"].lower().startswith("arch")})
    (out / "bacteria.txt").write_text("\n".join(bac_ids) + ("\n" if bac_ids else ""), encoding="utf-8")
    (out / "archaea.txt").write_text("\n".join(arc_ids) + ("\n" if arc_ids else ""), encoding="utf-8")

    # Place genomes from taxonomy split dirs
    bins_root = out / "bins"
    dst_bac = bins_root / "bacteria"
    dst_arc = bins_root / "archaea"
    _ensure_dir(dst_bac)
    _ensure_dir(dst_arc)

    bac_idx = _index_split_dir(bac_dir)
    arc_idx = _index_split_dir(arc_dir)

    missing_files: list[str] = []

    for gid in bac_ids:
        src = bac_idx.get(gid)
        if src is None:
            missing_files.append(gid)
            continue
        _place_file(src, dst_bac / src.name, mode=place_mode)

    for gid in arc_ids:
        src = arc_idx.get(gid)
        if src is None:
            missing_files.append(gid)
            continue
        _place_file(src, dst_arc / src.name, mode=place_mode)

    if missing_files:
        raise ValueError(
            "Filtered genome IDs could not be resolved to files in taxonomy split dirs. "
            f"First few: {missing_files[:10]}"
        )