from __future__ import annotations

from pathlib import Path
import json
import shutil
from collections import Counter

FASTA_SUFFIXES = (".fa", ".fna", ".fasta")

def _list_fastas(mags: Path) -> list[Path]:
    return sorted([p for p in mags.iterdir() if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES])

def _read_rename_map(fp: Path) -> dict[str, str]:
    """
    Read a 2-column TSV with no header: old_id<TAB>new_id
    """
    mapping: dict[str, str] = {}
    for i, line in enumerate(fp.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) != 2:
            raise ValueError(f"Invalid rename-map at line {i}: expected 2 tab-separated columns, got {len(parts)}")
        old, new = parts[0].strip(), parts[1].strip()
        if not old or not new:
            raise ValueError(f"Invalid rename-map at line {i}: empty old/new id")
        mapping[old] = new
    return mapping

def prep_step(*, mags: Path, out: Path, rename_map: Path | None, force: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)

    out_mags = out / "mags"
    out_mags.mkdir(parents=True, exist_ok=True)

    report_fp = out / "prep_report.json"
    mapping_fp = out / "rename_map.tsv"

    if (report_fp.exists() or mapping_fp.exists()) and not force:
        raise FileExistsError(f"{out} already contains prep outputs. Use --force to overwrite.")

    fasta_files = _list_fastas(mags)
    if not fasta_files:
        raise ValueError(f"No FASTA files found in: {mags}")

    # Default genome IDs: input filename stem
    input_ids = [p.stem for p in fasta_files]
    dupes = [k for k, v in Counter(input_ids).items() if v > 1]
    if dupes:
        raise ValueError(f"Duplicate genome IDs (filename stems) found: {sorted(dupes)}")

    # Optional renaming
    user_map: dict[str, str] = _read_rename_map(rename_map) if rename_map else {}
    # Create final mapping (identity unless renamed)
    final_map: dict[str, str] = {}
    for gid in input_ids:
        final_map[gid] = user_map.get(gid, gid)

    # Ensure renamed IDs are unique
    new_ids = list(final_map.values())
    dupes_new = [k for k, v in Counter(new_ids).items() if v > 1]
    if dupes_new:
        raise ValueError(f"Renaming produced duplicate IDs: {sorted(dupes_new)}")

    # Copy files -> out/mags with .fna suffix
    copied: list[dict[str, str]] = []
    for p in fasta_files:
        old_id = p.stem
        new_id = final_map[old_id]
        dest = out_mags / f"{new_id}.fna"
        if dest.exists() and not force:
            raise FileExistsError(f"Destination exists: {dest} (use --force)")
        shutil.copy2(p, dest)
        copied.append({"src": str(p), "dest": str(dest), "old_id": old_id, "new_id": new_id})

    # Write mapping file (always)
    lines = [f"{old}\t{new}" for old, new in sorted(final_map.items())]
    mapping_fp.write_text("\n".join(lines) + "\n")

    report = {
        "mags_dir": str(mags),
        "out_dir": str(out),
        "n_fastas_in": len(fasta_files),
        "n_fastas_out": len(copied),
        "rename_map_provided": bool(rename_map),
        "outputs": {
            "mags_dir": str(out_mags),
            "rename_map_tsv": str(mapping_fp),
            "prep_report_json": str(report_fp),
        },
        "copies": copied,
    }
    report_fp.write_text(json.dumps(report, indent=2))