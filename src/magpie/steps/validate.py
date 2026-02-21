from __future__ import annotations
from pathlib import Path
import json
from collections import Counter


FASTA_SUFFIXES = (".fa", ".fna", ".fasta")

def validate_step(*, mags: Path, out: Path, force: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    report_fp = out / "report.json"
    if report_fp.exists() and not force:
        raise FileExistsError(f"{report_fp} exists. Use --force to overwrite.")

    fasta_files = sorted([p for p in mags.iterdir() if p.is_file() and p.suffix.lower() in FASTA_SUFFIXES])
    if not fasta_files:
        raise ValueError(f"No FASTA files found in: {mags}")

    # Minimal checks: unique basenames, non-empty files
    names = [p.stem for p in fasta_files]
    dupes = sorted([n for n, c in Counter(names).items() if c > 1])
    empty = [str(p) for p in fasta_files if p.stat().st_size == 0]

    report = {
        "mags_dir": str(mags),
        "n_fastas": len(fasta_files),
        "duplicates": dupes,
        "empty_files": empty,
    }
    report_fp.write_text(json.dumps(report, indent=2))