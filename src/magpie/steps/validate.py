# src/magpie/steps/validate.py
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

from ..util.deps import check_gtdbtk

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def validate_step(*, mags: Path, out: Path, force: bool, require_gtdbtk: bool = False) -> None:
    out.mkdir(parents=True, exist_ok=True)
    report_fp = out / "report.json"
    if report_fp.exists() and not force:
        raise FileExistsError(f"{report_fp} exists. Use --force to overwrite.")

    # Accept plain + gz FASTA
    def is_fasta(p: Path) -> bool:
        name = p.name.lower()
        return any(name.endswith(s) for s in FASTA_SUFFIXES)

    fasta_files = sorted([p for p in mags.iterdir() if p.is_file() and is_fasta(p)])
    if not fasta_files:
        raise ValueError(f"No FASTA files found in: {mags}")

    # Genome IDs are derived from basename (same logic as earlier; quick/cheap check)
    # Note: For *.fa.gz, .stem removes only .gz; we normalise by stripping second suffix too.
    def genome_id(p: Path) -> str:
        name = p.name
        lower = name.lower()
        if lower.endswith(".gz"):
            base = name[:-3]  # drop ".gz"
            for suff in (".fa", ".fna", ".fasta"):
                if base.lower().endswith(suff):
                    return base[: -len(suff)]
            return Path(base).stem
        return p.stem

    names = [genome_id(p) for p in fasta_files]
    dupes = sorted([n for n, c in Counter(names).items() if c > 1])
    empty = [str(p) for p in fasta_files if p.stat().st_size == 0]
    suffix_counts = dict(Counter([(".gz" if p.name.lower().endswith(".gz") else p.suffix.lower()) for p in fasta_files]))

    # Dependency checks (GTDB-Tk only for now)
    gtdb = check_gtdbtk()
    dependencies = {
        "gtdbtk": {
            "found": gtdb.found,
            "path": gtdb.path,
            "version": gtdb.version,
            "error": gtdb.error,
        }
    }

    if require_gtdbtk and not gtdb.found:
        raise RuntimeError(
            "GTDB-Tk (gtdbtk) not found on PATH. "
            "Install it (e.g., conda/mamba) and re-run. "
            "Example: `conda install -c conda-forge -c bioconda gtdbtk`."
        )

    report = {
        "mags_dir": str(mags),
        "n_fastas": len(fasta_files),
        "fastas": [str(p) for p in fasta_files],
        "suffix_counts": suffix_counts,
        "duplicates": dupes,
        "empty_files": empty,
        "dependencies": dependencies,
        "require_gtdbtk": require_gtdbtk,
    }
    report_fp.write_text(json.dumps(report, indent=2), encoding="utf-8")