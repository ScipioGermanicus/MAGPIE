from __future__ import annotations

import json
from pathlib import Path
from typing import Final

from ..util.shell import run_cmd

# Latest GTDB-Tk versions typically place summaries in out_dir/classify/
SUMMARY_FILES: Final[tuple[str, str]] = (
    "gtdbtk.bac120.summary.tsv",
    "gtdbtk.ar53.summary.tsv",
)


def _find_classify_dir(out_dir: Path) -> Path:
    """
    Return the directory that contains GTDB-Tk summary files.
    Prefer out_dir/classify if present.
    """
    classify = out_dir / "classify"
    if all((classify / fn).exists() for fn in SUMMARY_FILES):
        return classify
    if all((out_dir / fn).exists() for fn in SUMMARY_FILES):
        return out_dir
    # If neither has both, still return classify if it exists (helpful for error msgs)
    return classify if classify.exists() else out_dir


def _assert_classify_dir_has_summaries(classify_dir: Path) -> None:
    missing = [fn for fn in SUMMARY_FILES if not (classify_dir / fn).exists()]
    if missing:
        raise FileNotFoundError(
            "GTDB-Tk summary files not found. Expected the GTDB-Tk 'classify' directory to contain:\n"
            f"  - {SUMMARY_FILES[0]}\n"
            f"  - {SUMMARY_FILES[1]}\n"
            f"Missing: {missing}\n"
            f"Given directory: {classify_dir}"
        )


def gtdbtk_step(
    *,
    mags_dir: Path,
    out: Path,
    gtdb_classify: Path | None,
    force: bool,
    cpus: int = 8,
) -> Path:
    """
    Either:
      - reuse an existing GTDB-Tk classify dir (gtdb_classify), OR
      - run gtdbtk classify_wf into `out`

    Returns:
      Path to the classify directory containing summary files.
    """
    out.mkdir(parents=True, exist_ok=True)
    report_fp = out / "gtdbtk_report.json"
    if report_fp.exists() and not force:
        raise FileExistsError(f"{report_fp} exists. Use --force to overwrite.")

    # Reuse mode: user points directly at the classify folder
    if gtdb_classify is not None:
        classify_dir = gtdb_classify.resolve()
        _assert_classify_dir_has_summaries(classify_dir)
        report = {
            "mode": "reuse",
            "mags_dir": str(mags_dir),
            "out_dir": str(out),
            "classify_dir": str(classify_dir),
            "summaries": list(SUMMARY_FILES),
        }
        report_fp.write_text(json.dumps(report, indent=2), encoding="utf-8")
        return classify_dir

    # Run mode
    mags_dir = mags_dir.resolve()
    out = out.resolve()

    # Basic input presence check
    fastas = [p for p in mags_dir.iterdir() if p.is_file() and p.name.lower().endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"))]
    if not fastas:
        raise ValueError(f"No MAG FASTA files found in: {mags_dir}")

    cmd = [
        "gtdbtk",
        "classify_wf",
        "--genome_dir",
        str(mags_dir),
        "--out_dir",
        str(out),
        "--cpus",
        str(cpus),
    ]

    run_cmd(cmd)

    classify_dir = _find_classify_dir(out)
    _assert_classify_dir_has_summaries(classify_dir)

    report = {
        "mode": "run",
        "mags_dir": str(mags_dir),
        "out_dir": str(out),
        "classify_dir": str(classify_dir),
        "cpus": cpus,
        "cmd": cmd,
        "summaries": list(SUMMARY_FILES),
    }
    report_fp.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return classify_dir