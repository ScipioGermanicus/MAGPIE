# src/magpie/steps/prep.py
from __future__ import annotations
from pathlib import Path

def prep_step(*, mags: Path, out: Path, rename_map: Path | None, force: bool) -> None:
    out.mkdir(parents=True, exist_ok=True)
    # TODO: implement prep logic
    (out / "PREP_NOT_IMPLEMENTED.txt").write_text(
        "prep_step is not implemented yet.\n"
        f"mags={mags}\nrename_map={rename_map}\nforce={force}\n"
    )