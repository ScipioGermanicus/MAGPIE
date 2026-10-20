"""
EC trait tables from eggNOG-mapper annotations (script 25_make_ec_table.py).

Same layout as the KO tables; only complete EC numbers (e.g. 1.1.1.1) are kept, as in
PICRUSt2's own ec.txt.gz files.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from .ko_table import trait_table_step

_EC_RE = re.compile(r"\d+\.\d+\.\d+\.\d+")


def normalise_ec(value: str) -> Optional[str]:
    v = value.strip()
    if v.lower().startswith("ec:"):
        v = v[3:].strip()
    return v if _EC_RE.fullmatch(v) else None


def ec_table_step(*, package_ref_dir: Path, eggnog_dir: Path, out: Path, keep_constant: bool) -> None:
    trait_table_step(package_ref_dir=package_ref_dir, eggnog_dir=eggnog_dir, out=out, column="EC",
                     table_name="ec.txt.gz", normalise=normalise_ec, keep_constant=keep_constant)
