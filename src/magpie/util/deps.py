# src/magpie/util/deps.py
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence

# Arguments that make a tool print its version (or at least exit quickly).
_VERSION_ARGS = {
    "checkm": ["-h"], "ssu-align": ["-h"], "esl-reformat": ["-h"], "hmmbuild": ["-h"],
    "prodigal": ["-v"], "diamond": ["version"],
}


@dataclass(frozen=True)
class DepCheck:
    exe: str
    purpose: str
    found: bool
    path: Optional[str]
    version: Optional[str]


def _version(exe: str, args: Sequence[str]) -> Optional[str]:
    try:
        p = subprocess.run([exe, *args], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
    except Exception:
        return None
    txt = (p.stdout or "").strip() or (p.stderr or "").strip()
    for line in txt.splitlines():
        if line.strip() and not line.startswith("#"):
            return line.strip()
    return txt.splitlines()[0].strip() if txt else None


def check_tool(exe: str, purpose: str) -> DepCheck:
    path = shutil.which(exe)
    if not path:
        return DepCheck(exe, purpose, False, None, None)
    name = exe.rsplit("/", 1)[-1]
    return DepCheck(exe, purpose, True, path, _version(exe, _VERSION_ARGS.get(name, ["--version"])))


def check_tools(specs: Iterable[tuple[str, str]]) -> List[DepCheck]:
    """specs: (executable, purpose)."""
    return [check_tool(exe, purpose) for exe, purpose in specs]
