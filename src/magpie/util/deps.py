# src/magpie/util/deps.py
from __future__ import annotations

from dataclasses import dataclass
import shutil
import subprocess
from typing import Optional


@dataclass(frozen=True)
class DepCheck:
    found: bool
    path: Optional[str]
    version: Optional[str]
    error: Optional[str]


def _run_version(cmd: list[str]) -> tuple[Optional[str], Optional[str]]:
    try:
        p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        return None, str(e)
    out = (p.stdout or "").strip()
    err = (p.stderr or "").strip()
    # Some tools write version to stderr; prefer stdout but fall back.
    txt = out if out else err
    if p.returncode != 0 and not txt:
        return None, f"non-zero exit ({p.returncode})"
    return txt.splitlines()[0] if txt else None, None


def check_gtdbtk() -> DepCheck:
    exe = "gtdbtk"
    path = shutil.which(exe)
    if not path:
        return DepCheck(found=False, path=None, version=None, error="not found on PATH")
    version, err = _run_version([exe, "--version"])
    return DepCheck(found=True, path=path, version=version, error=err)


def check_checkm() -> DepCheck:
    exe = "checkm"
    path = shutil.which(exe)
    if not path:
        return DepCheck(found=False, path=None, version=None, error="not found on PATH")
    # CheckM supports '--version' in many installs; if it fails, we still mark found=True.
    version, err = _run_version([exe, "--version"])
    return DepCheck(found=True, path=path, version=version, error=err)