from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class ToolStatus:
    name: str
    found: bool
    path: str | None
    version: str | None
    error: str | None


def _run_version(cmd: list[str]) -> tuple[str | None, str | None]:
    """
    Try to run a version command. Returns (version, error).
    """
    try:
        p = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:  # e.g., permission issues
        return None, f"{type(e).__name__}: {e}"

    out = (p.stdout or "").strip()
    err = (p.stderr or "").strip()

    if p.returncode != 0:
        # Some tools print version to stderr even on non-zero; still capture it.
        msg = err or out or f"non-zero exit code {p.returncode}"
        return None, msg

    # GTDB-Tk often prints version to stdout, but be tolerant.
    return (out or err or None), None


def check_gtdbtk() -> ToolStatus:
    path = shutil.which("gtdbtk")
    if not path:
        return ToolStatus(
            name="gtdbtk",
            found=False,
            path=None,
            version=None,
            error="Not found on PATH.",
        )

    # Best-effort version detection
    # GTDB-Tk generally supports `gtdbtk --version`.
    version, verr = _run_version(["gtdbtk", "--version"])

    # If version is very verbose, keep first line
    if version:
        version = version.splitlines()[0].strip()

    return ToolStatus(
        name="gtdbtk",
        found=True,
        path=path,
        version=version,
        error=verr,
    )