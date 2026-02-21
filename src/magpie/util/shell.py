from __future__ import annotations
import subprocess
from pathlib import Path
from typing import Sequence, Mapping, Optional

class ShellError(RuntimeError):
    pass

def run(cmd: Sequence[str], *, cwd: Path | None = None, env: Mapping[str, str] | None = None) -> None:
    p = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        env=dict(env) if env else None,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise ShellError(
            "Command failed:\n"
            f"  cmd: {' '.join(cmd)}\n"
            f"  exit: {p.returncode}\n"
            f"  stdout:\n{p.stdout}\n"
            f"  stderr:\n{p.stderr}\n"
        )