from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass
class ShellError(RuntimeError):
    cmd: list[str]
    returncode: int
    stdout: str
    stderr: str

    def __str__(self) -> str:
        return (
            f"Command failed (exit {self.returncode}): {' '.join(self.cmd)}\n"
            f"--- stdout ---\n{self.stdout}\n"
            f"--- stderr ---\n{self.stderr}\n"
        )


def run_cmd(
    cmd: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> None:
    p = subprocess.run(
        list(cmd),
        cwd=str(cwd) if cwd else None,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if p.returncode != 0:
        raise ShellError(list(cmd), p.returncode, p.stdout, p.stderr)