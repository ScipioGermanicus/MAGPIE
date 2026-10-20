from __future__ import annotations

import logging
import os
import shlex
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

_LOGGER = logging.getLogger("magpie.shell")

# Every external command is appended here (set by the pipeline), so a run can be audited
# and reproduced by hand.
_COMMAND_LOG: Path | None = None
_LOCK = threading.Lock()


def set_command_log(path: Path | None) -> None:
    global _COMMAND_LOG
    _COMMAND_LOG = path
    if path is not None:
        path.parent.mkdir(parents=True, exist_ok=True)


@dataclass
class ShellError(RuntimeError):
    cmd: list[str]
    returncode: int
    stdout: str
    stderr: str

    def __str__(self) -> str:
        def tail(s: str, n: int = 30) -> str:
            return "\n".join(s.strip().splitlines()[-n:])

        return (
            f"Command failed (exit {self.returncode}): {shlex.join(self.cmd)}\n"
            f"--- stdout (tail) ---\n{tail(self.stdout)}\n"
            f"--- stderr (tail) ---\n{tail(self.stderr)}\n"
        )


def run_cmd(
    cmd: Sequence[str],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    stdout_path: Path | None = None,
    log_path: Path | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run an external command (no shell) and raise ShellError if it fails.

    stdout_path: write the tool's stdout to this file.
    log_path:    append the command and its output to this file.
    env:         extra environment variables merged over os.environ.
    """
    cmd = [str(c) for c in cmd]
    line = shlex.join(cmd) + (f" > {shlex.quote(str(stdout_path))}" if stdout_path else "")
    if cwd is not None:
        line = f"(cd {shlex.quote(str(cwd))} && {line})"
    _LOGGER.debug("Running: %s", line)
    if _COMMAND_LOG is not None:
        with _LOCK, _COMMAND_LOG.open("a", encoding="utf-8") as fh:
            fh.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n{line}\n")

    full_env = None
    if env:
        full_env = dict(os.environ)
        full_env.update(env)

    if stdout_path is not None:
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("w", encoding="utf-8") as fout:
            p = subprocess.run(cmd, cwd=cwd, env=full_env, stdout=fout, stderr=subprocess.PIPE, text=True)
        out = ""
    else:
        p = subprocess.run(cmd, cwd=cwd, env=full_env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out = p.stdout or ""

    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"$ {shlex.join(cmd)}\n{out}{p.stderr or ''}\n")

    if check and p.returncode != 0:
        raise ShellError(cmd, p.returncode, out, p.stderr or "")
    return p


def run(cmd: Sequence[str], *, cwd: Path | None = None, env: Mapping[str, str] | None = None) -> None:
    """Compatibility wrapper used by MAGPIE steps."""
    run_cmd(cmd, cwd=cwd, env=env)
