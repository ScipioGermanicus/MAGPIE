from __future__ import annotations

import logging
from pathlib import Path

_FORMAT = "%(asctime)s | %(levelname)-7s | %(name)s | %(message)s"


def configure_console(level: int) -> None:
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)
    for h in list(root.handlers):
        if getattr(h, "_magpie_console", False):
            root.removeHandler(h)
    handler = logging.StreamHandler()
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(_FORMAT, datefmt="%H:%M:%S"))
    handler._magpie_console = True  # type: ignore[attr-defined]
    root.addHandler(handler)


def add_file_log(path: Path) -> None:
    """Keep a DEBUG-level log of every run next to its outputs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    for h in root.handlers:
        if isinstance(h, logging.FileHandler) and Path(h.baseFilename) == path.resolve():
            return
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setLevel(logging.DEBUG)
    handler.setFormatter(logging.Formatter(_FORMAT, datefmt="%Y-%m-%d %H:%M:%S"))
    root.addHandler(handler)
