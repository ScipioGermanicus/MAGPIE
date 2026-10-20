"""
Optimise the GTR+G model on the reference tree (script 17_run_raxml.sh; wiki: "Run raxml"):

  raxml-ng --evaluate --msa <domain>_raxml-check.raxml.reduced.phy --tree <domain>_16S_iqtree.treefile
           --model GTR+G --prefix <domain>_raxml --threads N
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..util.io import write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.raxml_evaluate")


def raxml_evaluate_step(*, raxml_check_dir: Path, iqtree_dir: Path, out: Path, raxml_ng_bin: str, threads: int) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        phy = raxml_check_dir / f"{domain}_raxml-check.raxml.reduced.phy"
        tree = iqtree_dir / f"{domain}_16S_iqtree.treefile"
        if not (phy.exists() and tree.exists()):
            report[domain] = {"status": "no_input"}
            continue
        prefix = out / f"{domain}_raxml"
        run_cmd([raxml_ng_bin, "--evaluate", "--msa", str(phy), "--tree", str(tree), "--model", "GTR+G",
                 "--prefix", str(prefix), "--threads", str(threads)], log_path=out / f"{domain}_raxml.stdout.log")
        model = Path(f"{prefix}.raxml.bestModel").read_text(encoding="utf-8").strip()
        _LOGGER.info("%s: %s", domain, model)
        report[domain] = {"status": "ok", "model": model}
    write_json(out / "report.json", report)
