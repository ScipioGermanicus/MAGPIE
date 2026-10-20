"""
Reference tree from the 16S alignment with IQ-TREE (script 16_build_trees.sh):

  iqtree -s <domain>_raxml-check.raxml.reduced.phy -m GTR+G -bb 1000 -nt N -seed 12345 -redo
         -pre <domain>_16S_iqtree

(The PICRUSt2 wiki prunes the GTDB genome tree instead; MAGs are not in that tree.)
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..util.io import write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.iqtree")


def iqtree_step(*, raxml_check_dir: Path, out: Path, iqtree_bin: str, threads: int, bootstrap: int, seed: int) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        phy = raxml_check_dir / f"{domain}_raxml-check.raxml.reduced.phy"
        if not phy.exists():
            report[domain] = {"status": "no_input"}
            continue
        prefix = out / f"{domain}_16S_iqtree"
        cmd = [iqtree_bin, "-s", str(phy), "-m", "GTR+G"]
        if bootstrap > 0:
            cmd += ["-bb", str(bootstrap)]
        cmd += ["-nt", str(threads), "-seed", str(seed), "-redo", "-pre", str(prefix)]
        _LOGGER.info("%s: IQ-TREE GTR+G%s", domain, f" with {bootstrap} ultrafast bootstraps" if bootstrap else "")
        run_cmd(cmd, log_path=out / f"{domain}_iqtree.stdout.log")
        report[domain] = {"status": "ok", "treefile": f"{prefix}.treefile"}
    write_json(out / "report.json", report)
