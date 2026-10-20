"""
Build the 16S HMMs (script 18, second part; wiki: "Run the HMMs"):

  hmmbuild --cpu N <domain>_raxml-check.raxml.reduced_dna.hmm <domain>_raxml-check.raxml.reduced_dna.sto
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..util.io import write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.hmm_build")


def hmm_build_step(*, hmm_prep_dir: Path, out: Path, hmmbuild_bin: str, cpus: int) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        sto = hmm_prep_dir / f"{domain}_raxml-check.raxml.reduced_dna.sto"
        if not sto.exists():
            report[domain] = {"status": "no_input"}
            continue
        hmm = out / f"{domain}_raxml-check.raxml.reduced_dna.hmm"
        run_cmd([hmmbuild_bin, "--cpu", str(cpus), str(hmm), str(sto)], log_path=out / f"{domain}_hmmbuild.log")
        _LOGGER.info("%s: HMM written to %s", domain, hmm.name)
        report[domain] = {"status": "ok"}
    write_json(out / "report.json", report)
