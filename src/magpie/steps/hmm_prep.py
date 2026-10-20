"""
Reformat the reference alignment for HMM building (script 18, first part; wiki: "Resave the phylip
files and reformat for the hmm files"):

  <domain>_raxml-check.raxml.reduced.phy -> <domain>_raxml-check.raxml.reduced.fna
  esl-reformat -d -o <...>.reduced_dna.fna afa <...>.reduced.fna
  esl-reformat -o <...>.reduced_dna.sto stockholm <...>.reduced_dna.fna
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..util.io import read_phylip, write_fasta, write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.hmm_prep")


def hmm_prep_step(*, raxml_check_dir: Path, out: Path, esl_reformat_bin: str) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        phy = raxml_check_dir / f"{domain}_raxml-check.raxml.reduced.phy"
        if not phy.exists():
            report[domain] = {"status": "no_input"}
            continue
        base = f"{domain}_raxml-check.raxml.reduced"
        records = read_phylip(phy)
        write_fasta(records, out / f"{base}.fna", width=60)  # as Biopython's AlignIO in script 18
        run_cmd([esl_reformat_bin, "-d", "-o", str(out / f"{base}_dna.fna"), "afa", str(out / f"{base}.fna")])
        run_cmd([esl_reformat_bin, "-o", str(out / f"{base}_dna.sto"), "stockholm", str(out / f"{base}_dna.fna")])
        _LOGGER.info("%s: %d sequences converted to DNA FASTA and Stockholm", domain, len(records))
        report[domain] = {"status": "ok", "n_taxa": len(records)}
    write_json(out / "report.json", report)
