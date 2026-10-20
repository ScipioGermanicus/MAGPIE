# src/magpie/steps/validate.py
"""Check the inputs and that every external tool needed by the selected steps is available."""
from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import List

from ..pipeline import MagpieError
from ..util.deps import check_tools
from ..util.io import open_text
from .prep import plan_genomes

_LOGGER = logging.getLogger("magpie.validate")


def required_tools(cfg, steps: List[str]) -> List[tuple[str, str]]:
    need = []
    if "gtdbtk" in steps and cfg.gtdb_classify is None:
        need.append((cfg.gtdbtk_bin, "GTDB-Tk classification (or give --gtdb-classify)"))
    if "qc" in steps and not (cfg.checkm_qa or cfg.checkm_qa_bacteria or cfg.checkm_qa_archaea):
        need.append((cfg.checkm_bin, "CheckM quality estimates (or give --checkm-qa)"))
    if "barrnap" in steps:
        need.append((cfg.barrnap_bin, "16S rRNA gene prediction"))
    if "rrna" in steps:
        need.append((cfg.vsearch_bin, "16S clustering"))
    if "align" in steps:
        need += [(cfg.ssu_align_bin, "16S alignment"), (cfg.esl_reformat_bin, "alignment conversion")]
    if "raxml_check" in steps or "raxml_evaluate" in steps:
        need.append((cfg.raxml_ng_bin, "alignment check / model evaluation"))
    if "iqtree" in steps:
        need.append((cfg.iqtree_bin, "16S reference tree"))
    if "hmm_prep" in steps:
        need.append((cfg.esl_reformat_bin, "alignment conversion"))
    if "hmm_build" in steps:
        need.append((cfg.hmmbuild_bin, "16S HMM"))
    if "eggnog" in steps and cfg.runs_emapper:
        need.append((cfg.emapper_bin, "functional annotation"))
    seen, out = set(), []
    for exe, purpose in need:
        if exe not in seen:
            seen.add(exe)
            out.append((exe, purpose))
    return out


def _barrnap_problem(exe: str, version: str | None) -> str | None:
    m = re.search(r"barrnap\s+(\d+)\.", version or "")
    if m and int(m.group(1)) >= 1:
        return (f"{version} is not supported: barrnap 1.x has no --reject option, which the PICRUSt2 "
                f"workflow uses (--reject 0.8). Install barrnap 0.9 (conda install -c bioconda barrnap=0.9).")
    return None


def validate_step(*, cfg, out: Path, steps: List[str]) -> None:
    problems: List[str] = []
    genomes = []
    if "prep" in steps:
        if cfg.mags is None:
            raise MagpieError("--mags is required to run the prep step")
        plan = plan_genomes(cfg.mags, cfg.rename_map, cfg.sequential_ids, cfg.out_prefix, cfg.pad)
        for g in plan:
            with open_text(g.path) as fh:
                first = next((ln for ln in fh if ln.strip()), "")
            if not first.startswith(">"):
                problems.append(f"{g.path} is empty or not FASTA")
        genomes = [{"file": str(g.path), "genome_id": g.genome_id} for g in plan]

    for opt, val in (("--gtdb-classify", cfg.gtdb_classify), ("--checkm-qa", cfg.checkm_qa),
                     ("--checkm-qa-bacteria", cfg.checkm_qa_bacteria), ("--checkm-qa-archaea", cfg.checkm_qa_archaea),
                     ("--eggnog-existing-dir", cfg.eggnog_existing_dir), ("--eggnog-data-dir", cfg.eggnog_data_dir)):
        if val is not None and not Path(val).exists():
            problems.append(f"{opt} does not exist: {val}")
    if "gtdbtk" in steps and cfg.gtdb_classify is None and not os.environ.get("GTDBTK_DATA_PATH"):
        problems.append("GTDB-Tk would be run but GTDBTK_DATA_PATH is not set (or give --gtdb-classify)")
    if "eggnog" in steps and cfg.eggnog_existing_dir is None and cfg.eggnog_data_dir is None:
        problems.append("Annotation needs --eggnog-existing-dir (existing annotations) and/or --eggnog-data-dir "
                        "(to run eggNOG-mapper); or stop before annotation with --to package_ref")

    checks = check_tools(required_tools(cfg, steps))
    for c in checks:
        if not c.found:
            problems.append(f"{c.exe} not found on PATH (needed for {c.purpose})")
        elif c.exe.endswith("barrnap"):
            msg = _barrnap_problem(c.exe, c.version)
            if msg:
                problems.append(msg)

    (out / "report.json").write_text(json.dumps({
        "steps": steps, "genomes": genomes,
        "tools": [c.__dict__ for c in checks], "problems": problems}, indent=2), encoding="utf-8")
    if problems:
        raise MagpieError("Validation failed:\n" + "\n".join(f"  * {p}" for p in problems))
    _LOGGER.info("Inputs and %d tools OK%s", len(checks), f" ({len(genomes)} genomes)" if genomes else "")
