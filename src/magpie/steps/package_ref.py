"""
Assemble the reference files (scripts 19_raxml_info_files.py and 20_reformat.sh; wiki: "Make the
raxml_info files", "Copy and rename the resulting files", "Filter 16S copy files").

  bac_ref/bac_ref.fna         <- 13_hmm_prep/bacteria_raxml-check.raxml.reduced_dna.fna
  bac_ref/bac_ref.hmm         <- 14_hmm_build/bacteria_raxml-check.raxml.reduced_dna.hmm
  bac_ref/bac_ref.tre         <- 11_iqtree/bacteria_16S_iqtree.treefile
  bac_ref/bac_ref.model       <- 12_raxml_evaluate/bacteria_raxml.raxml.bestModel
  bac_ref/bac_ref.raxml_info  <- RAxML 7-style info file made from bacteria_raxml.raxml.log
  bacteria_16S_copies.txt     <- 16S copies of the reference genomes (assembly, 16S_rRNA_Count; max 10)
(and the same for archaea / arc_ref)
"""
from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path
from typing import Dict, List

from ..pipeline import MagpieError
from ..util.io import ensure_dir, read_fasta_ids, write_id_list, write_json, write_tsv

_LOGGER = logging.getLogger("magpie.package_ref")

REF_NAME = {"bacteria": "bac_ref", "archaea": "arc_ref"}


def parse_raxml_ng_log(log_path: Path) -> Dict[str, str]:
    text = log_path.read_text(encoding="utf-8", errors="replace")

    def find(pattern: str, what: str) -> re.Match:
        m = re.search(pattern, text, flags=re.MULTILINE)
        if m is None:
            raise MagpieError(f"Could not find {what} in {log_path}")
        return m

    called = find(r"^(RAxML-NG was called at .+ as follows:)\s*\n\s*\n?(.+)$", "the RAxML-NG command line")
    patterns = re.search(r"Alignment comprises\s+\d+\s+partitions?\s+and\s+(\d+)\s+patterns", text) or \
        find(r"Alignment sites / patterns:\s*\d+\s*/\s*(\d+)", "the number of alignment patterns")
    base = find(r"Base frequencies \(\w+\):\s*([0-9.eE+\- ]+)$", "base frequencies").group(1).split()
    rates = find(r"Substitution rates \(\w+\):\s*([0-9.eE+\- ]+)$", "substitution rates").group(1).split()
    loglik = find(r"Final LogLikelihood:\s+(-?[0-9.eE+\-]+)", "the final log-likelihood").group(1)
    elapsed = re.search(r"Elapsed time:\s+([0-9.eE+\-]+)\s+seconds", text)
    return {"called_line": called.group(1).strip(), "cmd_line": called.group(2).strip(),
            "patterns": patterns.group(1), "base_freqs": " ".join(base), "subs_rates": " ".join(rates),
            "final_ll": loglik, "elapsed": elapsed.group(1) if elapsed else "0.0"}


def write_raxml_info(log_path: Path, out_path: Path) -> None:
    """RAxML 7.x-style info file for SEPP/pplacer, in the layout of script 19 (wiki template)."""
    p = parse_raxml_ng_log(log_path)
    txt = [
        "This is RAxML version 7.7.2 released by Alexandros Stamatakis on July 31 2013.",
        "",
        "This is a RAxML_info file from an --evaluate run, manually reformatted",
        "",
        "Partition: 0",
        f"Alignment Patterns: {p['patterns']}",
        "Name: No Name Provided",
        "DataType: DNA",
        "Substitution Matrix: GTR",
        "",
        p["called_line"],
        "",
        p["cmd_line"],
        "",
        f"Base frequencies: {p['base_freqs']}",
        "",
        f"Inference[0]: Time {p['elapsed']} CAT-based likelihood -0000, best rearrangement setting 5",
        f"alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: {p['subs_rates']}",
        "",
        "",
        "NOT conducting any final model optimizations on all 1 trees under CAT-based",
        "model ....",
        "",
        f"Final GAMMA  likelihood: {p['final_ll']}",
        "",
    ]
    out_path.write_text("\n".join(txt), encoding="utf-8")


def _read_copies(path: Path) -> Dict[str, int]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[1].strip().isdigit():
            out[parts[0].strip()] = int(parts[1])
    return out


def package_ref_step(*, rrna_dir: Path, iqtree_dir: Path, raxml_evaluate_dir: Path, hmm_prep_dir: Path,
                     hmm_build_dir: Path, out: Path, max_16s_copies: int) -> None:
    built: List[str] = []
    report = {}
    for domain, ref in REF_NAME.items():
        src = {
            "fna": hmm_prep_dir / f"{domain}_raxml-check.raxml.reduced_dna.fna",
            "hmm": hmm_build_dir / f"{domain}_raxml-check.raxml.reduced_dna.hmm",
            "tre": iqtree_dir / f"{domain}_16S_iqtree.treefile",
            "model": raxml_evaluate_dir / f"{domain}_raxml.raxml.bestModel",
            "log": raxml_evaluate_dir / f"{domain}_raxml.raxml.log",
        }
        if not all(p.exists() for p in src.values()):
            _LOGGER.warning("%s: no %s is built (inputs missing; see the warnings of earlier steps)", domain, ref)
            report[domain] = {"status": "not_built"}
            continue
        ref_dir = ensure_dir(out / ref)
        for ext in ("fna", "hmm", "tre", "model"):
            shutil.copyfile(src[ext], ref_dir / f"{ref}.{ext}")
        write_raxml_info(src["log"], ref_dir / f"{ref}.raxml_info")

        ids = read_fasta_ids(ref_dir / f"{ref}.fna")
        copies = _read_copies(rrna_dir / f"{domain}_16S_copies.txt")
        missing = [g for g in ids if g not in copies]
        if missing:
            raise MagpieError(f"{domain}: {len(missing)} reference genomes have no 16S copy number (e.g. {missing[:3]})")
        capped = sum(1 for g in ids if copies[g] > max_16s_copies)
        write_tsv(out / f"{domain}_16S_copies.txt", ["assembly", "16S_rRNA_Count"],
                  [[g, min(copies[g], max_16s_copies)] for g in ids])
        built.append(domain)
        _LOGGER.info("%s: %s with %d reference genomes (%d copy numbers capped at %d)", domain, ref, len(ids),
                     capped, max_16s_copies)
        report[domain] = {"status": "built", "n_genomes": len(ids), "capped": capped}
    if not built:
        raise MagpieError("No reference could be built for either domain")
    write_id_list(out / "domains_built.txt", built)
    write_json(out / "report.json", report)
