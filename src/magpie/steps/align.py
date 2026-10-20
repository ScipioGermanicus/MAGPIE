"""
Align the 16S centroids with SSU-ALIGN (script 11_align_16S.sh; wiki: "Align sequences"):

  ssu-align --rfonly <domain>_16S_centroids.fasta <domain>_16S_centroids_ssu_align
  esl-reformat -o <domain>_16S_centroids_ssu_align.fna afa \
      <domain>_16S_centroids_ssu_align/<domain>_16S_centroids_ssu_align.<domain>.stk

ssu-align assigns each sequence to its best-matching domain model; as in the wiki, only the
sequences aligned to the domain's own model are kept.
"""
from __future__ import annotations

import logging
import os
import shutil
from pathlib import Path

from ..util.io import has_fasta_records, read_fasta, read_fasta_ids, write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.align")


def ssualign_env(ssu_align_bin: str) -> dict | None:
    """Set SSUALIGNDIR from the conda install (as script 11 does) if it is not set already."""
    if os.environ.get("SSUALIGNDIR"):
        return None
    exe = shutil.which(ssu_align_bin)
    if exe is None:
        return None
    share = Path(exe).resolve().parent.parent / "share"
    for cand in [share / "ssu-align-0.1.1", share / "ssu-align", *sorted(share.glob("ssu-align*"))]:
        if cand.is_dir():
            return {"SSUALIGNDIR": str(cand)}
    return None


def align_step(*, rrna_dir: Path, out: Path, ssu_align_bin: str, esl_reformat_bin: str) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        centroids = rrna_dir / f"{domain}_16S_centroids.fasta"
        if not has_fasta_records(centroids):
            _LOGGER.warning("%s: no 16S centroids to align", domain)
            report[domain] = {"status": "no_input"}
            continue
        name = f"{domain}_16S_centroids_ssu_align"
        # ssu-align creates <name>/ in the working directory and does not accept a path for it.
        run_cmd([ssu_align_bin, "--rfonly", str(centroids.resolve()), name], cwd=out,
                env=ssualign_env(ssu_align_bin), log_path=out / f"{name}.log")
        stk = out / name / f"{name}.{domain}.stk"
        if not stk.exists() or stk.stat().st_size == 0:
            _LOGGER.warning("%s: ssu-align assigned none of the sequences to the %s model", domain, domain)
            report[domain] = {"status": "no_sequences_for_domain_model"}
            continue
        afa = out / f"{name}.fna"
        run_cmd([esl_reformat_bin, "-o", str(afa), "afa", str(stk)])

        recs = read_fasta(afa)
        dropped = sorted(set(read_fasta_ids(centroids)) - {h.split()[0] for h, _ in recs})
        if dropped:
            _LOGGER.warning("%s: %d centroids fit another domain's model better and are dropped (e.g. %s)",
                            domain, len(dropped), ", ".join(dropped[:5]))
        _LOGGER.info("%s: %d sequences aligned, %d columns", domain, len(recs), len(recs[0][1]))
        report[domain] = {"status": "aligned", "n_aligned": len(recs), "alignment_length": len(recs[0][1]),
                          "dropped_other_domain": dropped}
    write_json(out / "report.json", report)
