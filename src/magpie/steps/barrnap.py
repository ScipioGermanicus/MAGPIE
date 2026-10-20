"""
Predict rRNA genes with barrnap (script 06; command as in the PICRUSt2 wiki, "Run barrnap"):

  barrnap -q -k bac|arc <genome>.fna --outseq <domain>/<genome>.fna --reject 0.8 --threads N

Outputs per domain: <genome>.fna (all rRNA hits, from --outseq) and <genome>.gff.
barrnap 0.9 reads uncompressed FASTA only, so each genome is decompressed to a temporary file.
"""
from __future__ import annotations

import gzip
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import List

from ..util.io import ensure_dir, write_json
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.barrnap")

KINGDOM = {"bacteria": "bac", "archaea": "arc"}


def _barrnap_one(gid: str, genome_gz: Path, out_dir: Path, kingdom: str, barrnap_bin: str, reject: float,
                 threads: int) -> None:
    tmp = out_dir / f"{gid}.tmp.fna"
    with gzip.open(genome_gz, "rb") as fin, tmp.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
    try:
        run_cmd([barrnap_bin, "-q", "-k", kingdom, str(tmp), "--outseq", str(out_dir / f"{gid}.fna"),
                 "--reject", str(reject), "--threads", str(threads)], stdout_path=out_dir / f"{gid}.gff")
    finally:
        tmp.unlink(missing_ok=True)
        tmp.with_name(tmp.name + ".fai").unlink(missing_ok=True)
    if not (out_dir / f"{gid}.fna").exists():  # barrnap writes nothing when there are no hits
        (out_dir / f"{gid}.fna").write_text("", encoding="utf-8")


def barrnap_step(*, qc_dir: Path, out: Path, barrnap_bin: str, reject: float, cpus: int) -> None:
    report = {}
    for domain, kingdom in KINGDOM.items():
        ids: List[str] = [x for x in (qc_dir / f"{domain}.txt").read_text().split() if x]
        dom_out = ensure_dir(out / domain)
        if not ids:
            report[domain] = 0
            continue
        jobs = max(1, min(cpus, len(ids)))
        threads = max(1, cpus // jobs)
        _LOGGER.info("barrnap on %d %s genomes (%d parallel jobs)", len(ids), domain, jobs)
        bins = qc_dir / "bins" / domain
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(lambda g: _barrnap_one(g, bins / f"{g}_genomic.fna.gz", dom_out, kingdom, barrnap_bin,
                                                 reject, threads), ids))
        report[domain] = len(ids)
    write_json(out / "report.json", {"genomes": report, "reject": reject})
