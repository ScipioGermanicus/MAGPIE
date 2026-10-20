"""
raxml-ng --check (script 14; wiki: "Run the raxml check"):

  raxml-ng --check --msa <domain>_16S_centroids_ssu_align_best.fna --model GTR+G
           --prefix <domain>_raxml-check --threads N

RAxML-NG writes <prefix>.raxml.reduced.phy when it removes identical sequences or gap-only
columns. As in the wiki and script 14, the archaeal alignment is always converted to PHYLIP as-is
("Convert archaea alignment to phylip so we have the same for both"), replacing any reduced file;
for bacteria the reduced file is used, or the alignment is converted if RAxML-NG wrote none.
Domains with fewer than --min-genomes sequences are skipped (RAxML-NG/IQ-TREE need >= 4).
"""
from __future__ import annotations

import logging
from pathlib import Path

from ..util.io import read_fasta, read_phylip, write_json, write_phylip
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.raxml_check")


def raxml_check_step(*, choose_best_dir: Path, out: Path, raxml_ng_bin: str, threads: int, min_genomes: int) -> None:
    report = {}
    for domain in ("bacteria", "archaea"):
        msa = choose_best_dir / f"{domain}_16S_centroids_ssu_align_best.fna"
        if not msa.exists():
            report[domain] = {"status": "no_input"}
            continue
        records = [(h.split()[0], s) for h, s in read_fasta(msa)]
        if len(records) < min_genomes:
            _LOGGER.warning("%s: only %d genomes with a unique 16S gene; at least %d are needed, so no %s "
                            "reference is built (PICRUSt2's default one can be used instead)",
                            domain, len(records), min_genomes, domain)
            report[domain] = {"status": "too_few_genomes", "n": len(records)}
            continue

        prefix = f"{domain}_raxml-check"
        run_cmd([raxml_ng_bin, "--check", "--msa", str(msa.resolve()), "--model", "GTR+G", "--prefix", prefix,
                 "--threads", str(threads)], cwd=out, log_path=out / f"{prefix}.stdout.log")
        reduced = out / f"{prefix}.raxml.reduced.phy"
        source = "raxml-ng"
        if domain == "archaea" or not reduced.exists():
            write_phylip(records, reduced)
            source = "converted"
        kept = read_phylip(reduced)
        removed = sorted({g for g, _ in records} - {g for g, _ in kept})
        if removed:
            _LOGGER.info("%s: raxml-ng removed %d sequences identical to another after alignment", domain, len(removed))
        if len(kept) < min_genomes:
            _LOGGER.warning("%s: only %d sequences left after raxml-ng --check; no %s reference is built",
                            domain, len(kept), domain)
            reduced.unlink()
            report[domain] = {"status": "too_few_genomes", "n": len(kept)}
            continue
        _LOGGER.info("%s: %d sequences x %d columns (%s)", domain, len(kept), len(kept[0][1]), source)
        report[domain] = {"status": "ok", "phylip": source, "n_in": len(records), "n_out": len(kept),
                          "removed_identical": removed}
    write_json(out / "report.json", report)
