"""
Split genomes into bacteria and archaea from the GTDB-Tk classification
(script 03_domain_classification.py).

Outputs: domain_map.tsv, bacteria.txt, archaea.txt and genomes_to_search_barrnap/{bacteria,archaea}/
(symlinks to the formatted genomes; use --move-tax-split to copy instead).
"""
from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path

from ..pipeline import MagpieError
from ..util.io import ensure_dir, write_id_list, write_tsv

_LOGGER = logging.getLogger("magpie.taxonomy")


def infer_domain(classification: str) -> str | None:
    if classification.startswith("d__Bacteria"):
        return "bacteria"
    if classification.startswith("d__Archaea"):
        return "archaea"
    return None


def taxonomy_step(*, prep_dir: Path, gtdbtk_dir: Path, out: Path, copy_files: bool = False) -> None:
    ids = {"bacteria": [], "archaea": []}
    rows, unknown = [], []
    with (gtdbtk_dir / "classification.tsv").open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            dom = infer_domain(row["classification"])
            if dom is None:
                unknown.append(row["genome_id"])
                continue
            ids[dom].append(row["genome_id"])
            rows.append([row["genome_id"], dom.capitalize(), row["classification"]])
    if unknown:
        _LOGGER.warning("%d genomes have no GTDB domain and are excluded (e.g. %s)", len(unknown), ", ".join(unknown[:5]))
    if not rows:
        raise MagpieError("No genomes were classified as Bacteria or Archaea")

    write_tsv(out / "domain_map.tsv", ["genome_id", "domain", "classification"], sorted(rows))
    for dom, lst in ids.items():
        lst.sort()
        write_id_list(out / f"{dom}.txt", lst)
        gdir = ensure_dir(out / "genomes_to_search_barrnap" / dom)
        for gid in lst:
            src = (prep_dir / "mags" / f"{gid}_genomic.fna.gz").resolve()
            dst = gdir / src.name
            if copy_files:
                shutil.copyfile(src, dst)
            else:
                dst.symlink_to(src)
    _LOGGER.info("Domains: %d bacteria, %d archaea", len(ids["bacteria"]), len(ids["archaea"]))
