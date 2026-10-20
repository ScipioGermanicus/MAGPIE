"""
Format the input MAGs (script 01_formating_genomes.sh).

Every genome is written to <out>/mags/<genome_id>_genomic.fna.gz and recorded in id_map.tsv
(original_filename, original_id, new_id). Genome IDs are the file names without FASTA extensions,
or MAG0001, MAG0002, ... with --sequential-ids (files sorted by name, as in script 01), or taken
from --rename-map. Characters outside [A-Za-z0-9._-] are replaced by '_', because IDs become
FASTA, PHYLIP and Newick labels.
"""
from __future__ import annotations

import csv
import json
import logging
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

from ..pipeline import MagpieError
from ..util.io import copy_or_decompress, ensure_dir, iter_fasta, list_fasta_files, open_text, sanitise_id, strip_fasta_suffix

_LOGGER = logging.getLogger("magpie.prep")


@dataclass(frozen=True)
class GenomeInput:
    path: Path
    original_id: str
    genome_id: str


def _read_rename_map(path: Path) -> Dict[str, str]:
    """Two columns: original file name or ID, new ID. A header line is optional."""
    mapping: Dict[str, str] = {}
    with open_text(path) as fh:
        for i, row in enumerate(csv.reader(fh, delimiter="\t")):
            if not row or not row[0].strip() or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise MagpieError(f"--rename-map line {i + 1} needs two tab-separated columns: {path}")
            old, new = row[0].strip(), row[1].strip()
            if i == 0 and old.lower() in {"original_filename", "original", "original_id", "old_id", "raw_id"}:
                continue
            mapping[old] = new
    return mapping


def plan_genomes(mags: Path, rename_map: Path | None, sequential_ids: bool, out_prefix: str, pad: int) -> List[GenomeInput]:
    """Decide the genome ID of every input file (shared by validate and prep)."""
    files = list_fasta_files(mags)
    if not files:
        raise MagpieError(f"No MAG FASTA files (.fa/.fna/.fasta, optionally .gz) found in {mags}")
    if rename_map is not None and sequential_ids:
        raise MagpieError("--rename-map cannot be combined with --sequential-ids")

    rename = _read_rename_map(rename_map) if rename_map is not None else {}
    plan = []
    for i, path in enumerate(files, start=1):
        original_id = strip_fasta_suffix(path.name)
        if sequential_ids:
            new_id = f"{out_prefix}{i:0{pad}d}"
        else:
            new_id = rename.get(path.name, rename.get(original_id, original_id))
        plan.append(GenomeInput(path, original_id, sanitise_id(new_id)))

    dupes = sorted(g for g, c in Counter(p.genome_id for p in plan).items() if c > 1)
    if dupes:
        raise MagpieError(f"{len(dupes)} genome IDs are not unique (e.g. {', '.join(dupes[:5])}); "
                          "rename the files or use --rename-map / --sequential-ids")
    if rename:
        unmapped = [p.path.name for p in plan if p.path.name not in rename and p.original_id not in rename]
        if unmapped:
            _LOGGER.warning("%d input files are not in --rename-map and keep their own name (e.g. %s)",
                            len(unmapped), ", ".join(unmapped[:3]))
    return plan


def prep_step(*, mags: Path, out: Path, rename_map: Path | None, sequential_ids: bool, out_prefix: str, pad: int) -> None:
    mags_out = ensure_dir(out / "mags")
    plan = plan_genomes(mags, rename_map, sequential_ids, out_prefix, pad)

    rows = []
    for g in plan:
        dst = mags_out / f"{g.genome_id}_genomic.fna.gz"
        copy_or_decompress(g.path, dst)
        contigs, total = set(), 0
        for header, seq in iter_fasta(dst):
            cid = header.split()[0]
            if cid in contigs:
                raise MagpieError(f"Duplicate contig name '{cid}' in {g.path}")
            contigs.add(cid)
            total += len(seq)
        if total == 0:
            raise MagpieError(f"Genome file contains no sequence: {g.path}")
        rows.append((g.path.name, g.original_id, g.genome_id, len(contigs), total))

    with (out / "id_map.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["original_filename", "original_id", "new_id", "n_contigs", "length"])
        w.writerows(rows)

    sanitised = [r for r in rows if not sequential_ids and rename_map is None and r[1] != r[2]]
    if sanitised:
        _LOGGER.warning("%d genome IDs contained characters outside [A-Za-z0-9._-] and were changed "
                        "(e.g. %s -> %s); see id_map.tsv", len(sanitised), sanitised[0][1], sanitised[0][2])
    (out / "prep_report.json").write_text(json.dumps({"mags_dir": str(mags), "n_genomes": len(rows)}, indent=2),
                                          encoding="utf-8")
    _LOGGER.info("Formatted %d genomes into %s", len(rows), mags_out)
