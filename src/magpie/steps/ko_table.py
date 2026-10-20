"""
KO trait tables from eggNOG-mapper annotations (script 24_make_ko_table.py).

Rows are the genomes of bac_ref.fna / arc_ref.fna ("assembly" column), columns are KOs
("ko:K00001"), values are the number of genes annotated with the KO (a gene with several KOs
counts towards each). Gene families with the same count in every genome are removed unless
--keep-constant-traits is given, because PICRUSt2's default max-parsimony HSP method fails on
them; they are listed in <ref>/<table>_constant_removed.tsv.
"""
from __future__ import annotations

import gzip
import logging
import re
from collections import Counter
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from ..pipeline import MagpieError
from ..util.io import ensure_dir, open_text, read_fasta_ids, read_id_list, write_json

_LOGGER = logging.getLogger("magpie.traits")

REF_NAME = {"bacteria": "bac_ref", "archaea": "arc_ref"}
SUFFIX = ".emapper.annotations"


def normalise_ko(value: str) -> Optional[str]:
    v = value.strip()
    if v.startswith("ko:"):
        v = v[3:]
    return f"ko:{v}" if re.fullmatch(r"K\d{5}", v) else None


def count_annotations(path: Path, column: str, normalise: Callable[[str], Optional[str]]) -> Counter:
    counts: Counter = Counter()
    idx = None
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                cols = line[1:].rstrip("\n").split("\t")
                if column in cols:
                    idx = cols.index(column)
                continue
            if idx is None:
                raise MagpieError(f"No header with a '{column}' column before the data in {path}")
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= idx:
                continue
            values = {normalise(v) for v in re.split(r"[,|]", fields[idx]) if v.strip() not in ("", "-")}
            values.discard(None)
            counts.update(values)
    if idx is None:
        raise MagpieError(f"No header with a '{column}' column in {path}")
    return counts


def write_trait_table(path: Path, genomes: List[str], counts: Dict[str, Counter],
                      keep_constant: bool) -> Tuple[int, List[Tuple[str, int]]]:
    families = sorted(set().union(*(counts[g].keys() for g in genomes)))
    kept, constant = [], []
    for fam in families:
        values = {counts[g].get(fam, 0) for g in genomes}
        if len(values) == 1 and not keep_constant:
            constant.append((fam, values.pop()))
        else:
            kept.append(fam)
    with gzip.open(path, "wt", encoding="utf-8") as out:
        out.write("\t".join(["assembly"] + kept) + "\n")
        for g in genomes:
            out.write("\t".join([g] + [str(counts[g].get(f, 0)) for f in kept]) + "\n")
    return len(kept), constant


def trait_table_step(*, package_ref_dir: Path, eggnog_dir: Path, out: Path, column: str, table_name: str,
                     normalise: Callable[[str], Optional[str]], keep_constant: bool) -> None:
    report = {}
    for domain in read_id_list(package_ref_dir / "domains_built.txt"):
        ref = REF_NAME[domain]
        genomes = sorted(read_fasta_ids(package_ref_dir / ref / f"{ref}.fna"))
        counts = {g: count_annotations(eggnog_dir / "annotations" / f"{g}{SUFFIX}", column, normalise) for g in genomes}
        empty = [g for g in genomes if not counts[g]]
        if empty:
            _LOGGER.warning("%s: %d genomes have no %s annotation at all (e.g. %s)", domain, len(empty), column,
                            ", ".join(empty[:3]))
        ref_out = ensure_dir(out / ref)
        n_kept, constant = write_trait_table(ref_out / table_name, genomes, counts, keep_constant)
        if constant:
            with (ref_out / table_name.replace(".txt.gz", "_constant_removed.tsv")).open("w", encoding="utf-8") as fh:
                fh.write("gene_family\tcount_in_every_genome\n")
                fh.writelines(f"{f}\t{v}\n" for f, v in constant)
        if n_kept == 0:
            raise MagpieError(f"{domain}: no {column} gene family varies between the reference genomes; "
                              "PICRUSt2 cannot use an empty table (use more genomes or --keep-constant-traits)")
        _LOGGER.info("%s %s: %d genomes x %d gene families%s", domain, table_name, len(genomes), n_kept,
                     f" ({len(constant)} constant ones removed)" if constant else "")
        report[domain] = {"genomes": len(genomes), "families": n_kept, "constant_removed": len(constant),
                          "genomes_without_annotation": empty}
    write_json(out / "report.json", report)


def ko_table_step(*, package_ref_dir: Path, eggnog_dir: Path, out: Path, keep_constant: bool) -> None:
    trait_table_step(package_ref_dir=package_ref_dir, eggnog_dir=eggnog_dir, out=out, column="KEGG_ko",
                     table_name="ko.txt.gz", normalise=normalise_ko, keep_constant=keep_constant)
