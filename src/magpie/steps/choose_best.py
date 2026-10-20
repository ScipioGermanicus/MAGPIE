"""
Choose the best genome per 100%-identity 16S cluster (scripts 12/13; logic of the PICRUSt2 wiki,
"Now look at choosing best genome for each cluster").

Cluster membership is read from the vsearch .uc file exactly as in the wiki: rows whose last
field is not '*' pair a member with its centroid. For each cluster:
  * if the completeness values differ: the genome with the highest completeness
    (the first one in cluster order on a tie for the maximum),
  * else if the contamination values differ: the genome with the lowest contamination,
  * else: a random genome (seeded with --seed for reproducibility).
The aligned centroid sequence is renamed to the chosen genome.

Outputs: <domain>_16S_centroids_ssu_align_best.fna, <domain>_16S_clusters_processed.txt,
<domain>_metadata_clusters_ssu_align_centroids.csv
"""
from __future__ import annotations

import csv
import logging
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

from ..util.io import read_fasta, write_fasta, write_json

_LOGGER = logging.getLogger("magpie.choose_best")


def parse_clusters(uc_path: Path, genes_in_alignment: set) -> Dict[str, List[str]]:
    clusters: Dict[str, List[str]] = defaultdict(list)
    with uc_path.open("r", encoding="utf-8") as fh:
        for raw in fh:
            row = raw.rstrip("\n").split("\t")
            if len(row) < 2 or row[-1] == "*":
                continue
            a, b = row[-2], row[-1]
            if a in genes_in_alignment:
                clusters[a].append(b)
            elif b in genes_in_alignment:
                clusters[b].append(a)
            # neither: the centroid fitted another domain's model in ssu-align (dropped, as in the wiki)
    return clusters


def choose_best(members: List[str], quality: Dict[str, Dict[str, float]], rng: random.Random) -> str:
    """The wiki's selection rule."""
    comp = [quality[g]["completeness"] for g in members]
    cont = [quality[g]["contamination"] for g in members]
    if len(set(comp)) == 1:
        if len(set(cont)) == 1:
            return members[rng.randrange(len(members))]
        lowest = min(cont)
        return members[cont.index(lowest)]
    highest = max(comp)
    return members[comp.index(highest)]


def read_quality(checkm_filtered: Path) -> Dict[str, Dict[str, object]]:
    out: Dict[str, Dict[str, object]] = {}
    with checkm_filtered.open("r", encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["genome_id"]] = {**row, "completeness": float(row["checkm_completeness"]),
                                     "contamination": float(row["checkm_contamination"])}
    return out


def choose_best_step(*, qc_dir: Path, rrna_dir: Path, align_dir: Path, out: Path, seed: int) -> None:
    quality = read_quality(qc_dir / "checkm_filtered.tsv")
    report = {}
    for domain in ("bacteria", "archaea"):
        aligned = align_dir / f"{domain}_16S_centroids_ssu_align.fna"
        if not aligned.exists():
            report[domain] = {"status": "no_alignment"}
            continue
        records = [(h.split()[0], s) for h, s in read_fasta(aligned)]
        genes = {gid for gid, _ in records}
        clusters = parse_clusters(rrna_dir / f"{domain}_16S_clusters.uc", genes)

        best: Dict[str, str] = {}
        for centroid, others in clusters.items():
            members = [centroid] + others
            best[centroid] = choose_best(members, quality, random.Random(f"{seed}:{domain}:{centroid}"))

        renamed = [(best.get(gid, gid), seq) for gid, seq in records]
        write_fasta(renamed, out / f"{domain}_16S_centroids_ssu_align_best.fna", width=60)
        with (out / f"{domain}_16S_clusters_processed.txt").open("w", encoding="utf-8") as fh:
            fh.write("Centroid\tBest\tAll genomes\n")
            for centroid, others in clusters.items():
                fh.write(f"{centroid}\t{best[centroid]}\t{','.join([centroid] + others)}\n")

        ids = [gid for gid, _ in renamed]
        cols = ["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"]
        with (out / f"{domain}_metadata_clusters_ssu_align_centroids.csv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(cols)
            for gid in ids:
                w.writerow([quality[gid].get(c, "") for c in cols])

        n_changed = sum(1 for c, b in best.items() if c != b)
        _LOGGER.info("%s: %d genomes in the alignment, %d multi-genome clusters, %d centroids replaced by a better genome",
                     domain, len(ids), len(clusters), n_changed)
        report[domain] = {"n_genomes": len(ids), "n_multi_genome_clusters": len(clusters), "n_replaced": n_changed}
    write_json(out / "report.json", report)
