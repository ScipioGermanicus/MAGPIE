"""
16S copies, one 16S per genome, and clustering of identical 16S genes
(scripts 07-10; wiki: "Count copies per genome", "Now cluster the ones with multiple 16S copies",
"Make files with single 16S for each genome", "Cluster these single 16S genes for all genomes").

  1. Count records with '16S' in the barrnap --outseq header per genome -> <domain>_16S_copies.txt
     (all genomes, including 0 copies); write <domain>_16S_single/ and <domain>_16S_multiple/.
  2. vsearch --cluster_fast <multiple>/<g>.fna --id 0.9 --centroids <domain>_16S_clustered/<g>.fna
  3. One sequence per genome: a random centroid (wiki; seeded with --seed) or the longest
     (--multi-copy-pick longest, as in the earlier MAGPIE scripts) -> <domain>_16S_genes.fasta
  4. vsearch --cluster_fast <domain>_16S_genes.fasta --id 1 --centroids <domain>_16S_centroids.fasta
             --uc <domain>_16S_clusters.uc
"""
from __future__ import annotations

import logging
import random
import statistics
from pathlib import Path
from typing import Dict, List, Tuple

from ..util.io import ensure_dir, list_fasta_files, read_fasta, strip_fasta_suffix, write_fasta, write_json, write_tsv
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.rrna")


def _vsearch(vsearch_bin: str, fasta: Path, centroids: Path, cluster_id: float, threads: int, uc: Path | None = None) -> None:
    cmd = [vsearch_bin, "--cluster_fast", str(fasta), "--id", format(cluster_id, "g"),
           "--centroids", str(centroids), "--threads", str(threads), "--quiet"]
    if uc is not None:
        cmd += ["--uc", str(uc)]
    run_cmd(cmd)


def rrna_step(*, barrnap_dir: Path, out: Path, vsearch_bin: str, cpus: int, multi_cluster_id: float,
              final_cluster_id: float, multi_copy_pick: str, seed: int) -> None:
    if multi_copy_pick not in ("random", "longest"):
        raise ValueError("multi_copy_pick must be 'random' or 'longest'")
    report: Dict[str, dict] = {}

    for domain in ("bacteria", "archaea"):
        single = ensure_dir(out / f"{domain}_16S_single")
        multiple = ensure_dir(out / f"{domain}_16S_multiple")
        clustered = ensure_dir(out / f"{domain}_16S_clustered")

        # 1. copies per genome
        copies: List[Tuple[str, int]] = []
        for fp in list_fasta_files(barrnap_dir / domain):
            gid = strip_fasta_suffix(fp.name)
            recs = [(h, s) for h, s in read_fasta(fp) if "16S" in h.split()[0]]
            copies.append((gid, len(recs)))
            if len(recs) == 1:
                write_fasta(recs, single / f"{gid}.fna")
            elif len(recs) > 1:
                write_fasta(recs, multiple / f"{gid}.fna")
        with (out / f"{domain}_16S_copies.txt").open("w", encoding="utf-8") as fh:
            for gid, n in copies:
                fh.write(f"{gid}\t{n}\n")
        pos = [n for _, n in copies if n > 0]
        if pos:
            _LOGGER.info("%s: %d of %d genomes have a 16S gene (mean %.2f, median %g, max %d copies)", domain,
                         len(pos), len(copies), statistics.mean(pos), statistics.median(pos), max(pos))
        else:
            _LOGGER.warning("%s: no 16S genes found", domain)

        # 2. cluster multi-copy genomes at 90%
        for fp in list_fasta_files(multiple):
            _vsearch(vsearch_bin, fp, clustered / fp.name, multi_cluster_id, cpus)

        # 3. one 16S per genome
        genes: List[Tuple[str, str]] = []
        sources = []
        for fp in list_fasta_files(clustered):
            gid = strip_fasta_suffix(fp.name)
            recs = read_fasta(fp)
            if len(recs) > 1 and multi_copy_pick == "random":
                seq = recs[random.Random(f"{seed}:{gid}").randrange(len(recs))][1]
            else:
                seq = max(recs, key=lambda r: len(r[1]))[1]
            genes.append((gid, seq))
            sources.append((gid, "clustered", len(recs)))
        for fp in list_fasta_files(single):
            gid = strip_fasta_suffix(fp.name)
            genes.append((gid, read_fasta(fp)[0][1]))
            sources.append((gid, "single", 1))
        genes.sort(key=lambda r: r[0])
        sources.sort(key=lambda r: r[0])
        genes_fa = out / f"{domain}_16S_genes.fasta"
        write_fasta(genes, genes_fa, width=60)
        write_tsv(out / f"{domain}_16S_genes.map.tsv", ["genome_id", "source", "n_centroids"], sources)

        # 4. cluster identical 16S genes across genomes
        centroids = out / f"{domain}_16S_centroids.fasta"
        uc = out / f"{domain}_16S_clusters.uc"
        if genes:
            _vsearch(vsearch_bin, genes_fa, centroids, final_cluster_id, cpus, uc=uc)
            n_cent = len(read_fasta(centroids))
        else:
            centroids.write_text("", encoding="utf-8")
            uc.write_text("", encoding="utf-8")
            n_cent = 0
        lengths = [len(s) for _, s in genes]
        report[domain] = {"genomes": len(copies), "with_16S": len(pos), "multi_copy": sum(1 for n in pos if n > 1),
                          "divergent_copies": sum(1 for s in sources if s[2] > 1), "centroids": n_cent,
                          "length_min": min(lengths, default=0), "length_max": max(lengths, default=0)}
        _LOGGER.info("%s: %d genomes with 16S -> %d unique 16S clusters", domain, len(genes), n_cent)

    write_json(out / "report.json", {"multi_copy_pick": multi_copy_pick, "seed": seed, "domains": report})
