"""
eggNOG-mapper annotation of the reference genomes (scripts 22/23; command as in the PICRUSt2 wiki,
"Annotate genomes"):

  emapper.py -m diamond --itype genome --genepred prodigal --data_dir <eggnog_data>
             -i <genome>.fna -o <genome> --output_dir <out>/emapper --cpu N --override

Existing annotations are reused with --eggnog-existing-dir: files named <name>.emapper.annotations
(optionally .gz) where <name> is the MAGPIE genome ID, the original file name or its stem (e.g.
BC13.bin.1.143.emapper.annotations or MAG1988.emapper.annotations). eggNOG-mapper is run only
for reference genomes without an existing annotation, and only if --eggnog-data-dir is given.

Only genomes in bac_ref.fna / arc_ref.fna are annotated. The step keeps its directory between
runs: genomes whose annotation finished earlier are not annotated again.

Output: annotations/<genome_id>.emapper.annotations for every reference genome.
"""
from __future__ import annotations

import gzip
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Optional

from ..pipeline import MagpieError
from ..util.ids import IdLookup
from ..util.io import ensure_dir, read_fasta_ids, read_id_list, write_json, write_tsv
from ..util.shell import ShellError, run_cmd

_LOGGER = logging.getLogger("magpie.eggnog")

SUFFIX = ".emapper.annotations"
REF_NAME = {"bacteria": "bac_ref", "archaea": "arc_ref"}


def find_existing_annotations(ann_dir: Path, lookup: IdLookup) -> Dict[str, Path]:
    """MAGPIE genome ID -> existing annotation file in ann_dir."""
    found: Dict[str, Path] = {}
    for p in sorted(ann_dir.iterdir()):
        name = p.name
        if not p.is_file() or not (name.endswith(SUFFIX) or name.endswith(SUFFIX + ".gz")):
            continue
        gid = lookup.resolve(name[: name.index(SUFFIX)])
        if gid is None:
            continue
        if gid in found:
            raise MagpieError(f"Two eggNOG annotation files for genome {gid} in {ann_dir}: {found[gid].name}, {name}")
        found[gid] = p
    return found


def annotation_complete(path: Path, require_footer: bool) -> bool:
    """A file written by eggNOG-mapper has a '#query' header and ends with '## Total time'."""
    if not path.exists() or path.stat().st_size == 0:
        return False
    opener = gzip.open if path.name.endswith(".gz") else open
    header = footer = False
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#query"):
                header = True
            elif line.startswith("## Total time"):
                footer = True
    return header and (footer or not require_footer)


def _emapper_one(gid: str, genome_gz: Path, emapper_dir: Path, emapper_bin: str, data_dir: Path, threads: int) -> None:
    tmp = ensure_dir(emapper_dir / "tmp" / gid)
    fna = tmp / f"{gid}.fna"
    with gzip.open(genome_gz, "rb") as fin, fna.open("wb") as fout:
        shutil.copyfileobj(fin, fout)
    try:
        run_cmd([emapper_bin, "-m", "diamond", "--itype", "genome", "--genepred", "prodigal",
                 "--data_dir", str(data_dir), "-i", str(fna), "-o", gid, "--output_dir", str(emapper_dir),
                 "--temp_dir", str(tmp), "--cpu", str(threads), "--override"],
                log_path=emapper_dir / "logs" / f"{gid}.log")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if not annotation_complete(emapper_dir / f"{gid}{SUFFIX}", require_footer=True):
        raise MagpieError(f"eggNOG-mapper did not produce a complete annotation for {gid}")


def eggnog_step(*, prep_dir: Path, package_ref_dir: Path, out: Path, eggnog_existing_dir: Optional[Path],
                eggnog_data_dir: Optional[Path], emapper_bin: str, cpus: int, jobs: int) -> None:
    refs = {d: read_fasta_ids(package_ref_dir / REF_NAME[d] / f"{REF_NAME[d]}.fna")
            for d in read_id_list(package_ref_dir / "domains_built.txt")}
    domain_of = {g: d for d, ids in refs.items() for g in ids}
    shutil.rmtree(out / "annotations", ignore_errors=True)
    ann_dir = ensure_dir(out / "annotations")
    emapper_dir = ensure_dir(out / "emapper")

    existing = find_existing_annotations(eggnog_existing_dir, IdLookup(prep_dir / "id_map.tsv")) \
        if eggnog_existing_dir is not None else {}
    status: Dict[str, str] = {}
    todo: List[str] = []
    for gid in sorted(domain_of):
        if gid in existing:
            src = existing[gid]
            if not annotation_complete(src, require_footer=False):
                raise MagpieError(f"{src} has no '#query' header line; is it an eggNOG-mapper annotation file?")
            dst = ann_dir / f"{gid}{SUFFIX}"
            if src.name.endswith(".gz"):
                with gzip.open(src, "rb") as fin, dst.open("wb") as fout:
                    shutil.copyfileobj(fin, fout)
            else:
                dst.symlink_to(src.resolve())
            status[gid] = "existing"
        elif annotation_complete(emapper_dir / f"{gid}{SUFFIX}", require_footer=True):
            status[gid] = "annotated_earlier"
        else:
            todo.append(gid)

    if todo:
        if eggnog_data_dir is None:
            raise MagpieError(f"{len(todo)} reference genomes have no eggNOG annotation (e.g. {', '.join(todo[:5])}). "
                              "Add them to --eggnog-existing-dir, give --eggnog-data-dir to run eggNOG-mapper, or "
                              "rerun from the qc step with --eggnog-existing-dir so that MAGs without an annotation "
                              "are excluded before the reference is built (magpie run --from qc ...).")
        jobs = max(1, min(jobs, len(todo)))
        threads = max(1, cpus // jobs)
        _LOGGER.info("Annotating %d genomes with eggNOG-mapper (%d jobs x %d threads)", len(todo), jobs, threads)
        failed = []
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futs = {pool.submit(_emapper_one, g, prep_dir / "mags" / f"{g}_genomic.fna.gz", emapper_dir,
                                emapper_bin, eggnog_data_dir, threads): g for g in todo}
            for i, fut in enumerate(as_completed(futs), start=1):
                g = futs[fut]
                try:
                    fut.result()
                    status[g] = "annotated"
                    _LOGGER.info("eggNOG-mapper %d/%d: %s", i, len(todo), g)
                except (ShellError, MagpieError) as e:
                    failed.append(g)
                    _LOGGER.error("eggNOG-mapper failed for %s: %s", g, e)
        if failed:
            raise MagpieError(f"eggNOG-mapper failed for {len(failed)} genomes; rerun to retry them "
                              f"(finished genomes are kept). Logs: {emapper_dir / 'logs'}")

    for gid, st in status.items():
        if st in ("annotated", "annotated_earlier"):
            (ann_dir / f"{gid}{SUFFIX}").symlink_to((emapper_dir / f"{gid}{SUFFIX}").resolve())
    write_tsv(out / "annotation_sources.tsv", ["genome_id", "domain", "source"],
              [[g, domain_of[g], status[g]] for g in sorted(status)])
    counts: Dict[str, int] = {}
    for st in status.values():
        counts[st] = counts.get(st, 0) + 1
    write_json(out / "report.json", {"n_reference_genomes": len(domain_of), "sources": counts})
    _LOGGER.info("eggNOG annotations ready for %d reference genomes (%s)", len(status),
                 ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())))
