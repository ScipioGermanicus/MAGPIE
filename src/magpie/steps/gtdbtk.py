"""
GTDB-Tk classification (script 02_infer_taxonomy.sh), or reuse of existing results.

Run mode:  gtdbtk classify_wf --genome_dir <prep>/mags --extension fna.gz --out_dir <out>/gtdbtk_out
           --cpus N [--pplacer_cpus M] [--skip_ani_screen]   (the flag exists in GTDB-Tk 2.2-2.6 only)
Reuse mode (--gtdb-classify DIR): DIR contains gtdbtk.bac120.summary.tsv and/or
           gtdbtk.ar53.summary.tsv (a GTDB-Tk out dir or its classify/ subdir). Genome names may be
           MAGPIE IDs, original file names or file stems; they are mapped to MAGPIE IDs.

Output: classification.tsv (genome_id, marker_set, classification) and copies of the summaries.
"""
from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path
from typing import Final, List, Optional

from ..pipeline import MagpieError
from ..util.ids import IdLookup
from ..util.io import ensure_dir, open_text, write_tsv
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.gtdbtk")

MARKER_SETS: Final = ("bac120", "ar53")


def find_summary(gtdb_dir: Path, marker: str) -> Optional[Path]:
    for d in (gtdb_dir, gtdb_dir / "classify"):
        p = d / f"gtdbtk.{marker}.summary.tsv"
        if p.is_file():
            return p
    return None


def _run_gtdbtk(*, mags_dir: Path, out_dir: Path, gtdbtk_bin: str, cpus: int, pplacer_cpus: Optional[int]) -> Path:
    help_text = run_cmd([gtdbtk_bin, "classify_wf", "-h"], check=False)
    flags = (help_text.stdout or "") + (help_text.stderr or "")
    cmd = [gtdbtk_bin, "classify_wf", "--genome_dir", str(mags_dir), "--extension", "fna.gz",
           "--out_dir", str(out_dir), "--cpus", str(cpus), "--tmpdir", str(ensure_dir(out_dir.parent / "tmp"))]
    if pplacer_cpus:
        cmd += ["--pplacer_cpus", str(pplacer_cpus)]
    if "--skip_ani_screen" in flags:
        cmd.append("--skip_ani_screen")
    _LOGGER.info("Running GTDB-Tk classify_wf (needs GTDBTK_DATA_PATH and >=128 GB RAM)")
    run_cmd(cmd, log_path=out_dir.parent / "gtdbtk.log")
    return out_dir


def gtdbtk_step(*, prep_dir: Path, out: Path, gtdb_classify: Optional[Path], gtdbtk_bin: str, cpus: int,
                pplacer_cpus: Optional[int]) -> None:
    if gtdb_classify is not None:
        src = gtdb_classify
        _LOGGER.info("Reusing GTDB-Tk results from %s", src)
    else:
        src = _run_gtdbtk(mags_dir=prep_dir / "mags", out_dir=out / "gtdbtk_out", gtdbtk_bin=gtdbtk_bin,
                          cpus=cpus, pplacer_cpus=pplacer_cpus)

    lookup = IdLookup(prep_dir / "id_map.tsv")
    wanted = set(lookup.ids)
    rows: List[List[str]] = []
    seen = set()
    n_other = 0
    found = False
    for marker in MARKER_SETS:
        summary = find_summary(src, marker)
        if summary is None:
            continue
        found = True
        shutil.copyfile(summary, out / summary.name)
        with open_text(summary) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames or not {"user_genome", "classification"} <= set(reader.fieldnames):
                raise MagpieError(f"{summary} lacks user_genome/classification columns")
            for row in reader:
                gid = lookup.resolve(row["user_genome"])
                if gid is None:
                    n_other += 1  # e.g. genomes of a larger GTDB-Tk run that are not among our MAGs
                    continue
                if gid not in seen:
                    seen.add(gid)
                    rows.append([gid, marker, row["classification"].strip()])
    if not found:
        raise MagpieError(f"No gtdbtk.bac120.summary.tsv / gtdbtk.ar53.summary.tsv found in {src} or {src}/classify")
    if n_other:
        _LOGGER.info("Ignored %d GTDB-Tk rows for genomes that are not among the input MAGs", n_other)
    missing = sorted(wanted - seen)
    if missing:
        _LOGGER.warning("%d input genomes are not in the GTDB-Tk summaries and are excluded (e.g. %s)",
                        len(missing), ", ".join(missing[:5]))
    if not rows:
        raise MagpieError("None of the input genomes were found in the GTDB-Tk summaries. Check that GTDB-Tk was "
                          "run on these MAGs (names must match file names or the IDs in id_map.tsv).")
    write_tsv(out / "classification.tsv", ["genome_id", "marker_set", "classification"], sorted(rows))
    _LOGGER.info("GTDB-Tk classification available for %d of %d genomes", len(rows), len(wanted))
