"""
Genome quality filter (scripts 04_quality_control.sh, 04b_extract.sh, 05_quality_filtering.py;
wiki: "Filter to include only genomes >90% completion and <10% redundancy").

CheckM values are either computed here (checkm lineage_wf + checkm qa, per domain as in script 04)
or read from existing CheckM output (--checkm-qa, or --checkm-qa-bacteria/--checkm-qa-archaea).
Accepted formats: `checkm qa -o 2` fixed-width tables, tab-separated QA tables (--tab_table or the
merged table of script 04b) and minimal 3-column tables (genome_id, checkm_completeness,
checkm_contamination). Genome names may be MAGPIE IDs, original file names or stems.

The domain of each genome comes from the GTDB-Tk-based taxonomy step. When eggNOG annotations are
only reused (--eggnog-existing-dir without --eggnog-data-dir), genomes without an annotation are
excluded here, so that every reference genome has functional traits.
"""
from __future__ import annotations

import gzip
import logging
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from ..pipeline import MagpieError
from ..util.ids import IdLookup
from ..util.io import ensure_dir, open_text, write_id_list, write_json, write_tsv
from ..util.shell import run_cmd

_LOGGER = logging.getLogger("magpie.qc")

_SEP_RE = re.compile(r"^\s*-{5,}\s*$")
_SPACES_RE = re.compile(r"\s{2,}")


@dataclass(frozen=True)
class CheckMRow:
    name: str
    marker_lineage: str
    completeness: float
    contamination: float


def _col(cols: List[str], *names: str) -> Optional[int]:
    low = [c.strip().lower() for c in cols]
    for n in names:
        if n in low:
            return low.index(n)
    return None


def parse_checkm_table(path: Path) -> List[CheckMRow]:
    """Parse CheckM QA output in any of the formats used by the MAGPIE scripts."""
    with open_text(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() and not _SEP_RE.match(ln)]
    if not lines:
        raise MagpieError(f"CheckM table is empty: {path}")
    if "\t" in lines[0]:
        split = lambda s: s.split("\t")  # noqa: E731
    else:
        split = lambda s: _SPACES_RE.split(s.strip())  # noqa: E731

    header = split(lines[0])
    i_id = _col(header, "bin id", "genome_id", "name", "genome", "user_genome")
    i_comp = _col(header, "completeness", "checkm_completeness")
    i_cont = _col(header, "contamination", "checkm_contamination")
    i_lin = _col(header, "marker lineage", "marker_lineage")
    if None in (i_id, i_comp, i_cont):
        raise MagpieError(f"No genome ID / Completeness / Contamination columns in {path} (header: {header})")

    rows, seen = [], set()
    for ln in lines[1:]:
        parts = split(ln)
        if len(parts) <= max(i_id, i_comp, i_cont):
            continue
        try:
            row = CheckMRow(parts[i_id].strip(), parts[i_lin].strip() if i_lin is not None else "",
                            float(parts[i_comp]), float(parts[i_cont]))
        except ValueError:
            continue
        if row.name not in seen:
            seen.add(row.name)
            rows.append(row)
    if not rows:
        raise MagpieError(f"No rows could be parsed from CheckM table {path}")
    return rows


def _run_checkm(*, genome_dirs: Dict[str, Path], out: Path, checkm_bin: str, cpus: int) -> List[Path]:
    """checkm lineage_wf + qa per domain on unzipped genomes (script 04)."""
    tables = []
    for dom, src in genome_dirs.items():
        files = sorted(src.glob("*.fna.gz"))
        if not files:
            continue
        work = ensure_dir(out / "checkm" / "unzipped" / dom)
        for f in files:
            with gzip.open(f, "rb") as fin, (work / f.name[:-3]).open("wb") as fout:
                shutil.copyfileobj(fin, fout)
        res = out / "checkm" / dom
        run_cmd([checkm_bin, "lineage_wf", "-x", "fna", "-t", str(cpus), str(work), str(res)],
                log_path=out / "checkm" / f"{dom}.log")
        table = res / "checkm_qa.tsv"
        run_cmd([checkm_bin, "qa", str(res / "lineage.ms"), str(res), "-o", "2", "--tab_table", "-f", str(table)],
                log_path=out / "checkm" / f"{dom}.log")
        shutil.rmtree(work)
        tables.append(table)
    return tables


def qc_step(*, prep_dir: Path, tax_dir: Path, out: Path, checkm_tables: List[Path], checkm_bin: str, cpus: int,
            completeness_min: float, contamination_max: float,
            annotated: Optional[set] = None) -> None:
    domains = {d: [ln for ln in (tax_dir / f"{d}.txt").read_text().split() if ln] for d in ("bacteria", "archaea")}
    domain_of = {g: d for d, ids in domains.items() for g in ids}

    if not checkm_tables:
        _LOGGER.info("Running CheckM lineage_wf on %d genomes", len(domain_of))
        checkm_tables = _run_checkm(genome_dirs={d: tax_dir / "genomes_to_search_barrnap" / d for d in domains},
                                    out=out, checkm_bin=checkm_bin, cpus=cpus)

    lookup = IdLookup(prep_dir / "id_map.tsv")
    quality: Dict[str, CheckMRow] = {}
    n_other = 0
    for table in checkm_tables:
        for row in parse_checkm_table(table):
            gid = lookup.resolve(row.name)
            if gid is None:
                n_other += 1
            else:
                quality.setdefault(gid, row)
    if n_other:
        _LOGGER.info("Ignored %d CheckM rows for genomes that are not among the input MAGs", n_other)
    missing = sorted(g for g in domain_of if g not in quality)
    if missing:
        raise MagpieError(f"{len(missing)} genomes have no CheckM values (e.g. {', '.join(missing[:5])}). "
                          f"Tables used: {', '.join(map(str, checkm_tables))}")

    all_rows, kept, no_annotation = [], {"bacteria": [], "archaea": []}, []
    for gid in sorted(domain_of):
        q, dom = quality[gid], domain_of[gid]
        status = "pass"
        if not (q.completeness >= completeness_min and q.contamination <= contamination_max):
            status = "fail_quality"
        elif annotated is not None and gid not in annotated:
            status = "no_annotation"
            no_annotation.append(gid)
        all_rows.append([gid, q.marker_lineage, f"{q.completeness:.2f}", f"{q.contamination:.2f}",
                         dom.capitalize(), status])
        if status == "pass":
            kept[dom].append(gid)

    header = ["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"]
    write_tsv(out / "checkm_clean_all.tsv", header + ["status"], all_rows)
    write_tsv(out / "checkm_filtered.tsv", header, [r[:5] for r in all_rows if r[5] == "pass"])
    for dom, ids in kept.items():
        write_id_list(out / f"{dom}.txt", ids)
        bins = ensure_dir(out / "bins" / dom)
        for gid in ids:
            (bins / f"{gid}_genomic.fna.gz").symlink_to((prep_dir / "mags" / f"{gid}_genomic.fna.gz").resolve())

    n_fail = sum(1 for r in all_rows if r[5] == "fail_quality")
    _LOGGER.info("QC: %d bacteria and %d archaea pass (completeness >= %g, contamination <= %g); %d fail",
                 len(kept["bacteria"]), len(kept["archaea"]), completeness_min, contamination_max, n_fail)
    if no_annotation:
        _LOGGER.warning("%d genomes passed QC but have no eggNOG annotation and are excluded (e.g. %s)",
                        len(no_annotation), ", ".join(no_annotation[:5]))
    write_json(out / "report.json", {"checkm_tables": [str(t) for t in checkm_tables],
                                     "completeness_min": completeness_min, "contamination_max": contamination_max,
                                     "kept": {d: len(v) for d, v in kept.items()}, "failed_quality": n_fail,
                                     "excluded_no_annotation": no_annotation})
    if not any(kept.values()):
        raise MagpieError("No genomes passed quality control")
