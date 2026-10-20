"""
Assemble and check the final PICRUSt2 reference database (layout of gtdb_r226_picrust_ref):

  picrust2_db/bac_ref/bac_ref.{fna,hmm,tre,model,raxml_info}, ko.txt.gz, ec.txt.gz
  picrust2_db/arc_ref/...                     (only if an archaeal reference could be built)
  picrust2_db/bacteria_16S_copies.txt, archaea_16S_copies.txt
  picrust2_db/run_picrust2.sh                 PICRUSt2 (v2.6+) command using this database
"""
from __future__ import annotations

import logging
import re
import shutil
from pathlib import Path
from typing import List

from ..pipeline import MagpieError
from ..util.io import ensure_dir, open_text, read_fasta_ids, read_id_list, write_json

_LOGGER = logging.getLogger("magpie.picrust2_db")

REF_NAME = {"bacteria": "bac_ref", "archaea": "arc_ref"}
TABLES = ("ec.txt.gz", "ko.txt.gz")  # EC first: PICRUSt2 infers pathways from the first custom table


def tree_tips(path: Path) -> List[str]:
    """Leaf labels of a Newick tree (labels directly after '(' or ',')."""
    text = path.read_text(encoding="utf-8")
    return re.findall(r"[(,]\s*'?([^():,;']+?)'?\s*(?=[:),])", text)


def _table_rows(path: Path) -> List[str]:
    with open_text(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header[0] != "assembly":
            raise MagpieError(f"First column of {path} must be 'assembly'")
        return [ln.split("\t", 1)[0] for ln in fh if ln.strip()]


def check_db(db: Path) -> List[str]:
    """Return a list of problems with a database directory (empty if consistent)."""
    problems = []
    for domain, ref in REF_NAME.items():
        rdir = db / ref
        if not rdir.exists():
            continue
        files = {e: rdir / f"{ref}.{e}" for e in ("fna", "hmm", "tre", "model", "raxml_info")}
        absent = [str(p) for p in files.values() if not p.exists() or p.stat().st_size == 0]
        if absent:
            problems += [f"missing or empty: {p}" for p in absent]
            continue
        ids = read_fasta_ids(files["fna"])
        if len(set(ids)) != len(ids):
            problems.append(f"{ref}: duplicate IDs in {ref}.fna")
        if sorted(tree_tips(files["tre"])) != sorted(ids):
            problems.append(f"{ref}: tree tips differ from the sequences in {ref}.fna")
        nseq = re.search(r"^NSEQ\s+(\d+)", files["hmm"].read_text(encoding="utf-8", errors="replace"), re.M)
        if not nseq or int(nseq.group(1)) != len(ids):
            problems.append(f"{ref}: HMM NSEQ differs from the number of reference sequences")
        for table in [db / f"{domain}_16S_copies.txt"] + [rdir / t for t in TABLES]:
            if not table.exists():
                problems.append(f"missing: {table}")
            elif sorted(_table_rows(table)) != sorted(ids):
                problems.append(f"{ref}: genomes in {table.name} differ from {ref}.fna")
    return problems


def _run_script(db: Path, built: List[str]) -> str:
    lines = [
        "#!/usr/bin/env bash",
        "# Run PICRUSt2 (v2.6+) with this MAGPIE database. Run inside the PICRUSt2 conda environment:",
        "#   bash run_picrust2.sh <study_seqs.fna> <study_table.biom|.tsv> <output_dir> [threads] [extra options]",
        "set -euo pipefail",
        'DB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'SEQS="$1"; TABLE="$2"; OUT="$3"; THREADS="${4:-1}"; shift $(( $# < 4 ? $# : 4 ))',
    ]
    if len(built) < 2:
        lines.append('DEFAULTS="$(python -c \'import os, picrust2; '
                     'print(os.path.join(os.path.dirname(picrust2.__file__), "default_files"))\')"')
    cmd = ['picrust2_pipeline.py -s "$SEQS" -i "$TABLE" -o "$OUT" -p "$THREADS" \\']
    for n, domain in ((1, "bacteria"), (2, "archaea")):
        ref = REF_NAME[domain]
        if domain in built:
            tables = ",".join(f"$DB/{ref}/{t}" for t in TABLES)
            cmd.append(f'  --ref_dir{n} "$DB/{ref}" --marker_gene_table_ref{n} "$DB/{domain}_16S_copies.txt" \\')
        else:
            # Not enough genomes for this domain: PICRUSt2's default reference (and tables) are used.
            tables = ",".join(f"$DEFAULTS/{domain}/{t}" for t in TABLES)
        cmd.append(f'  --custom_trait_tables_ref{n} "{tables}" \\')
    # Custom tables are named after their file (ec.txt.gz -> "ec.txt"); pathway inference needs that name.
    cmd += ['  --reaction_func ec.txt \\', '  "$@"']
    return "\n".join(lines + cmd) + "\n"


def picrust2_db_step(*, package_ref_dir: Path, ko_dir: Path, ec_dir: Path, db: Path) -> None:
    built = read_id_list(package_ref_dir / "domains_built.txt")
    if db.exists():
        shutil.rmtree(db)
    for domain in built:
        ref = REF_NAME[domain]
        shutil.copytree(package_ref_dir / ref, db / ref)
        shutil.copyfile(package_ref_dir / f"{domain}_16S_copies.txt", db / f"{domain}_16S_copies.txt")
        for table_dir, name in ((ko_dir, "ko.txt.gz"), (ec_dir, "ec.txt.gz")):
            shutil.copyfile(table_dir / ref / name, db / ref / name)
            for extra in (table_dir / ref).glob("*_constant_removed.tsv"):
                shutil.copyfile(extra, db / ref / extra.name)
    ensure_dir(db)
    (db / "run_picrust2.sh").write_text(_run_script(db, built), encoding="utf-8")
    (db / "run_picrust2.sh").chmod(0o755)

    problems = check_db(db)
    write_json(db / "magpie_db_report.json", {"domains": built, "problems": problems})
    if problems:
        raise MagpieError("The assembled database is inconsistent:\n" + "\n".join(f"  * {p}" for p in problems))
    missing = [d for d in REF_NAME if d not in built]
    _LOGGER.info("PICRUSt2 database: %s (%s)%s", db, ", ".join(REF_NAME[d] for d in built),
                 f"; PICRUSt2's default reference will be used for {', '.join(missing)}" if missing else "")
