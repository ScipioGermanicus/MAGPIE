from __future__ import annotations

import csv
import gzip
import logging
import shutil
from pathlib import Path

from ..util.shell import run_cmd as shell_run
from ..util.deps import check_checkm

_LOGGER = logging.getLogger("magpie.checkm")

FASTA_SUFFIXES = (".fa", ".fna", ".fasta")


def _ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def _iter_fasta_files(d: Path) -> list[Path]:
    files: list[Path] = []
    for p in d.iterdir():
        if not p.is_file():
            continue
        name = p.name.lower()
        if any(name.endswith(sfx) for sfx in FASTA_SUFFIXES) or any(name.endswith(sfx + ".gz") for sfx in FASTA_SUFFIXES):
            files.append(p)
    return sorted(files)


def _gunzip_to(src_gz: Path, dst: Path) -> None:
    with gzip.open(src_gz, "rb") as fin, open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _stage_fastas(src_dir: Path, dst_dir: Path) -> str:
    """
    Stage FASTAs into dst_dir, normalising to .fna.
    The output file stem is preserved (including dots), which is crucial for ID matching.
    """
    _ensure_dir(dst_dir)
    files = _iter_fasta_files(src_dir)
    if not files:
        raise FileNotFoundError(f"No FASTA(.gz) files found in: {src_dir}")

    for p in files:
        name = p.name
        if name.lower().endswith(".gz"):
            name = name[:-3]

        stem = name
        for sfx in FASTA_SUFFIXES:
            if stem.lower().endswith(sfx):
                stem = stem[: -len(sfx)]
                break

        out_fp = dst_dir / f"{stem}.fna"
        if p.name.lower().endswith(".gz"):
            _gunzip_to(p, out_fp)
        else:
            shutil.copy2(p, out_fp)

    return "fna"


def _write_min_tsv(checkm_qa_tsv: Path, out_min_tsv: Path) -> None:
    with open(checkm_qa_tsv, "r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        raise ValueError(f"Empty CheckM QA file: {checkm_qa_tsv}")

    header = rows[0]
    idx = {h: i for i, h in enumerate(header)}
    needed = ["Bin Id", "Completeness", "Contamination"]
    missing = [c for c in needed if c not in idx]
    if missing:
        raise ValueError(f"CheckM QA TSV missing columns {missing}. Seen: {header}")

    _ensure_dir(out_min_tsv.parent)
    with open(out_min_tsv, "w", encoding="utf-8", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["genome_id", "checkm_completeness", "checkm_contamination"])
        for r in rows[1:]:
            if not r or len(r) <= max(idx.values()):
                continue
            gid = r[idx["Bin Id"]].strip()
            cpl = r[idx["Completeness"]].strip()
            cnt = r[idx["Contamination"]].strip()
            if gid:
                w.writerow([gid, cpl, cnt])


def checkm_step(
    *,
    bacteria_dir: Path,
    archaea_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    checkm_bin: str = "checkm",
) -> Path:
    out = out.resolve()
    merged_min = out / "merged" / "checkm_results.min.tsv"
    if merged_min.exists() and not force:
        _LOGGER.info("Reusing existing merged CheckM TSV: %s", merged_min)
        return merged_min

    dep = check_checkm()
    if not dep.found:
        raise RuntimeError(
            "CheckM (checkm) not found on PATH. "
            "Install it (e.g., conda/mamba) and re-run. "
            "Example: `conda install -c conda-forge -c bioconda checkm-genome`."
        )

    _ensure_dir(out / "bacteria")
    _ensure_dir(out / "archaea")
    _ensure_dir(out / "merged")
    _ensure_dir(out / "staged" / "bacteria")
    _ensure_dir(out / "staged" / "archaea")

    ext = _stage_fastas(bacteria_dir, out / "staged" / "bacteria")
    _stage_fastas(archaea_dir, out / "staged" / "archaea")

    shell_run([checkm_bin, "lineage_wf", "-x", ext, "-t", str(cpus), str(out / "staged" / "bacteria"), str(out / "bacteria")])
    shell_run([checkm_bin, "lineage_wf", "-x", ext, "-t", str(cpus), str(out / "staged" / "archaea"), str(out / "archaea")])

    bac_qa = out / "bacteria" / "checkm_qa.tsv"
    arc_qa = out / "archaea" / "checkm_qa.tsv"
    shell_run([checkm_bin, "qa", str(out / "bacteria" / "lineage.ms"), str(out / "bacteria"), "-o", "2", "-f", str(bac_qa)])
    shell_run([checkm_bin, "qa", str(out / "archaea" / "lineage.ms"), str(out / "archaea"), "-o", "2", "-f", str(arc_qa)])

    bac_min = out / "bacteria" / "checkm_results.min.tsv"
    arc_min = out / "archaea" / "checkm_results.min.tsv"
    _write_min_tsv(bac_qa, bac_min)
    _write_min_tsv(arc_qa, arc_min)

    with open(merged_min, "w", encoding="utf-8", newline="") as out_fh:
        w = csv.writer(out_fh, delimiter="\t")
        w.writerow(["genome_id", "checkm_completeness", "checkm_contamination"])
        for fp in (bac_min, arc_min):
            with open(fp, "r", encoding="utf-8", newline="") as fh:
                r = csv.reader(fh, delimiter="\t")
                next(r, None)  # skip header
                for row in r:
                    if row:
                        w.writerow(row)

    return merged_min