from __future__ import annotations

import csv
import json
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

_LOGGER = logging.getLogger("magpie.eggnog")

FASTA_EXTS = (".fa", ".fna", ".fasta")
GZ_EXT = ".gz"


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _strip_fasta_suffix(name: str) -> str:
    s = name.strip()
    if s.lower().endswith(GZ_EXT):
        s = s[:-len(GZ_EXT)]
    for ext in FASTA_EXTS:
        if s.lower().endswith(ext):
            s = s[:-len(ext)]
            break
    return s


def _iter_genome_files(dir_: Path) -> Iterable[Path]:
    if not dir_.exists():
        return
    for p in sorted(dir_.iterdir(), key=lambda x: x.name):
        if not p.is_file():
            continue
        low = p.name.lower()
        if any(low.endswith(ext) for ext in FASTA_EXTS) or any(low.endswith(ext + GZ_EXT) for ext in FASTA_EXTS):
            yield p


def _run_command(cmd: list[str]) -> None:
    _LOGGER.debug("Running command: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"Command failed:\n"
            f"  {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\n"
            f"stderr:\n{proc.stderr}"
        )


def _find_annotation_files(dir_: Path) -> Dict[str, Path]:
    """
    Return mapping:
      genome_id -> <dir>/<genome_id>.emapper.annotations
    """
    out: Dict[str, Path] = {}
    if not dir_.exists():
        return out

    for p in sorted(dir_.iterdir(), key=lambda x: x.name):
        if not p.is_file():
            continue
        if not p.name.endswith(".emapper.annotations"):
            continue
        gid = p.name[: -len(".emapper.annotations")]
        if gid in out:
            raise ValueError(f"Duplicate eggNOG annotation for genome_id '{gid}' in {dir_}")
        out[gid] = p
    return out


def _copy_existing_annotations(
    *,
    genome_ids: List[str],
    existing_dir: Path,
    out_annotations_dir: Path,
) -> Tuple[int, List[str]]:
    found = _find_annotation_files(existing_dir)
    copied = 0
    missing: List[str] = []

    _ensure_dir(out_annotations_dir)

    for gid in genome_ids:
        src = found.get(gid)
        if src is None:
            missing.append(gid)
            continue
        dst = out_annotations_dir / src.name
        shutil.copy2(src, dst)
        copied += 1

    return copied, missing


def _run_prodigal(
    *,
    genome_fp: Path,
    proteins_out: Path,
    prodigal_bin: str,
) -> None:
    cmd = [
        prodigal_bin,
        "-i",
        str(genome_fp),
        "-a",
        str(proteins_out),
        "-p",
        "single",
    ]
    _run_command(cmd)


def _run_emapper(
    *,
    proteins_faa: Path,
    sample_id: str,
    output_dir: Path,
    emapper_bin: str,
    eggnog_data_dir: Path,
    cpus: int,
) -> None:
    cmd = [
        emapper_bin,
        "-m",
        "diamond",
        "--data_dir",
        str(eggnog_data_dir),
        "-i",
        str(proteins_faa),
        "-o",
        sample_id,
        "--output_dir",
        str(output_dir),
        "--cpu",
        str(cpus),
    ]
    _run_command(cmd)


def eggnog_step(
    *,
    prep_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    eggnog_existing_dir: Path | None,
    eggnog_data_dir: Path | None,
    prodigal_bin: str,
    emapper_bin: str,
) -> None:
    """
    Either reuse existing eggNOG annotations or run Prodigal + eggNOG-mapper
    on prepared MAGs under prep_dir/mags.

    Outputs:
      out/proteins/*.faa
      out/annotations/*.emapper.annotations
      out/summary.tsv
      out/report.json
    """
    mags_dir = prep_dir / "mags"
    if not mags_dir.exists():
        raise FileNotFoundError(f"Prepared MAG directory not found: {mags_dir}")

    proteins_dir = out / "proteins"
    annotations_dir = out / "annotations"

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"eggNOG outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)
    _ensure_dir(proteins_dir)
    _ensure_dir(annotations_dir)

    genome_files = list(_iter_genome_files(mags_dir))
    if not genome_files:
        raise ValueError(f"No prepared MAG FASTA files found in: {mags_dir}")

    genome_ids = [_strip_fasta_suffix(p.name) for p in genome_files]

    rows: List[List[str]] = []
    report = {
        "inputs": {
            "prep_dir": str(prep_dir),
            "mags_dir": str(mags_dir),
            "eggnog_existing_dir": str(eggnog_existing_dir) if eggnog_existing_dir is not None else None,
            "eggnog_data_dir": str(eggnog_data_dir) if eggnog_data_dir is not None else None,
            "prodigal_bin": prodigal_bin,
            "emapper_bin": emapper_bin,
            "cpus": cpus,
        },
        "mode": "",
        "counts": {
            "total_genomes": len(genome_ids),
            "completed": 0,
            "missing_existing_annotations": 0,
        },
        "outputs": {
            "proteins_dir": str(proteins_dir),
            "annotations_dir": str(annotations_dir),
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    if eggnog_existing_dir is not None:
        report["mode"] = "reuse_existing_annotations"
        copied, missing = _copy_existing_annotations(
            genome_ids=genome_ids,
            existing_dir=eggnog_existing_dir,
            out_annotations_dir=annotations_dir,
        )

        report["counts"]["completed"] = copied
        report["counts"]["missing_existing_annotations"] = len(missing)

        for gid in genome_ids:
            ann = annotations_dir / f"{gid}.emapper.annotations"
            status = "reused" if ann.exists() else "missing_existing_annotation"
            rows.append([
                gid,
                status,
                "",
                str(ann if ann.exists() else ""),
            ])

    else:
        if eggnog_data_dir is None:
            raise ValueError(
                "When not using --eggnog-existing-dir, you must provide --eggnog-data-dir "
                "so MAGPIE can run eggNOG-mapper."
            )

        report["mode"] = "run_prodigal_and_emapper"

        completed = 0
        for genome_fp in genome_files:
            gid = _strip_fasta_suffix(genome_fp.name)
            proteins_out = proteins_dir / f"{gid}.faa"
            annotation_out = annotations_dir / f"{gid}.emapper.annotations"

            _run_prodigal(
                genome_fp=genome_fp,
                proteins_out=proteins_out,
                prodigal_bin=prodigal_bin,
            )

            _run_emapper(
                proteins_faa=proteins_out,
                sample_id=gid,
                output_dir=annotations_dir,
                emapper_bin=emapper_bin,
                eggnog_data_dir=eggnog_data_dir,
                cpus=cpus,
            )

            if not annotation_out.exists():
                raise FileNotFoundError(
                    f"eggNOG-mapper did not produce expected annotation file: {annotation_out}"
                )

            completed += 1
            rows.append([
                gid,
                "annotated",
                str(proteins_out),
                str(annotation_out),
            ])

        report["counts"]["completed"] = completed

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome_id", "status", "proteins_faa", "annotation_file"])
        w.writerows(rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("eggNOG step complete. Completed=%d", report["counts"]["completed"])