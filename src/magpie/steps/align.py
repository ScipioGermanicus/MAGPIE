from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List

_LOGGER = logging.getLogger("magpie.align")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _read_fasta_lengths(path: Path) -> List[int]:
    lengths: List[int] = []
    current_len = 0
    seen_header = False

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header:
                    lengths.append(current_len)
                seen_header = True
                current_len = 0
            else:
                if not seen_header:
                    raise ValueError(f"Malformed FASTA without header in {path}")
                current_len += len(line)

    if seen_header:
        lengths.append(current_len)

    return lengths


def _has_fasta_records(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                return True
    return False


def _run_command(cmd: list[str], stdout_path: Path | None = None) -> None:
    _LOGGER.debug("Running command: %s", " ".join(cmd))

    if stdout_path is None:
        proc = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    else:
        with stdout_path.open("w", encoding="utf-8") as fout:
            proc = subprocess.run(
                cmd,
                stdout=fout,
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


def _align_one_domain(
    *,
    domain: str,
    centroid_fasta: Path,
    model_path: Path,
    out_dir: Path,
    cmalign_bin: str,
    esl_reformat_bin: str,
    cpus: int,
    mxsize: int,
) -> Dict[str, object]:
    prefix = f"{domain}_16S_centroids"
    stk_out = out_dir / f"{prefix}.stk"
    afa_out = out_dir / f"{prefix}_ssu_align.fna"

    if not _has_fasta_records(centroid_fasta):
        _LOGGER.warning("Skipping %s alignment: no sequences found in %s", domain, centroid_fasta)
        return {
            "domain": domain,
            "input_fasta": str(centroid_fasta),
            "model": str(model_path),
            "stockholm_out": str(stk_out),
            "aligned_fasta_out": str(afa_out),
            "status": "skipped_empty_input",
            "n_sequences": 0,
            "min_len": 0,
            "max_len": 0,
            "mean_len": 0.0,
            "median_len": 0.0,
        }

    _LOGGER.info("Aligning %s centroids with model: %s", domain, model_path)

    _run_command(
        [
            cmalign_bin,
            "--dna",
            "--matchonly",
            "--mxsize",
            str(mxsize),
            "--cpu",
            str(cpus),
            str(model_path),
            str(centroid_fasta),
        ],
        stdout_path=stk_out,
    )

    _run_command(
        [
            esl_reformat_bin,
            "-o",
            str(afa_out),
            "afa",
            str(stk_out),
        ]
    )

    lengths = _read_fasta_lengths(afa_out)
    n = len(lengths)

    if n:
        sorted_lengths = sorted(lengths)
        mid = n // 2
        if n % 2 == 1:
            median_len = float(sorted_lengths[mid])
        else:
            median_len = (sorted_lengths[mid - 1] + sorted_lengths[mid]) / 2.0
    else:
        median_len = 0.0

    return {
        "domain": domain,
        "input_fasta": str(centroid_fasta),
        "model": str(model_path),
        "stockholm_out": str(stk_out),
        "aligned_fasta_out": str(afa_out),
        "status": "aligned",
        "n_sequences": n,
        "min_len": min(lengths) if lengths else 0,
        "max_len": max(lengths) if lengths else 0,
        "mean_len": (sum(lengths) / n) if n else 0.0,
        "median_len": median_len,
    }


def align_step(
    *,
    rrna_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    cmalign_bin: str,
    esl_reformat_bin: str,
    ssu_models_dir: Path,
    mxsize_archaea: int = 4096,
    mxsize_bacteria: int = 8192,
) -> None:
    """
    Align 16S centroid sequences from rrna outputs using cmalign and convert to aligned FASTA.

    Expected inputs:
      rrna_dir/archaea/archaea_16S_centroids.fasta
      rrna_dir/bacteria/bacteria_16S_centroids.fasta

    Required models:
      ssu_models_dir/archaea.cm
      ssu_models_dir/bacteria.cm
    """
    arch_centroids = rrna_dir / "archaea" / "archaea_16S_centroids.fasta"
    bact_centroids = rrna_dir / "bacteria" / "bacteria_16S_centroids.fasta"

    arch_model = ssu_models_dir / "archaea.cm"
    bact_model = ssu_models_dir / "bacteria.cm"

    missing = [str(p) for p in [arch_centroids, bact_centroids, arch_model, bact_model] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required align inputs/models:\n" + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [
        out / "archaea_16S_centroids.stk",
        out / "archaea_16S_centroids_ssu_align.fna",
        out / "bacteria_16S_centroids.stk",
        out / "bacteria_16S_centroids_ssu_align.fna",
        out / "summary.tsv",
        out / "report.json",
    ]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"Alignment outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    rows: List[List[str]] = []
    report = {
        "inputs": {
            "rrna_dir": str(rrna_dir),
            "ssu_models_dir": str(ssu_models_dir),
            "cmalign_bin": cmalign_bin,
            "esl_reformat_bin": esl_reformat_bin,
            "cpus": cpus,
            "mxsize_archaea": mxsize_archaea,
            "mxsize_bacteria": mxsize_bacteria,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    archaea_stats = _align_one_domain(
        domain="archaea",
        centroid_fasta=arch_centroids,
        model_path=arch_model,
        out_dir=out,
        cmalign_bin=cmalign_bin,
        esl_reformat_bin=esl_reformat_bin,
        cpus=cpus,
        mxsize=mxsize_archaea,
    )

    bacteria_stats = _align_one_domain(
        domain="bacteria",
        centroid_fasta=bact_centroids,
        model_path=bact_model,
        out_dir=out,
        cmalign_bin=cmalign_bin,
        esl_reformat_bin=esl_reformat_bin,
        cpus=cpus,
        mxsize=mxsize_bacteria,
    )

    for stats in (archaea_stats, bacteria_stats):
        report["domains"][stats["domain"]] = stats
        rows.append([
            str(stats["domain"]),
            str(stats["status"]),
            str(stats["n_sequences"]),
            str(stats["min_len"]),
            str(stats["max_len"]),
            f"{float(stats['mean_len']):.3f}",
            f"{float(stats['median_len']):.3f}",
            str(stats["input_fasta"]),
            str(stats["model"]),
            str(stats["stockholm_out"]),
            str(stats["aligned_fasta_out"]),
        ])

    with (out / "summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_sequences",
            "min_len",
            "max_len",
            "mean_len",
            "median_len",
            "input_fasta",
            "model",
            "stockholm_out",
            "aligned_fasta_out",
        ])
        w.writerows(rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info(
        "Alignment complete. Archaea status=%s n=%d; Bacteria status=%s n=%d",
        str(archaea_stats["status"]),
        int(archaea_stats["n_sequences"]),
        str(bacteria_stats["status"]),
        int(bacteria_stats["n_sequences"]),
    )