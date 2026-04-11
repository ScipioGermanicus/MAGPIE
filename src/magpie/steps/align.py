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


def _run_ssu_align(
    *,
    input_fasta: Path,
    output_dir: Path,
    ssu_align_bin: str,
    force: bool,
) -> None:
    cmd = [ssu_align_bin]
    if force:
        cmd.append("-f")
    cmd.extend(["--rfonly", str(input_fasta), str(output_dir)])
    _run_command(cmd)


def _run_esl_reformat(
    *,
    stk_path: Path,
    out_fasta: Path,
    esl_reformat_bin: str,
) -> None:
    cmd = [
        esl_reformat_bin,
        "-o",
        str(out_fasta),
        "afa",
        str(stk_path),
    ]
    _run_command(cmd)


def _align_one_domain(
    *,
    domain: str,
    input_fasta: Path,
    out: Path,
    ssu_align_bin: str,
    esl_reformat_bin: str,
    force: bool,
) -> Dict[str, object]:
    out_prefix_dir = out / f"{domain}_16S_centroids_ssu_align"
    out_fasta = out / f"{domain}_16S_centroids_ssu_align.fna"
    stk_path = out_prefix_dir / f"{domain}_16S_centroids_ssu_align.{domain}.stk"

    if not _has_fasta_records(input_fasta):
        _LOGGER.warning("Skipping %s alignment: no sequences found in %s", domain, input_fasta)
        return {
            "domain": domain,
            "input_fasta": str(input_fasta),
            "ssu_align_dir": str(out_prefix_dir),
            "stockholm_out": str(stk_path),
            "aligned_fasta_out": str(out_fasta),
            "status": "skipped_empty_input",
            "n_sequences": 0,
            "min_len": 0,
            "max_len": 0,
            "mean_len": 0.0,
            "median_len": 0.0,
        }

    _LOGGER.info("Aligning %s 16S centroids with ssu-align", domain)

    _run_ssu_align(
        input_fasta=input_fasta,
        output_dir=out_prefix_dir,
        ssu_align_bin=ssu_align_bin,
        force=force,
    )

    if not stk_path.exists():
        raise FileNotFoundError(
            f"Expected Stockholm alignment not found for {domain}: {stk_path}"
        )

    _run_esl_reformat(
        stk_path=stk_path,
        out_fasta=out_fasta,
        esl_reformat_bin=esl_reformat_bin,
    )

    lengths = _read_fasta_lengths(out_fasta)
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
        "input_fasta": str(input_fasta),
        "ssu_align_dir": str(out_prefix_dir),
        "stockholm_out": str(stk_path),
        "aligned_fasta_out": str(out_fasta),
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
    force: bool,
    ssu_align_bin: str,
    esl_reformat_bin: str,
) -> None:
    """
    Align 16S centroid sequences from rrna outputs using ssu-align and convert to aligned FASTA.

    Expected inputs:
      rrna_dir/archaea/archaea_16S_centroids.fasta
      rrna_dir/bacteria/bacteria_16S_centroids.fasta
    """
    arch_centroids = rrna_dir / "archaea" / "archaea_16S_centroids.fasta"
    bact_centroids = rrna_dir / "bacteria" / "bacteria_16S_centroids.fasta"

    missing = [str(p) for p in [arch_centroids, bact_centroids] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required align inputs:\n" + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [
        out / "archaea_16S_centroids_ssu_align.fna",
        out / "bacteria_16S_centroids_ssu_align.fna",
        out / "summary.tsv",
        out / "report.json",
    ]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"Alignment outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "rrna_dir": str(rrna_dir),
            "ssu_align_bin": ssu_align_bin,
            "esl_reformat_bin": esl_reformat_bin,
            "force": force,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    archaea_stats = _align_one_domain(
        domain="archaea",
        input_fasta=arch_centroids,
        out=out,
        ssu_align_bin=ssu_align_bin,
        esl_reformat_bin=esl_reformat_bin,
        force=force,
    )

    bacteria_stats = _align_one_domain(
        domain="bacteria",
        input_fasta=bact_centroids,
        out=out,
        ssu_align_bin=ssu_align_bin,
        esl_reformat_bin=esl_reformat_bin,
        force=force,
    )

    rows = []
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
            str(stats["ssu_align_dir"]),
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
            "ssu_align_dir",
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