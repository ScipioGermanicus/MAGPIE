from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

_LOGGER = logging.getLogger("magpie.raxml_check")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _has_fasta_records(path: Path) -> bool:
    if not path.exists():
        return False
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                return True
    return False


def _read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    header: str | None = None
    seq_parts: List[str] = []

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:].strip().split()[0]
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA without header in {path}")
                seq_parts.append(line)

    if header is not None:
        records.append((header, "".join(seq_parts)))

    return records


def _write_relaxed_phylip(records: List[Tuple[str, str]], out_path: Path) -> None:
    if not records:
        out_path.write_text("", encoding="utf-8")
        return

    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(
            f"Cannot write PHYLIP: sequences have unequal lengths in alignment for {out_path}"
        )

    aln_len = next(iter(lengths))

    with out_path.open("w", encoding="utf-8") as f:
        f.write(f"{len(records)} {aln_len}\n")
        for rec_id, seq in records:
            f.write(f"{rec_id} {seq}\n")


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


def _run_raxml_check(
    *,
    msa: Path,
    prefix: Path,
    raxml_ng_bin: str,
    threads: int,
) -> None:
    cmd = [
        raxml_ng_bin,
        "--check",
        "--msa",
        str(msa),
        "--model",
        "GTR+G",
        "--prefix",
        str(prefix),
        "--threads",
        str(threads),
    ]
    _run_command(cmd)


def _check_one_domain(
    *,
    domain: str,
    msa_fasta: Path,
    out_dir: Path,
    raxml_ng_bin: str,
    threads: int,
) -> Dict[str, object]:
    prefix = out_dir / f"{domain}_raxml-check"
    reduced_phy = out_dir / f"{domain}_raxml-check.raxml.reduced.phy"
    log_path = out_dir / f"{domain}_raxml-check.raxml.log"

    if not _has_fasta_records(msa_fasta):
        _LOGGER.warning("Skipping %s raxml-check: no sequences found in %s", domain, msa_fasta)
        return {
            "domain": domain,
            "status": "skipped_empty_input",
            "msa_fasta": str(msa_fasta),
            "prefix": str(prefix),
            "reduced_phy": str(reduced_phy),
            "log_path": str(log_path),
            "n_sequences": 0,
            "alignment_length": 0,
            "phy_source": "none",
        }

    records = _read_fasta(msa_fasta)
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Alignment for {domain} has unequal sequence lengths: {msa_fasta}")
    aln_len = next(iter(lengths))

    _LOGGER.info("Running raxml-ng --check for %s", domain)
    _run_raxml_check(
        msa=msa_fasta,
        prefix=prefix,
        raxml_ng_bin=raxml_ng_bin,
        threads=threads,
    )

    phy_source = "raxml_reduced"
    if not reduced_phy.exists():
        _LOGGER.warning(
            "RAxML-NG did not produce reduced PHYLIP for %s; writing PHYLIP directly from aligned FASTA.",
            domain,
        )
        _write_relaxed_phylip(records, reduced_phy)
        phy_source = "manual_from_fasta"

    return {
        "domain": domain,
        "status": "processed",
        "msa_fasta": str(msa_fasta),
        "prefix": str(prefix),
        "reduced_phy": str(reduced_phy),
        "log_path": str(log_path),
        "n_sequences": len(records),
        "alignment_length": aln_len,
        "phy_source": phy_source,
    }


def raxml_check_step(
    *,
    choose_best_dir: Path,
    out: Path,
    raxml_ng_bin: str,
    threads: int,
    force: bool,
) -> None:
    """
    Run raxml-ng --check on the best aligned 16S FASTAs and ensure a PHYLIP output exists.

    Inputs:
      choose_best_dir/archaea/archaea_16S_centroids_ssu_align_best.fna
      choose_best_dir/bacteria/bacteria_16S_centroids_ssu_align_best.fna
    """
    arch_msa = choose_best_dir / "archaea" / "archaea_16S_centroids_ssu_align_best.fna"
    bact_msa = choose_best_dir / "bacteria" / "bacteria_16S_centroids_ssu_align_best.fna"

    missing = [str(p) for p in [arch_msa, bact_msa] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required choose-best FASTAs for raxml-check:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"RAxML-check outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "choose_best_dir": str(choose_best_dir),
            "raxml_ng_bin": raxml_ng_bin,
            "threads": threads,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, msa in (("archaea", arch_msa), ("bacteria", bact_msa)):
        stats = _check_one_domain(
            domain=domain,
            msa_fasta=msa,
            out_dir=out,
            raxml_ng_bin=raxml_ng_bin,
            threads=threads,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats["n_sequences"]),
            str(stats["alignment_length"]),
            str(stats["phy_source"]),
            str(stats["msa_fasta"]),
            str(stats["reduced_phy"]),
            str(stats["log_path"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_sequences",
            "alignment_length",
            "phy_source",
            "msa_fasta",
            "reduced_phy",
            "log_path",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("RAxML-check complete.")