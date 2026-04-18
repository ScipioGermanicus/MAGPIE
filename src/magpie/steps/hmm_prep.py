from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

_LOGGER = logging.getLogger("magpie.hmm_prep")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


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


def _read_relaxed_phylip(path: Path) -> Tuple[List[Tuple[str, str]], int]:
    """
    Read relaxed PHYLIP written by RAxML-NG.
    Returns:
      records: list[(name, seq)]
      aln_len: alignment length from header
    """
    if not path.exists():
        return ([], 0)

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().strip()
        if not header:
            return ([], 0)

        parts = header.split()
        if len(parts) < 2:
            raise ValueError(f"Could not parse PHYLIP header from: {path}")

        try:
            n_taxa = int(parts[0])
            aln_len = int(parts[1])
        except ValueError as e:
            raise ValueError(f"Invalid PHYLIP header in {path}: {header}") from e

        records: List[Tuple[str, str]] = []
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                continue

            name, seq = parts[0], parts[1].replace(" ", "")
            records.append((name, seq))

    if len(records) != n_taxa:
        _LOGGER.warning(
            "PHYLIP header for %s says %d taxa, but %d records were parsed.",
            path,
            n_taxa,
            len(records),
        )

    lengths = {len(seq) for _, seq in records}
    if records and len(lengths) != 1:
        raise ValueError(f"Unequal alignment lengths found in PHYLIP: {path}")

    return (records, aln_len)


def _write_fasta(records: List[Tuple[str, str]], out_path: Path) -> None:
    with out_path.open("w", encoding="utf-8") as f:
        for name, seq in records:
            f.write(f">{name}\n")
            f.write(f"{seq}\n")


def _prep_one_domain(
    *,
    domain: str,
    phy_path: Path,
    out_dir: Path,
    esl_reformat_bin: str,
) -> Dict[str, object]:
    fasta_out = out_dir / f"{domain}_raxml-check.raxml.reduced.fna"
    dna_fasta_out = out_dir / f"{domain}_raxml-check.raxml.reduced_dna.fna"
    stockholm_out = out_dir / f"{domain}_raxml-check.raxml.reduced_dna.sto"

    records, aln_len = _read_relaxed_phylip(phy_path)

    if not records:
        _LOGGER.warning("Skipping %s hmm-prep: PHYLIP missing or empty: %s", domain, phy_path)
        return {
            "domain": domain,
            "status": "skipped_empty_input",
            "input_phy": str(phy_path),
            "fasta_out": str(fasta_out),
            "dna_fasta_out": str(dna_fasta_out),
            "stockholm_out": str(stockholm_out),
            "n_taxa": 0,
            "alignment_length": 0,
        }

    _write_fasta(records, fasta_out)

    _run_command([
        esl_reformat_bin,
        "-d",
        "-o",
        str(dna_fasta_out),
        "afa",
        str(fasta_out),
    ])

    _run_command([
        esl_reformat_bin,
        "-o",
        str(stockholm_out),
        "stockholm",
        str(dna_fasta_out),
    ])

    return {
        "domain": domain,
        "status": "processed",
        "input_phy": str(phy_path),
        "fasta_out": str(fasta_out),
        "dna_fasta_out": str(dna_fasta_out),
        "stockholm_out": str(stockholm_out),
        "n_taxa": len(records),
        "alignment_length": aln_len,
    }


def hmm_prep_step(
    *,
    raxml_check_dir: Path,
    out: Path,
    esl_reformat_bin: str,
    force: bool,
) -> None:
    """
    Convert reduced PHYLIP alignments to FASTA, DNA FASTA, and Stockholm for HMM preparation.

    Inputs:
      raxml_check_dir/archaea_raxml-check.raxml.reduced.phy
      raxml_check_dir/bacteria_raxml-check.raxml.reduced.phy
    """
    arch_phy = raxml_check_dir / "archaea_raxml-check.raxml.reduced.phy"
    bact_phy = raxml_check_dir / "bacteria_raxml-check.raxml.reduced.phy"

    missing = [str(p) for p in [arch_phy, bact_phy] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required RAxML-check PHYLIP inputs for hmm-prep:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"HMM-prep outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "raxml_check_dir": str(raxml_check_dir),
            "esl_reformat_bin": esl_reformat_bin,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, phy in (("archaea", arch_phy), ("bacteria", bact_phy)):
        stats = _prep_one_domain(
            domain=domain,
            phy_path=phy,
            out_dir=out,
            esl_reformat_bin=esl_reformat_bin,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats["n_taxa"]),
            str(stats["alignment_length"]),
            str(stats["input_phy"]),
            str(stats["fasta_out"]),
            str(stats["dna_fasta_out"]),
            str(stats["stockholm_out"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_taxa",
            "alignment_length",
            "input_phy",
            "fasta_out",
            "dna_fasta_out",
            "stockholm_out",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("HMM-prep step complete.")