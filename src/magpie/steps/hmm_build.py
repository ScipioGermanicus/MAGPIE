from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List

_LOGGER = logging.getLogger("magpie.hmm_build")


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


def _stockholm_has_sequences(path: Path) -> bool:
    if not path.exists():
        return False

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line == "//":
                continue
            return True

    return False


def _build_one_domain(
    *,
    domain: str,
    stockholm_in: Path,
    out_dir: Path,
    hmmbuild_bin: str,
    cpus: int,
) -> Dict[str, object]:
    hmm_out = out_dir / f"{domain}_raxml-check.raxml.reduced_dna.hmm"

    if not stockholm_in.exists():
        _LOGGER.warning("Skipping %s hmm-build: Stockholm input missing: %s", domain, stockholm_in)
        return {
            "domain": domain,
            "status": "skipped_missing_input",
            "stockholm_in": str(stockholm_in),
            "hmm_out": str(hmm_out),
        }

    if not _stockholm_has_sequences(stockholm_in):
        _LOGGER.warning("Skipping %s hmm-build: Stockholm input empty: %s", domain, stockholm_in)
        return {
            "domain": domain,
            "status": "skipped_empty_input",
            "stockholm_in": str(stockholm_in),
            "hmm_out": str(hmm_out),
        }

    _LOGGER.info("Running hmmbuild for %s", domain)
    _run_command([
        hmmbuild_bin,
        "--cpu",
        str(cpus),
        str(hmm_out),
        str(stockholm_in),
    ])

    return {
        "domain": domain,
        "status": "processed",
        "stockholm_in": str(stockholm_in),
        "hmm_out": str(hmm_out),
    }


def hmm_build_step(
    *,
    hmm_prep_dir: Path,
    out: Path,
    hmmbuild_bin: str,
    cpus: int,
    force: bool,
) -> None:
    """
    Build domain-specific HMMs from Stockholm alignments produced by hmm-prep.

    Inputs:
      hmm_prep_dir/archaea_raxml-check.raxml.reduced_dna.sto
      hmm_prep_dir/bacteria_raxml-check.raxml.reduced_dna.sto
    """
    arch_sto = hmm_prep_dir / "archaea_raxml-check.raxml.reduced_dna.sto"
    bact_sto = hmm_prep_dir / "bacteria_raxml-check.raxml.reduced_dna.sto"

    missing = [str(p) for p in [arch_sto, bact_sto] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required HMM-prep Stockholm inputs for hmm-build:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"HMM-build outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "hmm_prep_dir": str(hmm_prep_dir),
            "hmmbuild_bin": hmmbuild_bin,
            "cpus": cpus,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, sto in (("archaea", arch_sto), ("bacteria", bact_sto)):
        stats = _build_one_domain(
            domain=domain,
            stockholm_in=sto,
            out_dir=out,
            hmmbuild_bin=hmmbuild_bin,
            cpus=cpus,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats["stockholm_in"]),
            str(stats["hmm_out"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "stockholm_in",
            "hmm_out",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("HMM-build step complete.")