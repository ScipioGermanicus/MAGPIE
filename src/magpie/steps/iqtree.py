from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List

_LOGGER = logging.getLogger("magpie.iqtree")


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


def _phylip_has_enough_taxa(path: Path) -> tuple[int, int]:
    """
    Read the first line of a PHYLIP alignment and return:
      (n_taxa, alignment_length)
    """
    if not path.exists():
        return (0, 0)

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        first = fh.readline().strip()

    if not first:
        return (0, 0)

    parts = first.split()
    if len(parts) < 2:
        raise ValueError(f"Could not parse PHYLIP header from: {path}")

    try:
        n_taxa = int(parts[0])
        aln_len = int(parts[1])
    except ValueError as e:
        raise ValueError(f"Invalid PHYLIP header in {path}: {first}") from e

    return (n_taxa, aln_len)


def _run_iqtree(
    *,
    msa: Path,
    prefix: Path,
    iqtree_bin: str,
    threads: int,
    bootstrap: int,
    seed: int,
) -> None:
    cmd = [
        iqtree_bin,
        "-s",
        str(msa),
        "-m",
        "GTR+G",
        "-bb",
        str(bootstrap),
        "-nt",
        str(threads),
        "-seed",
        str(seed),
        "-redo",
        "-pre",
        str(prefix),
    ]
    _run_command(cmd)


def _build_one_domain(
    *,
    domain: str,
    msa_phy: Path,
    out_dir: Path,
    iqtree_bin: str,
    threads: int,
    bootstrap: int,
    seed: int,
) -> Dict[str, object]:
    prefix = out_dir / f"{domain}_16S_iqtree"
    treefile = out_dir / f"{domain}_16S_iqtree.treefile"
    log = out_dir / f"{domain}_16S_iqtree.log"
    iqtree_report = out_dir / f"{domain}_16S_iqtree.iqtree"
    contree = out_dir / f"{domain}_16S_iqtree.contree"

    n_taxa, aln_len = _phylip_has_enough_taxa(msa_phy)

    if n_taxa == 0:
        _LOGGER.warning("Skipping %s IQ-TREE: PHYLIP file missing or empty: %s", domain, msa_phy)
        return {
            "domain": domain,
            "status": "skipped_empty_input",
            "msa_phy": str(msa_phy),
            "prefix": str(prefix),
            "treefile": str(treefile),
            "log": str(log),
            "iqtree_report": str(iqtree_report),
            "contree": str(contree),
            "n_taxa": 0,
            "alignment_length": 0,
        }

    if n_taxa < 4:
        _LOGGER.warning(
            "Skipping %s IQ-TREE: fewer than 4 taxa in PHYLIP (%d taxa).",
            domain,
            n_taxa,
        )
        return {
            "domain": domain,
            "status": "skipped_too_few_taxa",
            "msa_phy": str(msa_phy),
            "prefix": str(prefix),
            "treefile": str(treefile),
            "log": str(log),
            "iqtree_report": str(iqtree_report),
            "contree": str(contree),
            "n_taxa": n_taxa,
            "alignment_length": aln_len,
        }

    _LOGGER.info("Running IQ-TREE for %s", domain)
    _run_iqtree(
        msa=msa_phy,
        prefix=prefix,
        iqtree_bin=iqtree_bin,
        threads=threads,
        bootstrap=bootstrap,
        seed=seed,
    )

    return {
        "domain": domain,
        "status": "processed",
        "msa_phy": str(msa_phy),
        "prefix": str(prefix),
        "treefile": str(treefile),
        "log": str(log),
        "iqtree_report": str(iqtree_report),
        "contree": str(contree),
        "n_taxa": n_taxa,
        "alignment_length": aln_len,
    }


def iqtree_step(
    *,
    raxml_check_dir: Path,
    out: Path,
    iqtree_bin: str,
    threads: int,
    bootstrap: int,
    seed: int,
    force: bool,
) -> None:
    """
    Build archaeal and bacterial 16S trees from RAxML-check reduced PHYLIP alignments.

    Inputs:
      raxml_check_dir/archaea_raxml-check.raxml.reduced.phy
      raxml_check_dir/bacteria_raxml-check.raxml.reduced.phy
    """
    arch_phy = raxml_check_dir / "archaea_raxml-check.raxml.reduced.phy"
    bact_phy = raxml_check_dir / "bacteria_raxml-check.raxml.reduced.phy"

    missing = [str(p) for p in [arch_phy, bact_phy] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required RAxML-check PHYLIP inputs for IQ-TREE:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"IQ-TREE outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "raxml_check_dir": str(raxml_check_dir),
            "iqtree_bin": iqtree_bin,
            "threads": threads,
            "bootstrap": bootstrap,
            "seed": seed,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, msa_phy in (("archaea", arch_phy), ("bacteria", bact_phy)):
        stats = _build_one_domain(
            domain=domain,
            msa_phy=msa_phy,
            out_dir=out,
            iqtree_bin=iqtree_bin,
            threads=threads,
            bootstrap=bootstrap,
            seed=seed,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats["n_taxa"]),
            str(stats["alignment_length"]),
            str(stats["msa_phy"]),
            str(stats["treefile"]),
            str(stats["contree"]),
            str(stats["log"]),
            str(stats["iqtree_report"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_taxa",
            "alignment_length",
            "msa_phy",
            "treefile",
            "contree",
            "log",
            "iqtree_report",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("IQ-TREE step complete.")