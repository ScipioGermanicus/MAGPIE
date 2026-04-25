from __future__ import annotations

import csv
import json
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

_LOGGER = logging.getLogger("magpie.raxml_evaluate")


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


def _phylip_has_enough_taxa(path: Path) -> Tuple[int, int]:
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


def _run_raxml_evaluate(
    *,
    msa: Path,
    tree: Path,
    prefix: Path,
    raxml_ng_bin: str,
    threads: int,
    force: bool,
) -> None:
    cmd = [
        raxml_ng_bin,
        "--evaluate",
        "--msa",
        str(msa),
        "--tree",
        str(tree),
        "--model",
        "GTR+G",
        "--prefix",
        str(prefix),
        "--threads",
        str(threads),
    ]

    if force:
        cmd.append("--redo")

    _run_command(cmd)


def _evaluate_one_domain(
    *,
    domain: str,
    msa_phy: Path,
    start_tree: Path,
    out_dir: Path,
    raxml_ng_bin: str,
    threads: int,
    force: bool,
) -> Dict[str, object]:
    prefix = out_dir / f"{domain}_raxml"
    best_tree = out_dir / f"{domain}_raxml.raxml.bestTree"
    best_model = out_dir / f"{domain}_raxml.raxml.bestModel"
    best_log = out_dir / f"{domain}_raxml.raxml.log"

    n_taxa, aln_len = _phylip_has_enough_taxa(msa_phy)

    if n_taxa == 0:
        _LOGGER.warning("Skipping %s RAxML evaluate: PHYLIP file missing or empty: %s", domain, msa_phy)
        return {
            "domain": domain,
            "status": "skipped_empty_input",
            "msa_phy": str(msa_phy),
            "start_tree": str(start_tree),
            "prefix": str(prefix),
            "best_tree": str(best_tree),
            "best_model": str(best_model),
            "log": str(best_log),
            "n_taxa": 0,
            "alignment_length": 0,
        }

    if n_taxa < 4:
        _LOGGER.warning(
            "Skipping %s RAxML evaluate: fewer than 4 taxa in PHYLIP (%d taxa).",
            domain,
            n_taxa,
        )
        return {
            "domain": domain,
            "status": "skipped_too_few_taxa",
            "msa_phy": str(msa_phy),
            "start_tree": str(start_tree),
            "prefix": str(prefix),
            "best_tree": str(best_tree),
            "best_model": str(best_model),
            "log": str(best_log),
            "n_taxa": n_taxa,
            "alignment_length": aln_len,
        }

    if not start_tree.exists():
        raise FileNotFoundError(f"Missing start tree for {domain}: {start_tree}")

    _LOGGER.info("Running RAxML-NG --evaluate for %s", domain)
    _run_raxml_evaluate(
        msa=msa_phy,
        tree=start_tree,
        prefix=prefix,
        raxml_ng_bin=raxml_ng_bin,
        threads=threads,
        force=force,
    )

    return {
        "domain": domain,
        "status": "processed",
        "msa_phy": str(msa_phy),
        "start_tree": str(start_tree),
        "prefix": str(prefix),
        "best_tree": str(best_tree),
        "best_model": str(best_model),
        "log": str(best_log),
        "n_taxa": n_taxa,
        "alignment_length": aln_len,
    }


def raxml_evaluate_step(
    *,
    raxml_check_dir: Path,
    iqtree_dir: Path,
    out: Path,
    raxml_ng_bin: str,
    threads: int,
    force: bool,
) -> None:
    """
    Evaluate domain-specific trees with RAxML-NG using reduced PHYLIP alignments
    and IQ-TREE starting trees.

    Inputs:
      raxml_check_dir/archaea_raxml-check.raxml.reduced.phy
      raxml_check_dir/bacteria_raxml-check.raxml.reduced.phy
      iqtree_dir/archaea_16S_iqtree.treefile
      iqtree_dir/bacteria_16S_iqtree.treefile
    """
    arch_phy = raxml_check_dir / "archaea_raxml-check.raxml.reduced.phy"
    bact_phy = raxml_check_dir / "bacteria_raxml-check.raxml.reduced.phy"

    arch_tree = iqtree_dir / "archaea_16S_iqtree.treefile"
    bact_tree = iqtree_dir / "bacteria_16S_iqtree.treefile"

    missing = [str(p) for p in [arch_phy, bact_phy] if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required RAxML-check PHYLIP inputs for raxml-evaluate:\n"
            + "\n".join(f"  - {x}" for x in missing)
        )

    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"RAxML-evaluate outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "raxml_check_dir": str(raxml_check_dir),
            "iqtree_dir": str(iqtree_dir),
            "raxml_ng_bin": raxml_ng_bin,
            "threads": threads,
            "force": force,
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain, msa_phy, start_tree in (
        ("archaea", arch_phy, arch_tree),
        ("bacteria", bact_phy, bact_tree),
    ):
        stats = _evaluate_one_domain(
            domain=domain,
            msa_phy=msa_phy,
            start_tree=start_tree,
            out_dir=out,
            raxml_ng_bin=raxml_ng_bin,
            threads=threads,
            force=force,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats["n_taxa"]),
            str(stats["alignment_length"]),
            str(stats["msa_phy"]),
            str(stats["start_tree"]),
            str(stats["best_tree"]),
            str(stats["best_model"]),
            str(stats["log"]),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "n_taxa",
            "alignment_length",
            "msa_phy",
            "start_tree",
            "best_tree",
            "best_model",
            "log",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("RAxML-evaluate step complete.")