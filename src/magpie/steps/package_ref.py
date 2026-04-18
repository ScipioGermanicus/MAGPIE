from __future__ import annotations

import csv
import json
import logging
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

_LOGGER = logging.getLogger("magpie.package_ref")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _read_fasta_ids(path: Path) -> List[str]:
    ids: List[str] = []
    if not path.exists():
        return ids

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if raw.startswith(">"):
                ids.append(raw[1:].strip().split()[0])
    return ids


def _read_copy_table(path: Path) -> Dict[str, int]:
    out: Dict[str, int] = {}
    if not path.exists():
        return out

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            gid = parts[0].strip()
            try:
                count = int(float(parts[1].strip()))
            except ValueError:
                continue
            out[gid] = count
    return out


def _write_filtered_copy_table(
    *,
    included_genomes: List[str],
    copy_table_in: Path,
    out_path: Path,
) -> Tuple[int, int]:
    """
    Write filtered copy table with columns:
      assembly    16S_rRNA_Count
    Counts >10 are capped to 10.
    Returns:
      (n_written, n_capped)
    """
    copies = _read_copy_table(copy_table_in)

    n_written = 0
    n_capped = 0

    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["assembly", "16S_rRNA_Count"])

        for gid in included_genomes:
            if gid not in copies:
                continue
            count = copies[gid]
            if count > 10:
                count = 10
                n_capped += 1
            w.writerow([gid, str(count)])
            n_written += 1

    return (n_written, n_capped)


def _copy_required(src: Path, dst: Path) -> None:
    if not src.exists():
        raise FileNotFoundError(f"Required input file missing: {src}")
    shutil.copy2(src, dst)


def _parse_raxml_ng_log(log_path: Path) -> Dict[str, str]:
    """
    Parse a RAxML-NG --evaluate log into fields needed for a RAxML 7.x-style raxml_info.
    """
    patterns = None
    base_freqs = None
    subs_rates = None
    final_ll = None
    elapsed = None
    called_line = None
    cmd_line = None

    with log_path.open("r", encoding="utf-8", errors="replace") as fh:
        lines = fh.readlines()

    for i, line in enumerate(lines):
        if "RAxML-NG was called at" in line:
            called_line = line.rstrip("\n")
            for j in range(i + 1, len(lines)):
                cmd = lines[j].strip()
                if cmd:
                    cmd_line = cmd
                    break
            break

    for raw in lines:
        line = raw.strip()

        m = re.search(r"Alignment comprises\s+(\d+)\s+partitions\s+and\s+(\d+)\s+patterns", line)
        if m:
            patterns = m.group(2)
            continue

        if line.startswith("Base frequencies (ML):"):
            base_freqs = line.split(":", 1)[1].strip().split()
            continue

        if line.startswith("Substitution rates (ML):"):
            subs_rates = line.split(":", 1)[1].strip().split()
            continue

        m2 = re.search(r"Final LogLikelihood:\s+(-?[0-9.eE+-]+)", line)
        if m2:
            final_ll = m2.group(1)
            continue

        m3 = re.search(r"Elapsed time:\s+([0-9.eE+-]+)\s+seconds", line)
        if m3:
            elapsed = m3.group(1)
            continue

    if patterns is None:
        raise RuntimeError(f"Could not parse alignment patterns from {log_path}")
    if base_freqs is None:
        raise RuntimeError(f"Could not parse base frequencies from {log_path}")
    if subs_rates is None:
        raise RuntimeError(f"Could not parse substitution rates from {log_path}")
    if final_ll is None:
        raise RuntimeError(f"Could not parse final log-likelihood from {log_path}")
    if elapsed is None:
        elapsed = "0.0"
    if called_line is None or cmd_line is None:
        raise RuntimeError(f"Could not find 'RAxML-NG was called at' section in {log_path}")

    return {
        "patterns": patterns,
        "base_freqs": " ".join(base_freqs),
        "subs_rates": " ".join(subs_rates),
        "final_ll": final_ll,
        "elapsed": elapsed,
        "called_line": called_line,
        "cmd_line": cmd_line,
    }


def _write_raxml_info(log_path: Path, out_path: Path) -> None:
    parsed = _parse_raxml_ng_log(log_path)

    txt: List[str] = []
    txt.append("This is RAxML version 7.7.2 released by Alexandros Stamatakis on July 31 2013.")
    txt.append("")
    txt.append("This is a RAxML_info file from an --evaluate run, manually reformatted")
    txt.append("")
    txt.append("Partition: 0")
    txt.append(f"Alignment Patterns: {parsed['patterns']}")
    txt.append("Name: No Name Provided")
    txt.append("DataType: DNA")
    txt.append("Substitution Matrix: GTR")
    txt.append("")
    txt.append(parsed["called_line"])
    txt.append("")
    txt.append(parsed["cmd_line"])
    txt.append("")
    txt.append(f"Base frequencies: {parsed['base_freqs']}")
    txt.append("")
    txt.append(f"Inference[0]: Time {parsed['elapsed']} CAT-based likelihood -0000, best rearrangement setting 5")
    txt.append(f"alpha[0]: 1.000000 rates[0] ac ag at cg ct gt: {parsed['subs_rates']}")
    txt.append("")
    txt.append("")
    txt.append("NOT conducting any final model optimizations on all 1 trees under CAT-based")
    txt.append("model ....")
    txt.append("")
    txt.append(f"Final GAMMA  likelihood: {parsed['final_ll']}")
    txt.append("")

    out_path.write_text("\n".join(txt), encoding="utf-8")


def _package_one_domain(
    *,
    domain: str,
    hmm_prep_dir: Path,
    hmm_build_dir: Path,
    iqtree_dir: Path,
    raxml_evaluate_dir: Path,
    rrna_dir: Path,
    out_root: Path,
) -> Dict[str, object]:
    ref_name = "bac_ref" if domain == "bacteria" else "arc_ref"
    ref_dir = out_root / ref_name
    _ensure_dir(ref_dir)

    prefix = "bac_ref" if domain == "bacteria" else "arc_ref"

    fna_src = hmm_prep_dir / f"{domain}_raxml-check.raxml.reduced_dna.fna"
    hmm_src = hmm_build_dir / f"{domain}_raxml-check.raxml.reduced_dna.hmm"
    tree_src = iqtree_dir / f"{domain}_16S_iqtree.treefile"
    model_src = raxml_evaluate_dir / f"{domain}_raxml.raxml.bestModel"
    log_src = raxml_evaluate_dir / f"{domain}_raxml.raxml.log"
    copies_src = rrna_dir / f"{domain}_16S_copies.txt"

    fna_dst = ref_dir / f"{prefix}.fna"
    hmm_dst = ref_dir / f"{prefix}.hmm"
    tree_dst = ref_dir / f"{prefix}.tre"
    model_dst = ref_dir / f"{prefix}.model"
    info_dst = ref_dir / f"{prefix}.raxml_info"

    if not fna_src.exists():
        _LOGGER.warning("Skipping %s package-ref: missing FASTA input %s", domain, fna_src)
        return {
            "domain": domain,
            "status": "skipped_missing_input",
            "ref_dir": str(ref_dir),
        }

    _copy_required(fna_src, fna_dst)
    _copy_required(hmm_src, hmm_dst)
    _copy_required(tree_src, tree_dst)
    _copy_required(model_src, model_dst)
    _write_raxml_info(log_src, info_dst)

    fasta_ids = _read_fasta_ids(fna_dst)
    copy_table_dst = out_root / f"{domain}_16S_copies.txt"
    n_written, n_capped = _write_filtered_copy_table(
        included_genomes=fasta_ids,
        copy_table_in=copies_src,
        out_path=copy_table_dst,
    )

    return {
        "domain": domain,
        "status": "processed",
        "ref_dir": str(ref_dir),
        "fna": str(fna_dst),
        "hmm": str(hmm_dst),
        "tree": str(tree_dst),
        "model": str(model_dst),
        "raxml_info": str(info_dst),
        "copies_table": str(copy_table_dst),
        "n_fasta_ids": len(fasta_ids),
        "n_copy_rows_written": n_written,
        "n_copy_rows_capped": n_capped,
    }


def package_ref_step(
    *,
    rrna_dir: Path,
    iqtree_dir: Path,
    raxml_evaluate_dir: Path,
    hmm_prep_dir: Path,
    hmm_build_dir: Path,
    out: Path,
    force: bool,
) -> None:
    """
    Package PICRUSt2-style reference folders and auxiliary copy-count tables.

    Outputs:
      out/bac_ref/*
      out/arc_ref/*
      out/bacteria_16S_copies.txt
      out/archaea_16S_copies.txt
    """
    key_outputs = [out / "summary.tsv", out / "report.json"]
    if any(p.exists() for p in key_outputs) and not force:
        raise FileExistsError(f"package-ref outputs already exist under {out}. Use --force to overwrite.")

    _ensure_dir(out)

    report = {
        "inputs": {
            "rrna_dir": str(rrna_dir),
            "iqtree_dir": str(iqtree_dir),
            "raxml_evaluate_dir": str(raxml_evaluate_dir),
            "hmm_prep_dir": str(hmm_prep_dir),
            "hmm_build_dir": str(hmm_build_dir),
        },
        "domains": {},
        "outputs": {
            "summary_tsv": str(out / "summary.tsv"),
            "report_json": str(out / "report.json"),
            "bac_ref": str(out / "bac_ref"),
            "arc_ref": str(out / "arc_ref"),
        },
    }

    summary_rows: List[List[str]] = []

    for domain in ("bacteria", "archaea"):
        stats = _package_one_domain(
            domain=domain,
            hmm_prep_dir=hmm_prep_dir,
            hmm_build_dir=hmm_build_dir,
            iqtree_dir=iqtree_dir,
            raxml_evaluate_dir=raxml_evaluate_dir,
            rrna_dir=rrna_dir,
            out_root=out,
        )

        report["domains"][domain] = stats
        summary_rows.append([
            domain,
            str(stats["status"]),
            str(stats.get("ref_dir", "")),
            str(stats.get("fna", "")),
            str(stats.get("hmm", "")),
            str(stats.get("tree", "")),
            str(stats.get("model", "")),
            str(stats.get("raxml_info", "")),
            str(stats.get("copies_table", "")),
            str(stats.get("n_fasta_ids", 0)),
            str(stats.get("n_copy_rows_written", 0)),
            str(stats.get("n_copy_rows_capped", 0)),
        ])

    with (out / "summary.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "domain",
            "status",
            "ref_dir",
            "fna",
            "hmm",
            "tree",
            "model",
            "raxml_info",
            "copies_table",
            "n_fasta_ids",
            "n_copy_rows_written",
            "n_copy_rows_capped",
        ])
        w.writerows(summary_rows)

    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    _LOGGER.info("package-ref step complete.")