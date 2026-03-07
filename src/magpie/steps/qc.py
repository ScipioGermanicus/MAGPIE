from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import gzip
import json
import re
import shutil
from typing import Dict, Iterable, List, Optional, Tuple

from magpie.util.shell import run as shell_run


SEP_RE = re.compile(r"^\s*-{5,}\s*$")
SPLIT_RE = re.compile(r"\s{2,}")  # CheckM QA tables align with >=2 spaces
FASTA_EXTS = (".fa", ".fna", ".fasta")
GZ_EXT = ".gz"


@dataclass(frozen=True)
class CheckMRow:
    genome_id: str
    marker_lineage: str
    completeness: float
    contamination: float


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
    for p in dir_.iterdir():
        if not p.is_file():
            continue
        low = p.name.lower()
        if any(low.endswith(ext) for ext in FASTA_EXTS) or any(low.endswith(ext + GZ_EXT) for ext in FASTA_EXTS):
            yield p


def _read_domain_map(domain_map_tsv: Path) -> Dict[str, str]:
    """
    domain_map.tsv written by taxonomy step.
    Expected columns: genome_id \t domain
    """
    domain_by_id: Dict[str, str] = {}
    with domain_map_tsv.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        genome_key = "genome_id" if "genome_id" in fields else ("Genome" if "Genome" in fields else None)
        domain_key = "domain" if "domain" in fields else ("Domain" if "Domain" in fields else None)

        if genome_key is None or domain_key is None:
            raise ValueError(f"domain_map.tsv missing expected columns. Found: {fields}")

        for row in reader:
            gid = (row.get(genome_key) or "").strip()
            dom = (row.get(domain_key) or "").strip()
            if gid:
                domain_by_id[gid] = dom
    return domain_by_id


def parse_checkm_qa_table(tsv_like: Path) -> List[CheckMRow]:
    """
    Parse CheckM QA (-o 2) wide table (space-aligned).
    Works with the style:
    Bin Id, Marker lineage, ..., Completeness, Contamination, ...
    """
    header_cols: Optional[List[str]] = None
    idx_comp = idx_cont = None
    rows: List[CheckMRow] = []

    with tsv_like.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            s = line.strip()
            if not s:
                continue
            if SEP_RE.match(s):
                continue

            if s.startswith("Bin Id"):
                header_cols = SPLIT_RE.split(s)
                try:
                    idx_comp = header_cols.index("Completeness")
                    idx_cont = header_cols.index("Contamination")
                except ValueError as e:
                    raise ValueError(f"Could not locate Completeness/Contamination in header: {header_cols}") from e
                continue

            if header_cols is None or idx_comp is None or idx_cont is None:
                continue

            parts = SPLIT_RE.split(s)
            if len(parts) <= max(idx_comp, idx_cont):
                continue

            genome_id = parts[0].strip()
            marker_lineage = parts[1].strip() if len(parts) > 1 else ""

            try:
                completeness = float(parts[idx_comp])
                contamination = float(parts[idx_cont])
            except ValueError:
                continue

            rows.append(
                CheckMRow(
                    genome_id=genome_id,
                    marker_lineage=marker_lineage,
                    completeness=completeness,
                    contamination=contamination,
                )
            )

    seen = set()
    dedup: List[CheckMRow] = []
    for r in rows:
        if r.genome_id in seen:
            continue
        seen.add(r.genome_id)
        dedup.append(r)

    if not dedup:
        raise ValueError(f"No rows parsed from CheckM QA table: {tsv_like}")

    return dedup


def parse_legacy_minimal_tsv(path: Path) -> List[CheckMRow]:
    """
    Legacy minimal TSV with columns:
      genome_id, checkm_completeness, checkm_contamination
    """
    out: List[CheckMRow] = []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []

        gid_key = "genome_id" if "genome_id" in fields else None
        c_key = "checkm_completeness" if "checkm_completeness" in fields else None
        t_key = "checkm_contamination" if "checkm_contamination" in fields else None

        if not (gid_key and c_key and t_key):
            raise ValueError(f"Legacy CheckM results TSV missing required columns. Found: {fields}")

        for row in reader:
            gid = (row.get(gid_key) or "").strip()
            if not gid:
                continue
            try:
                cpl = float((row.get(c_key) or "").strip())
                cnt = float((row.get(t_key) or "").strip())
            except ValueError:
                continue

            out.append(
                CheckMRow(
                    genome_id=gid,
                    marker_lineage="",
                    completeness=cpl,
                    contamination=cnt,
                )
            )

    seen = set()
    dedup: List[CheckMRow] = []
    for r in out:
        if r.genome_id in seen:
            continue
        seen.add(r.genome_id)
        dedup.append(r)

    if not dedup:
        raise ValueError(f"No rows parsed from legacy CheckM TSV: {path}")

    return dedup


def write_tsv(path: Path, header: List[str], rows: Iterable[List[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def place_bins(
    kept_ids: List[str],
    src_dir: Path,
    dst_dir: Path,
    mode: str = "symlink",
) -> Tuple[int, List[str]]:
    """
    Place genomes by exact filename-stem match inside src_dir.
    """
    dst_dir.mkdir(parents=True, exist_ok=True)

    index: Dict[str, Path] = {}
    for p in _iter_genome_files(src_dir):
        gid = _strip_fasta_suffix(p.name)
        if gid in index:
            raise ValueError(f"Ambiguous genome_id '{gid}' in {src_dir}: {index[gid].name} and {p.name}")
        index[gid] = p

    placed = 0
    missing: List[str] = []

    for gid in kept_ids:
        src = index.get(gid)
        if src is None:
            missing.append(gid)
            continue

        dst = dst_dir / src.name
        if dst.exists():
            continue

        try:
            if mode == "symlink":
                dst.symlink_to(src)
            elif mode == "hardlink":
                dst.hardlink_to(src)
            elif mode == "copy":
                shutil.copy2(src, dst)
            else:
                raise ValueError(f"Unknown placement mode: {mode}")
        except OSError:
            if mode in ("symlink", "hardlink"):
                shutil.copy2(src, dst)
            else:
                raise

        placed += 1

    return placed, missing


def _gunzip_or_copy_dir(src: Path, dst: Path) -> None:
    _ensure_dir(dst)
    any_seen = False

    for p in _iter_genome_files(src):
        any_seen = True
        if p.name.lower().endswith(GZ_EXT):
            out_fp = dst / p.name[:-len(GZ_EXT)]
            with gzip.open(p, "rb") as fin, out_fp.open("wb") as fout:
                shutil.copyfileobj(fin, fout)
        else:
            shutil.copy2(p, dst / p.name)

    if not any_seen:
        raise ValueError(f"No FASTA files found in: {src}")


def _run_checkm(
    *,
    bacteria_dir: Path,
    archaea_dir: Path,
    out: Path,
    cpus: int,
    checkm_bin: str,
) -> Tuple[Optional[Path], Optional[Path]]:
    """
    Run CheckM on taxonomy-split genomes and return QA file paths.
    """
    checkm_root = out / "checkm_run"
    work_bac = checkm_root / "unzipped" / "bacteria"
    work_arc = checkm_root / "unzipped" / "archaea"
    _ensure_dir(work_bac)
    _ensure_dir(work_arc)

    qa_bac: Optional[Path] = None
    qa_arc: Optional[Path] = None

    bac_files = list(_iter_genome_files(bacteria_dir)) if bacteria_dir.exists() else []
    arc_files = list(_iter_genome_files(archaea_dir)) if archaea_dir.exists() else []

    if bac_files:
        _gunzip_or_copy_dir(bacteria_dir, work_bac)
        out_bac = checkm_root / "bacteria"
        _ensure_dir(out_bac)
        ext = "fna" if list(work_bac.glob("*.fna")) else "fa"
        shell_run([checkm_bin, "lineage_wf", "-x", ext, "-t", str(cpus), str(work_bac), str(out_bac)])
        qa_bac = out_bac / "checkm_qa.tsv"
        shell_run([checkm_bin, "qa", str(out_bac / "lineage.ms"), str(out_bac), "-o", "2", "-f", str(qa_bac)])

    if arc_files:
        _gunzip_or_copy_dir(archaea_dir, work_arc)
        out_arc = checkm_root / "archaea"
        _ensure_dir(out_arc)
        ext = "fna" if list(work_arc.glob("*.fna")) else "fa"
        shell_run([checkm_bin, "lineage_wf", "-x", ext, "-t", str(cpus), str(work_arc), str(out_arc)])
        qa_arc = out_arc / "checkm_qa.tsv"
        shell_run([checkm_bin, "qa", str(out_arc / "lineage.ms"), str(out_arc), "-o", "2", "-f", str(qa_arc)])

    if qa_bac is None and qa_arc is None:
        raise ValueError("No genomes found in bacteria/archaea dirs; cannot run CheckM.")

    return qa_bac, qa_arc

def _resolve_split_dirs(tax_dir: Path) -> Tuple[Path, Path]:
    """
    Resolve taxonomy split genome directories using a small set of explicit layouts.
    No fuzzy searching.
    """
    candidates = [
        (tax_dir / "bacteria", tax_dir / "archaea"),
        (tax_dir / "split" / "bacteria", tax_dir / "split" / "archaea"),
        (tax_dir / "genomes" / "bacteria", tax_dir / "genomes" / "archaea"),
        (tax_dir / "genomes_to_search_barrnap" / "bacteria", tax_dir / "genomes_to_search_barrnap" / "archaea"),
    ]

    for bac, arc in candidates:
        if bac.exists() and arc.exists():
            return bac, arc

    checked = "\n".join(
        f"  - {bac}\n  - {arc}"
        for bac, arc in candidates
    )
    raise FileNotFoundError(
        "Could not find taxonomy split genome directories. Checked:\n"
        f"{checked}"
    )

def qc_step(
    *,
    tax_dir: Path,
    out: Path,
    cpus: int,
    force: bool,
    checkm_qa: Path | None,
    checkm_qa_bacteria: Path | None,
    checkm_qa_archaea: Path | None,
    checkm_results: Path | None,
    completeness_min: float,
    contamination_max: float,
    place_mode: str,
    checkm_bin: str,
) -> None:
    """
    QC step used by MAGPIE CLI.
    Expected taxonomy outputs under tax_dir:
      - domain_map.tsv
      - bacteria/archaea OR split/bacteria/archaea OR similar explicit layouts
    """
    out.mkdir(parents=True, exist_ok=True)

    out_all = out / "checkm_clean_all.tsv"
    out_filt = out / "checkm_filtered.tsv"
    if (out_all.exists() or out_filt.exists()) and not force:
        raise FileExistsError(f"QC outputs exist under {out}. Use --force to overwrite.")

    domain_map = tax_dir / "domain_map.tsv"
    bacteria_dir, archaea_dir = _resolve_split_dirs(tax_dir)

    domain_by_id = _read_domain_map(domain_map)

    merged_rows: List[CheckMRow] = []

    if checkm_results is not None:
        merged_rows.extend(parse_legacy_minimal_tsv(checkm_results))
    elif checkm_qa is not None:
        merged_rows.extend(parse_checkm_qa_table(checkm_qa))
    elif checkm_qa_bacteria is not None or checkm_qa_archaea is not None:
        if checkm_qa_bacteria is not None:
            merged_rows.extend(parse_checkm_qa_table(checkm_qa_bacteria))
        if checkm_qa_archaea is not None:
            merged_rows.extend(parse_checkm_qa_table(checkm_qa_archaea))
    else:
        qa_bac, qa_arc = _run_checkm(
            bacteria_dir=bacteria_dir,
            archaea_dir=archaea_dir,
            out=out,
            cpus=cpus,
            checkm_bin=checkm_bin,
        )
        if qa_bac is not None:
            merged_rows.extend(parse_checkm_qa_table(qa_bac))
        if qa_arc is not None:
            merged_rows.extend(parse_checkm_qa_table(qa_arc))

    seen = set()
    dedup_rows: List[CheckMRow] = []
    for r in merged_rows:
        if r.genome_id in seen:
            continue
        seen.add(r.genome_id)
        dedup_rows.append(r)

    if not dedup_rows:
        raise ValueError("No CheckM rows available after parsing inputs.")

    merged = []
    missing_in_domain = []
    for r in dedup_rows:
        dom = domain_by_id.get(r.genome_id)
        if dom is None:
            missing_in_domain.append(r.genome_id)
            continue
        merged.append((r, dom))

    if missing_in_domain:
        missing_preview = ", ".join(missing_in_domain[:10])
        raise ValueError(
            f"{len(missing_in_domain)} CheckM IDs not present in domain_map.tsv "
            f"(first 10: {missing_preview}). This usually means an ID mismatch between "
            f"prepared genomes and CheckM input."
        )

    write_tsv(
        out / "checkm_results.min.tsv",
        header=["genome_id", "checkm_completeness", "checkm_contamination"],
        rows=[
            [r.genome_id, f"{r.completeness:.2f}", f"{r.contamination:.2f}"]
            for (r, _) in merged
        ],
    )

    write_tsv(
        out_all,
        header=["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"],
        rows=[
            [r.genome_id, r.marker_lineage, f"{r.completeness:.2f}", f"{r.contamination:.2f}", dom]
            for (r, dom) in merged
        ],
    )

    kept = [
        (r, dom) for (r, dom) in merged
        if r.completeness >= completeness_min and r.contamination <= contamination_max
    ]

    write_tsv(
        out_filt,
        header=["genome_id", "marker_lineage", "checkm_completeness", "checkm_contamination", "domain"],
        rows=[
            [r.genome_id, r.marker_lineage, f"{r.completeness:.2f}", f"{r.contamination:.2f}", dom]
            for (r, dom) in kept
        ],
    )

    bac_ids = sorted({r.genome_id for (r, dom) in kept if dom.lower().startswith("bact")})
    arc_ids = sorted({r.genome_id for (r, dom) in kept if dom.lower().startswith("arch")})

    (out / "bacteria.txt").write_text("\n".join(bac_ids) + ("\n" if bac_ids else ""), encoding="utf-8")
    (out / "archaea.txt").write_text("\n".join(arc_ids) + ("\n" if arc_ids else ""), encoding="utf-8")

    bins_root = out / "bins"
    placed_bac, miss_bac = place_bins(bac_ids, bacteria_dir, bins_root / "bacteria", mode=place_mode)
    placed_arc, miss_arc = place_bins(arc_ids, archaea_dir, bins_root / "archaea", mode=place_mode)

    report = {
        "inputs": {
            "tax_dir": str(tax_dir),
            "domain_map": str(domain_map),
            "bacteria_dir": str(bacteria_dir),
            "archaea_dir": str(archaea_dir),
            "checkm_qa": str(checkm_qa) if checkm_qa is not None else None,
            "checkm_qa_bacteria": str(checkm_qa_bacteria) if checkm_qa_bacteria is not None else None,
            "checkm_qa_archaea": str(checkm_qa_archaea) if checkm_qa_archaea is not None else None,
            "checkm_results": str(checkm_results) if checkm_results is not None else None,
        },
        "thresholds": {
            "completeness_min": completeness_min,
            "contamination_max": contamination_max,
        },
        "kept": {
            "bacteria": len(bac_ids),
            "archaea": len(arc_ids),
            "total": len(kept),
        },
        "placed": {
            "bacteria": placed_bac,
            "archaea": placed_arc,
        },
        "missing_files": {
            "bacteria": miss_bac[:20],
            "archaea": miss_arc[:20],
        },
        "placement_mode": place_mode,
    }
    (out / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")