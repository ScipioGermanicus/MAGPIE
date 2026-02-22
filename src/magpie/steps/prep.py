from __future__ import annotations

from pathlib import Path
import gzip
import json
import shutil
from typing import Iterable

FASTA_SUFFIXES = (".fa", ".fna", ".fasta")
FASTA_GZ_SUFFIXES = (".fa.gz", ".fna.gz", ".fasta.gz")


def _is_fasta(p: Path) -> bool:
    s = p.name.lower()
    return s.endswith(FASTA_SUFFIXES) or s.endswith(FASTA_GZ_SUFFIXES)


def _iter_mag_fastas(mags_dir: Path) -> list[Path]:
    files = [p for p in mags_dir.iterdir() if p.is_file() and _is_fasta(p)]
    return sorted(files, key=lambda x: x.name.lower())


def _decompress_gz_to_bytes(src_gz: Path) -> bytes:
    with gzip.open(src_gz, "rb") as f:
        return f.read()


def _read_bytes_from_plain(src: Path) -> bytes:
    return src.read_bytes()


def _write_gz_bytes(dst_gz: Path, data: bytes) -> None:
    dst_gz.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dst_gz, "wb") as f:
        f.write(data)


def _as_fna_bytes(src: Path) -> bytes:
    """
    Read a FASTA file (.fa/.fna/.fasta or gzipped) and return raw bytes.

    v0.1: no header rewriting; just bytes-in/bytes-out.
    """
    name = src.name.lower()
    if name.endswith(".gz"):
        return _decompress_gz_to_bytes(src)
    return _read_bytes_from_plain(src)


def prep_step(
    *,
    mags: Path,
    out: Path,
    rename_map: Path | None,
    force: bool,
    # New knobs for your formatting step:
    sequential_ids: bool = False,
    out_prefix: str = "MAG",
    pad: int = 4,
) -> None:
    """
    Prepare MAGs for downstream steps.

    v0.1 behaviour:
      - copies MAG FASTAs into out/mags/
      - optionally assigns sequential IDs (MAG0001...) and writes gzipped files
      - writes a mapping file (id_map.tsv) and a prep_report.json
    """
    out.mkdir(parents=True, exist_ok=True)

    out_mags = out / "mags"
    out_mags.mkdir(parents=True, exist_ok=True)

    report_fp = out / "prep_report.json"
    map_fp = out / "id_map.tsv"

    if (report_fp.exists() or map_fp.exists()) and not force:
        raise FileExistsError(f"{out} already contains prep outputs. Use --force to overwrite.")

    files = _iter_mag_fastas(mags)
    if not files:
        raise ValueError(f"No MAG FASTA files found in: {mags}")

    # For v0.1: ignore rename_map unless you later want it for non-sequential mode
    # (keeping signature for compatibility with your CLI).
    if rename_map is not None and sequential_ids:
        raise ValueError("--rename-map cannot be used together with --sequential-ids in v0.1")

    # Write mapping header
    map_lines = ["original_filename\tnew_id"]

    copies: list[dict[str, str]] = []
    i = 1
    for src in files:
        original = src.name

        if sequential_ids:
            new_id = f"{out_prefix}{i:0{pad}d}"
            dst = out_mags / f"{new_id}_genomic.fna.gz"
            data = _as_fna_bytes(src)
            if dst.exists() and not force:
                raise FileExistsError(f"Destination exists: {dst} (use --force)")
            _write_gz_bytes(dst, data)
        else:
            # Preserve ID (stem) and standardise suffix to .fna (not gz) for simplicity
            # You can change this to always .fna.gz if you prefer.
            new_id = src.stem  # note: for .fa.gz, stem becomes ".fa" -> not ideal; handle below
            name_lower = src.name.lower()
            if name_lower.endswith(".fa.gz") or name_lower.endswith(".fna.gz") or name_lower.endswith(".fasta.gz"):
                # For *.gz, Path.stem strips only ".gz", leaving ".fa" etc. behind.
                # Normalize by removing the second suffix as well.
                new_id = src.name[: -len(".gz")]
                for suff in (".fa", ".fna", ".fasta"):
                    if new_id.lower().endswith(suff):
                        new_id = new_id[: -len(suff)]
                        break

            dst = out_mags / f"{new_id}.fna"
            if dst.exists() and not force:
                raise FileExistsError(f"Destination exists: {dst} (use --force)")
            # Just copy bytes (decompress if needed)
            data = _as_fna_bytes(src)
            dst.write_bytes(data)

        map_lines.append(f"{original}\t{new_id}")
        copies.append({"src": str(src), "dest": str(dst), "new_id": new_id})
        i += 1

    map_fp.write_text("\n".join(map_lines) + "\n", encoding="utf-8")

    report = {
        "mags_dir": str(mags),
        "out_dir": str(out),
        "n_inputs": len(files),
        "sequential_ids": sequential_ids,
        "out_prefix": out_prefix if sequential_ids else None,
        "pad": pad if sequential_ids else None,
        "outputs": {
            "mags_dir": str(out_mags),
            "id_map_tsv": str(map_fp),
            "prep_report_json": str(report_fp),
        },
        "copies": copies,
    }
    report_fp.write_text(json.dumps(report, indent=2), encoding="utf-8")