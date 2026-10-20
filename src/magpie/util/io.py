from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
from pathlib import Path
from typing import IO, Dict, Iterable, Iterator, List, Sequence, Tuple

FASTA_EXTS: Tuple[str, ...] = (".fa", ".fna", ".fasta")
GZ_EXT = ".gz"

# Genome IDs end up as FASTA headers, PHYLIP names and Newick tip labels, so they must not
# contain whitespace or characters with a special meaning in those formats.
_UNSAFE_CHARS_RE = re.compile(r"[^A-Za-z0-9._\-]")


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def open_text(path: Path, mode: str = "rt") -> IO[str]:
    """Open a plain or gzipped text file."""
    if path.name.lower().endswith(GZ_EXT):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")  # type: ignore[return-value]
    return path.open(mode.replace("t", ""), encoding="utf-8", errors="replace")


def is_fasta_name(name: str) -> bool:
    low = name.lower()
    if low.endswith(GZ_EXT):
        low = low[: -len(GZ_EXT)]
    return low.endswith(FASTA_EXTS)


def strip_fasta_suffix(name: str) -> str:
    """
    'bin.1.fa.gz' -> 'bin.1'; 'MAG0001_genomic.fna.gz' -> 'MAG0001'.
    Dots inside IDs (e.g. BC13.bin.1.603) are preserved. A trailing '_genomic' (NCBI/GTDB and
    MAGPIE file naming) is not part of the genome ID.
    """
    s = name.strip()
    if s.lower().endswith(GZ_EXT):
        s = s[: -len(GZ_EXT)]
    for ext in FASTA_EXTS:
        if s.lower().endswith(ext):
            s = s[: -len(ext)]
            break
    if s.endswith("_genomic"):
        s = s[: -len("_genomic")]
    return s


def list_fasta_files(dir_: Path) -> List[Path]:
    if not dir_.exists():
        return []
    return sorted(
        (p for p in dir_.iterdir() if p.is_file() and not p.name.startswith(".") and is_fasta_name(p.name)),
        key=lambda p: p.name,
    )


def sanitise_id(genome_id: str) -> str:
    return _UNSAFE_CHARS_RE.sub("_", genome_id.strip())


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """Yield (header without '>', sequence) from a plain or gzipped FASTA."""
    header: str | None = None
    parts: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Malformed FASTA (sequence before first header) in {path}")
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    return list(iter_fasta(path))


def read_fasta_ids(path: Path) -> List[str]:
    ids: List[str] = []
    if not path.exists():
        return ids
    with open_text(path) as fh:
        for raw in fh:
            if raw.startswith(">"):
                ids.append(raw[1:].strip().split()[0])
    return ids


def has_fasta_records(path: Path) -> bool:
    if not path.exists():
        return False
    with open_text(path) as fh:
        return any(raw.startswith(">") for raw in fh)


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 0) -> int:
    """Write records; width > 0 wraps sequence lines (60 = Biopython's FASTA writer)."""
    n = 0
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as out:
        for header, seq in records:
            out.write(f">{header}\n")
            if width > 0:
                for i in range(0, len(seq), width):
                    out.write(seq[i: i + width] + "\n")
            else:
                out.write(seq + "\n")
            n += 1
    return n


def copy_or_decompress(src: Path, dst: Path) -> None:
    """Copy src to dst, decompressing if src is gzipped and dst is not (and vice versa)."""
    ensure_dir(dst.parent)
    src_gz = src.name.lower().endswith(GZ_EXT)
    dst_gz = dst.name.lower().endswith(GZ_EXT)
    if src_gz and not dst_gz:
        with gzip.open(src, "rb") as fin, dst.open("wb") as fout:
            shutil.copyfileobj(fin, fout)
    elif dst_gz and not src_gz:
        with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src, dst)


def gzip_file(src: Path, dst: Path) -> None:
    ensure_dir(dst.parent)
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(list(header))
        for r in rows:
            w.writerow([str(x) for x in r])


def read_tsv_dicts(path: Path) -> List[Dict[str, str]]:
    with open_text(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_json(path: Path, obj: object) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(obj, indent=2, default=str) + "\n", encoding="utf-8")


def read_id_list(path: Path) -> List[str]:
    if not path.exists():
        return []
    return [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]


def write_id_list(path: Path, ids: Iterable[str]) -> None:
    ids = list(ids)
    ensure_dir(path.parent)
    path.write_text("\n".join(ids) + ("\n" if ids else ""), encoding="utf-8")


def read_phylip(path: Path) -> List[Tuple[str, str]]:
    """Read a sequential (relaxed) PHYLIP alignment as written by RAxML-NG."""
    records: List[Tuple[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().split()
        if len(header) < 2:
            return []
        n_taxa, aln_len = int(header[0]), int(header[1])
        for raw in fh:
            parts = raw.split(maxsplit=1)
            if len(parts) < 2:
                continue
            records.append((parts[0], parts[1].replace(" ", "").strip()))
    if len(records) != n_taxa:
        raise ValueError(f"PHYLIP header in {path} declares {n_taxa} taxa but {len(records)} were read")
    if any(len(seq) != aln_len for _, seq in records):
        raise ValueError(f"PHYLIP sequences in {path} do not all have the declared length {aln_len}")
    return records


def write_phylip(records: Sequence[Tuple[str, str]], path: Path) -> None:
    lengths = {len(seq) for _, seq in records}
    if len(lengths) != 1:
        raise ValueError(f"Cannot write PHYLIP {path}: sequences have unequal lengths {sorted(lengths)[:5]}")
    ensure_dir(path.parent)
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"{len(records)} {lengths.pop()}\n")
        for name, seq in records:
            fh.write(f"{name} {seq}\n")
