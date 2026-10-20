"""
Resolve genome names used by external tools/files to MAGPIE genome IDs.

GTDB-Tk summaries, CheckM tables and eggNOG annotation files may refer to a genome by its
original file name, its original file stem, or the MAGPIE ID (with or without '_genomic').
prep writes id_map.tsv; IdLookup maps all of these names to the MAGPIE ID.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Iterable, List, Optional

from .io import strip_fasta_suffix


class IdLookup:
    def __init__(self, id_map_tsv: Path):
        self._map: Dict[str, str] = {}
        self.ids: List[str] = []
        with id_map_tsv.open("r", encoding="utf-8", newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                gid = row["new_id"]
                self.ids.append(gid)
                for key in (row["original_filename"], row.get("original_id") or "", gid,
                            strip_fasta_suffix(row["original_filename"])):
                    if key:
                        self._map[key] = gid

    def resolve(self, name: str) -> Optional[str]:
        name = name.strip()
        for cand in (name, strip_fasta_suffix(name)):
            if cand in self._map:
                return self._map[cand]
        return None

    def resolve_all(self, names: Iterable[str]) -> Dict[str, str]:
        """name -> genome ID for every name that resolves."""
        out = {}
        for n in names:
            gid = self.resolve(n)
            if gid is not None:
                out[n] = gid
        return out
