"""
Install a MAGPIE database as the default reference of a PICRUSt2 (v2.6+) conda environment.

PICRUSt2 2.6 reads <env>/lib/python*/site-packages/picrust2/default_files/:
  bacteria/bac_ref/bac_ref.{fna,hmm,model,raxml_info,tre}
  bacteria/16S.txt.gz, bacteria/ko.txt.gz, bacteria/ec.txt.gz (+ other optional trait tables)
  archaea/arc_ref/arc_ref.*, archaea/16S.txt.gz, ...

Only the domains present in the MAGPIE database are replaced. Every replaced file is backed up
first. Default trait tables that MAGPIE does not build (go, pfam, cazy, ...) and the metadata
table describe the old reference genomes, so they are moved to the backup for the replaced domains.
"""
from __future__ import annotations

import gzip
import json
import logging
import shutil
import time
from pathlib import Path
from typing import Dict, List

from ..pipeline import MagpieError
from ..util.io import ensure_dir
from .picrust2_db import check_db

_LOGGER = logging.getLogger("magpie.install_picrust2_ref")

REF_NAME = {"bacteria": "bac_ref", "archaea": "arc_ref"}


def find_default_files(picrust2_env: Path) -> Path:
    found = {p.resolve() for p in picrust2_env.glob("lib/python*/site-packages/picrust2/default_files")
             if (p / "bacteria").is_dir()}
    if not found:
        raise MagpieError(f"No picrust2/default_files found under {picrust2_env}/lib/python*/site-packages")
    if len(found) > 1:
        raise MagpieError("Several picrust2/default_files directories found:\n" + "\n".join(map(str, found)))
    return found.pop()


def _gzip_to(src: Path, dst: Path) -> None:
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def install_picrust2_ref_step(*, db: Path, picrust2_env: Path, out: Path, backup: bool = True,
                              dry_run: bool = False) -> List[Dict[str, str]]:
    problems = check_db(db)
    if problems:
        raise MagpieError("The database is not consistent:\n" + "\n".join(problems))
    domains = [d for d, ref in REF_NAME.items() if (db / ref).is_dir()]
    if not domains:
        raise MagpieError(f"No bac_ref/ or arc_ref/ found in {db}")
    default_files = find_default_files(picrust2_env)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    backup_root = out / f"backup_picrust2_default_files_{stamp}"
    ops: List[Dict[str, str]] = []

    for domain in domains:
        ref = REF_NAME[domain]
        target = default_files / domain
        new: Dict[str, Path] = {f"{ref}/{ref}.{e}": db / ref / f"{ref}.{e}" for e in ("fna", "hmm", "model", "raxml_info", "tre")}
        new["16S.txt.gz"] = db / f"{domain}_16S_copies.txt"
        for table in sorted((db / ref).glob("*.txt.gz")):
            new[table.name] = table
        stale = [p for p in sorted(target.glob("*.txt.gz")) + sorted(target.glob("*_metadata.csv.gz"))
                 if p.name not in new]
        stale += [p for p in sorted((target / ref).glob("*")) if f"{ref}/{p.name}" not in new]

        for path in [target / rel for rel in new] + stale:
            if backup and path.exists():
                dst = backup_root / domain / path.relative_to(target)
                ops.append({"action": "backup", "file": str(path), "to": str(dst)})
                if not dry_run:
                    ensure_dir(dst.parent)
                    shutil.copy2(path, dst)
        for path in stale:
            ops.append({"action": "remove", "file": str(path)})
            if not dry_run:
                path.unlink()
        for rel, src in new.items():
            dst = target / rel
            ops.append({"action": "install", "file": str(dst), "from": str(src)})
            if not dry_run:
                ensure_dir(dst.parent)
                if rel == "16S.txt.gz":
                    _gzip_to(src, dst)
                else:
                    shutil.copyfile(src, dst)

    ensure_dir(out)
    (out / f"install_manifest_{stamp}.json").write_text(json.dumps(
        {"db": str(db), "default_files": str(default_files), "domains": domains, "backup": str(backup_root) if backup else None,
         "dry_run": dry_run, "operations": ops}, indent=2), encoding="utf-8")
    _LOGGER.info("%s %s into %s%s", "Would install" if dry_run else "Installed", ", ".join(REF_NAME[d] for d in domains),
                 default_files, f" (backup: {backup_root})" if backup and not dry_run else "")
    return ops
