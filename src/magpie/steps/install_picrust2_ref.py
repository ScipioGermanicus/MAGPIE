from __future__ import annotations

import csv
import gzip
import json
import logging
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

_LOGGER = logging.getLogger("magpie.install_picrust2_ref")


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def _find_picrust2_default_files(picrust2_env: Path) -> Path:
    """
    Find <conda-env>/lib/pythonX/site-packages/picrust2/default_files.
    """
    if not picrust2_env.exists():
        raise FileNotFoundError(f"PICRUSt2 environment path does not exist: {picrust2_env}")

    candidates = sorted(picrust2_env.glob("lib/python*/site-packages/picrust2/default_files"))

    if not candidates:
        raise FileNotFoundError(
            f"Could not find picrust2/default_files under: {picrust2_env}\n"
            f"Expected something like:\n"
            f"  {picrust2_env}/lib/python3.*/site-packages/picrust2/default_files"
        )

    if len(candidates) > 1:
        joined = "\n".join(f"  - {p}" for p in candidates)
        raise ValueError(
            f"Multiple PICRUSt2 default_files directories found under {picrust2_env}:\n{joined}\n"
            f"Please remove ambiguity or use a cleaner PICRUSt2 environment."
        )

    return candidates[0]


def _validate_package_ref(package_ref_dir: Path) -> None:
    required = [
        package_ref_dir / "bac_ref" / "bac_ref.fna",
        package_ref_dir / "bac_ref" / "bac_ref.hmm",
        package_ref_dir / "bac_ref" / "bac_ref.tre",
        package_ref_dir / "bac_ref" / "bac_ref.model",
        package_ref_dir / "bac_ref" / "bac_ref.raxml_info",
        package_ref_dir / "bac_ref" / "ko.txt.gz",
        package_ref_dir / "bac_ref" / "ec.txt.gz",
        package_ref_dir / "arc_ref" / "arc_ref.fna",
        package_ref_dir / "arc_ref" / "arc_ref.hmm",
        package_ref_dir / "arc_ref" / "arc_ref.tre",
        package_ref_dir / "arc_ref" / "arc_ref.model",
        package_ref_dir / "arc_ref" / "arc_ref.raxml_info",
        package_ref_dir / "arc_ref" / "ko.txt.gz",
        package_ref_dir / "arc_ref" / "ec.txt.gz",
        package_ref_dir / "bacteria_16S_copies.txt",
        package_ref_dir / "archaea_16S_copies.txt",
    ]

    missing = [p for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError(
            "Package reference directory is incomplete. Missing:\n"
            + "\n".join(f"  - {p}" for p in missing)
        )


def _gzip_copy(src: Path, dst: Path) -> None:
    """
    Copy plain text src to gzipped dst.
    """
    _ensure_dir(dst.parent)
    with src.open("rb") as fin, gzip.open(dst, "wb") as fout:
        shutil.copyfileobj(fin, fout)


def _plain_copy(src: Path, dst: Path) -> None:
    _ensure_dir(dst.parent)
    shutil.copy2(src, dst)


def _backup_file(path: Path, backup_root: Path) -> Path | None:
    if not path.exists():
        return None

    rel = path.relative_to(path.parents[2]) if len(path.parents) >= 3 else Path(path.name)
    dst = backup_root / rel
    _ensure_dir(dst.parent)
    shutil.copy2(path, dst)
    return dst


def _target_dirs(default_files_dir: Path) -> Tuple[Path, Path]:
    """
    Resolve PICRUSt2 bacteria/archaea default directories.

    This assumes the default files are structured as:
      default_files/bacteria
      default_files/archaea
    """
    bac_dir = default_files_dir / "bacteria"
    arc_dir = default_files_dir / "archaea"

    if not bac_dir.exists():
        raise FileNotFoundError(f"Could not find PICRUSt2 bacterial default directory: {bac_dir}")
    if not arc_dir.exists():
        raise FileNotFoundError(f"Could not find PICRUSt2 archaeal default directory: {arc_dir}")

    return bac_dir, arc_dir


def _find_copy_count_target(domain_dir: Path, domain: str) -> Path:
    """
    Try to replace the existing PICRUSt2 16S copy-count file if present.

    Different installs/custom workflows may use slightly different names.
    """
    candidates = [
        domain_dir / "16S.txt",
        domain_dir / "16S.txt.gz",
        domain_dir / "16S_copies.txt",
        domain_dir / "16S_copies.txt.gz",
        domain_dir / f"{domain}_16S_copies.txt",
        domain_dir / f"{domain}_16S_copies.txt.gz",
    ]

    for p in candidates:
        if p.exists():
            return p

    # Stable fallback if no existing file is found.
    return domain_dir / "16S.txt"


def _install_file(
    *,
    src: Path,
    dst: Path,
    force: bool,
    backup_root: Path | None,
    gzip_output: bool = False,
) -> Dict[str, str]:
    if dst.exists() and not force:
        raise FileExistsError(
            f"Target file already exists: {dst}\n"
            f"Use --force to overwrite."
        )

    backup_path = None
    if backup_root is not None:
        backup_path = _backup_file(dst, backup_root)

    if gzip_output:
        _gzip_copy(src, dst)
    else:
        _plain_copy(src, dst)

    return {
        "source": str(src),
        "target": str(dst),
        "backup": str(backup_path) if backup_path is not None else "",
        "gzip_output": str(gzip_output),
    }


def install_picrust2_ref_step(
    *,
    package_ref_dir: Path,
    picrust2_env: Path,
    out: Path,
    force: bool,
    backup: bool,
) -> None:
    """
    Install a MAGPIE-produced package_ref directory into a selected PICRUSt2 conda environment.

    This modifies:
      <picrust2_env>/lib/python*/site-packages/picrust2/default_files

    It should only be called explicitly through:
      magpie install-picrust2-ref
    """
    package_ref_dir = package_ref_dir.resolve()
    picrust2_env = picrust2_env.resolve()
    out = out.resolve()

    _validate_package_ref(package_ref_dir)

    default_files_dir = _find_picrust2_default_files(picrust2_env)
    bac_target_dir, arc_target_dir = _target_dirs(default_files_dir)

    _ensure_dir(out)

    timestamp = time.strftime("%Y%m%d-%H%M%S")
    backup_root = out / f"backup_picrust2_default_files_{timestamp}" if backup else None
    if backup_root is not None:
        _ensure_dir(backup_root)

    operations: List[Dict[str, str]] = []

    # Bacteria reference files
    bac_src = package_ref_dir / "bac_ref"
    operations.append(_install_file(src=bac_src / "bac_ref.fna", dst=bac_target_dir / "bac_ref.fna.gz", force=force, backup_root=backup_root, gzip_output=True))
    operations.append(_install_file(src=bac_src / "bac_ref.hmm", dst=bac_target_dir / "bac_ref.hmm", force=force, backup_root=backup_root))
    operations.append(_install_file(src=bac_src / "bac_ref.tre", dst=bac_target_dir / "bac_ref.tre", force=force, backup_root=backup_root))
    operations.append(_install_file(src=bac_src / "bac_ref.model", dst=bac_target_dir / "bac_ref.model", force=force, backup_root=backup_root))
    operations.append(_install_file(src=bac_src / "bac_ref.raxml_info", dst=bac_target_dir / "bac_ref.raxml_info", force=force, backup_root=backup_root))
    operations.append(_install_file(src=bac_src / "ko.txt.gz", dst=bac_target_dir / "ko.txt.gz", force=force, backup_root=backup_root))
    operations.append(_install_file(src=bac_src / "ec.txt.gz", dst=bac_target_dir / "ec.txt.gz", force=force, backup_root=backup_root))

    bac_copy_target = _find_copy_count_target(bac_target_dir, "bacteria")
    operations.append(_install_file(
        src=package_ref_dir / "bacteria_16S_copies.txt",
        dst=bac_copy_target,
        force=force,
        backup_root=backup_root,
        gzip_output=bac_copy_target.name.endswith(".gz"),
    ))

    # Archaea reference files
    arc_src = package_ref_dir / "arc_ref"
    operations.append(_install_file(src=arc_src / "arc_ref.fna", dst=arc_target_dir / "arc_ref.fna.gz", force=force, backup_root=backup_root, gzip_output=True))
    operations.append(_install_file(src=arc_src / "arc_ref.hmm", dst=arc_target_dir / "arc_ref.hmm", force=force, backup_root=backup_root))
    operations.append(_install_file(src=arc_src / "arc_ref.tre", dst=arc_target_dir / "arc_ref.tre", force=force, backup_root=backup_root))
    operations.append(_install_file(src=arc_src / "arc_ref.model", dst=arc_target_dir / "arc_ref.model", force=force, backup_root=backup_root))
    operations.append(_install_file(src=arc_src / "arc_ref.raxml_info", dst=arc_target_dir / "arc_ref.raxml_info", force=force, backup_root=backup_root))
    operations.append(_install_file(src=arc_src / "ko.txt.gz", dst=arc_target_dir / "ko.txt.gz", force=force, backup_root=backup_root))
    operations.append(_install_file(src=arc_src / "ec.txt.gz", dst=arc_target_dir / "ec.txt.gz", force=force, backup_root=backup_root))

    arc_copy_target = _find_copy_count_target(arc_target_dir, "archaea")
    operations.append(_install_file(
        src=package_ref_dir / "archaea_16S_copies.txt",
        dst=arc_copy_target,
        force=force,
        backup_root=backup_root,
        gzip_output=arc_copy_target.name.endswith(".gz"),
    ))

    manifest = {
        "package_ref_dir": str(package_ref_dir),
        "picrust2_env": str(picrust2_env),
        "picrust2_default_files": str(default_files_dir),
        "backup": backup,
        "backup_root": str(backup_root) if backup_root is not None else None,
        "force": force,
        "operations": operations,
    }

    manifest_path = out / "install_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    with (out / "install_manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["source", "target", "backup", "gzip_output"], delimiter="\t")
        w.writeheader()
        w.writerows(operations)

    _LOGGER.info("Installed MAGPIE reference into PICRUSt2 default files: %s", default_files_dir)
    _LOGGER.info("Install manifest written to: %s", manifest_path)