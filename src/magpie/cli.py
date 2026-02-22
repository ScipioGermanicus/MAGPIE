from __future__ import annotations

import logging
from pathlib import Path

import typer

from .steps.validate import validate_step
from .steps.prep import prep_step
from .steps.gtdbtk import gtdbtk_step
from .steps.taxonomy import taxonomy_step
from .steps.qc import qc_step

app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="MAGPIE: build components of a custom PICRUSt2 reference database from MAGs.",
)

_LOGGER = logging.getLogger("magpie")


def _configure_logging(verbose: bool, debug: bool) -> None:
    level = logging.WARNING
    if verbose:
        level = logging.INFO
    if debug:
        level = logging.DEBUG

    root = logging.getLogger()
    if not root.handlers:
        logging.basicConfig(
            level=level,
            format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
        )
    else:
        root.setLevel(level)


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Enable INFO-level logging."),
    debug: bool = typer.Option(False, "--debug", help="Enable DEBUG-level logging (very noisy)."),
) -> None:
    _configure_logging(verbose=verbose, debug=debug)
    _LOGGER.debug("Logging configured (verbose=%s, debug=%s).", verbose, debug)


@app.command("run")
def cmd_run(
    mags: Path = typer.Option(..., "--mags", exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True),
    out: Path = typer.Option(..., "--out", file_okay=False, dir_okay=True, resolve_path=True),

    # prep options
    rename_map: Path | None = typer.Option(None, "--rename-map", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    sequential_ids: bool = typer.Option(False, "--sequential-ids/--no-sequential-ids"),
    out_prefix: str = typer.Option("MAG", "--out-prefix"),
    pad: int = typer.Option(4, "--pad", min=1, max=12),

    # gtdb options
    gtdb_classify: Path | None = typer.Option(
        None,
        "--gtdb-classify",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Path to an existing GTDB-Tk *classify/* directory. If provided, GTDB-Tk is skipped.",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for GTDB-Tk / CheckM (if run)."),

    # taxonomy options
    move_tax_split: bool = typer.Option(
        False,
        "--move-tax-split",
        help="Move genomes into bacteria/archaea split dirs (default is copy; safer).",
    ),

    # QC (CheckM + filtering) options
    checkm_results: Path | None = typer.Option(
        None,
        "--checkm-results",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Optional: path to existing merged minimal CheckM TSV (genome_id,checkm_completeness,checkm_contamination). If provided, CheckM is skipped.",
    ),
    completeness_min: float = typer.Option(90.0, "--completeness-min", min=0.0, max=100.0),
    contamination_max: float = typer.Option(10.0, "--contamination-max", min=0.0, max=100.0),
    place_mode: str = typer.Option("copy", "--place-mode", help="How to place filtered genomes: copy (default), symlink, hardlink."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """
    Run MAGPIE workflow (validate -> prep -> gtdbtk (optional) -> taxonomy -> qc).
    """
    validate_dir = out / "01_validate"
    prep_dir = out / "02_prep"
    gtdb_dir = out / "03_gtdbtk"
    tax_dir = out / "04_taxonomy"
    qc_dir = out / "05_qc"

    _LOGGER.info("Running MAGPIE workflow in: %s", out)

    require_gtdbtk = gtdb_classify is None
    require_checkm = checkm_results is None  # only require CheckM if we will run it

    _LOGGER.info("Step 1/5: validate -> %s", validate_dir)
    validate_step(
        mags=mags,
        out=validate_dir,
        force=force,
        require_gtdbtk=require_gtdbtk,
        require_checkm=require_checkm,
    )

    _LOGGER.info("Step 2/5: prep -> %s", prep_dir)
    prep_step(
        mags=mags,
        out=prep_dir,
        rename_map=rename_map,
        force=force,
        sequential_ids=sequential_ids,
        out_prefix=out_prefix,
        pad=pad,
    )

    prepared_mags_dir = prep_dir / "mags"

    _LOGGER.info("Step 3/5: gtdbtk (or reuse) -> %s", gtdb_dir)
    classify_dir = gtdbtk_step(
        mags_dir=prepared_mags_dir,
        out=gtdb_dir,
        gtdb_classify=gtdb_classify,
        cpus=cpus,
        force=force,
    )

    _LOGGER.info("Step 4/5: taxonomy -> %s", tax_dir)
    taxonomy_step(
        prep_dir=prep_dir,
        classify_dir=classify_dir,
        out=tax_dir,
        force=force,
        move_files=move_tax_split,
    )

    _LOGGER.info("Step 5/5: qc -> %s", qc_dir)
    qc_step(
        tax_dir=tax_dir,
        out=qc_dir,
        cpus=cpus,
        force=force,
        checkm_results=checkm_results,
        completeness_min=completeness_min,
        contamination_max=contamination_max,
        place_mode=place_mode,
    )

    _LOGGER.info("MAGPIE run complete.")


@app.command("tax")
def cmd_tax(
    prep_dir: Path = typer.Option(
        ...,
        "--prep-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE prep output directory (e.g. out/02_prep).",
    ),
    gtdb_classify: Path = typer.Option(
        ...,
        "--gtdb-classify",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="GTDB-Tk classify directory containing summary TSVs.",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for taxonomy artefacts and split genomes.",
    ),
    move_files: bool = typer.Option(False, "--move-files", help="Move genomes instead of copying."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Infer domain taxonomy from GTDB-Tk summaries and split prepared genomes."""
    taxonomy_step(prep_dir=prep_dir, classify_dir=gtdb_classify, out=out, force=force, move_files=move_files)


@app.command("qc")
def cmd_qc(
    tax_dir: Path = typer.Option(
        ...,
        "--tax-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE taxonomy output directory (e.g. out/04_taxonomy).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for QC artefacts and filtered genomes.",
    ),
    checkm_results: Path | None = typer.Option(
        None,
        "--checkm-results",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Optional: merged minimal CheckM TSV. If provided, CheckM is skipped.",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for CheckM (if run)."),
    completeness_min: float = typer.Option(90.0, "--completeness-min", min=0.0, max=100.0),
    contamination_max: float = typer.Option(10.0, "--contamination-max", min=0.0, max=100.0),
    place_mode: str = typer.Option("copy", "--place-mode", help="copy (default), symlink, hardlink."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Run CheckM (optional) and filter genomes by completeness/contamination."""
    qc_step(
        tax_dir=tax_dir,
        out=out,
        cpus=cpus,
        force=force,
        checkm_results=checkm_results,
        completeness_min=completeness_min,
        contamination_max=contamination_max,
        place_mode=place_mode,
    )