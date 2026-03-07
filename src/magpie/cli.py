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


def _validate_checkm_reuse_args(
    *,
    checkm_qa: Path | None,
    checkm_qa_bacteria: Path | None,
    checkm_qa_archaea: Path | None,
    checkm_results: Path | None,
) -> None:
    # Mutually exclusive sources:
    # - single merged QA
    # - two QA files (bac+arc, bac-only, arc-only)
    # - legacy minimal TSV
    if checkm_results is not None:
        if checkm_qa is not None or checkm_qa_bacteria is not None or checkm_qa_archaea is not None:
            raise typer.BadParameter("Use either --checkm-results (legacy) OR CheckM QA options, not both.")

    if checkm_qa is not None:
        if checkm_qa_bacteria is not None or checkm_qa_archaea is not None:
            raise typer.BadParameter("Use either --checkm-qa OR --checkm-qa-bacteria/--checkm-qa-archaea, not both.")

    if (checkm_qa_bacteria is None) ^ (checkm_qa_archaea is None):
        # Allow one of them (e.g., bacteria-only datasets).
        # QC will handle any resulting missing IDs explicitly.
        pass


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
        help="Path to an existing GTDB-Tk classify/ directory. If provided, GTDB-Tk is skipped.",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for GTDB-Tk / CheckM (if run)."),

    # taxonomy options
    move_tax_split: bool = typer.Option(
        False,
        "--move-tax-split",
        help="Move genomes into bacteria/archaea split dirs (default is copy; safer).",
    ),

    # qc options
    skip_qc: bool = typer.Option(False, "--skip-qc", help="Stop after taxonomy (do not run QC)."),

    # Standard CheckM QA inputs (preferred)
    checkm_qa: Path | None = typer.Option(
        None,
        "--checkm-qa",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help=(
            "Path to a standard CheckM QA table produced by: "
            "`checkm qa <lineage.ms> <outdir> -o 2 -f checkm_qa.tsv`. "
            "Use this if you already have a single merged QA file."
        ),
    ),
    checkm_qa_bacteria: Path | None = typer.Option(
        None,
        "--checkm-qa-bacteria",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to standard CheckM QA table for bacteria (checkm qa ... -o 2 -f checkm_qa.tsv).",
    ),
    checkm_qa_archaea: Path | None = typer.Option(
        None,
        "--checkm-qa-archaea",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to standard CheckM QA table for archaea (checkm qa ... -o 2 -f checkm_qa.tsv).",
    ),

    # Legacy support
    checkm_results: Path | None = typer.Option(
        None,
        "--checkm-results",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help=(
            "LEGACY: Path to minimal CheckM TSV with columns: "
            "genome_id, checkm_completeness, checkm_contamination."
        ),
    ),

    completeness_min: float = typer.Option(90.0, "--completeness-min"),
    contamination_max: float = typer.Option(10.0, "--contamination-max"),
    place_mode: str = typer.Option("copy", "--place-mode", help="copy|symlink|hardlink"),
    checkm_bin: str = typer.Option("checkm", "--checkm-bin", help="CheckM executable name or path."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """
    Run MAGPIE workflow (validate -> prep -> gtdbtk (optional) -> taxonomy -> qc).
    """
    _validate_checkm_reuse_args(
        checkm_qa=checkm_qa,
        checkm_qa_bacteria=checkm_qa_bacteria,
        checkm_qa_archaea=checkm_qa_archaea,
        checkm_results=checkm_results,
    )

    validate_dir = out / "01_validate"
    prep_dir = out / "02_prep"
    gtdb_dir = out / "03_gtdbtk"
    tax_dir = out / "04_taxonomy"
    qc_dir = out / "05_qc"

    _LOGGER.info("Running MAGPIE workflow in: %s", out)

    require_gtdbtk = gtdb_classify is None

    # Require CheckM only if QC will run AND no reuse source was provided
    reuse_provided = (
        (checkm_qa is not None)
        or (checkm_qa_bacteria is not None)
        or (checkm_qa_archaea is not None)
        or (checkm_results is not None)
    )
    require_checkm = (not skip_qc) and (not reuse_provided)

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

    if skip_qc:
        _LOGGER.warning("Skipping QC (--skip-qc). Pipeline stops after taxonomy.")
        return

    _LOGGER.info("Step 5/5: qc -> %s", qc_dir)
    qc_step(
        tax_dir=tax_dir,
        out=qc_dir,
        cpus=cpus,
        force=force,
        checkm_qa=checkm_qa,
        checkm_qa_bacteria=checkm_qa_bacteria,
        checkm_qa_archaea=checkm_qa_archaea,
        checkm_results=checkm_results,
        completeness_min=completeness_min,
        contamination_max=contamination_max,
        place_mode=place_mode,
        checkm_bin=checkm_bin,
    )

    _LOGGER.info("MAGPIE run complete.")


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
        help="Output directory for QC artefacts (e.g. out/05_qc).",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for CheckM (if run)."),

    checkm_qa: Path | None = typer.Option(
        None,
        "--checkm-qa",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to a single merged standard CheckM QA table (checkm qa ... -o 2 -f checkm_qa.tsv).",
    ),
    checkm_qa_bacteria: Path | None = typer.Option(
        None,
        "--checkm-qa-bacteria",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to standard CheckM QA table for bacteria.",
    ),
    checkm_qa_archaea: Path | None = typer.Option(
        None,
        "--checkm-qa-archaea",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to standard CheckM QA table for archaea.",
    ),

    checkm_results: Path | None = typer.Option(
        None,
        "--checkm-results",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="LEGACY: minimal CheckM TSV (genome_id, checkm_completeness, checkm_contamination).",
    ),

    completeness_min: float = typer.Option(90.0, "--completeness-min"),
    contamination_max: float = typer.Option(10.0, "--contamination-max"),
    place_mode: str = typer.Option("copy", "--place-mode", help="copy|symlink|hardlink"),
    checkm_bin: str = typer.Option("checkm", "--checkm-bin"),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Run QC on MAGPIE taxonomy outputs (optionally reusing CheckM results)."""
    _validate_checkm_reuse_args(
        checkm_qa=checkm_qa,
        checkm_qa_bacteria=checkm_qa_bacteria,
        checkm_qa_archaea=checkm_qa_archaea,
        checkm_results=checkm_results,
    )

    qc_step(
        tax_dir=tax_dir,
        out=out,
        cpus=cpus,
        force=force,
        checkm_qa=checkm_qa,
        checkm_qa_bacteria=checkm_qa_bacteria,
        checkm_qa_archaea=checkm_qa_archaea,
        checkm_results=checkm_results,
        completeness_min=completeness_min,
        contamination_max=contamination_max,
        place_mode=place_mode,
        checkm_bin=checkm_bin,
    )


if __name__ == "__main__":
    app()