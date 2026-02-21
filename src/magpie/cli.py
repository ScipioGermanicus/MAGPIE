from __future__ import annotations

import logging
from pathlib import Path

import typer

from .steps.validate import validate_step
from .steps.prep import prep_step

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


@app.command("validate")
def cmd_validate(
    mags: Path = typer.Option(
        ...,
        "--mags",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing MAG FASTA files (.fa/.fna/.fasta).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Directory to write validation outputs (e.g., report.json).",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs in --out if present."),
) -> None:
    _LOGGER.info("Running validate: mags=%s out=%s", mags, out)
    validate_step(mags=mags, out=out, force=force)
    _LOGGER.info("Validate complete.")


@app.command("prep")
def cmd_prep(
    mags: Path = typer.Option(
        ...,
        "--mags",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing MAG FASTA files to prepare.",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Directory to write prepared MAGs and prep artefacts.",
    ),
    rename_map: Path | None = typer.Option(
        None,
        "--rename-map",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Optional TSV mapping old_id -> new_id (two columns, no header). "
             "Not compatible with sequential renaming in v0.1.",
    ),
    sequential_ids: bool = typer.Option(
        True,
        "--sequential-ids/--no-sequential-ids",
        help="Rename genomes to sequential IDs (default: enabled). Use --no-sequential-ids to preserve IDs.",
    ),
    out_prefix: str = typer.Option(
        "MAG",
        "--out-prefix",
        help="Prefix for sequential IDs (used when sequential renaming is enabled).",
    ),
    pad: int = typer.Option(
        4,
        "--pad",
        min=1,
        max=12,
        help="Zero-padding width for sequential IDs (used when sequential renaming is enabled).",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs in --out if present."),
) -> None:
    _LOGGER.info("Running prep: mags=%s out=%s", mags, out)
    prep_step(
        mags=mags,
        out=out,
        rename_map=rename_map,
        force=force,
        sequential_ids=sequential_ids,
        out_prefix=out_prefix,
        pad=pad,
    )
    _LOGGER.info("Prep complete.")


@app.command("run")
def cmd_run(
    mags: Path = typer.Option(
        ...,
        "--mags",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing MAG FASTA files.",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Working directory for MAGPIE outputs.",
    ),
    rename_map: Path | None = typer.Option(
        None,
        "--rename-map",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Optional TSV mapping old_id -> new_id (two columns, no header). "
             "Not compatible with sequential renaming in v0.1.",
    ),
    sequential_ids: bool = typer.Option(
        True,
        "--sequential-ids/--no-sequential-ids",
        help="Rename genomes to sequential IDs (default: enabled). Use --no-sequential-ids to preserve IDs.",
    ),
    out_prefix: str = typer.Option(
        "MAG",
        "--out-prefix",
        help="Prefix for sequential IDs (used when sequential renaming is enabled).",
    ),
    pad: int = typer.Option(
        4,
        "--pad",
        min=1,
        max=12,
        help="Zero-padding width for sequential IDs (used when sequential renaming is enabled).",
    ),
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs for steps in --out."),
) -> None:
    validate_dir = out / "01_validate"
    prep_dir = out / "02_prep"

    _LOGGER.info("Running MAGPIE v0.1 workflow in: %s", out)

    _LOGGER.info("Step 1/2: validate -> %s", validate_dir)
    validate_step(mags=mags, out=validate_dir, force=force)

    _LOGGER.info("Step 2/2: prep -> %s", prep_dir)
    prep_step(
        mags=mags,
        out=prep_dir,
        rename_map=rename_map,
        force=force,
        sequential_ids=sequential_ids,
        out_prefix=out_prefix,
        pad=pad,
    )

    _LOGGER.info("MAGPIE run complete.")