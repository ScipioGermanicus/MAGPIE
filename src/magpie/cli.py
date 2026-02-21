from __future__ import annotations

from pathlib import Path
import typer

from .steps.validate import validate_step
from .steps.prep import prep_step
from .steps.annotate import annotate_step
from .steps.build_db import build_db_step
from .pipeline import pipeline_run

app = typer.Typer(no_args_is_help=True, add_completion=False)

@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", help="Verbose logging."),
    debug: bool = typer.Option(False, "--debug", help="Debug logging."),
):
    # configure logging here if you want
    pass

@app.command("validate")
def cmd_validate(
    mags: Path = typer.Option(..., "--mags", exists=True, readable=True, help="Directory with MAG FASTAs."),
    out: Path = typer.Option(..., "--out", help="Output directory for validation report."),
    force: bool = typer.Option(False, "--force", help="Overwrite existing outputs."),
):
    validate_step(mags=mags, out=out, force=force)

@app.command("prep")
def cmd_prep(
    mags: Path = typer.Option(..., "--mags", exists=True, readable=True),
    out: Path = typer.Option(..., "--out"),
    rename_map: Path | None = typer.Option(None, "--rename-map", exists=True, readable=True),
    force: bool = typer.Option(False, "--force"),
):
    prep_step(mags=mags, out=out, rename_map=rename_map, force=force)

@app.command("annotate")
def cmd_annotate(
    mags: Path = typer.Option(..., "--mags", exists=True, readable=True),
    out: Path = typer.Option(..., "--out"),
    threads: int = typer.Option(8, "--threads"),
    force: bool = typer.Option(False, "--force"),
):
    annotate_step(mags=mags, out=out, threads=threads, force=force)

@app.command("build-db")
def cmd_build_db(
    annot: Path = typer.Option(..., "--annot", exists=True, readable=True),
    tree: Path = typer.Option(..., "--tree", exists=True, readable=True),
    out: Path = typer.Option(..., "--out"),
    force: bool = typer.Option(False, "--force"),
):
    build_db_step(annot=annot, tree=tree, out=out, force=force)

@app.command("pipeline")
def cmd_pipeline(
    mags: Path = typer.Option(..., "--mags", exists=True, readable=True),
    out: Path = typer.Option(..., "--out", help="Working directory."),
    db_out: Path = typer.Option(..., "--db-out", help="Final database output directory."),
    threads: int = typer.Option(8, "--threads"),
    force: bool = typer.Option(False, "--force"),
):
    pipeline_run(mags=mags, out=out, db_out=db_out, threads=threads, force=force)

if __name__ == "__main__":
    app()