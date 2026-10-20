from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import List, Optional

import typer

from . import __version__
from .pipeline import STEP_DIRS, STEP_NAMES, Config, MagpieError, load_saved_config, run_pipeline, status
from .util.logging import configure_console

app = typer.Typer(
    no_args_is_help=True,
    add_completion=False,
    help="MAGPIE: build a custom PICRUSt2 reference database from MAGs.",
)

_LOGGER = logging.getLogger("magpie")

DIR_IN = dict(exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True)
FILE_IN = dict(exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True)
P_IO, P_STEPS, P_TAX, P_QC, P_16S, P_TREE, P_FUNC, P_BIN = (
    "Input / output", "Step selection", "Taxonomy (GTDB-Tk)", "Genome quality (CheckM)", "16S rRNA genes",
    "Alignment, trees and HMMs", "Functional annotation (eggNOG)", "Executables")


def _version(value: bool) -> None:
    if value:
        typer.echo(f"magpie {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Show DEBUG messages on the console."),
    quiet: bool = typer.Option(False, "--quiet", "-q", help="Only show warnings and errors."),
    debug: bool = typer.Option(False, "--debug", hidden=True),
    version: bool = typer.Option(False, "--version", callback=_version, is_eager=True, help="Print version and exit."),
) -> None:
    configure_console(logging.DEBUG if (verbose or debug) else logging.WARNING if quiet else logging.INFO)


def _run_or_exit(fn) -> None:
    try:
        fn()
    except MagpieError as e:
        _LOGGER.error("%s", e)
        sys.exit(1)
    except KeyboardInterrupt:
        _LOGGER.error("Interrupted. Rerun the same command to resume.")
        sys.exit(130)


def _steps_option(value: Optional[str]) -> Optional[List[str]]:
    if value is None:
        return None
    return [s.strip() for s in value.replace(" ", ",").split(",") if s.strip()]


@app.command("run")
def cmd_run(
    ctx: typer.Context,
    out: Path = typer.Option(..., "--out", file_okay=False, resolve_path=True, rich_help_panel=P_IO,
                             help="Output directory. Rerunning with the same directory resumes the workflow."),
    mags: Optional[Path] = typer.Option(None, "--mags", **DIR_IN, rich_help_panel=P_IO,
                                        help="Directory with MAG FASTA files (.fa/.fna/.fasta, optionally .gz)."),
    rename_map: Optional[Path] = typer.Option(None, "--rename-map", **FILE_IN, rich_help_panel=P_IO,
                                              help="TSV: original file name (or ID) <tab> new genome ID."),
    sequential_ids: bool = typer.Option(False, "--sequential-ids/--no-sequential-ids", rich_help_panel=P_IO,
                                        help="Rename genomes to MAG0001, MAG0002, ... (files sorted by name)."),
    out_prefix: str = typer.Option("MAG", "--out-prefix", rich_help_panel=P_IO),
    pad: int = typer.Option(4, "--pad", min=1, max=12, rich_help_panel=P_IO),
    cpus: int = typer.Option(8, "--cpus", min=1, rich_help_panel=P_IO, help="CPUs for external tools."),
    seed: int = typer.Option(100, "--seed", rich_help_panel=P_IO,
                             help="Seed for the random choices of the workflow (16S copy, tied genomes)."),

    only: Optional[str] = typer.Option(None, "--only", rich_help_panel=P_STEPS,
                                       help="Run only these steps (comma-separated), e.g. --only barrnap,rrna."),
    start: Optional[str] = typer.Option(None, "--from", rich_help_panel=P_STEPS,
                                        help="(Re)run from this step onwards."),
    stop: Optional[str] = typer.Option(None, "--to", rich_help_panel=P_STEPS,
                                       help="Stop after this step, e.g. --to package_ref to skip annotation."),
    force: bool = typer.Option(False, "--force", rich_help_panel=P_STEPS,
                               help="Recompute steps that are already done (eggNOG annotations are kept)."),
    reset_options: bool = typer.Option(False, "--reset-options", rich_help_panel=P_STEPS,
                                       help="Do not reuse the options saved by the previous run in --out."),

    gtdb_classify: Optional[Path] = typer.Option(
        None, "--gtdb-classify", **DIR_IN, rich_help_panel=P_TAX,
        help="Existing GTDB-Tk results (dir with gtdbtk.bac120/ar53.summary.tsv, or its parent). Skips GTDB-Tk."),
    pplacer_cpus: Optional[int] = typer.Option(None, "--pplacer-cpus", min=1, rich_help_panel=P_TAX),

    checkm_qa: Optional[Path] = typer.Option(None, "--checkm-qa", **FILE_IN, rich_help_panel=P_QC,
                                             help="Existing CheckM QA table for all genomes. Skips CheckM."),
    checkm_qa_bacteria: Optional[Path] = typer.Option(None, "--checkm-qa-bacteria", **FILE_IN, rich_help_panel=P_QC),
    checkm_qa_archaea: Optional[Path] = typer.Option(None, "--checkm-qa-archaea", **FILE_IN, rich_help_panel=P_QC),
    completeness_min: float = typer.Option(90.0, "--completeness-min", rich_help_panel=P_QC),
    contamination_max: float = typer.Option(10.0, "--contamination-max", rich_help_panel=P_QC),

    barrnap_reject: float = typer.Option(0.8, "--barrnap-reject", min=0.0, max=1.0, rich_help_panel=P_16S),
    multi_cluster_id: float = typer.Option(0.90, "--multi-cluster-id", min=0.0, max=1.0, rich_help_panel=P_16S),
    final_cluster_id: float = typer.Option(1.0, "--final-cluster-id", min=0.0, max=1.0, rich_help_panel=P_16S),
    multi_copy_pick: str = typer.Option("random", "--multi-copy-pick", rich_help_panel=P_16S,
                                        help="16S kept for genomes with divergent copies: random (wiki) or longest."),
    max_16s_copies: int = typer.Option(10, "--max-16s-copies", min=1, rich_help_panel=P_16S),

    iqtree_bootstrap: int = typer.Option(1000, "--iqtree-bootstrap", min=0, rich_help_panel=P_TREE),
    iqtree_seed: int = typer.Option(12345, "--iqtree-seed", rich_help_panel=P_TREE),
    min_genomes: int = typer.Option(4, "--min-genomes", min=4, rich_help_panel=P_TREE,
                                    help="Minimum genomes needed to build a domain's reference."),

    eggnog_existing_dir: Optional[Path] = typer.Option(
        None, "--eggnog-existing-dir", **DIR_IN, rich_help_panel=P_FUNC,
        help="Existing <genome>.emapper.annotations files. Without --eggnog-data-dir, MAGs without one are excluded."),
    eggnog_data_dir: Optional[Path] = typer.Option(None, "--eggnog-data-dir", **DIR_IN, rich_help_panel=P_FUNC,
                                                   help="eggNOG-mapper data dir; enables running eggNOG-mapper."),
    eggnog_jobs: int = typer.Option(1, "--eggnog-jobs", min=1, rich_help_panel=P_FUNC),
    keep_constant_traits: bool = typer.Option(
        False, "--keep-constant-traits", rich_help_panel=P_FUNC,
        help="Keep gene families with the same count in every genome (PICRUSt2's default mp method fails on them)."),

    gtdbtk_bin: str = typer.Option("gtdbtk", "--gtdbtk-bin", rich_help_panel=P_BIN),
    checkm_bin: str = typer.Option("checkm", "--checkm-bin", rich_help_panel=P_BIN),
    barrnap_bin: str = typer.Option("barrnap", "--barrnap-bin", rich_help_panel=P_BIN),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", rich_help_panel=P_BIN),
    ssu_align_bin: str = typer.Option("ssu-align", "--ssu-align-bin", rich_help_panel=P_BIN),
    esl_reformat_bin: str = typer.Option("esl-reformat", "--esl-reformat-bin", rich_help_panel=P_BIN),
    raxml_ng_bin: str = typer.Option("raxml-ng", "--raxml-ng-bin", rich_help_panel=P_BIN),
    iqtree_bin: str = typer.Option("iqtree", "--iqtree-bin", rich_help_panel=P_BIN),
    hmmbuild_bin: str = typer.Option("hmmbuild", "--hmmbuild-bin", rich_help_panel=P_BIN),
    emapper_bin: str = typer.Option("emapper.py", "--emapper-bin", rich_help_panel=P_BIN),
) -> None:
    """
    Run the MAGPIE workflow (or some of its steps) and write OUT/picrust2_db.

    Options that are not given are taken from the previous run in OUT (OUT/magpie_config.json),
    so a partial run such as `magpie run --out OUT --only iqtree` needs no other options.
    """
    params = dict(ctx.params)
    for k in ("only", "start", "stop", "force", "reset_options"):
        params.pop(k)
    saved = {} if reset_options else load_saved_config(out)
    for name, value in saved.items():
        if name in params and _is_default(ctx, name) and value is not None:
            params[name] = value
    for name in _PATH_OPTIONS:
        if params.get(name) is not None:
            params[name] = Path(params[name])
    cfg = Config(force=force, **params)
    _run_or_exit(lambda: run_pipeline(cfg, only=_steps_option(only), start=start, stop=stop))


def _is_default(ctx: typer.Context, name: str) -> bool:
    source = ctx.get_parameter_source(name)
    return source is None or getattr(source, "name", "") == "DEFAULT"


_PATH_OPTIONS = {"out", "mags", "rename_map", "gtdb_classify", "checkm_qa", "checkm_qa_bacteria", "checkm_qa_archaea",
                 "eggnog_existing_dir", "eggnog_data_dir"}


@app.command("steps")
def cmd_steps() -> None:
    """List the workflow steps (names for --only / --from / --to)."""
    for name in STEP_NAMES:
        typer.echo(f"{STEP_DIRS[name]:18s} {name}")


@app.command("status")
def cmd_status(out: Path = typer.Option(..., "--out", **DIR_IN)) -> None:
    """Show which steps are done in an output directory."""
    for dirname, name, state in status(out):
        typer.echo(f"{dirname:18s} {name:15s} {state}")


@app.command("check-db")
def cmd_check_db(db: Path = typer.Option(..., "--db", **DIR_IN, help="A MAGPIE picrust2_db directory.")) -> None:
    """Check that FASTA, tree, HMM, 16S copies and trait tables of a database agree."""
    from .steps.picrust2_db import check_db

    problems = check_db(db)
    if problems:
        typer.echo("\n".join(problems), err=True)
        raise typer.Exit(1)
    typer.echo(f"OK: {db}")


@app.command("install-picrust2-ref")
def cmd_install_picrust2_ref(
    db: Path = typer.Option(..., "--db", "--package-ref-dir", **DIR_IN, help="MAGPIE output OUT/picrust2_db."),
    picrust2_env: Path = typer.Option(..., "--picrust2-env", **DIR_IN,
                                      help="PICRUSt2 conda environment whose default files are replaced."),
    out: Path = typer.Option(..., "--out", file_okay=False, resolve_path=True,
                             help="Directory for the install manifest and the backup of the replaced files."),
    backup: bool = typer.Option(True, "--backup/--no-backup"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Only report what would be changed."),
) -> None:
    """Install a MAGPIE database as the default reference of a PICRUSt2 (v2.6+) environment."""
    from .steps.install_picrust2_ref import install_picrust2_ref_step

    _run_or_exit(lambda: install_picrust2_ref_step(db=db, picrust2_env=picrust2_env, out=out, backup=backup,
                                                   dry_run=dry_run))


if __name__ == "__main__":
    app()
