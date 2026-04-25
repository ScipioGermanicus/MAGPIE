from __future__ import annotations

import logging
from pathlib import Path

import typer

from .steps.validate import validate_step
from .steps.prep import prep_step
from .steps.gtdbtk import gtdbtk_step
from .steps.taxonomy import taxonomy_step
from .steps.qc import qc_step
from .steps.barrnap import barrnap_step
from .steps.rrna import rrna_step
from .steps.align import align_step
from .steps.choose_best import choose_best_step
from .steps.raxml_check import raxml_check_step
from .steps.iqtree import iqtree_step
from .steps.raxml_evaluate import raxml_evaluate_step
from .steps.hmm_prep import hmm_prep_step
from .steps.hmm_build import hmm_build_step
from .steps.package_ref import package_ref_step
from .steps.ko_table import ko_table_step

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
    if checkm_results is not None:
        if checkm_qa is not None or checkm_qa_bacteria is not None or checkm_qa_archaea is not None:
            raise typer.BadParameter("Use either --checkm-results (legacy) OR CheckM QA options, not both.")

    if checkm_qa is not None:
        if checkm_qa_bacteria is not None or checkm_qa_archaea is not None:
            raise typer.BadParameter("Use either --checkm-qa OR --checkm-qa-bacteria/--checkm-qa-archaea, not both.")

    if (checkm_qa_bacteria is None) ^ (checkm_qa_archaea is None):
        raise typer.BadParameter("Provide both --checkm-qa-bacteria and --checkm-qa-archaea together.")


@app.command("run")
def cmd_run(
    mags: Path = typer.Option(..., "--mags", exists=True, file_okay=False, dir_okay=True, readable=True, resolve_path=True),
    out: Path = typer.Option(..., "--out", file_okay=False, dir_okay=True, resolve_path=True),

    rename_map: Path | None = typer.Option(None, "--rename-map", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    sequential_ids: bool = typer.Option(False, "--sequential-ids/--no-sequential-ids"),
    out_prefix: str = typer.Option("MAG", "--out-prefix"),
    pad: int = typer.Option(4, "--pad", min=1, max=12),

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
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for external tools."),

    move_tax_split: bool = typer.Option(
        False,
        "--move-tax-split",
        help="Move genomes into bacteria/archaea split dirs (default is copy; safer).",
    ),

    skip_qc: bool = typer.Option(False, "--skip-qc", help="Stop after taxonomy."),

    checkm_qa: Path | None = typer.Option(
        None,
        "--checkm-qa",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
        help="Path to a single merged CheckM QA table.",
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
        help="LEGACY: minimal CheckM TSV.",
    ),
    completeness_min: float = typer.Option(90.0, "--completeness-min"),
    contamination_max: float = typer.Option(10.0, "--contamination-max"),
    place_mode: str = typer.Option("copy", "--place-mode", help="copy|symlink|hardlink"),
    checkm_bin: str = typer.Option("checkm", "--checkm-bin"),

    skip_barrnap: bool = typer.Option(False, "--skip-barrnap", help="Stop after QC."),
    barrnap_bin: str = typer.Option("barrnap", "--barrnap-bin"),
    barrnap_reject: float = typer.Option(0.8, "--barrnap-reject", min=0.0, max=1.0),

    skip_rrna: bool = typer.Option(False, "--skip-rrna", help="Stop after Barrnap."),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin"),
    multi_cluster_id: float = typer.Option(0.90, "--multi-cluster-id", min=0.0, max=1.0),
    final_cluster_id: float = typer.Option(1.0, "--final-cluster-id", min=0.0, max=1.0),

    skip_align: bool = typer.Option(False, "--skip-align", help="Stop after rrna."),
    cmalign_bin: str = typer.Option("cmalign", "--cmalign-bin"),
    esl_reformat_bin: str = typer.Option("esl-reformat", "--esl-reformat-bin"),
    ssu_models_dir: Path = typer.Option(
        ...,
        "--ssu-models-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing archaea.cm and bacteria.cm covariance models.",
    ),
    cmalign_mxsize_archaea: int = typer.Option(4096, "--cmalign-mxsize-archaea", min=1),
    cmalign_mxsize_bacteria: int = typer.Option(8192, "--cmalign-mxsize-bacteria", min=1),

    skip_choose_best: bool = typer.Option(False, "--skip-choose-best", help="Stop after align."),

    skip_raxml_check: bool = typer.Option(False, "--skip-raxml-check", help="Stop after choose-best."),
    raxml_ng_bin: str = typer.Option("raxml-ng", "--raxml-ng-bin"),

    skip_iqtree: bool = typer.Option(False, "--skip-iqtree", help="Stop after raxml-check."),
    iqtree_bin: str = typer.Option("iqtree", "--iqtree-bin"),
    iqtree_bootstrap: int = typer.Option(1000, "--iqtree-bootstrap", min=0),
    iqtree_seed: int = typer.Option(12345, "--iqtree-seed"),

    skip_raxml_evaluate: bool = typer.Option(False, "--skip-raxml-evaluate", help="Stop after IQ-TREE."),

    skip_hmm_prep: bool = typer.Option(False, "--skip-hmm-prep", help="Stop after RAxML-evaluate."),

    skip_hmm_build: bool = typer.Option(False, "--skip-hmm-build", help="Stop after HMM-prep."),
    hmmbuild_bin: str = typer.Option("hmmbuild", "--hmmbuild-bin"),

    skip_package_ref: bool = typer.Option(False, "--skip-package-ref", help="Stop after HMM-build."),

    skip_ko_table: bool = typer.Option(False, "--skip-ko-table", help="Stop after package-ref."),
    eggnog_existing_dir: Path | None = typer.Option(
        None,
        "--eggnog-existing-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing existing eggNOG files named <genome_id>.emapper.annotations.",
    ),
    allow_missing_eggnog: bool = typer.Option(
        False,
        "--allow-missing-eggnog",
        help="Allow genomes without matching eggNOG annotations and write zero-filled KO rows.",
    ),

    force: bool = typer.Option(False, "--force"),
) -> None:
    """
    Run MAGPIE workflow through final PICRUSt2 marker-reference packaging and KO table construction.
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
    barrnap_dir = out / "06_barrnap"
    rrna_dir = out / "07_rrna"
    align_dir = out / "08_align"
    choose_best_dir = out / "09_choose_best"
    raxml_dir = out / "10_raxml_check"
    iqtree_dir = out / "11_iqtree"
    raxml_eval_dir = out / "12_raxml_evaluate"
    hmm_prep_dir = out / "13_hmm_prep"
    hmm_build_dir = out / "14_hmm_build"
    package_ref_dir = out / "15_package_ref"
    ko_table_dir = out / "16_ko_table"

    require_gtdbtk = gtdb_classify is None
    reuse_provided = (
        checkm_qa is not None
        or checkm_qa_bacteria is not None
        or checkm_qa_archaea is not None
        or checkm_results is not None
    )
    require_checkm = (not skip_qc) and (not reuse_provided)

    _LOGGER.info("Step 1/16: validate -> %s", validate_dir)
    validate_step(
        mags=mags,
        out=validate_dir,
        force=force,
        require_gtdbtk=require_gtdbtk,
        require_checkm=require_checkm,
    )

    _LOGGER.info("Step 2/16: prep -> %s", prep_dir)
    prep_step(
        mags=mags,
        out=prep_dir,
        rename_map=rename_map,
        force=force,
        sequential_ids=sequential_ids,
        out_prefix=out_prefix,
        pad=pad,
    )

    _LOGGER.info("Step 3/16: gtdbtk -> %s", gtdb_dir)
    classify_dir = gtdbtk_step(
        mags_dir=prep_dir / "mags",
        out=gtdb_dir,
        gtdb_classify=gtdb_classify,
        cpus=cpus,
        force=force,
    )

    _LOGGER.info("Step 4/16: taxonomy -> %s", tax_dir)
    taxonomy_step(
        prep_dir=prep_dir,
        classify_dir=classify_dir,
        out=tax_dir,
        force=force,
        move_files=move_tax_split,
    )

    if skip_qc:
        _LOGGER.warning("Skipping QC (--skip-qc).")
        return

    _LOGGER.info("Step 5/16: qc -> %s", qc_dir)
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

    if skip_barrnap:
        _LOGGER.warning("Skipping Barrnap (--skip-barrnap).")
        return

    _LOGGER.info("Step 6/16: barrnap -> %s", barrnap_dir)
    barrnap_step(
        qc_dir=qc_dir,
        out=barrnap_dir,
        cpus=cpus,
        force=force,
        barrnap_bin=barrnap_bin,
        reject=barrnap_reject,
    )

    if skip_rrna:
        _LOGGER.warning("Skipping rrna (--skip-rrna).")
        return

    _LOGGER.info("Step 7/16: rrna -> %s", rrna_dir)
    rrna_step(
        barrnap_dir=barrnap_dir,
        out=rrna_dir,
        cpus=cpus,
        force=force,
        vsearch_bin=vsearch_bin,
        multi_cluster_id=multi_cluster_id,
        final_cluster_id=final_cluster_id,
    )

    if skip_align:
        _LOGGER.warning("Skipping align (--skip-align).")
        return

    _LOGGER.info("Step 8/16: align -> %s", align_dir)
    align_step(
        rrna_dir=rrna_dir,
        out=align_dir,
        cpus=cpus,
        force=force,
        cmalign_bin=cmalign_bin,
        esl_reformat_bin=esl_reformat_bin,
        ssu_models_dir=ssu_models_dir,
        mxsize_archaea=cmalign_mxsize_archaea,
        mxsize_bacteria=cmalign_mxsize_bacteria,
    )

    if skip_choose_best:
        _LOGGER.warning("Skipping choose-best (--skip-choose-best).")
        return

    _LOGGER.info("Step 9/16: choose-best -> %s", choose_best_dir)
    choose_best_step(
        prep_dir=prep_dir,
        qc_dir=qc_dir,
        rrna_dir=rrna_dir,
        align_dir=align_dir,
        out=choose_best_dir,
        force=force,
    )

    if skip_raxml_check:
        _LOGGER.warning("Skipping raxml-check (--skip-raxml-check).")
        return

    _LOGGER.info("Step 10/16: raxml-check -> %s", raxml_dir)
    raxml_check_step(
        choose_best_dir=choose_best_dir,
        out=raxml_dir,
        raxml_ng_bin=raxml_ng_bin,
        threads=cpus,
        force=force,
    )

    if skip_iqtree:
        _LOGGER.warning("Skipping iqtree (--skip-iqtree).")
        return

    _LOGGER.info("Step 11/16: iqtree -> %s", iqtree_dir)
    iqtree_step(
        raxml_check_dir=raxml_dir,
        out=iqtree_dir,
        iqtree_bin=iqtree_bin,
        threads=cpus,
        bootstrap=iqtree_bootstrap,
        seed=iqtree_seed,
        force=force,
    )

    if skip_raxml_evaluate:
        _LOGGER.warning("Skipping raxml-evaluate (--skip-raxml-evaluate).")
        return

    _LOGGER.info("Step 12/16: raxml-evaluate -> %s", raxml_eval_dir)
    raxml_evaluate_step(
        raxml_check_dir=raxml_dir,
        iqtree_dir=iqtree_dir,
        out=raxml_eval_dir,
        raxml_ng_bin=raxml_ng_bin,
        threads=cpus,
        force=force,
    )

    if skip_hmm_prep:
        _LOGGER.warning("Skipping hmm-prep (--skip-hmm-prep).")
        return

    _LOGGER.info("Step 13/16: hmm-prep -> %s", hmm_prep_dir)
    hmm_prep_step(
        raxml_check_dir=raxml_dir,
        out=hmm_prep_dir,
        esl_reformat_bin=esl_reformat_bin,
        force=force,
    )

    if skip_hmm_build:
        _LOGGER.warning("Skipping hmm-build (--skip-hmm-build).")
        return

    _LOGGER.info("Step 14/16: hmm-build -> %s", hmm_build_dir)
    hmm_build_step(
        hmm_prep_dir=hmm_prep_dir,
        out=hmm_build_dir,
        hmmbuild_bin=hmmbuild_bin,
        cpus=cpus,
        force=force,
    )

    if skip_package_ref:
        _LOGGER.warning("Skipping package-ref (--skip-package-ref).")
        return

    _LOGGER.info("Step 15/16: package-ref -> %s", package_ref_dir)
    package_ref_step(
        rrna_dir=rrna_dir,
        iqtree_dir=iqtree_dir,
        raxml_evaluate_dir=raxml_eval_dir,
        hmm_prep_dir=hmm_prep_dir,
        hmm_build_dir=hmm_build_dir,
        out=package_ref_dir,
        force=force,
    )

    if skip_ko_table:
        _LOGGER.warning("Skipping ko-table (--skip-ko-table).")
        return

    if eggnog_existing_dir is None:
        raise typer.BadParameter(
            "To build KO tables, provide --eggnog-existing-dir or use --skip-ko-table."
        )

    _LOGGER.info("Step 16/16: ko-table -> %s", ko_table_dir)
    ko_table_step(
        package_ref_dir=package_ref_dir,
        eggnog_existing_dir=eggnog_existing_dir,
        out=ko_table_dir,
        force=force,
        allow_missing_eggnog=allow_missing_eggnog,
    )

    _LOGGER.info("MAGPIE run complete.")


@app.command("ko-table")
def cmd_ko_table(
    package_ref_dir: Path = typer.Option(
        ...,
        "--package-ref-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE package-ref output directory, e.g. out/15_package_ref.",
    ),
    eggnog_existing_dir: Path = typer.Option(
        ...,
        "--eggnog-existing-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="Directory containing existing eggNOG files named <genome_id>.emapper.annotations.",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for KO table reports, e.g. out/16_ko_table.",
    ),
    allow_missing_eggnog: bool = typer.Option(
        False,
        "--allow-missing-eggnog",
        help="Allow genomes without matching eggNOG annotations and write zero-filled KO rows.",
    ),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Build PICRUSt2-style ko.txt.gz tables from existing eggNOG annotations."""
    ko_table_step(
        package_ref_dir=package_ref_dir,
        eggnog_existing_dir=eggnog_existing_dir,
        out=out,
        force=force,
        allow_missing_eggnog=allow_missing_eggnog,
    )


if __name__ == "__main__":
    app()