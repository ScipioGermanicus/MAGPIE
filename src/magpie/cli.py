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
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for GTDB-Tk / CheckM / Barrnap / VSEARCH / alignment / tree / HMM tools."),

    # taxonomy options
    move_tax_split: bool = typer.Option(
        False,
        "--move-tax-split",
        help="Move genomes into bacteria/archaea split dirs (default is copy; safer).",
    ),

    # qc options
    skip_qc: bool = typer.Option(False, "--skip-qc", help="Stop after taxonomy (do not run QC)."),

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
    checkm_bin: str = typer.Option("checkm", "--checkm-bin", help="CheckM executable name or path."),

    # barrnap options
    skip_barrnap: bool = typer.Option(False, "--skip-barrnap", help="Stop after QC (do not run Barrnap)."),
    barrnap_bin: str = typer.Option("barrnap", "--barrnap-bin", help="Barrnap executable name or path."),
    barrnap_reject: float = typer.Option(0.8, "--barrnap-reject", min=0.0, max=1.0),

    # rrna options
    skip_rrna: bool = typer.Option(False, "--skip-rrna", help="Stop after Barrnap (do not run 16S representative selection/clustering)."),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin", help="VSEARCH executable name or path."),
    multi_cluster_id: float = typer.Option(0.90, "--multi-cluster-id", min=0.0, max=1.0),
    final_cluster_id: float = typer.Option(1.0, "--final-cluster-id", min=0.0, max=1.0),

    # align options
    skip_align: bool = typer.Option(False, "--skip-align", help="Stop after rrna (do not run alignment)."),
    cmalign_bin: str = typer.Option("cmalign", "--cmalign-bin", help="cmalign executable name or path."),
    esl_reformat_bin: str = typer.Option("esl-reformat", "--esl-reformat-bin", help="esl-reformat executable name or path."),
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

    # choose-best options
    skip_choose_best: bool = typer.Option(False, "--skip-choose-best", help="Stop after align (do not choose best genome per 16S cluster)."),

    # raxml-check options
    skip_raxml_check: bool = typer.Option(False, "--skip-raxml-check", help="Stop after choose-best (do not run raxml-ng --check)."),
    raxml_ng_bin: str = typer.Option("raxml-ng", "--raxml-ng-bin", help="raxml-ng executable name or path."),

    # iqtree options
    skip_iqtree: bool = typer.Option(False, "--skip-iqtree", help="Stop after raxml-check (do not build IQ-TREE trees)."),
    iqtree_bin: str = typer.Option("iqtree", "--iqtree-bin", help="IQ-TREE executable name or path."),
    iqtree_bootstrap: int = typer.Option(1000, "--iqtree-bootstrap", min=0, help="Ultrafast bootstrap replicates for IQ-TREE."),
    iqtree_seed: int = typer.Option(12345, "--iqtree-seed", help="Random seed for IQ-TREE."),

    # raxml-evaluate options
    skip_raxml_evaluate: bool = typer.Option(False, "--skip-raxml-evaluate", help="Stop after IQ-TREE (do not run RAxML-NG --evaluate)."),

    # hmm-prep options
    skip_hmm_prep: bool = typer.Option(False, "--skip-hmm-prep", help="Stop after RAxML-evaluate (do not convert reduced PHYLIP to FASTA/DNA FASTA/Stockholm)."),

    # hmm-build options
    skip_hmm_build: bool = typer.Option(False, "--skip-hmm-build", help="Stop after HMM-prep (do not run hmmbuild on Stockholm files)."),
    hmmbuild_bin: str = typer.Option("hmmbuild", "--hmmbuild-bin", help="hmmbuild executable name or path."),

    # package-ref options
    skip_package_ref: bool = typer.Option(False, "--skip-package-ref", help="Stop after HMM-build (do not package PICRUSt2-style reference folders)."),

    force: bool = typer.Option(False, "--force"),
) -> None:
    """
    Run MAGPIE workflow (validate -> prep -> gtdbtk (optional) -> taxonomy -> qc -> barrnap -> rrna -> align -> choose-best -> raxml-check -> iqtree -> raxml-evaluate -> hmm-prep -> hmm-build -> package-ref).
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

    _LOGGER.info("Running MAGPIE workflow in: %s", out)

    require_gtdbtk = gtdb_classify is None
    reuse_provided = (
        (checkm_qa is not None)
        or (checkm_qa_bacteria is not None)
        or (checkm_qa_archaea is not None)
        or (checkm_results is not None)
    )
    require_checkm = (not skip_qc) and (not reuse_provided)

    _LOGGER.info("Step 1/15: validate -> %s", validate_dir)
    validate_step(
        mags=mags,
        out=validate_dir,
        force=force,
        require_gtdbtk=require_gtdbtk,
        require_checkm=require_checkm,
    )

    _LOGGER.info("Step 2/15: prep -> %s", prep_dir)
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

    _LOGGER.info("Step 3/15: gtdbtk (or reuse) -> %s", gtdb_dir)
    classify_dir = gtdbtk_step(
        mags_dir=prepared_mags_dir,
        out=gtdb_dir,
        gtdb_classify=gtdb_classify,
        cpus=cpus,
        force=force,
    )

    _LOGGER.info("Step 4/15: taxonomy -> %s", tax_dir)
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

    _LOGGER.info("Step 5/15: qc -> %s", qc_dir)
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
        _LOGGER.warning("Skipping Barrnap (--skip-barrnap). Pipeline stops after QC.")
        return

    _LOGGER.info("Step 6/15: barrnap -> %s", barrnap_dir)
    barrnap_step(
        qc_dir=qc_dir,
        out=barrnap_dir,
        cpus=cpus,
        force=force,
        barrnap_bin=barrnap_bin,
        reject=barrnap_reject,
    )

    if skip_rrna:
        _LOGGER.warning("Skipping rrna (--skip-rrna). Pipeline stops after Barrnap.")
        return

    _LOGGER.info("Step 7/15: rrna -> %s", rrna_dir)
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
        _LOGGER.warning("Skipping align (--skip-align). Pipeline stops after rrna.")
        return

    _LOGGER.info("Step 8/15: align -> %s", align_dir)
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
        _LOGGER.warning("Skipping choose-best (--skip-choose-best). Pipeline stops after align.")
        return

    _LOGGER.info("Step 9/15: choose-best -> %s", choose_best_dir)
    choose_best_step(
        prep_dir=prep_dir,
        qc_dir=qc_dir,
        rrna_dir=rrna_dir,
        align_dir=align_dir,
        out=choose_best_dir,
        force=force,
    )

    if skip_raxml_check:
        _LOGGER.warning("Skipping raxml-check (--skip-raxml-check). Pipeline stops after choose-best.")
        return

    _LOGGER.info("Step 10/15: raxml-check -> %s", raxml_dir)
    raxml_check_step(
        choose_best_dir=choose_best_dir,
        out=raxml_dir,
        raxml_ng_bin=raxml_ng_bin,
        threads=cpus,
        force=force,
    )

    if skip_iqtree:
        _LOGGER.warning("Skipping iqtree (--skip-iqtree). Pipeline stops after raxml-check.")
        return

    _LOGGER.info("Step 11/15: iqtree -> %s", iqtree_dir)
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
        _LOGGER.warning("Skipping raxml-evaluate (--skip-raxml-evaluate). Pipeline stops after iqtree.")
        return

    _LOGGER.info("Step 12/15: raxml-evaluate -> %s", raxml_eval_dir)
    raxml_evaluate_step(
        raxml_check_dir=raxml_dir,
        iqtree_dir=iqtree_dir,
        out=raxml_eval_dir,
        raxml_ng_bin=raxml_ng_bin,
        threads=cpus,
        force=force,
    )

    if skip_hmm_prep:
        _LOGGER.warning("Skipping hmm-prep (--skip-hmm-prep). Pipeline stops after raxml-evaluate.")
        return

    _LOGGER.info("Step 13/15: hmm-prep -> %s", hmm_prep_dir)
    hmm_prep_step(
        raxml_check_dir=raxml_dir,
        out=hmm_prep_dir,
        esl_reformat_bin=esl_reformat_bin,
        force=force,
    )

    if skip_hmm_build:
        _LOGGER.warning("Skipping hmm-build (--skip-hmm-build). Pipeline stops after hmm-prep.")
        return

    _LOGGER.info("Step 14/15: hmm-build -> %s", hmm_build_dir)
    hmm_build_step(
        hmm_prep_dir=hmm_prep_dir,
        out=hmm_build_dir,
        hmmbuild_bin=hmmbuild_bin,
        cpus=cpus,
        force=force,
    )

    if skip_package_ref:
        _LOGGER.warning("Skipping package-ref (--skip-package-ref). Pipeline stops after hmm-build.")
        return

    _LOGGER.info("Step 15/15: package-ref -> %s", package_ref_dir)
    package_ref_step(
        rrna_dir=rrna_dir,
        iqtree_dir=iqtree_dir,
        raxml_evaluate_dir=raxml_eval_dir,
        hmm_prep_dir=hmm_prep_dir,
        hmm_build_dir=hmm_build_dir,
        out=package_ref_dir,
        force=force,
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
    checkm_qa: Path | None = typer.Option(None, "--checkm-qa", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    checkm_qa_bacteria: Path | None = typer.Option(None, "--checkm-qa-bacteria", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    checkm_qa_archaea: Path | None = typer.Option(None, "--checkm-qa-archaea", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    checkm_results: Path | None = typer.Option(None, "--checkm-results", exists=True, file_okay=True, dir_okay=False, readable=True, resolve_path=True),
    completeness_min: float = typer.Option(90.0, "--completeness-min"),
    contamination_max: float = typer.Option(10.0, "--contamination-max"),
    place_mode: str = typer.Option("copy", "--place-mode", help="copy|symlink|hardlink"),
    checkm_bin: str = typer.Option("checkm", "--checkm-bin"),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Run QC on MAGPIE taxonomy outputs."""
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


@app.command("barrnap")
def cmd_barrnap(
    qc_dir: Path = typer.Option(
        ...,
        "--qc-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE QC output directory (e.g. out/05_qc).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for Barrnap artefacts (e.g. out/06_barrnap).",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="Threads for Barrnap."),
    barrnap_bin: str = typer.Option("barrnap", "--barrnap-bin"),
    barrnap_reject: float = typer.Option(0.8, "--barrnap-reject", min=0.0, max=1.0),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Predict 16S genes from QC-filtered MAG bins using Barrnap."""
    barrnap_step(
        qc_dir=qc_dir,
        out=out,
        cpus=cpus,
        force=force,
        barrnap_bin=barrnap_bin,
        reject=barrnap_reject,
    )


@app.command("rrna")
def cmd_rrna(
    barrnap_dir: Path = typer.Option(
        ...,
        "--barrnap-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE Barrnap output directory (e.g. out/06_barrnap).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for 16S copy counting / clustering artefacts (e.g. out/07_rrna).",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="Threads for VSEARCH."),
    vsearch_bin: str = typer.Option("vsearch", "--vsearch-bin"),
    multi_cluster_id: float = typer.Option(0.90, "--multi-cluster-id", min=0.0, max=1.0),
    final_cluster_id: float = typer.Option(1.0, "--final-cluster-id", min=0.0, max=1.0),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Count, cluster, and select representative 16S genes from Barrnap outputs."""
    rrna_step(
        barrnap_dir=barrnap_dir,
        out=out,
        cpus=cpus,
        force=force,
        vsearch_bin=vsearch_bin,
        multi_cluster_id=multi_cluster_id,
        final_cluster_id=final_cluster_id,
    )


@app.command("align")
def cmd_align(
    rrna_dir: Path = typer.Option(
        ...,
        "--rrna-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE rrna output directory (e.g. out/07_rrna).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for aligned 16S centroids (e.g. out/08_align).",
    ),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for cmalign."),
    cmalign_bin: str = typer.Option("cmalign", "--cmalign-bin", help="cmalign executable name or path."),
    esl_reformat_bin: str = typer.Option("esl-reformat", "--esl-reformat-bin", help="esl-reformat executable name or path."),
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
    mxsize_archaea: int = typer.Option(4096, "--mxsize-archaea", min=1),
    mxsize_bacteria: int = typer.Option(8192, "--mxsize-bacteria", min=1),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Align domain-specific 16S centroid FASTAs using cmalign and convert to aligned FASTA."""
    align_step(
        rrna_dir=rrna_dir,
        out=out,
        cpus=cpus,
        force=force,
        cmalign_bin=cmalign_bin,
        esl_reformat_bin=esl_reformat_bin,
        ssu_models_dir=ssu_models_dir,
        mxsize_archaea=mxsize_archaea,
        mxsize_bacteria=mxsize_bacteria,
    )


@app.command("choose-best")
def cmd_choose_best(
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
    qc_dir: Path = typer.Option(
        ...,
        "--qc-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE QC output directory (e.g. out/05_qc).",
    ),
    rrna_dir: Path = typer.Option(
        ...,
        "--rrna-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE rrna output directory (e.g. out/07_rrna).",
    ),
    align_dir: Path = typer.Option(
        ...,
        "--align-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE align output directory (e.g. out/08_align).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for choose-best artefacts (e.g. out/09_choose_best).",
    ),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Choose the best genome per 16S cluster for archaea and bacteria."""
    choose_best_step(
        prep_dir=prep_dir,
        qc_dir=qc_dir,
        rrna_dir=rrna_dir,
        align_dir=align_dir,
        out=out,
        force=force,
    )


@app.command("raxml-check")
def cmd_raxml_check(
    choose_best_dir: Path = typer.Option(
        ...,
        "--choose-best-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE choose-best output directory (e.g. out/09_choose_best).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for RAxML-check artefacts (e.g. out/10_raxml_check).",
    ),
    raxml_ng_bin: str = typer.Option("raxml-ng", "--raxml-ng-bin", help="raxml-ng executable name or path."),
    threads: int = typer.Option(8, "--threads", min=1, help="Threads for raxml-ng."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Run raxml-ng --check on best aligned 16S FASTAs and ensure reduced PHYLIP outputs exist."""
    raxml_check_step(
        choose_best_dir=choose_best_dir,
        out=out,
        raxml_ng_bin=raxml_ng_bin,
        threads=threads,
        force=force,
    )


@app.command("iqtree")
def cmd_iqtree(
    raxml_check_dir: Path = typer.Option(
        ...,
        "--raxml-check-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE raxml-check output directory (e.g. out/10_raxml_check).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for IQ-TREE artefacts (e.g. out/11_iqtree).",
    ),
    iqtree_bin: str = typer.Option("iqtree", "--iqtree-bin", help="IQ-TREE executable name or path."),
    threads: int = typer.Option(8, "--threads", min=1, help="Threads for IQ-TREE."),
    bootstrap: int = typer.Option(1000, "--bootstrap", min=0, help="Ultrafast bootstrap replicates for IQ-TREE."),
    seed: int = typer.Option(12345, "--seed", help="Random seed for IQ-TREE."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Build archaeal and bacterial 16S trees from reduced PHYLIP alignments using IQ-TREE."""
    iqtree_step(
        raxml_check_dir=raxml_check_dir,
        out=out,
        iqtree_bin=iqtree_bin,
        threads=threads,
        bootstrap=bootstrap,
        seed=seed,
        force=force,
    )


@app.command("raxml-evaluate")
def cmd_raxml_evaluate(
    raxml_check_dir: Path = typer.Option(
        ...,
        "--raxml-check-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE raxml-check output directory (e.g. out/10_raxml_check).",
    ),
    iqtree_dir: Path = typer.Option(
        ...,
        "--iqtree-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE IQ-TREE output directory (e.g. out/11_iqtree).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for RAxML-evaluate artefacts (e.g. out/12_raxml_evaluate).",
    ),
    raxml_ng_bin: str = typer.Option("raxml-ng", "--raxml-ng-bin", help="raxml-ng executable name or path."),
    threads: int = typer.Option(8, "--threads", min=1, help="Threads for raxml-ng."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Evaluate IQ-TREE starting trees with RAxML-NG using reduced PHYLIP alignments."""
    raxml_evaluate_step(
        raxml_check_dir=raxml_check_dir,
        iqtree_dir=iqtree_dir,
        out=out,
        raxml_ng_bin=raxml_ng_bin,
        threads=threads,
        force=force,
    )


@app.command("hmm-prep")
def cmd_hmm_prep(
    raxml_check_dir: Path = typer.Option(
        ...,
        "--raxml-check-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE raxml-check output directory (e.g. out/10_raxml_check).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for HMM-prep artefacts (e.g. out/13_hmm_prep).",
    ),
    esl_reformat_bin: str = typer.Option("esl-reformat", "--esl-reformat-bin", help="esl-reformat executable name or path."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Convert reduced PHYLIP alignments to FASTA, DNA FASTA, and Stockholm for HMM preparation."""
    hmm_prep_step(
        raxml_check_dir=raxml_check_dir,
        out=out,
        esl_reformat_bin=esl_reformat_bin,
        force=force,
    )


@app.command("hmm-build")
def cmd_hmm_build(
    hmm_prep_dir: Path = typer.Option(
        ...,
        "--hmm-prep-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE HMM-prep output directory (e.g. out/13_hmm_prep).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for HMM-build artefacts (e.g. out/14_hmm_build).",
    ),
    hmmbuild_bin: str = typer.Option("hmmbuild", "--hmmbuild-bin", help="hmmbuild executable name or path."),
    cpus: int = typer.Option(8, "--cpus", min=1, help="CPUs for hmmbuild."),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Build archaeal and bacterial 16S HMMs from Stockholm alignments."""
    hmm_build_step(
        hmm_prep_dir=hmm_prep_dir,
        out=out,
        hmmbuild_bin=hmmbuild_bin,
        cpus=cpus,
        force=force,
    )


@app.command("package-ref")
def cmd_package_ref(
    rrna_dir: Path = typer.Option(
        ...,
        "--rrna-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE rrna output directory (e.g. out/07_rrna).",
    ),
    iqtree_dir: Path = typer.Option(
        ...,
        "--iqtree-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE IQ-TREE output directory (e.g. out/11_iqtree).",
    ),
    raxml_evaluate_dir: Path = typer.Option(
        ...,
        "--raxml-evaluate-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE RAxML-evaluate output directory (e.g. out/12_raxml_evaluate).",
    ),
    hmm_prep_dir: Path = typer.Option(
        ...,
        "--hmm-prep-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE HMM-prep output directory (e.g. out/13_hmm_prep).",
    ),
    hmm_build_dir: Path = typer.Option(
        ...,
        "--hmm-build-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        resolve_path=True,
        help="MAGPIE HMM-build output directory (e.g. out/14_hmm_build).",
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
        help="Output directory for packaged PICRUSt2-style reference files (e.g. out/15_package_ref).",
    ),
    force: bool = typer.Option(False, "--force"),
) -> None:
    """Package PICRUSt2-style reference folders and filtered 16S copy-count tables."""
    package_ref_step(
        rrna_dir=rrna_dir,
        iqtree_dir=iqtree_dir,
        raxml_evaluate_dir=raxml_evaluate_dir,
        hmm_prep_dir=hmm_prep_dir,
        hmm_build_dir=hmm_build_dir,
        out=out,
        force=force,
    )


if __name__ == "__main__":
    app()