"""Wire each workflow step to its input/output directories and options."""
from __future__ import annotations

from .pipeline import Context, MagpieError, Step
from .steps.align import align_step
from .steps.barrnap import barrnap_step
from .steps.choose_best import choose_best_step
from .steps.ec_table import ec_table_step
from .steps.eggnog import eggnog_step, find_existing_annotations
from .steps.gtdbtk import gtdbtk_step
from .steps.hmm_build import hmm_build_step
from .steps.hmm_prep import hmm_prep_step
from .steps.iqtree import iqtree_step
from .steps.ko_table import ko_table_step
from .steps.package_ref import package_ref_step
from .steps.picrust2_db import picrust2_db_step
from .steps.prep import prep_step
from .steps.qc import qc_step
from .steps.raxml_check import raxml_check_step
from .steps.raxml_evaluate import raxml_evaluate_step
from .steps.rrna import rrna_step
from .steps.taxonomy import taxonomy_step
from .steps.validate import validate_step
from .util.ids import IdLookup


def _validate(ctx: Context) -> None:
    validate_step(cfg=ctx.cfg, out=ctx.d("validate"), steps=ctx.selected)


def _prep(ctx: Context) -> None:
    c = ctx.cfg
    if c.mags is None:
        raise MagpieError("--mags is required for the prep step")
    prep_step(mags=c.mags, out=ctx.d("prep"), rename_map=c.rename_map, sequential_ids=c.sequential_ids,
              out_prefix=c.out_prefix, pad=c.pad)


def _gtdbtk(ctx: Context) -> None:
    c = ctx.cfg
    gtdbtk_step(prep_dir=ctx.d("prep"), out=ctx.d("gtdbtk"), gtdb_classify=c.gtdb_classify, gtdbtk_bin=c.gtdbtk_bin,
                cpus=c.cpus, pplacer_cpus=c.pplacer_cpus)


def _taxonomy(ctx: Context) -> None:
    taxonomy_step(prep_dir=ctx.d("prep"), gtdbtk_dir=ctx.d("gtdbtk"), out=ctx.d("taxonomy"))


def _qc(ctx: Context) -> None:
    c = ctx.cfg
    tables = [p for p in (c.checkm_qa, c.checkm_qa_bacteria, c.checkm_qa_archaea) if p is not None]
    annotated = None
    if c.eggnog_existing_dir is not None and not c.runs_emapper:
        # Annotations are only reused: genomes without one cannot get trait values, so exclude them now.
        annotated = set(find_existing_annotations(c.eggnog_existing_dir, IdLookup(ctx.d("prep") / "id_map.tsv")))
    qc_step(prep_dir=ctx.d("prep"), tax_dir=ctx.d("taxonomy"), out=ctx.d("qc"), checkm_tables=tables,
            checkm_bin=c.checkm_bin, cpus=c.cpus, completeness_min=c.completeness_min,
            contamination_max=c.contamination_max, annotated=annotated)


def _barrnap(ctx: Context) -> None:
    c = ctx.cfg
    barrnap_step(qc_dir=ctx.d("qc"), out=ctx.d("barrnap"), barrnap_bin=c.barrnap_bin, reject=c.barrnap_reject,
                 cpus=c.cpus)


def _rrna(ctx: Context) -> None:
    c = ctx.cfg
    rrna_step(barrnap_dir=ctx.d("barrnap"), out=ctx.d("rrna"), vsearch_bin=c.vsearch_bin, cpus=c.cpus,
              multi_cluster_id=c.multi_cluster_id, final_cluster_id=c.final_cluster_id,
              multi_copy_pick=c.multi_copy_pick, seed=c.seed)


def _align(ctx: Context) -> None:
    c = ctx.cfg
    align_step(rrna_dir=ctx.d("rrna"), out=ctx.d("align"), ssu_align_bin=c.ssu_align_bin,
               esl_reformat_bin=c.esl_reformat_bin)


def _choose_best(ctx: Context) -> None:
    choose_best_step(qc_dir=ctx.d("qc"), rrna_dir=ctx.d("rrna"), align_dir=ctx.d("align"), out=ctx.d("choose_best"),
                     seed=ctx.cfg.seed)


def _raxml_check(ctx: Context) -> None:
    c = ctx.cfg
    raxml_check_step(choose_best_dir=ctx.d("choose_best"), out=ctx.d("raxml_check"), raxml_ng_bin=c.raxml_ng_bin,
                     threads=c.cpus, min_genomes=c.min_genomes)


def _iqtree(ctx: Context) -> None:
    c = ctx.cfg
    iqtree_step(raxml_check_dir=ctx.d("raxml_check"), out=ctx.d("iqtree"), iqtree_bin=c.iqtree_bin, threads=c.cpus,
                bootstrap=c.iqtree_bootstrap, seed=c.iqtree_seed)


def _raxml_evaluate(ctx: Context) -> None:
    c = ctx.cfg
    raxml_evaluate_step(raxml_check_dir=ctx.d("raxml_check"), iqtree_dir=ctx.d("iqtree"), out=ctx.d("raxml_evaluate"),
                        raxml_ng_bin=c.raxml_ng_bin, threads=c.cpus)


def _hmm_prep(ctx: Context) -> None:
    hmm_prep_step(raxml_check_dir=ctx.d("raxml_check"), out=ctx.d("hmm_prep"), esl_reformat_bin=ctx.cfg.esl_reformat_bin)


def _hmm_build(ctx: Context) -> None:
    hmm_build_step(hmm_prep_dir=ctx.d("hmm_prep"), out=ctx.d("hmm_build"), hmmbuild_bin=ctx.cfg.hmmbuild_bin,
                   cpus=ctx.cfg.cpus)


def _package_ref(ctx: Context) -> None:
    package_ref_step(rrna_dir=ctx.d("rrna"), iqtree_dir=ctx.d("iqtree"), raxml_evaluate_dir=ctx.d("raxml_evaluate"),
                     hmm_prep_dir=ctx.d("hmm_prep"), hmm_build_dir=ctx.d("hmm_build"), out=ctx.d("package_ref"),
                     max_16s_copies=ctx.cfg.max_16s_copies)


def _eggnog(ctx: Context) -> None:
    c = ctx.cfg
    eggnog_step(prep_dir=ctx.d("prep"), package_ref_dir=ctx.d("package_ref"), out=ctx.d("eggnog"),
                eggnog_existing_dir=c.eggnog_existing_dir, eggnog_data_dir=c.eggnog_data_dir,
                emapper_bin=c.emapper_bin, cpus=c.cpus, jobs=c.eggnog_jobs)


def _ko_table(ctx: Context) -> None:
    ko_table_step(package_ref_dir=ctx.d("package_ref"), eggnog_dir=ctx.d("eggnog"), out=ctx.d("ko_table"),
                  keep_constant=ctx.cfg.keep_constant_traits)


def _ec_table(ctx: Context) -> None:
    ec_table_step(package_ref_dir=ctx.d("package_ref"), eggnog_dir=ctx.d("eggnog"), out=ctx.d("ec_table"),
                  keep_constant=ctx.cfg.keep_constant_traits)


def _picrust2_db(ctx: Context) -> None:
    picrust2_db_step(package_ref_dir=ctx.d("package_ref"), ko_dir=ctx.d("ko_table"), ec_dir=ctx.d("ec_table"),
                     db=ctx.db_dir)


STEPS = [
    Step("validate", _validate, "check inputs and external tools"),
    Step("prep", _prep, "format MAGs and assign genome IDs"),
    Step("gtdbtk", _gtdbtk, "GTDB-Tk classification (run or reuse)"),
    Step("taxonomy", _taxonomy, "split genomes into bacteria/archaea"),
    Step("qc", _qc, "CheckM completeness/contamination filter"),
    Step("barrnap", _barrnap, "predict 16S rRNA genes"),
    Step("rrna", _rrna, "count 16S copies, one 16S per genome, cluster at 100%"),
    Step("align", _align, "align 16S centroids with ssu-align"),
    Step("choose_best", _choose_best, "choose the best genome per 16S cluster"),
    Step("raxml_check", _raxml_check, "raxml-ng --check"),
    Step("iqtree", _iqtree, "IQ-TREE 16S reference tree"),
    Step("raxml_evaluate", _raxml_evaluate, "raxml-ng --evaluate (model parameters)"),
    Step("hmm_prep", _hmm_prep, "reformat alignment for hmmbuild"),
    Step("hmm_build", _hmm_build, "build 16S HMMs"),
    Step("package_ref", _package_ref, "assemble bac_ref/arc_ref and 16S copy tables"),
    Step("eggnog", _eggnog, "eggNOG annotations (reuse or run)", resumable=True),
    Step("ko_table", _ko_table, "KO trait tables"),
    Step("ec_table", _ec_table, "EC trait tables"),
    Step("picrust2_db", _picrust2_db, "assemble and check the final PICRUSt2 database"),
]
