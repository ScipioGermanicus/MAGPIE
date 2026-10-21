![Version](https://img.shields.io/badge/version-0.2.0-blue)

# MAG-based PICRUSt2 Integrated Extension (MAGPIE)

MAGPIE builds a custom PICRUSt2 reference database (16S reference alignment, tree, HMM, model,
16S copy numbers and KO/EC trait tables) from a set of MAGs. It follows the PICRUSt2 developers'
instructions for [creating the PICRUSt2 database](https://github.com/picrust/picrust2/wiki/Creating-the-PICRUSt2-database),
adapted to MAGs, in one command. Every step can also be run on its own, and the slow tools
(GTDB-Tk, CheckM, eggNOG-mapper) can be skipped by pointing MAGPIE to results you already have.

**Important note**: The current version of MAGPIE does not support integration of GTDB-Tk, checkM, and EGGNOG. Taxonomic- and functional annotation, as well as MAG quality files need to be computed separately. MAGPIE is instead pointed to existing files.

## Installation

```bash
conda create -n magpie -c conda-forge -c bioconda python=3.11 "barrnap=0.9" vsearch ssu-align \
    raxml-ng iqtree hmmer
# optional, only if MAGPIE should run them itself:
#   gtdbtk checkm-genome eggnog-mapper
conda activate magpie
pip install .
magpie --help
```

barrnap must be version 0.9: barrnap 1.x has no `--reject` option, which the workflow uses.

## Quick start

```bash
magpie run --mags MAGs/ --out magpie_out --sequential-ids --cpus 16 \
    --gtdb-classify gtdbtk_out/classify \
    --checkm-qa-bacteria checkm/bacteria/checkm_qa.tsv \
    --checkm-qa-archaea checkm/archaea/checkm_qa.tsv \
    --eggnog-existing-dir eggnog_annotations/
```

Here GTDB-Tk, CheckM and eggNOG-mapper results already exist. Leave out `--gtdb-classify` /
`--checkm-qa*` to let MAGPIE run GTDB-Tk / CheckM, and use `--eggnog-data-dir` to let it run
eggNOG-mapper.

The database is written to `magpie_out/picrust2_db/`. Run PICRUSt2 (v2.6+) with it:

```bash
conda activate picrust2
bash magpie_out/picrust2_db/run_picrust2.sh study_seqs.fna study_table.biom picrust2_out 8
```

or install it as the default reference of a PICRUSt2 environment (the replaced files are backed up):

```bash
magpie install-picrust2-ref --db magpie_out/picrust2_db --picrust2-env $CONDA_PREFIX --out install_backup
```

## Workflow

| Step | Directory | What it does | Original script |
|---|---|---|---|
| validate | 01_validate | check inputs and tools | |
| prep | 02_prep | copy MAGs to `<ID>_genomic.fna.gz`, `id_map.tsv` | 01 |
| gtdbtk | 03_gtdbtk | `gtdbtk classify_wf`, or reuse `--gtdb-classify` | 02 |
| taxonomy | 04_taxonomy | bacteria/archaea split from the GTDB classification | 03 |
| qc | 05_qc | CheckM `lineage_wf` (or reuse `--checkm-qa*`); keep >= 90% completeness, <= 10% contamination | 04, 04b, 05 |
| barrnap | 06_barrnap | `barrnap -q -k bac/arc --outseq ... --reject 0.8` | 06 |
| rrna | 07_rrna | 16S copies per genome; vsearch 0.9 within genomes; one 16S per genome; vsearch 1.0 across genomes | 07-10 |
| align | 08_align | `ssu-align --rfonly` + `esl-reformat afa` | 11 |
| choose_best | 09_choose_best | best genome (CheckM) per identical-16S cluster | 12, 13 |
| raxml_check | 10_raxml_check | `raxml-ng --check`, PHYLIP alignment | 14 |
| iqtree | 11_iqtree | 16S reference tree, `iqtree -m GTR+G -bb 1000` | 16 |
| raxml_evaluate | 12_raxml_evaluate | `raxml-ng --evaluate` model parameters | 17 |
| hmm_prep, hmm_build | 13_, 14_ | DNA FASTA/Stockholm, `hmmbuild` | 18 |
| package_ref | 15_package_ref | `bac_ref/`, `arc_ref/`, `.raxml_info`, 16S copy tables | 19, 20 |
| eggnog | 16_eggnog | reuse annotations, or `emapper.py --itype genome` | 22, 23 |
| ko_table, ec_table | 17_, 18_ | `ko.txt.gz`, `ec.txt.gz` trait tables | 24, 25 |
| picrust2_db | `picrust2_db/` | final database, consistency check, `run_picrust2.sh` | |

`magpie steps` lists the steps; `magpie status --out DIR` shows which are done.

### Running steps separately, resuming

* Rerunning the same `magpie run` command skips completed steps, so an interrupted run resumes.
* `--only STEP[,STEP]` runs single steps, `--from STEP` reruns from a step on, `--to STEP` stops after
  a step (e.g. `--to package_ref` builds the reference without annotation).
* Options that are not repeated are taken from the previous run in the same output directory
  (`OUT/magpie_config.json`), so `magpie run --out OUT --only iqtree` is enough.
* When a step is rerun, all later steps are marked as outdated and are recomputed by the next run.
* `--force` recomputes everything (eggNOG annotations already made are kept).

### Reusing existing results

| Instead of running | Give | Accepted files |
|---|---|---|
| GTDB-Tk | `--gtdb-classify DIR` | `gtdbtk.bac120.summary.tsv`, `gtdbtk.ar53.summary.tsv` (in DIR or DIR/classify) |
| CheckM | `--checkm-qa FILE` or `--checkm-qa-bacteria/--checkm-qa-archaea` | `checkm qa -o 2` output (fixed-width or `--tab_table`), or `genome_id, checkm_completeness, checkm_contamination` |
| eggNOG-mapper | `--eggnog-existing-dir DIR` | `<name>.emapper.annotations[.gz]` |

Genomes may be named by their original file name, file stem, or MAGPIE ID in all of these.
With `--eggnog-existing-dir` alone, MAGs without an annotation are excluded at the QC step so that
every reference genome has traits; add `--eggnog-data-dir` to annotate the missing ones instead.


## Dependencies

MAGPIE is a Python package that calls external tools. All of them are available from
[bioconda](https://bioconda.github.io/); the Installation section below creates a conda
environment with them. The "Tested" column lists the versions MAGPIE was tested with.

### Python

| Package | Tested | Notes |
|---|---|---|
| Python | 3.11 | 3.10 or newer |
| [typer](https://typer.tiangolo.com/) | 0.27 | installed automatically by `pip install .` |
| pytest | | only to run the unit tests in `tests/` |

### Always required

| Tool | Conda package | Tested | Used in step |
|---|---|---|---|
| barrnap | `barrnap=0.9` | 0.9 | barrnap. Must be 0.9: barrnap 1.x has no `--reject` option. |
| VSEARCH | `vsearch` | 2.32.0 | rrna |
| SSU-ALIGN | `ssu-align` | 0.1.1 | align |
| esl-reformat (Easel, installed with HMMER) | `hmmer` | 0.49 | align, hmm_prep |
| RAxML-NG | `raxml-ng` | 2.0.3 | raxml_check, raxml_evaluate |
| IQ-TREE | `iqtree` | 3.1.4 | iqtree. IQ-TREE 2 or 3; the command uses `-bb`, `-nt` and `-pre`. |
| hmmbuild (HMMER) | `hmmer` | 3.4 | hmm_build |

### Only needed if MAGPIE runs these tools itself

These steps can be skipped by giving existing results (see
[Reusing existing results](#reusing-existing-results)).

| Tool | Conda package | Tested | Needed when | Reference data |
|---|---|---|---|---|
| GTDB-Tk | `gtdbtk` | 2.7.2 | `--gtdb-classify` is not given | GTDB-Tk reference package (about 140 GB, set `GTDBTK_DATA_PATH`); classify_wf needs about 128 GB RAM |
| CheckM | `checkm-genome` | 1.2.5 | `--checkm-qa*` is not given | CheckM data (about 1.4 GB, set with `checkm data setRoot`); lineage_wf needs about 40 GB RAM |
| eggNOG-mapper, with Prodigal and DIAMOND | `eggnog-mapper` | 2.1.13 | `--eggnog-data-dir` is given | eggNOG data (about 50 GB, `download_eggnog_data.py`) |

Prodigal (2.6.3) and DIAMOND (2.2.8) are installed with the `eggnog-mapper` package.

### Using the database

PICRUSt2 2.6 or newer (tested with 2.6.3) is needed to run `run_picrust2.sh` and for
`magpie install-picrust2-ref`. Install it in its own conda environment
(`conda create -n picrust2 -c conda-forge -c bioconda picrust2`).

The first step of `magpie run` (validate) checks that every tool needed by your run is on your
`PATH` and stops with a list of missing tools if not.


### Notes on the method

* The reference tree is an IQ-TREE tree of the 16S alignment (the PICRUSt2 developers prune the
  GTDB genome tree, which does not contain new MAGs).
* Where a genome has divergent 16S copies, one is chosen at random (`--multi-copy-pick random`,
  as in the PICRUSt2 instructions; `longest` reproduces the earlier MAGPIE scripts). Random
  choices are seeded (`--seed`).
* The best genome of a 16S cluster is chosen as in the PICRUSt2 instructions: highest
  completeness; if all are equal, lowest contamination; if those are equal too, at random.
* Gene families with the same count in every reference genome are removed from the trait tables
  (listed in `*_constant_removed.tsv`), because PICRUSt2's default HSP method (`mp`) fails on them.
  Use `--keep-constant-traits` if you run PICRUSt2 with `-m emp_prob` or `-m pic`.
* A domain needs at least 4 reference genomes. If e.g. no archaeal reference can be built,
  `run_picrust2.sh` uses PICRUSt2's default archaeal reference for that domain.
