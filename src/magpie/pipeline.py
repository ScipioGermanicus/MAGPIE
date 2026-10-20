"""
MAGPIE workflow orchestration: step registry, checkpoints and run configuration.

Each step writes into its own numbered directory under the output directory and marks
completion with a `.magpie_done` file. Re-running `magpie run` skips completed steps, so an
interrupted run resumes where it stopped. Steps can also be run on their own (`--only`) or as a
range (`--from` / `--to`) against an existing output directory; re-running a step marks all later
steps as outdated, so they are recomputed on the next run.
"""
from __future__ import annotations

import json
import logging
import shutil
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

from . import __version__
from .util.io import read_id_list, write_json
from .util.logging import add_file_log
from .util.shell import set_command_log

_LOGGER = logging.getLogger("magpie")

DOMAINS = ("bacteria", "archaea")
DONE_MARKER = ".magpie_done"


class MagpieError(RuntimeError):
    """An error whose message is meant for the user (printed without a traceback)."""


@dataclass
class Config:
    out: Path
    mags: Optional[Path] = None
    cpus: int = 8
    seed: int = 100
    force: bool = False

    # 02 prep
    rename_map: Optional[Path] = None
    sequential_ids: bool = False
    out_prefix: str = "MAG"
    pad: int = 4

    # 03 gtdbtk
    gtdb_classify: Optional[Path] = None
    gtdbtk_bin: str = "gtdbtk"
    pplacer_cpus: Optional[int] = None

    # 05 qc
    checkm_qa: Optional[Path] = None
    checkm_qa_bacteria: Optional[Path] = None
    checkm_qa_archaea: Optional[Path] = None
    checkm_bin: str = "checkm"
    completeness_min: float = 90.0
    contamination_max: float = 10.0

    # 06-07 16S genes
    barrnap_bin: str = "barrnap"
    barrnap_reject: float = 0.8
    vsearch_bin: str = "vsearch"
    multi_cluster_id: float = 0.90
    final_cluster_id: float = 1.0
    multi_copy_pick: str = "random"  # random (wiki) | longest (previous MAGPIE scripts)

    # 08 align
    ssu_align_bin: str = "ssu-align"
    esl_reformat_bin: str = "esl-reformat"

    # 10-14 trees, models, HMMs
    raxml_ng_bin: str = "raxml-ng"
    iqtree_bin: str = "iqtree"
    iqtree_bootstrap: int = 1000
    iqtree_seed: int = 12345
    hmmbuild_bin: str = "hmmbuild"
    min_genomes: int = 4

    # 16 eggnog, 17-18 trait tables
    eggnog_existing_dir: Optional[Path] = None
    eggnog_data_dir: Optional[Path] = None
    emapper_bin: str = "emapper.py"
    eggnog_jobs: int = 1
    keep_constant_traits: bool = False
    max_16s_copies: int = 10

    def to_dict(self) -> dict:
        return {k: (str(v) if isinstance(v, Path) else v) for k, v in asdict(self).items()}

    @property
    def runs_emapper(self) -> bool:
        """eggNOG-mapper is run for genomes without an existing annotation only if data is given."""
        return self.eggnog_data_dir is not None


@dataclass
class Context:
    cfg: Config
    selected: List[str] = field(default_factory=lambda: list(STEP_NAMES))

    def d(self, step: str) -> Path:
        """Output directory of a step."""
        return self.cfg.out / STEP_DIRS[step]

    @property
    def db_dir(self) -> Path:
        return self.cfg.out / "picrust2_db"

    def domain_ids(self, step: str, domain: str) -> List[str]:
        return read_id_list(self.d(step) / f"{domain}.txt")


@dataclass(frozen=True)
class Step:
    name: str
    func: Callable[[Context], None]
    help: str
    # Keep the step directory between runs so that expensive per-genome work can resume.
    resumable: bool = False


STEP_NAMES: List[str] = [
    "validate", "prep", "gtdbtk", "taxonomy", "qc", "barrnap", "rrna", "align", "choose_best",
    "raxml_check", "iqtree", "raxml_evaluate", "hmm_prep", "hmm_build", "package_ref", "eggnog",
    "ko_table", "ec_table", "picrust2_db",
]
STEP_DIRS: Dict[str, str] = {name: f"{i:02d}_{name}" for i, name in enumerate(STEP_NAMES, start=1)}


def steps() -> List[Step]:
    from . import stages

    s = stages.STEPS
    assert [x.name for x in s] == STEP_NAMES
    return s


def _select(only: Optional[List[str]], start: Optional[str], stop: Optional[str]) -> List[str]:
    for name in (only or []) + [x for x in (start, stop) if x]:
        if name not in STEP_NAMES:
            raise MagpieError(f"Unknown step '{name}'. Steps: {', '.join(STEP_NAMES)}")
    if only:
        return [s for s in STEP_NAMES if s in only]
    i = STEP_NAMES.index(start) if start else 0
    j = STEP_NAMES.index(stop) if stop else len(STEP_NAMES) - 1
    if i > j:
        raise MagpieError(f"--from {start} comes after --to {stop}")
    return STEP_NAMES[i: j + 1]


def _save_config(cfg: Config) -> None:
    """Remember the options of the last run, so later partial runs can reuse them."""
    write_json(cfg.out / "magpie_config.json", {"magpie_version": __version__, **cfg.to_dict()})


def load_saved_config(out: Path) -> dict:
    fp = out / "magpie_config.json"
    if not fp.exists():
        return {}
    data = json.loads(fp.read_text(encoding="utf-8"))
    data.pop("magpie_version", None)
    return data


def run_pipeline(cfg: Config, only: Optional[List[str]] = None, start: Optional[str] = None,
                 stop: Optional[str] = None) -> None:
    registry = {s.name: s for s in steps()}
    selected = _select(only, start, stop)
    explicit = bool(only or start)

    cfg.out.mkdir(parents=True, exist_ok=True)
    add_file_log(cfg.out / "magpie.log")
    set_command_log(cfg.out / "commands.log")
    _save_config(cfg)
    ctx = Context(cfg, selected)
    t_start = time.time()

    for name in selected:
        step = registry[name]
        sdir = ctx.d(name)
        idx = STEP_NAMES.index(name) + 1
        done = (sdir / DONE_MARKER).exists()
        # Explicitly requested steps (--only/--from) and everything after a recomputed step rerun.
        if done and not cfg.force and not explicit:
            _LOGGER.info("[%02d] %-15s already done, skipping (%s)", idx, name, sdir)
            continue

        missing_prereq = [p for p in STEP_NAMES[: idx - 1]
                          if p not in selected and not (ctx.d(p) / DONE_MARKER).exists()]
        if missing_prereq:
            raise MagpieError(f"Step '{name}' needs the output of earlier steps that have not been run in "
                              f"{cfg.out}: {', '.join(missing_prereq)}")

        if sdir.exists() and not step.resumable:
            shutil.rmtree(sdir)
        sdir.mkdir(parents=True, exist_ok=True)
        (sdir / DONE_MARKER).unlink(missing_ok=True)
        # Outputs of later steps were derived from the old version of this step.
        for later in STEP_NAMES[idx:]:
            (ctx.d(later) / DONE_MARKER).unlink(missing_ok=True)
        explicit = True  # everything after a (re)computed step is recomputed too

        _LOGGER.info("[%02d] %-15s %s", idx, name, step.help)
        t0 = time.time()
        step.func(ctx)
        (sdir / DONE_MARKER).write_text(
            json.dumps({"step": name, "version": __version__, "seconds": round(time.time() - t0, 1)}) + "\n",
            encoding="utf-8")
        _LOGGER.info("[%02d] %-15s done in %.1fs", idx, name, time.time() - t0)

    if selected[-1] == STEP_NAMES[-1]:
        _LOGGER.info("MAGPIE finished in %.1f min. PICRUSt2 database: %s", (time.time() - t_start) / 60, ctx.db_dir)


def status(out: Path) -> List[tuple]:
    rows = []
    for name in STEP_NAMES:
        d = out / STEP_DIRS[name]
        marker = d / DONE_MARKER
        state = "done" if marker.exists() else ("incomplete" if d.exists() else "-")
        rows.append((STEP_DIRS[name], name, state))
    return rows


