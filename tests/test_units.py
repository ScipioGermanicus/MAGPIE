import gzip
import random
from collections import Counter
from pathlib import Path

import pytest

from magpie.pipeline import MagpieError, _select
from magpie.steps.choose_best import choose_best, parse_clusters
from magpie.steps.ec_table import normalise_ec
from magpie.steps.ko_table import count_annotations, normalise_ko, write_trait_table
from magpie.steps.package_ref import parse_raxml_ng_log
from magpie.steps.picrust2_db import tree_tips
from magpie.steps.qc import parse_checkm_table
from magpie.util.ids import IdLookup
from magpie.util.io import sanitise_id, strip_fasta_suffix


def test_strip_fasta_suffix_keeps_dotted_ids():
    assert strip_fasta_suffix("BC13.bin.1.603.fa") == "BC13.bin.1.603"
    assert strip_fasta_suffix("MAG0001_genomic.fna.gz") == "MAG0001"
    assert strip_fasta_suffix("x.fasta") == "x"
    assert sanitise_id("S3 bin(5)") == "S3_bin_5_"


def test_id_lookup(tmp_path: Path):
    m = tmp_path / "id_map.tsv"
    m.write_text("original_filename\toriginal_id\tnew_id\nBC13.bin.1.143.fa\tBC13.bin.1.143\tMAG1988\n")
    lk = IdLookup(m)
    for name in ("BC13.bin.1.143", "BC13.bin.1.143.fa", "MAG1988", "MAG1988_genomic", "MAG1988_genomic.fna.gz"):
        assert lk.resolve(name) == "MAG1988"
    assert lk.resolve("other") is None


CHECKM_FIXED = """\
------------------------------------------------------------------------------
  Bin Id                   Marker lineage         # genomes   # markers   # marker sets   Completeness   Contamination   Strain heterogeneity
------------------------------------------------------------------------------
  BC13.bin.1.22074     k__Bacteria (UID2328)         3167        126            75           98.56            0.00               0.00
  BC13.bin.1.422     o__Clostridiales (UID1120)      304         247           141           95.15            1.42              42.86
------------------------------------------------------------------------------
"""
CHECKM_TAB = "\tBin Id\tMarker lineage\t# genomes\tCompleteness\tContamination\n\tMAG2107\tk__Bacteria (UID209)\t5443\t100.00\t0.00\n"
CHECKM_MIN = "genome_id\tcheckm_completeness\tcheckm_contamination\nMAG2107\t100.00\t0.00\n"


@pytest.mark.parametrize("text,name,comp,cont", [
    (CHECKM_FIXED, "BC13.bin.1.422", 95.15, 1.42),
    (CHECKM_TAB, "MAG2107", 100.0, 0.0),
    (CHECKM_MIN, "MAG2107", 100.0, 0.0),
])
def test_parse_checkm_formats(tmp_path: Path, text, name, comp, cont):
    p = tmp_path / "qa.tsv"
    p.write_text(text)
    rows = {r.name: r for r in parse_checkm_table(p)}
    assert rows[name].completeness == comp and rows[name].contamination == cont


def test_choose_best_follows_wiki():
    q = {"a": {"completeness": 97.9, "contamination": 1.75}, "b": {"completeness": 97.9, "contamination": 1.05},
         "c": {"completeness": 95.0, "contamination": 0.0}}
    rng = random.Random(1)
    # completeness differs -> first genome with the highest completeness, contamination ignored
    assert choose_best(["c", "a", "b"], q, rng) == "a"
    q2 = {"a": {"completeness": 99.0, "contamination": 2.0}, "b": {"completeness": 99.0, "contamination": 1.0}}
    assert choose_best(["a", "b"], q2, rng) == "b"  # equal completeness -> lowest contamination
    q3 = {g: {"completeness": 99.0, "contamination": 0.0} for g in "abcd"}
    picks = {choose_best(list("abcd"), q3, random.Random(s)) for s in range(50)}
    assert len(picks) > 1  # full tie -> random


def test_parse_uc(tmp_path: Path):
    uc = tmp_path / "x.uc"
    uc.write_text("S\t0\t1500\t*\t*\t*\t*\t*\tA\t*\n"
                  "H\t0\t1500\t100.0\t+\t0\t0\t1500M\tB\tA\n"
                  "S\t1\t1400\t*\t*\t*\t*\t*\tC\t*\n"
                  "H\t2\t1400\t100.0\t+\t0\t0\t1400M\tE\tD\n"
                  "C\t0\t2\t100.0\t*\t*\t*\t*\tA\t*\n")
    assert dict(parse_clusters(uc, {"A", "C"})) == {"A": ["B"]}  # D not in alignment -> ignored


def test_trait_parsing_and_constant_columns(tmp_path: Path):
    ann = tmp_path / "g.emapper.annotations"
    ann.write_text("## emapper\n#query\tEC\tKEGG_ko\n"
                   "g_1\t1.1.1.1,2.7.-.-\tko:K00001,ko:K00002\n"
                   "g_2\t-\tko:K00001\n"
                   "## 2 queries scanned\n")
    ko = count_annotations(ann, "KEGG_ko", normalise_ko)
    ec = count_annotations(ann, "EC", normalise_ec)
    assert ko == Counter({"ko:K00001": 2, "ko:K00002": 1})
    assert ec == Counter({"1.1.1.1": 1})  # incomplete EC numbers are dropped

    counts = {"g1": Counter({"ko:K1": 1, "ko:K2": 1}), "g2": Counter({"ko:K1": 1, "ko:K2": 3})}
    out = tmp_path / "ko.txt.gz"
    n, constant = write_trait_table(out, ["g1", "g2"], counts, keep_constant=False)
    assert n == 1 and constant == [("ko:K1", 1)]
    with gzip.open(out, "rt") as fh:
        assert fh.read().splitlines() == ["assembly\tko:K2", "g1\t1", "g2\t3"]
    n, constant = write_trait_table(out, ["g1", "g2"], counts, keep_constant=True)
    assert n == 2 and constant == []


RAXML_LOG = """\
RAxML-NG was called at 01-Apr-2026 13:11:42 as follows:

raxml-ng --evaluate --msa a.phy --tree t.tre --model GTR+G --prefix x --threads 4

Alignment comprises 1 partitions and 557 patterns
   Base frequencies (ML): 0.220090 0.249290 0.290844 0.239776
   Substitution rates (ML): 1.222614 4.378139 1.877652 0.999680 5.988560 1.000000
Final LogLikelihood: -282040.011309
Elapsed time: 16.044 seconds
"""


def test_parse_raxml_log(tmp_path: Path):
    p = tmp_path / "x.raxml.log"
    p.write_text(RAXML_LOG)
    r = parse_raxml_ng_log(p)
    assert r["patterns"] == "557" and r["final_ll"] == "-282040.011309" and r["elapsed"] == "16.044"
    assert r["cmd_line"].startswith("raxml-ng --evaluate")
    assert r["subs_rates"].split()[-1] == "1.000000"


def test_tree_tips(tmp_path: Path):
    t = tmp_path / "t.tre"
    t.write_text("((MAG1:0.1,MAG2:0.2)100:0.3,MAG3:0.1,(BC13.bin.1.603:0.1,MAG4:0.1)95:0.2);\n")
    assert sorted(tree_tips(t)) == ["BC13.bin.1.603", "MAG1", "MAG2", "MAG3", "MAG4"]


def test_step_selection():
    assert _select(["rrna", "align"], None, None) == ["rrna", "align"]
    assert _select(None, "ko_table", None) == ["ko_table", "ec_table", "picrust2_db"]
    assert _select(None, None, "prep") == ["validate", "prep"]
    with pytest.raises(MagpieError):
        _select(["nope"], None, None)
