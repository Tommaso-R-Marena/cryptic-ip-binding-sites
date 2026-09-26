"""scripts/coevolution.py (docs/COEVOLUTION_PLAN.md) on synthetic proteomes."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def coev():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("coevolution", ROOT / "scripts" / "coevolution.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["coevolution"] = module
    spec.loader.exec_module(module)
    return module


def _proteome(organism, n, score, plddt=90.0, length=400, basic=0.12, start=0, seed=0):
    """One synthetic proteome; ``score`` may be a scalar or a callable of the index."""
    rng = np.random.default_rng(seed)
    scores = [score(i) if callable(score) else score + rng.normal(0, 0.01) for i in range(n)]
    return pd.DataFrame({
        "uniprot_id": [f"{organism[:3]}{start + i}" for i in range(n)],
        "organism_key": organism, "combined": scores, "rule_score": scores,
        "plddt_mean": plddt, "hull_depth": 10.0, "annotated": False, "seen": False,
        "cluster": [f"{organism[:3]}c{start + i}" for i in range(n)],
        "sequence_length": length, "basic_fraction": basic, "eligible": plddt >= coev_min(),
    })


def coev_min():
    return 70.0


def _normal(organism, n, mean, sd=1.0, plddt=90.0, length=400, basic=0.12, start=0, seed=1, rule=None):
    """A proteome whose score is continuous, so a percentile anchor cuts it proportionally."""
    rng = np.random.default_rng(seed)
    scores = rng.normal(mean, sd, n)
    return pd.DataFrame({
        "uniprot_id": [f"{organism[:3]}{start + i}" for i in range(n)],
        "organism_key": organism, "combined": scores,
        "rule_score": scores if rule is None else rng.normal(rule, sd, n),
        "plddt_mean": plddt, "hull_depth": 10.0, "annotated": False, "seen": False,
        "cluster": [f"{organism[:3]}c{start + i}" for i in range(n)],
        "sequence_length": length, "basic_fraction": basic, "eligible": plddt >= coev_min(),
    })


def _table(*frames):
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------- basics
def test_basic_fraction_counts_krh(coev):
    assert coev.basic_fraction("KRHAAAAAAA") == pytest.approx(0.3)
    assert np.isnan(coev.basic_fraction(""))


def test_build_table_refuses_a_table_without_the_screen_columns(coev):
    thin = pd.DataFrame({"uniprot_id": ["A"], "combined": [0.5]})
    with pytest.raises(ValueError, match="proteins_combined"):
        coev.build_table(thin, {}, pd.Series(dtype=float))


def test_build_table_adds_length_and_composition(coev):
    proteins = pd.DataFrame({"uniprot_id": ["A", "B"], "combined": [0.5, 0.2], "organism_key": "human",
                             "plddt_mean": [90.0, 50.0], "annotated": False, "seen": False,
                             "cluster": ["c1", "c2"]})
    out = coev.build_table(proteins, {"A": "KKKAAAAAAA"}, pd.Series({"A": 0.7}))
    row = out.set_index("uniprot_id").loc["A"]
    assert row["sequence_length"] == 10 and row["basic_fraction"] == pytest.approx(0.3)
    assert row["rule_score"] == pytest.approx(0.7) and bool(row["eligible"]) is True
    assert bool(out.set_index("uniprot_id").loc["B", "eligible"]) is False   # below the pLDDT floor


# ------------------------------------------------------------ threshold
def test_the_threshold_is_anchored_on_the_pooled_background(coev):
    table = _table(_proteome("human", 100, lambda i: i / 100.0),
                   _proteome("dictyostelium", 100, lambda i: i / 100.0, start=1000))
    t = coev.threshold_at_fpr(table, "combined", 0.10)
    part = coev.hits(table, "combined", t)
    assert 0.08 <= part["hit"].mean() <= 0.12          # about a tenth of the pooled background


def test_the_anchor_does_not_force_the_organisms_to_be_equal(coev):
    """Pooling the anchor must leave a real rate difference visible, not normalise it away."""
    table = _table(_normal("human", 400, 0.0, seed=1),
                   _normal("yeast", 400, 0.0, start=1000, seed=2),
                   _normal("dictyostelium", 400, 2.0, start=2000, seed=3))
    t = coev.threshold_at_fpr(table, "combined", 0.10)
    rates = coev.hits(table, "combined", t).groupby("organism_key")["hit"].mean()
    assert rates["dictyostelium"] > 4 * max(rates["human"], rates["yeast"])


def test_a_background_without_scores_refuses_to_set_a_threshold(coev):
    table = _table(_proteome("human", 5, 0.5))
    table["combined"] = np.nan
    with pytest.raises(ValueError, match="no pooled background"):
        coev.threshold_at_fpr(table, "combined", 0.01)


# ------------------------------------------------- the matching, which is the point
def test_a_real_difference_survives_the_matching(coev):
    """Same pLDDT, length and composition in both arms: the matching must not erase a true effect."""
    table = _table(_proteome("human", 150, 0.10),
                   _proteome("dictyostelium", 150, 0.90, start=1000))
    part = coev.hits(table, "combined", 0.5)
    r = coev.matched_difference(part, "dictyostelium", n_bootstrap=300)
    assert r["decision"] == "higher in dictyostelium"
    assert r["estimate"]["point"] == pytest.approx(1.0, abs=0.05)


def test_a_confounded_difference_is_removed_by_the_matching(coev):
    """The test the whole design exists for.

    Dictyostelium scores higher, but only because its pLDDT is higher, and score tracks
    pLDDT identically in both organisms. The raw comparison sees a difference; the matched
    one must not, because within a pLDDT bin the two organisms behave the same.
    """
    hi = _proteome("dictyostelium", 150, 0.90, plddt=95.0, start=1000)
    lo = _proteome("human", 150, 0.10, plddt=75.0)
    # the same relationship in both: high pLDDT -> high score. No organism effect at all.
    mixed = _table(hi, lo,
                   _proteome("dictyostelium", 150, 0.10, plddt=75.0, start=3000),
                   _proteome("human", 150, 0.90, plddt=95.0, start=4000))
    part = coev.hits(mixed, "combined", 0.5)
    raw = coev.raw_difference(part, "dictyostelium", n_bootstrap=300)
    matched = coev.matched_difference(part, "dictyostelium", n_bootstrap=300)
    assert abs(matched["estimate"]["point"]) < 0.02
    assert matched["decision"] == "no material difference"
    assert abs(raw["estimate"]["point"]) < 0.02    # balanced by construction here too


def test_matching_ignores_bins_only_one_organism_occupies(coev):
    """A bin with no comparator carries no information and must not contribute."""
    table = _table(_proteome("dictyostelium", 60, 0.9, plddt=99.0, start=1000),   # unmatched band
                   _proteome("dictyostelium", 60, 0.1, plddt=80.0, start=2000),
                   _proteome("human", 60, 0.1, plddt=80.0))
    part = coev.hits(table, "combined", 0.5)
    r = coev.matched_difference(part, "dictyostelium", n_bootstrap=300)
    assert r["bins"] == 1                                   # only the pLDDT-80 band is shared
    assert abs(r["estimate"]["point"]) < 0.05               # and inside it the arms agree


def test_few_clusters_are_not_evaluable(coev):
    table = _table(_proteome("dictyostelium", 2, 0.9, start=1000), _proteome("human", 2, 0.1))
    part = coev.hits(table, "combined", 0.5)
    assert coev.matched_difference(part, "dictyostelium", n_bootstrap=100)["decision"].startswith("not evaluable")


def test_no_shared_bin_is_not_evaluable_rather_than_zero(coev):
    table = _table(_proteome("dictyostelium", 40, 0.9, plddt=99.0, start=1000),
                   _proteome("human", 40, 0.9, plddt=75.0))
    part = coev.hits(table, "combined", 0.5)
    r = coev.matched_difference(part, "dictyostelium", n_bootstrap=100)
    assert r["decision"].startswith("not evaluable") and r["bins"] == 0


def test_decide_labels(coev):
    assert coev.decide({"point": 0.5, "low": 0.1, "high": 0.9, "p_value": 0.0}, 30) == "higher in dictyostelium"
    assert coev.decide({"point": -0.5, "low": -0.9, "high": -0.1, "p_value": 0.0}, 30) == "lower in dictyostelium"
    assert coev.decide({"point": 0.0, "low": -0.001, "high": 0.001, "p_value": 1.0}, 30) == "no material difference"
    assert coev.decide({"point": 0.0, "low": -0.5, "high": 0.5, "p_value": 1.0}, 30) == "inconclusive"
    assert coev.decide(None, 30).startswith("not evaluable")


# ------------------------------------------------------------- robustness
def test_disagreement_between_the_two_scores_is_flagged(coev):
    assert coev._agree({"point": 0.4}, {"point": 0.3}) is True
    assert coev._agree({"point": 0.4}, {"point": -0.3}) is False
    assert coev._agree({"point": float("nan")}, {"point": 0.3}) is None


def test_report_marks_the_primary_not_robust_when_the_rule_disagrees(coev):
    """combined says Dictyostelium is higher; the rule, which carries no training bias, says lower.

    The plan makes this void the claim rather than report the favourable arm.
    """
    r = coev.build(_table(_normal("dictyostelium", 400, 2.0, start=1000, seed=3, rule=-2.0),
                          _normal("human", 400, 0.0, seed=1, rule=0.0),
                          _normal("yeast", 400, 0.0, start=2000, seed=2, rule=0.0)),
                   n_bootstrap=200)
    assert r["C2_matched"]["decision"] == "higher in dictyostelium"
    assert r["C4_rule_only"]["decision"] == "lower in dictyostelium"
    assert r["robust"] is False


def test_report_agrees_when_both_scores_point_the_same_way(coev):
    r = coev.build(_table(_normal("dictyostelium", 400, 2.0, start=1000, seed=3, rule=2.0),
                          _normal("human", 400, 0.0, seed=1, rule=0.0),
                          _normal("yeast", 400, 0.0, start=2000, seed=2, rule=0.0)),
                   n_bootstrap=200)
    assert r["robust"] is True
    assert r["holm"]["C2"] <= 1.0 and "C4" in r["holm"]


# ----------------------------------------------------------------- outputs
def test_report_writes_its_outputs(coev, tmp_path):
    table = _table(_proteome("human", 80, lambda i: i / 80.0),
                   _proteome("yeast", 80, lambda i: i / 80.0, start=1000),
                   _proteome("dictyostelium", 80, lambda i: i / 80.0, start=2000))
    table.to_csv(tmp_path / "t.csv", index=False)
    assert coev.main(["report", "--table", str(tmp_path / "t.csv"),
                      "--out-dir", str(tmp_path / "out"), "--n-bootstrap", "100"]) == 0
    text = (tmp_path / "out" / "COEVOLUTION.md").read_text()
    assert "C2 (primary, matched)" in text and "C1 (raw, not the test)" in text
    payload = json.loads((tmp_path / "out" / "coevolution.json").read_text())
    assert payload["focus"] == "dictyostelium" and "sensitivity_fpr" in payload
