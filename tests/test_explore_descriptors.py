"""The pre-registered descriptor exploration (docs/EXPLORATION_PLAN.md)."""

from __future__ import annotations

import itertools
import json

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.benchmark import protocol
from scripts import benchmark
from scripts import explore_descriptors as ex
from tests.test_benchmark_script import _inputs


@pytest.fixture(scope="module")
def table(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("explore")
    pockets, entries = _inputs(tmp)
    frame = pd.read_csv(pockets)
    frame["n_basic_nitrogens"] = frame["n_basic_nitrogens"] + 3.0 * frame["label"]  # a real signal, prior +1
    frame.to_csv(pockets, index=False)
    out = tmp / "table.csv.gz"
    assert benchmark.main(["prepare", "--pockets-csv", str(pockets), "--entry-csv", str(entries),
                           "--output", str(out), "--summary-json", str(tmp / "s.json")]) == 0
    return out


def test_every_descriptor_has_a_declared_direction():
    assert set(ex.DIRECTIONS) == set(protocol.BENCHMARK_FEATURES)
    assert set(ex.DIRECTIONS.values()) <= {-1, 0, 1}
    for parts in ex.COMPOSITES.values():
        assert set(parts) <= set(ex.DIRECTIONS)


def test_missing_values_rank_least_binding_like():
    s = ex.signed(pd.Series([1.0, np.nan, 3.0]), -1)
    assert np.isfinite(s).all() and s.argmin() == 1 and s[0] > s[2]


def _brute_force(labels, scores, k):
    """Average over every tie-breaking order of P(positive in top k)."""
    labels, scores = np.asarray(labels), np.asarray(scores)
    hits, total = 0, 0
    for perm in itertools.permutations(range(len(labels))):
        order = sorted(perm, key=lambda i: -scores[i])  # stable: ties keep perm order
        hits += int(labels[order[:k]].any())
        total += 1
    return hits / total


@pytest.mark.parametrize("k", [1, 2, 3])
@pytest.mark.parametrize("seed", range(6))
def test_hit_probability_with_ties_matches_enumeration(k, seed):
    rng = np.random.default_rng(seed)
    n = 6
    labels = (rng.random(n) < 0.3).astype(int)
    labels[rng.integers(n)] = 1
    scores = rng.integers(0, 3, size=n).astype(float)
    table = ex.recovery_table(labels, scores, np.array(["S"] * n), np.array(["G"] * n), k)
    assert table["hit"].iloc[0] == pytest.approx(_brute_force(labels, scores, k))
    # A constant score gives exactly the chance rate.
    flat = ex.recovery_table(labels, np.zeros(n), np.array(["S"] * n), np.array(["G"] * n), k)
    assert flat["hit"].iloc[0] == pytest.approx(flat["chance"].iloc[0])


def test_group_mean_gives_groups_equal_weight():
    values = np.array([1.0, 1.0, 1.0, 0.0])
    groups = np.array(["A", "A", "A", "B"])
    assert ex.group_mean_statistic(values, groups)(np.ones(4)) == pytest.approx(0.5)


def test_benjamini_hochberg():
    adjusted = ex.benjamini_hochberg({"a": 0.01, "b": 0.02, "c": 0.03, "d": 0.5})
    assert adjusted == pytest.approx({"a": 0.04, "b": 0.04, "c": 0.04, "d": 0.5})


def test_explore_never_reads_the_holdout(table, tmp_path):
    frame = pd.read_csv(table)
    assert frame["holdout"].any()
    flipped = frame.copy()
    hold = flipped["holdout"]
    for task in ("ip_site", "cryptic_ip_site"):
        col = f"label_{task}"
        flipped.loc[hold & (flipped[col] >= 0), col] = 1 - flipped.loc[hold & (flipped[col] >= 0), col]
    flipped.loc[hold, "enclosure"] = 1e6
    other = tmp_path / "flipped.csv.gz"
    flipped.to_csv(other, index=False)
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    common = ["--n-bootstrap", "30", "--out-dir", str(tmp_path / "out")]
    assert ex.main(["explore", "--table", str(table), "--ledger", str(a), *common]) == 0
    assert ex.main(["explore", "--table", str(other), "--ledger", str(b), *common]) == 0
    ra, rb = ex.read_ledger(a)[0]["tasks"], ex.read_ledger(b)[0]["tasks"]
    assert json.dumps(ra, sort_keys=True) == json.dumps(rb, sort_keys=True)


def test_selection_rule_and_ledger_discipline(table, tmp_path):
    ledger = tmp_path / "ledger.jsonl"
    common = ["--table", str(table), "--ledger", str(ledger), "--n-bootstrap", "200", "--out-dir", str(tmp_path)]
    with pytest.raises(ex.LedgerError):
        ex.main(["confirm", *common])  # nothing explored yet
    assert ex.main(["explore", *common]) == 0
    assert ex.main(["explore", *common]) == 0  # the same table again: re-rendered, not recomputed
    assert len(ex.read_ledger(ledger)) == 1
    record = ex.read_ledger(ledger)[0]
    ip = record["tasks"]["ip_site"]
    sel = ip["selection"]
    assert sel["family_size"] == 2 * (len(ex.DIRECTIONS) + len(ex.COMPOSITES))
    assert "n_basic_nitrogens" in sel["eligible"]
    assert len(sel["confirmatory_set"]) <= ex.MAX_CONFIRMATORY
    for name in sel["eligible"]:
        assert name not in ex.NOT_BLIND and ex.prior(name) != 0
    assert set(sel["vs_rule_score"]) == set(sel["eligible"])
    assert ip["results"]["n_basic_nitrogens"]["verdict"] == "matches prior"
    assert record["tasks"]["cryptic_ip_site"]["role"] == "descriptive"
    assert "p_value" not in json.dumps(record["tasks"]["cryptic_ip_site"])
    assert (tmp_path / "EXPLORATION.md").read_text().startswith("## Descriptor exploration")

    assert ex.main(["confirm", *common]) == 0
    confirm = ex.read_ledger(ledger)[-1]["confirm"]
    assert confirm["confirmatory_set"] == sel["confirmatory_set"]
    for r in confirm.get("results", {}).values():
        assert r["decision"] in ("confirmed", "not confirmed (underpowered)")
    with pytest.raises(ex.LedgerError):
        ex.main(["confirm", *common])  # once only
