"""Operating characteristics of the group-bootstrap decisions used by the docking studies.

The fast tests pin the behaviour that the decisions depend on. The ``slow`` tests
are the calibration runs quoted in docs/METHODS.md (run with ``-m slow``):

- With 30 or more homology groups, the one-sided false-positive rate of "interval
  above 0.5" for a group-bootstrap ROC-AUC is near its nominal 2.5 %.
- With 5–9 groups it is about three times nominal. The plans' 5-group minimum
  keeps an interval from being called evidence below 5 groups, but it does not
  make one with 5–9 groups exact.
"""

from __future__ import annotations

import numpy as np
import pytest

from cryptic_ip.docking import stats
from cryptic_ip.rescoring import crossfit as cf


def _null_auc_rate(n_groups: int, n_sims: int, n_bootstrap: int, seed0: int) -> float:
    hits = 0
    for s in range(n_sims):
        rng = np.random.default_rng(seed0 + s)
        n = max(42, 2 * n_groups)
        g = np.array([f"G{i % n_groups}" for i in range(n)])
        family = rng.normal(0, 0.7, n_groups)[np.arange(n) % n_groups]  # correlated within a family
        a, b = family + rng.normal(-7, 1, n), family + rng.normal(-7, 1, n)
        est = stats.auc_estimate(np.r_[np.ones(n), np.zeros(n)], -np.r_[a, b], np.r_[g, g],
                                 n_bootstrap=n_bootstrap)["roc_auc"]
        hits += est["low"] > 0.5
    return hits / n_sims


def _rerank_runs(seed: int, n_groups: int = 30, effect: float = 0.0):
    """Pose lists with a near-native pose near the top in 60 % of runs; E_el marks it with ``effect``."""
    rng = np.random.default_rng(seed)
    runs = []
    for g in range(n_groups):
        for c in range(2):
            for s in (1, 2, 3):
                vina, rmsd = np.sort(rng.normal(-6, 1, 20)), rng.uniform(3, 10, 20)
                k = int(rng.integers(0, 6))
                native = rng.random() < 0.6
                if native:
                    rmsd[k] = 1.0
                eel = rng.normal(0, 30, 20)
                if native and effect and rng.random() < 0.5:
                    eel[k] -= 3 * effect
                runs.append(cf.Run(f"G{g}C{c}", f"G{g}", s, vina, rmsd, eel))
    return runs


def test_rerank_null_never_improves():
    decisions = [cf.evaluate(_rerank_runs(s), n_bootstrap=100, n_permutations=20)["F1"]["decision"]
                 for s in range(6)]
    assert "improves" not in decisions


def test_rerank_detects_a_strong_partial_signal():
    decisions = [cf.evaluate(_rerank_runs(s, effect=40.0), n_bootstrap=100, n_permutations=20)["F1"]["decision"]
                 for s in range(4)]
    assert decisions.count("improves") == 4


def test_cross_fitting_is_not_optimistic_under_the_null():
    """Choosing w in-sample would inflate the gain; the out-of-fold gain must not exceed the in-sample best."""
    runs = _rerank_runs(3)
    res = cf.evaluate(runs, n_bootstrap=50, n_permutations=10)
    in_sample_best = max(res["success_by_w_in_sample"].values()) - res["success_by_w_in_sample"]["0.0"]
    assert res["F1"]["group"]["difference"]["point"] <= in_sample_best + 1e-12


@pytest.mark.slow
@pytest.mark.parametrize("n_groups", [30, 80])
def test_group_bootstrap_auc_is_calibrated_with_many_groups(n_groups):
    assert _null_auc_rate(n_groups, n_sims=300, n_bootstrap=400, seed0=1000 * n_groups) <= 0.05


@pytest.mark.slow
def test_group_bootstrap_auc_is_liberal_with_few_groups():
    """Documents the known small-sample anticonservatism (about 6-8 % against 2.5 % nominal)."""
    assert _null_auc_rate(9, n_sims=300, n_bootstrap=400, seed0=9000) > 0.03


@pytest.mark.slow
def test_rerank_power_and_size():
    null = [cf.evaluate(_rerank_runs(s, n_groups=40), n_bootstrap=200, n_permutations=40)["F1"]["decision"]
            for s in range(40)]
    assert null.count("improves") <= 2
    alt = [cf.evaluate(_rerank_runs(s, n_groups=40, effect=20.0), n_bootstrap=200,
                       n_permutations=40)["F1"]["decision"] for s in range(20)]
    assert alt.count("improves") >= 16
