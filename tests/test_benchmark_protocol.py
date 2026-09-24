"""The benchmark protocol's leak guards and statistics (docs/ANALYSIS_PLAN.md)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import average_precision_score, roc_auc_score

from cryptic_ip.analysis.features import FEATURE_NAMES
from cryptic_ip.benchmark import protocol


def _data(n_groups=40, per_group=15, signal=2.0, seed=0):
    """Pockets in groups; one descriptor carries the label, the rest are noise."""
    rng = np.random.default_rng(seed)
    groups = np.repeat([f"G{i:03d}" for i in range(n_groups)], per_group)
    y = (rng.random(len(groups)) < 0.2).astype(int)
    X = pd.DataFrame(rng.normal(size=(len(groups), len(protocol.BENCHMARK_FEATURES))),
                     columns=list(protocol.BENCHMARK_FEATURES))
    X["enclosure"] += signal * y
    rule = X["enclosure"].to_numpy() + rng.normal(scale=2.0, size=len(y))
    return X, y, groups, rule


SMALL = dict(n_outer=3, n_inner=2, n_draws=1)


class TestFeatureGuards:
    def test_b_factor_descriptors_never_reach_a_model(self):
        for name in ("plddt_mean", "plddt_min", "plddt_fraction_above_cutoff"):
            assert name in FEATURE_NAMES and name not in protocol.BENCHMARK_FEATURES
            with pytest.raises(protocol.LeakError):
                protocol.assert_allowed_features(["enclosure", name])

    def test_non_descriptor_columns_are_refused(self):
        with pytest.raises(protocol.LeakError):
            protocol.assert_allowed_features(["enclosure", "label"])

    def test_arms_differ_only_by_hull_depth(self):
        full, reduced = protocol.ARMS["full"], protocol.ARMS["no_hull_depth"]
        assert set(full) - set(reduced) == {"hull_depth"} and set(reduced) < set(full)
        assert len(full) == 39

    def test_run_cv_refuses_an_excluded_descriptor(self):
        X, y, groups, rule = _data()
        X["plddt_mean"] = 0.0
        with pytest.raises(protocol.LeakError):
            protocol.run_cv(X, y, groups, rule, **SMALL)


def test_disjointness_guard():
    protocol.assert_disjoint(["a", "b"], ["c"])
    with pytest.raises(protocol.LeakError):
        protocol.assert_disjoint(["a", "b"], ["b", "c"])


def test_candidates_are_reproducible_and_span_the_families():
    one, two = protocol.draw_candidates(4, seed=7), protocol.draw_candidates(4, seed=7)
    assert one == two
    assert {c.family for c in one} == set(protocol.FAMILIES)


class TestRunCv:
    def test_every_row_scored_once_and_signal_found(self):
        X, y, groups, rule = _data()
        out = protocol.run_cv(X, y, groups, rule, **SMALL)
        assert (out["fold"] >= 0).all() and out["score"].notna().all()
        assert roc_auc_score(y, out["score"]) > 0.8
        # Folds never share a group.
        for fold in out["fold"].unique():
            protocol.assert_disjoint(groups[out["fold"] != fold], groups[out["fold"] == fold])

    def test_permuted_labels_give_chance(self):
        X, y, groups, rule = _data(n_groups=60)
        y_perm = np.random.default_rng(1).permutation(y)
        out = protocol.run_cv(X, y_perm, groups, rule, **SMALL)
        assert abs(roc_auc_score(y_perm, out["score"]) - 0.5) < 0.12

    def test_outer_test_labels_cannot_influence_their_own_predictions(self, monkeypatch):
        """Flip every label in one outer test fold: its predictions, threshold and model must not move."""
        X, y, groups, rule = _data()
        fixed = protocol._group_splits(y, groups, 3, seed=0)
        real = protocol._group_splits

        def splits(y_, groups_, n_splits, seed):
            if len(y_) == len(y):
                return fixed
            return real(y_, groups_, n_splits, seed)

        monkeypatch.setattr(protocol, "_group_splits", splits)
        before = protocol.run_cv(X, y, groups, rule, **SMALL)
        test0 = fixed[0][1]
        y_flipped = y.copy()
        y_flipped[test0] = 1 - y_flipped[test0]
        after = protocol.run_cv(X, y_flipped, groups, rule, **SMALL)
        cols = ["score", "calibrated", "threshold", "candidate", "rule_threshold"]
        pd.testing.assert_frame_equal(before.iloc[test0][cols], after.iloc[test0][cols])

    def test_paired_arms_share_folds(self):
        X, y, groups, rule = _data()
        full = protocol.run_cv(X, y, groups, rule, fold_seed=3, **SMALL)
        reduced = protocol.run_cv(X[list(protocol.ARMS["no_hull_depth"])], y, groups, rule, fold_seed=3, **SMALL)
        assert (full["fold"] == reduced["fold"]).all()


class TestHoldout:
    def _entries(self):
        return pd.DataFrame({
            "pdb_id": ["1AAA", "2AAA", "1BBB", "1CCC", "1DDD", "2DDD"],
            "group": ["A", "A", "B", "C", "D", "D"],
            "release_date": ["2001-01-01", "2023-05-01", "2010-01-01", "2024-01-01", "2022-01-01", "2023-01-01"],
        })

    def test_latest_groups_by_first_release(self):
        held = protocol.temporal_holdout(self._entries(), group_column="group", fraction=0.5)
        # First releases: A 2001, B 2010, D 2022, C 2024 -> latest half is D and C.
        assert held == {"1CCC", "1DDD", "2DDD"}  # A's 2023 entry stays with its group

    def test_missing_date_is_an_error(self):
        entries = self._entries()
        entries.loc[0, "release_date"] = None
        with pytest.raises(ValueError):
            protocol.temporal_holdout(entries, group_column="group")

    def test_locked_model_refuses_overlapping_groups(self):
        X, y, groups, rule = _data()
        with pytest.raises(protocol.LeakError):
            protocol.run_locked(X, y, groups, X.iloc[:5], groups[:5], rule, n_inner=2, n_draws=1)


class TestRanked:
    @pytest.mark.parametrize("seed", range(5))
    def test_matches_scikit_learn_with_ties_and_weights(self, seed):
        rng = np.random.default_rng(seed)
        y = (rng.random(400) < 0.3).astype(int)
        s = np.round(rng.normal(size=400) + y, 1)  # rounding creates ties
        w = rng.integers(0, 4, size=400).astype(float)
        ranked = protocol.Ranked(y, s)
        keep = w > 0
        assert ranked.roc_auc(w) == pytest.approx(roc_auc_score(y[keep], s[keep], sample_weight=w[keep]))
        assert ranked.average_precision(w) == pytest.approx(
            average_precision_score(y[keep], s[keep], sample_weight=w[keep])
        )
        assert ranked.roc_auc(np.ones(400)) == pytest.approx(roc_auc_score(y, s))

    def test_single_class_is_nan(self):
        ranked = protocol.Ranked([1, 1, 1], [0.1, 0.2, 0.3])
        assert np.isnan(ranked.roc_auc(np.ones(3)))


class TestInference:
    def test_bootstrap_interval_and_p_value(self):
        groups = np.repeat(np.arange(50), 4)
        values = np.random.default_rng(0).normal(loc=1.0, size=200)
        est = protocol.bootstrap_statistic(lambda w: float(np.average(values, weights=w)), groups, n_bootstrap=500)
        assert est.low < 1.0 < est.high and est.p_value < 0.01
        null = protocol.bootstrap_statistic(lambda w: float(np.average(values - 1.0, weights=w)), groups, n_bootstrap=500)
        assert null.p_value > 0.05

    def test_holm(self):
        assert protocol.holm({"a": 0.01, "b": 0.04}) == pytest.approx({"a": 0.02, "b": 0.04})
        assert protocol.holm({"a": 0.03, "b": 0.02}) == pytest.approx({"b": 0.04, "a": 0.04})

    def _est(self, low, high, point=None):
        return protocol.Estimate(point if point is not None else (low + high) / 2, low, high, 0.001, 2000)

    def test_decisions(self):
        supported = {"sequence": self._est(0.01, 0.05), "strict": self._est(0.005, 0.04)}
        assert protocol.decide(supported, 0.01, None, 0) == "supported"
        # A powered holdout must agree.
        assert protocol.decide(supported, 0.01, self._est(-0.02, 0.03), 12) == "inconclusive"
        assert protocol.decide(supported, 0.01, self._est(-0.02, 0.03), 5) == "supported"  # underpowered
        # Both groupings are required.
        mixed = {"sequence": self._est(0.01, 0.05), "strict": self._est(-0.01, 0.04)}
        assert protocol.decide(mixed, 0.01, None, 0) == "inconclusive"
        assert protocol.decide(supported, 0.2, None, 0) == "inconclusive"
        refuted = {"sequence": self._est(-0.05, -0.01), "strict": self._est(-0.04, -0.002)}
        assert protocol.decide(refuted, 0.01, None, 0) == "refuted"
        null = {"sequence": self._est(-0.004, 0.006), "strict": self._est(-0.008, 0.009)}
        assert protocol.decide(null, 0.6, None, 0) == "refuted"
