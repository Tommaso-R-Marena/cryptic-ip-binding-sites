"""Tests for the machine-learning pocket classifier."""

from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.ml_classifier import (
    DEFAULT_FEATURE_COLUMNS,
    FEATURE_COLUMNS,
    CrypticSiteMLClassifier,
    MLPocketScorer,
    classification_metrics,
    compare_models,
    delong_roc_test,
    enrichment_factor,
    expected_calibration_error,
    grouped_bootstrap_ci,
    model_comparison_table,
    precision_at_k,
    select_threshold,
)

try:
    import joblib  # noqa: F401
    import sklearn  # noqa: F401

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

pytestmark = pytest.mark.skipif(
    not SKLEARN_AVAILABLE, reason="scikit-learn/joblib not installed"
)


def _build_training_data() -> Tuple[pd.DataFrame, np.ndarray]:
    """Create a small separable dataset in the legacy six-feature schema."""
    positives = [
        (18.0, 2.0, 7.5, 6, 650.0, 89.0),
        (16.5, 3.2, 6.8, 5, 540.0, 86.0),
        (17.2, 2.8, 6.9, 5, 610.0, 82.0),
        (19.1, 1.6, 8.0, 7, 700.0, 91.0),
        (15.8, 4.2, 5.9, 4, 500.0, 80.0),
    ]
    negatives = [
        (6.1, 65.0, 2.5, 2, 380.0, 88.0),
        (7.5, 55.0, 3.2, 3, 430.0, 84.0),
        (5.2, 70.0, 1.8, 1, 350.0, 76.0),
        (8.0, 45.0, 3.1, 2, 470.0, 79.0),
        (4.8, 72.0, 1.4, 1, 325.0, 74.0),
    ]
    rows = [dict(zip(FEATURE_COLUMNS, record)) for record in positives + negatives]
    labels = np.asarray([1] * len(positives) + [0] * len(negatives))
    return pd.DataFrame(rows), labels


def _legacy_classifier(**kwargs) -> CrypticSiteMLClassifier:
    """Build a classifier on the legacy schema with a small search budget."""
    params = {
        "model_type": "random_forest",
        "random_state": 7,
        "n_splits": 2,
        "inner_splits": 2,
        "n_search_iter": 2,
        "feature_names": FEATURE_COLUMNS,
    }
    params.update(kwargs)
    return CrypticSiteMLClassifier(**params)


def _grouped_dataset(
    n_groups: int = 24, seed: int = 0
) -> Tuple[pd.DataFrame, np.ndarray, List[str]]:
    """Build a grouped dataset in the full descriptor schema.

    Each group stands for a structure contributing several pockets, so grouped
    splitting has something real to protect against.
    """
    rng = np.random.default_rng(seed)
    rows: List[dict] = []
    labels: List[int] = []
    groups: List[str] = []

    for group_index in range(n_groups):
        positive_group = group_index % 3 == 0
        # A structure-level offset shared by all of a group's pockets: this is
        # exactly the nuisance signal that leaks when splits are not grouped.
        offset = rng.normal(0.0, 1.0)
        for pocket_index in range(rng.integers(2, 6)):
            is_site = positive_group and pocket_index == 0
            base = {name: float(rng.normal(0.0, 1.0)) for name in DEFAULT_FEATURE_COLUMNS}
            base["burial_depth"] = float(rng.normal(18.0 if is_site else 7.0, 2.0) + offset)
            base["enclosure"] = float(np.clip(rng.normal(0.95 if is_site else 0.5, 0.08), 0, 1))
            base["sasa_mean"] = float(abs(rng.normal(3.0 if is_site else 45.0, 8.0)))
            base["n_basic_residues"] = float(max(0, rng.normal(6.0 if is_site else 2.0, 1.2)))
            base["pocket_volume"] = float(abs(rng.normal(600.0 if is_site else 900.0, 150.0)))
            base["coulomb_potential_kt"] = float(rng.normal(6.0 if is_site else 0.5, 1.5))
            base["plddt_mean"] = float(np.clip(rng.normal(88.0, 5.0), 20, 100))
            rows.append(base)
            labels.append(int(is_site))
            groups.append(f"G{group_index:03d}")

    return pd.DataFrame(rows), np.asarray(labels, dtype=int), groups


# --------------------------------------------------------------- metric units


def test_expected_calibration_error_is_zero_for_perfect_calibration():
    y = np.array([0, 0, 1, 1])
    probabilities = np.array([0.0, 0.0, 1.0, 1.0])
    assert expected_calibration_error(y, probabilities) == pytest.approx(0.0, abs=1e-9)


def test_enrichment_factor_rewards_ranking_positives_first():
    y = np.array([1] * 5 + [0] * 95)
    perfect = np.linspace(1.0, 0.0, 100)
    assert enrichment_factor(y, perfect, fraction=0.05) == pytest.approx(20.0)
    # A ranking that puts every positive last cannot enrich.
    assert enrichment_factor(y, -perfect, fraction=0.05) == pytest.approx(0.0)


def test_precision_at_k_matches_manual_count():
    y = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    scores = np.linspace(1.0, 0.0, 10)
    assert precision_at_k(y, scores, fraction=0.2) == pytest.approx(1.0)


def test_delong_test_reports_no_difference_for_identical_scores():
    y = np.array([0, 0, 1, 1, 0, 1])
    scores = np.array([0.1, 0.2, 0.8, 0.9, 0.3, 0.7])
    difference, p_value = delong_roc_test(y, scores, scores)
    assert difference == pytest.approx(0.0)
    assert p_value == pytest.approx(1.0)


def test_delong_test_detects_a_clear_difference():
    rng = np.random.default_rng(0)
    y = np.array([0] * 100 + [1] * 100)
    good = np.concatenate([rng.normal(0, 1, 100), rng.normal(3, 1, 100)])
    useless = rng.normal(0, 1, 200)
    difference, p_value = delong_roc_test(y, good, useless)
    assert difference > 0.2
    assert p_value < 0.01


def test_grouped_bootstrap_ci_brackets_the_point_estimate():
    rng = np.random.default_rng(1)
    y = np.array([0, 1] * 40)
    probabilities = np.clip(y * 0.6 + rng.normal(0.2, 0.15, y.size), 0, 1)
    groups = [f"g{i // 4}" for i in range(y.size)]
    from sklearn.metrics import roc_auc_score

    point, low, high = grouped_bootstrap_ci(
        y, probabilities, groups, metric=lambda a, b: float(roc_auc_score(a, b)), n_bootstrap=200
    )
    assert low <= point <= high


def test_select_threshold_recovers_a_separating_cut():
    y = np.array([0] * 20 + [1] * 20)
    probabilities = np.array([0.1] * 20 + [0.9] * 20)
    threshold = select_threshold(y, probabilities)
    assert 0.1 < threshold <= 0.9


def test_classification_metrics_handles_single_class_without_raising():
    metrics = classification_metrics([0, 0, 0], [0.1, 0.2, 0.3])
    assert np.isnan(metrics["roc_auc"])
    assert metrics["n"] == 3.0


# ------------------------------------------------------------- model behaviour


def test_random_forest_training_and_curves():
    X, y = _build_training_data()
    clf = _legacy_classifier(random_state=7)

    result = clf.fit(X, y)
    curves = clf.compute_curves(X, y)

    assert clf.best_estimator_ is not None
    assert result.roc_auc >= 0.0
    assert result.pr_auc >= 0.0
    assert result.roc_auc_ci[0] <= result.roc_auc_ci[1] or np.isnan(result.roc_auc_ci[0])
    assert curves["fpr"].shape[0] > 1
    assert curves["precision"].shape[0] > 1
    # Curves must come from out-of-fold predictions, not from training rows.
    assert curves["source"][0] == "out_of_fold"


def test_ml_scorer_adapter_returns_probability():
    X, y = _build_training_data()
    clf = _legacy_classifier(random_state=11)
    clf.fit(X, y)

    scorer = MLPocketScorer(clf)
    score = scorer.calculate_composite_score(
        volume=600.0, depth=17.0, sasa=3.0, basic_count=5, potential=6.5, plddt_confidence=85.0
    )

    assert 0.0 <= score <= 1.0
    assert isinstance(scorer.classify_site(score), str)


def test_ml_scorer_adds_missing_columns_for_richer_schema():
    """A caller with only legacy columns can still score a full-schema model."""
    X, y, groups = _grouped_dataset(n_groups=12, seed=3)
    clf = CrypticSiteMLClassifier(
        model_type="random_forest", random_state=5, n_splits=2, inner_splits=2, n_search_iter=2
    )
    clf.fit(X, y, groups=groups)

    scorer = MLPocketScorer(clf)
    legacy_frame = pd.DataFrame([{"burial_depth": 18.0, "sasa_mean": 3.0}])
    scores = scorer.calculate_composite_scores(legacy_frame)
    assert scores.shape == (1,)
    assert 0.0 <= float(scores[0]) <= 1.0


def test_model_serialization_round_trip(tmp_path: Path):
    X, y = _build_training_data()
    clf = _legacy_classifier(random_state=3)
    clf.fit(X, y)

    model_path = tmp_path / "cryptic_rf.joblib"
    clf.save(str(model_path))

    loaded = CrypticSiteMLClassifier.load(str(model_path))
    assert np.allclose(clf.predict_proba(X), loaded.predict_proba(X))
    assert tuple(loaded.feature_names) == tuple(FEATURE_COLUMNS)
    assert loaded.decision_threshold_ == pytest.approx(clf.decision_threshold_)


def test_missing_feature_columns_raise_a_clear_error():
    X, y = _build_training_data()
    clf = _legacy_classifier()
    with pytest.raises(ValueError, match="Missing feature columns"):
        clf.fit(X.drop(columns=["sasa"]), y)


def test_non_binary_labels_are_rejected():
    X, y = _build_training_data()
    clf = _legacy_classifier()
    with pytest.raises(ValueError, match="Binary 0/1 labels"):
        clf.fit(X, np.asarray([0, 1, 2] + [0] * 7))


def test_grouped_training_keeps_structures_out_of_their_own_evaluation():
    """Every out-of-fold prediction must come from a fold that excluded its group."""
    X, y, groups = _grouped_dataset(n_groups=24, seed=11)
    clf = CrypticSiteMLClassifier(
        model_type="random_forest", random_state=17, n_splits=3, inner_splits=2, n_search_iter=3
    )
    result = clf.fit(X, y, groups=groups)

    group_arr = np.asarray(groups, dtype=object)
    for fold in result.folds:
        test_groups = set(group_arr[fold.test_indices].tolist())
        train_groups = set(group_arr[np.setdiff1d(np.arange(len(X)), fold.test_indices)].tolist())
        assert not (test_groups & train_groups), "grouped split leaked a group across the boundary"

    assert result.n_groups == 24
    assert result.oof_probabilities.size == len(X)
    assert 0.0 <= result.selected_threshold <= 1.0


def test_learned_model_beats_chance_on_a_separable_grouped_dataset():
    X, y, groups = _grouped_dataset(n_groups=30, seed=5)
    clf = CrypticSiteMLClassifier(
        model_type="random_forest", random_state=1, n_splits=3, inner_splits=2, n_search_iter=3
    )
    result = clf.fit(X, y, groups=groups)
    assert result.metrics["roc_auc"] > 0.75
    assert result.metrics["pr_auc"] > result.metrics["prevalence"]


def test_compare_models_produces_an_ordered_table():
    X, y, groups = _grouped_dataset(n_groups=20, seed=8)
    results = compare_models(
        X,
        y,
        groups,
        model_names=["logistic_regression", "random_forest"],
        n_splits=2,
        inner_splits=2,
        n_search_iter=2,
        n_bootstrap=50,
    )
    table = model_comparison_table(results)
    assert set(table["model"]) == {"logistic_regression", "random_forest"}
    assert table["pr_auc"].is_monotonic_decreasing


def test_xgboost_training_runs():
    pytest.importorskip("xgboost")
    X, y = _build_training_data()
    clf = _legacy_classifier(model_type="xgboost", random_state=5)
    result = clf.fit(X, y)
    assert result.pr_auc >= 0.0


# --------------------------------------------------- deployment safety guard


def test_analyzer_refuses_a_model_that_validated_at_chance(tmp_path: Path):
    """A model with no demonstrated skill must not be deployed silently."""
    import json

    from cryptic_ip.analysis.analyzer import ProteinAnalyzer
    from cryptic_ip.analysis.scorer import PocketScorer
    from cryptic_ip.testing.synthetic import (
        SyntheticStructureSpec,
        write_synthetic_structure,
    )

    X, y = _build_training_data()
    clf = _legacy_classifier(random_state=2)
    clf.fit(X, y)

    model_dir = tmp_path / "models"
    model_dir.mkdir()
    model_path = model_dir / "chance_model.pkl"
    clf.save(str(model_path))
    (model_dir / "chance_model.metadata.json").write_text(
        json.dumps(
            {
                "selected_model_type": "random_forest",
                "model_candidates": {"random_forest": {"metrics": {"roc_auc": 0.496}}},
            }
        ),
        encoding="utf-8",
    )

    structure = write_synthetic_structure(
        SyntheticStructureSpec(name="GUARD", seed=3), tmp_path / "structures"
    )
    analyzer = ProteinAnalyzer(
        str(structure), use_ml_model=True, model_path=str(model_path), skip_electrostatics=True
    )
    assert isinstance(analyzer.scorer, PocketScorer)


def test_analyzer_deploys_a_model_with_demonstrated_skill(tmp_path: Path):
    import json

    from cryptic_ip.analysis.analyzer import ProteinAnalyzer
    from cryptic_ip.testing.synthetic import (
        SyntheticStructureSpec,
        write_synthetic_structure,
    )

    X, y = _build_training_data()
    clf = _legacy_classifier(random_state=2)
    clf.fit(X, y)

    model_dir = tmp_path / "models"
    model_dir.mkdir()
    model_path = model_dir / "good_model.pkl"
    clf.save(str(model_path))
    (model_dir / "good_model.metadata.json").write_text(
        json.dumps(
            {
                "selected_model_type": "random_forest",
                "model_candidates": {"random_forest": {"metrics": {"roc_auc": 0.93}}},
            }
        ),
        encoding="utf-8",
    )

    structure = write_synthetic_structure(
        SyntheticStructureSpec(name="GOOD", seed=4), tmp_path / "structures"
    )
    analyzer = ProteinAnalyzer(
        str(structure), use_ml_model=True, model_path=str(model_path), skip_electrostatics=True
    )
    assert isinstance(analyzer.scorer, MLPocketScorer)
