"""Tests for machine-learning pocket classifier module."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.ml_classifier import (
    FEATURE_COLUMNS,
    CrypticSiteMLClassifier,
    MLPocketScorer,
)

from typing import Tuple


try:
    import sklearn  # noqa: F401
    import joblib  # noqa: F401

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


def _build_training_data() -> Tuple[pd.DataFrame, np.ndarray]:
    """Create a small but separable dataset for unit testing."""
    rows = []
    labels = []

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

    for record in positives:
        rows.append(dict(zip(FEATURE_COLUMNS, record)))
        labels.append(1)
    for record in negatives:
        rows.append(dict(zip(FEATURE_COLUMNS, record)))
        labels.append(0)

    return pd.DataFrame(rows), np.asarray(labels)


@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="scikit-learn/joblib not installed")
def test_random_forest_training_and_curves():
    """Random forest fit should produce trained estimator and valid metrics."""
    X, y = _build_training_data()
    clf = CrypticSiteMLClassifier(model_type="random_forest", random_state=7)

    result = clf.fit(X, y)
    curves = clf.compute_curves(X, y)

    assert clf.best_estimator_ is not None
    assert result.roc_auc >= 0.0
    assert result.pr_auc >= 0.0
    assert result.roc_auc_ci[0] <= result.roc_auc_ci[1]
    assert result.pr_auc_ci[0] <= result.pr_auc_ci[1]
    assert curves["fpr"].shape[0] > 1
    assert curves["precision"].shape[0] > 1


@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="scikit-learn/joblib not installed")
def test_ml_scorer_adapter_returns_probability():
    """ML scorer adapter should emit a bounded score and classification label."""
    X, y = _build_training_data()
    clf = CrypticSiteMLClassifier(model_type="random_forest", random_state=11)
    clf.fit(X, y)

    scorer = MLPocketScorer(clf)
    score = scorer.calculate_composite_score(
        volume=600.0,
        depth=17.0,
        sasa=3.0,
        basic_count=5,
        potential=6.5,
        plddt_confidence=85.0,
    )

    assert 0.0 <= score <= 1.0
    assert isinstance(scorer.classify_site(score), str)


@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="scikit-learn/joblib not installed")
def test_model_serialization_round_trip(tmp_path: Path):
    """Serialized models should preserve prediction behavior."""
    X, y = _build_training_data()
    clf = CrypticSiteMLClassifier(model_type="random_forest", random_state=3)
    clf.fit(X, y)

    model_path = tmp_path / "cryptic_rf.joblib"
    clf.save(str(model_path))

    loaded = CrypticSiteMLClassifier.load(str(model_path))
    original = clf.predict_proba(X)
    restored = loaded.predict_proba(X)

    assert np.allclose(original, restored)


try:
    import xgboost  # noqa: F401

    XGBOOST_AVAILABLE = True
except ImportError:
    XGBOOST_AVAILABLE = False


@pytest.mark.skipif(
    not (SKLEARN_AVAILABLE and XGBOOST_AVAILABLE), reason="ml dependencies not installed"
)
def test_xgboost_training_runs():
    """XGBoost backend should fit when dependency is available."""
    X, y = _build_training_data()
    clf = CrypticSiteMLClassifier(model_type="xgboost", random_state=5)
    result = clf.fit(X, y)
    assert result.pr_auc >= 0.0
