"""Guards for two upstream removals that broke the pipeline on current installs.

Both were live failures, not warnings: ``CalibratedClassifierCV(cv="prefit")`` was
removed in scikit-learn 1.8 and raises, and ``np.trapz`` was removed from the NumPy
namespace after 2.2. The pinned versions in requirements.txt still carry both, so CI
did not catch either; these tests fail on the pins only if the compatibility shims
are removed, and pass on every version from the pins forward.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _fitted_base(n: int = 200):
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(0)
    x = rng.normal(size=(n, 4))
    y = (x[:, 0] + rng.normal(scale=0.3, size=n) > 0).astype(int)
    return LogisticRegression().fit(x, y), x, y


def test_prefit_calibrator_does_not_refit_the_base_estimator():
    from cryptic_ip.analysis.ml_classifier import _prefit_calibrator

    base, x, y = _fitted_base()
    coefficients = base.coef_.copy()

    calibrated = _prefit_calibrator(base, method="sigmoid", n_samples=len(y))
    calibrated.fit(x, y)

    # The whole point of freezing: the calibrator must not refit the base.
    assert np.allclose(base.coef_, coefficients)
    probabilities = calibrated.predict_proba(x)
    assert probabilities.shape == (len(y), 2)
    assert np.all(np.isfinite(probabilities))
    assert np.all((probabilities >= 0) & (probabilities <= 1))


def test_prefit_calibrator_fits_exactly_one_calibrator_over_all_the_data():
    """``cv="prefit"`` fitted one calibrator on everything; a default cv would fit five."""
    from cryptic_ip.analysis.ml_classifier import _prefit_calibrator

    base, x, y = _fitted_base()
    calibrated = _prefit_calibrator(base, method="sigmoid", n_samples=len(y)).fit(x, y)
    assert len(calibrated.calibrated_classifiers_) == 1


def test_prefit_calibrator_survives_a_calibration_set_smaller_than_a_cv_fold():
    """The regression this guards: four rows is fewer than the default five folds.

    The synthetic integration fixture carves a four-row calibration set, so a frozen
    estimator left on the default ``cv`` raised ValueError and no model trained at all.
    """
    from cryptic_ip.analysis.ml_classifier import _prefit_calibrator

    base, x, _ = _fitted_base()
    tiny_x, tiny_y = x[:4], np.array([0, 1, 0, 1])
    calibrated = _prefit_calibrator(base, method="sigmoid", n_samples=4).fit(tiny_x, tiny_y)
    assert len(calibrated.calibrated_classifiers_) == 1
    assert np.all(np.isfinite(calibrated.predict_proba(tiny_x)))


def _figures_module():
    spec = importlib.util.spec_from_file_location(
        "generate_publication_figures", ROOT / "scripts" / "generate_publication_figures.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["generate_publication_figures"] = module
    spec.loader.exec_module(module)
    return module


def test_trapezoidal_auc_does_not_depend_on_the_removed_numpy_alias():
    figures = _figures_module()
    assert figures._auc([0.0, 0.0, 1.0], [0.0, 1.0, 1.0]) == pytest.approx(1.0)
    assert figures._auc([0.0, 1.0], [0.0, 1.0]) == pytest.approx(0.5)


def test_statistical_validation_shares_the_same_numpy_guard():
    from cryptic_ip.analysis.statistical_validation import _trapz

    assert _trapz(np.array([0.0, 1.0]), np.array([0.0, 1.0])) == pytest.approx(0.5)
