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


def test_frozen_lets_a_calibrator_accept_an_already_fitted_estimator():
    from sklearn.calibration import CalibratedClassifierCV
    from sklearn.linear_model import LogisticRegression

    from cryptic_ip.analysis.ml_classifier import _frozen, _frozen_kwargs

    rng = np.random.default_rng(0)
    x = rng.normal(size=(200, 4))
    y = (x[:, 0] + rng.normal(scale=0.3, size=200) > 0).astype(int)
    base = LogisticRegression().fit(x[:120], y[:120])
    coefficients = base.coef_.copy()

    calibrated = CalibratedClassifierCV(_frozen(base), method="sigmoid", **_frozen_kwargs())
    calibrated.fit(x[120:], y[120:])

    # The point of freezing: the calibrator must not refit the base estimator.
    assert np.allclose(base.coef_, coefficients)
    probabilities = calibrated.predict_proba(x)
    assert probabilities.shape == (200, 2)
    assert np.all(np.isfinite(probabilities))
    assert np.all((probabilities >= 0) & (probabilities <= 1))


def test_frozen_kwargs_only_asks_for_prefit_when_freezing_is_unavailable():
    from cryptic_ip.analysis.ml_classifier import _frozen_kwargs

    try:
        import sklearn.frozen  # noqa: F401
    except ImportError:  # scikit-learn < 1.6
        assert _frozen_kwargs() == {"cv": "prefit"}
    else:
        assert _frozen_kwargs() == {}


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
