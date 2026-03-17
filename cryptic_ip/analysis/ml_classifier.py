"""Machine-learning models for cryptic IP-binding pocket classification."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    import joblib
    from sklearn.base import BaseEstimator
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.impute import SimpleImputer
    from sklearn.metrics import (
        auc,
        average_precision_score,
        precision_recall_curve,
        roc_auc_score,
        roc_curve,
    )
    from sklearn.model_selection import GridSearchCV, StratifiedKFold, cross_val_predict
    from sklearn.pipeline import Pipeline

    SKLEARN_AVAILABLE = True
except ImportError:  # pragma: no cover - import guard for optional dependency
    joblib = None
    BaseEstimator = Any
    RandomForestClassifier = Any
    SimpleImputer = Any
    GridSearchCV = Any
    StratifiedKFold = Any
    cross_val_predict = Any
    Pipeline = Any
    SKLEARN_AVAILABLE = False


def _ensure_sklearn_installed() -> None:
    """Validate that required ML dependencies are installed."""
    if not SKLEARN_AVAILABLE:
        raise ImportError(
            "scikit-learn and joblib are required for ML classification. "
            "Install dependencies with `pip install scikit-learn joblib`."
        )


FEATURE_COLUMNS: Tuple[str, ...] = (
    "pocket_depth",
    "sasa",
    "electrostatic_potential",
    "n_basic_residues",
    "pocket_volume",
    "plddt_confidence",
)


@dataclass
class TrainingResults:
    """Container for model training artifacts and validation metrics."""

    best_params: Dict[str, Any]
    roc_auc: float
    pr_auc: float
    roc_auc_ci: Tuple[float, float]
    pr_auc_ci: Tuple[float, float]


class CrypticSiteMLClassifier:
    """Train, validate, and deploy tree-based classifiers for pocket classification.

    The classifier supports Random Forest and XGBoost backends with hyperparameter
    optimization, 5-fold stratified cross-validation, SHAP-based feature attribution,
    and model serialization.
    """

    def __init__(
        self,
        model_type: str = "random_forest",
        random_state: int = 42,
        n_splits: int = 5,
    ) -> None:
        _ensure_sklearn_installed()
        self.model_type = model_type
        self.random_state = random_state
        self.n_splits = n_splits
        self.pipeline: Pipeline = self._build_pipeline(model_type)
        self.best_estimator_: Optional[Pipeline] = None

    def _build_pipeline(self, model_type: str) -> Pipeline:
        if model_type == "random_forest":
            estimator: BaseEstimator = RandomForestClassifier(random_state=self.random_state)
        elif model_type == "xgboost":
            try:
                from xgboost import XGBClassifier
            except ImportError as exc:  # pragma: no cover - optional dependency
                raise ImportError(
                    "xgboost is required for model_type='xgboost'. "
                    "Install with `pip install xgboost`."
                ) from exc
            estimator = XGBClassifier(
                random_state=self.random_state,
                eval_metric="logloss",
                use_label_encoder=False,
            )
        else:
            raise ValueError("model_type must be 'random_forest' or 'xgboost'.")

        return Pipeline(
            steps=[
                ("imputer", SimpleImputer(strategy="median")),
                ("classifier", estimator),
            ]
        )

    def _param_grid(self) -> Dict[str, List[Any]]:
        if self.model_type == "random_forest":
            return {
                "classifier__n_estimators": [200, 500],
                "classifier__max_depth": [None, 8, 12],
                "classifier__min_samples_split": [2, 4],
                "classifier__min_samples_leaf": [1, 2],
                "classifier__class_weight": [None, "balanced"],
            }
        return {
            "classifier__n_estimators": [200, 500],
            "classifier__max_depth": [3, 6],
            "classifier__learning_rate": [0.03, 0.1],
            "classifier__subsample": [0.8, 1.0],
            "classifier__colsample_bytree": [0.8, 1.0],
        }

    def _validate_inputs(
        self, features: pd.DataFrame, labels: Iterable[int]
    ) -> Tuple[pd.DataFrame, np.ndarray]:
        missing = [col for col in FEATURE_COLUMNS if col not in features.columns]
        if missing:
            raise ValueError(f"Missing feature columns: {missing}")

        y = np.asarray(list(labels))
        if len(np.unique(y)) != 2:
            raise ValueError("Binary labels are required (0/1).")

        counts = np.bincount(y)
        if np.min(counts) < self.n_splits:
            raise ValueError(
                f"5-fold stratified CV requires at least {self.n_splits} examples in each class; "
                f"class counts are {counts.tolist()}."
            )

        return features.loc[:, FEATURE_COLUMNS].copy(), y

    def fit(self, features: pd.DataFrame, labels: Iterable[int]) -> TrainingResults:
        """Fit the classifier with hyperparameter tuning and CV diagnostics."""
        X, y = self._validate_inputs(features, labels)
        splitter = StratifiedKFold(
            n_splits=self.n_splits, shuffle=True, random_state=self.random_state
        )

        grid = GridSearchCV(
            estimator=self.pipeline,
            param_grid=self._param_grid(),
            scoring="average_precision",
            cv=splitter,
            n_jobs=-1,
            refit=True,
        )
        grid.fit(X, y)
        self.best_estimator_ = grid.best_estimator_

        probabilities = cross_val_predict(
            self.best_estimator_,
            X,
            y,
            cv=splitter,
            method="predict_proba",
            n_jobs=-1,
        )[:, 1]

        roc_auc_value = roc_auc_score(y, probabilities)
        pr_auc_value = average_precision_score(y, probabilities)
        roc_ci, pr_ci = self._bootstrap_auc_ci(y, probabilities)

        return TrainingResults(
            best_params=grid.best_params_,
            roc_auc=float(roc_auc_value),
            pr_auc=float(pr_auc_value),
            roc_auc_ci=roc_ci,
            pr_auc_ci=pr_ci,
        )

    def _bootstrap_auc_ci(
        self,
        labels: np.ndarray,
        probabilities: np.ndarray,
        n_bootstrap: int = 1000,
        alpha: float = 0.95,
    ) -> Tuple[Tuple[float, float], Tuple[float, float]]:
        rng = np.random.default_rng(self.random_state)
        roc_values: List[float] = []
        pr_values: List[float] = []

        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(labels), len(labels))
            y_boot = labels[idx]
            p_boot = probabilities[idx]
            if len(np.unique(y_boot)) < 2:
                continue
            roc_values.append(roc_auc_score(y_boot, p_boot))
            pr_values.append(average_precision_score(y_boot, p_boot))

        lower = (1 - alpha) / 2
        upper = 1 - lower
        return (
            (float(np.quantile(roc_values, lower)), float(np.quantile(roc_values, upper))),
            (float(np.quantile(pr_values, lower)), float(np.quantile(pr_values, upper))),
        )

    def compute_curves(
        self, features: pd.DataFrame, labels: Iterable[int]
    ) -> Dict[str, np.ndarray]:
        """Compute ROC and PR curve coordinates from fitted model predictions."""
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")

        X, y = self._validate_inputs(features, labels)
        probabilities = self.best_estimator_.predict_proba(X)[:, 1]

        fpr, tpr, _ = roc_curve(y, probabilities)
        precision, recall, _ = precision_recall_curve(y, probabilities)

        return {
            "fpr": fpr,
            "tpr": tpr,
            "roc_auc": np.array([auc(fpr, tpr)]),
            "precision": precision,
            "recall": recall,
            "pr_auc": np.array([auc(recall, precision)]),
        }

    def shap_values(self, features: pd.DataFrame) -> pd.DataFrame:
        """Compute mean absolute SHAP values for each feature."""
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")

        try:
            import shap
        except ImportError as exc:  # pragma: no cover - optional dependency
            raise ImportError("SHAP is required for feature importance analysis.") from exc

        X = features.loc[:, FEATURE_COLUMNS].copy()
        transformed = self.best_estimator_.named_steps["imputer"].transform(X)
        model = self.best_estimator_.named_steps["classifier"]

        explainer = shap.TreeExplainer(model)
        values = explainer.shap_values(transformed)
        if isinstance(values, list):
            values = values[1]

        importance = np.abs(values).mean(axis=0)
        return pd.DataFrame({"feature": FEATURE_COLUMNS, "mean_abs_shap": importance}).sort_values(
            "mean_abs_shap", ascending=False
        )

    def predict_proba(self, features: pd.DataFrame) -> np.ndarray:
        """Predict cryptic-site probabilities."""
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        X = features.loc[:, FEATURE_COLUMNS].copy()
        return self.best_estimator_.predict_proba(X)[:, 1]

    def save(self, path: str) -> None:
        """Serialize model artifact for deployment."""
        if self.best_estimator_ is None:
            raise RuntimeError("Model is not trained. Run fit() first.")
        payload = {
            "model_type": self.model_type,
            "random_state": self.random_state,
            "n_splits": self.n_splits,
            "model": self.best_estimator_,
        }
        joblib.dump(payload, path)

    @classmethod
    def load(cls, path: str) -> "CrypticSiteMLClassifier":
        """Load a serialized model artifact."""
        _ensure_sklearn_installed()
        payload = joblib.load(path)
        classifier = cls(
            model_type=payload["model_type"],
            random_state=payload["random_state"],
            n_splits=payload["n_splits"],
        )
        classifier.best_estimator_ = payload["model"]
        return classifier


class MLPocketScorer:
    """Drop-in scorer adapter that replaces threshold scoring with ML probabilities."""

    def __init__(self, classifier: CrypticSiteMLClassifier) -> None:
        if classifier.best_estimator_ is None:
            raise ValueError("Classifier must be trained before using MLPocketScorer.")
        self.classifier = classifier

    def calculate_composite_score(
        self,
        volume: float,
        depth: float,
        sasa: float,
        basic_count: int,
        potential: Optional[float] = None,
        plddt_confidence: float = np.nan,
    ) -> float:
        """Return probability estimate of a cryptic IP-binding site."""
        sample = pd.DataFrame(
            [
                {
                    "pocket_depth": depth,
                    "sasa": sasa,
                    "electrostatic_potential": potential,
                    "n_basic_residues": basic_count,
                    "pocket_volume": volume,
                    "plddt_confidence": plddt_confidence,
                }
            ]
        )
        return float(self.classifier.predict_proba(sample)[0])


    def calculate_composite_scores(self, samples: pd.DataFrame) -> np.ndarray:
        """Vectorized probability inference over many pockets."""
        return self.classifier.predict_proba(samples)

    def classify_site(self, score: float) -> str:
        """Classify model probability into confidence buckets."""
        if score >= 0.8:
            return "High confidence cryptic IP site"
        if score >= 0.6:
            return "Moderate confidence candidate"
        if score >= 0.4:
            return "Low confidence - manual inspection recommended"
        return "Unlikely cryptic IP site"
