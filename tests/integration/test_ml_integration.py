"""Integration tests for ML classifier handoff into batch processing."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cryptic_ip.analysis.ml_classifier import CrypticSiteMLClassifier, MLPocketScorer
from cryptic_ip.database.batch_processing import AnalysisCache, append_results_to_file

try:
    import sklearn  # noqa: F401

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


@pytest.mark.integration
@pytest.mark.skipif(not SKLEARN_AVAILABLE, reason="scikit-learn not installed")
def test_classifier_training_integrates_with_batch_outputs(tmp_path: Path, validation_feature_set):
    """Model probabilities should flow into batch result rows and cache/export outputs."""
    X, y = validation_feature_set

    clf = CrypticSiteMLClassifier(random_state=12)
    clf.fit(X, y)
    scorer = MLPocketScorer(clf)

    proteins = [
        {"uniprot_id": "Q9TEST1", "organism": "human", "feature_tuple": (620.0, 17.4, 3.0, 6, 7.1, 90.0)},
        {"uniprot_id": "Q9TEST2", "organism": "human", "feature_tuple": (370.0, 6.8, 55.0, 2, 2.0, 79.0)},
    ]

    rows = []
    for item in proteins:
        volume, depth, sasa, basic_count, potential, plddt = item["feature_tuple"]
        proba = scorer.calculate_composite_score(
            volume=volume,
            depth=depth,
            sasa=sasa,
            basic_count=basic_count,
            potential=potential,
            plddt_confidence=plddt,
        )
        rows.append(
            {
                "uniprot_id": item["uniprot_id"],
                "organism": item["organism"],
                "composite_score": proba,
                "prediction_label": scorer.classify_site(proba),
            }
        )

    out_csv = append_results_to_file(rows, tmp_path / "batch_predictions.csv")
    assert out_csv.exists()

    parsed = pd.read_csv(out_csv)
    assert parsed["composite_score"].between(0.0, 1.0).all()

    cache = AnalysisCache(tmp_path / "analysis.sqlite", pipeline_version="ml-v1", pipeline_params={"model": "rf"})
    for row in rows:
        cache.set_cached_result(row["uniprot_id"], row["organism"], row)

    restored = cache.get_cached_result("Q9TEST1", "human")
    assert restored is not None
    assert 0.0 <= restored["composite_score"] <= 1.0

    exported_json = cache.export_results(tmp_path / "analysis_cache.json", "json")
    cache.close()
    assert exported_json.exists()
