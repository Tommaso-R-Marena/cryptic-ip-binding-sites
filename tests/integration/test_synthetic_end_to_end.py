"""Offline end-to-end run: synthetic structures -> features -> labels -> model.

This exercises the complete pipeline, including the two command-line scripts, in
an environment with no access to the PDB. It is the check that the project is
actually runnable, as opposed to individually unit-tested.
"""

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]

pytestmark = pytest.mark.skipif(
    shutil.which("fpocket") is None, reason="fpocket not installed"
)


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory) -> dict:
    """Build structures, extract features and train a model once for the module."""
    from cryptic_ip.testing.synthetic import build_synthetic_benchmark

    work = tmp_path_factory.mktemp("synthetic_e2e")
    build_synthetic_benchmark(work / "structures", n_buried=6, n_surface=6, n_decoy=8)

    features_csv = work / "pocket_features.csv"
    summary_json = work / "labeling_summary.json"
    extract = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "extract_pocket_features.py"),
            "--structures-dir", str(work / "structures"),
            "--entry-csv", str(work / "absent.csv"),
            "--output-csv", str(features_csv),
            "--cache-dir", str(work / "cache"),
            "--summary-json", str(summary_json),
            "--jobs", "2",
            "--sasa-points", "128",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert extract.returncode == 0, extract.stderr[-3000:]

    train = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts" / "train_ml_classifier.py"),
            "--features-csv", str(features_csv),
            "--work-dir", str(work / "ml"),
            "--model-dir", str(work / "models"),
            "--model-name", "synthetic_model",
            "--models", "logistic_regression", "random_forest",
            "--n-splits", "2",
            "--inner-splits", "2",
            "--n-search-iter", "3",
            "--n-bootstrap", "50",
            "--no-figures",
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert train.returncode == 0, train.stderr[-3000:]

    return {
        "work": work,
        "features_csv": features_csv,
        "summary": json.loads(summary_json.read_text()),
    }


def test_feature_extraction_produces_both_classes(synthetic_run):
    frame = pd.read_csv(synthetic_run["features_csv"])
    assert len(frame) > 0
    assert set(frame["label"].unique()) >= {0, 1}
    assert frame["group_key"].nunique() > 1


def test_every_descriptor_column_is_present(synthetic_run):
    from cryptic_ip.analysis.features import FEATURE_NAMES

    frame = pd.read_csv(synthetic_run["features_csv"])
    assert set(FEATURE_NAMES) <= set(frame.columns)


def test_labels_carry_their_provenance(synthetic_run):
    """Each label must be explainable, not just a number."""
    frame = pd.read_csv(synthetic_run["features_csv"])
    assert "label_reason" in frame.columns
    assert "overlap_fraction" in frame.columns
    positives = frame[frame["label"] == 1]
    assert (positives["overlap_fraction"] >= 0.30).all()


def test_detector_recall_is_reported(synthetic_run):
    summary = synthetic_run["summary"]
    assert "site_recall" in summary
    assert summary["n_structures_with_site"] > 0
    # Every structure that yielded no matching pocket is accounted for.
    assert (
        summary["n_structures_with_site"] + summary["n_structures_with_ligand_but_no_site"]
        > 0
    )


def test_training_writes_model_card_and_metadata(synthetic_run):
    models = synthetic_run["work"] / "models"
    assert (models / "synthetic_model.pkl").exists()
    card = (models / "synthetic_model_model_card.md").read_text()
    assert "Intended use" in card
    assert "Known limitations" in card
    assert "grouped by protein" in card

    metadata = json.loads((models / "synthetic_model.metadata.json").read_text())
    assert metadata["split"].startswith("nested StratifiedGroupKFold")
    assert metadata["labeling"].startswith("atom-level ligand overlap")
    assert metadata["n_groups"] > 1
    assert 0.0 <= metadata["decision_threshold"] <= 1.0


def test_model_comparison_includes_the_rule_based_baseline(synthetic_run):
    comparison = pd.read_csv(synthetic_run["work"] / "ml" / "ml_vs_threshold_comparison.csv")
    assert set(comparison["method"]).issuperset({"Rule-based scorer"})
    assert comparison["roc_auc"].notna().all()


def test_saved_model_scores_new_pockets(synthetic_run):
    """The serialised artefact must be usable for inference, not just storable."""
    from cryptic_ip.analysis.ml_classifier import CrypticSiteMLClassifier

    model_path = synthetic_run["work"] / "models" / "synthetic_model.pkl"
    classifier = CrypticSiteMLClassifier.load(str(model_path))
    frame = pd.read_csv(synthetic_run["features_csv"])
    probabilities = classifier.predict_proba(frame)
    assert probabilities.shape == (len(frame),)
    assert ((probabilities >= 0.0) & (probabilities <= 1.0)).all()


def test_analyzer_uses_the_saved_model_end_to_end(synthetic_run):
    """A trained model must drop straight into single-structure analysis."""
    from cryptic_ip.analysis.analyzer import ProteinAnalyzer

    structure = next((synthetic_run["work"] / "structures").glob("SYN_BURIED_*.pdb"))
    model_path = synthetic_run["work"] / "models" / "synthetic_model.pkl"
    analyzer = ProteinAnalyzer(
        str(structure),
        use_ml_model=True,
        model_path=str(model_path),
        skip_electrostatics=True,
    )
    scored = analyzer.run_pipeline(include_electrostatics=False)

    assert len(scored) > 0
    assert "composite_score" in scored.columns
    assert scored["composite_score"].between(0.0, 1.0).all()
    # The descriptor suite travels with the results, not just the score.
    assert "burial_depth" in scored.columns
    assert "enclosure" in scored.columns
    analyzer.cleanup()
