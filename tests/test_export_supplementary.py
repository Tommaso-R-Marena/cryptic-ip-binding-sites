"""Tests for supplementary table export."""

import json
from pathlib import Path

import pandas as pd
import pytest

from scripts.export_supplementary_tables import export_supplementary


@pytest.fixture
def publication_tree(tmp_path: Path) -> dict:
    pub = tmp_path / "publication"
    (pub / "validation").mkdir(parents=True)
    (pub / "ml_training").mkdir()
    (pub / "comparative").mkdir()
    (pub / "gallery").mkdir()

    dataset = tmp_path / "dataset.csv"
    dataset.write_text("pdb_id,label\n1ZY7,1\n", encoding="utf-8")

    pd.DataFrame({"pdb_id": ["1ZY7"], "composite_score": [0.8]}).to_csv(
        pub / "validation" / "control_benchmark.csv", index=False
    )
    pd.DataFrame({"method": ["threshold"], "roc_auc": [0.55]}).to_csv(
        pub / "ml_training" / "ml_vs_threshold_comparison.csv", index=False
    )
    pd.DataFrame({"organism": ["yeast"], "hit_rate": [0.02]}).to_csv(
        pub / "comparative" / "hit_rates.csv", index=False
    )
    pd.DataFrame({"pdb_id": ["1ZY7"], "rank": [1]}).to_csv(
        pub / "gallery" / "gallery_inputs.csv", index=False
    )
    pd.DataFrame({"pdb_id": ["1ZY7"], "label": [1], "volume": [500]}).to_csv(
        pub / "ml_training" / "validation_pocket_features.csv", index=False
    )

    yeast = tmp_path / "yeast"
    yeast.mkdir()
    (yeast / "yeast_pilot_summary.json").write_text('{"n_structures": 10}', encoding="utf-8")
    pd.DataFrame({"uniprot_id": ["P12345"], "composite_score": [0.8]}).to_csv(
        yeast / "yeast_pilot_hits.csv", index=False
    )

    return {"pub": pub, "dataset": dataset, "yeast": yeast}


def test_export_supplementary_writes_manifest_and_index(publication_tree: dict):
    out = export_supplementary(
        publication_tree["pub"],
        publication_tree["dataset"],
        publication_tree["yeast"],
    )

    assert (out / "S1_validation_dataset.csv").exists()
    assert (out / "S2_control_benchmark.csv").exists()
    assert (out / "S6_yeast_pilot_summary.json").exists()
    assert (out / "SUPPLEMENTARY_INDEX.md").exists()

    manifest = json.loads((out / "export_manifest.json").read_text(encoding="utf-8"))
    table_names = {entry["table"] for entry in manifest}
    assert "S2_control_benchmark.csv" in table_names
    assert "S3_ml_feature_summary.json" in table_names

    index = (out / "SUPPLEMENTARY_INDEX.md").read_text(encoding="utf-8")
    assert "S1_validation_dataset.csv" in index
    assert "included" in index
