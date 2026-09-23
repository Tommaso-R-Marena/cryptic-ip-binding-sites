"""Tests for the publication figure generator (headless, synthetic inputs)."""

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
STYLE = ROOT / "scripts" / "nature_publication.mplstyle"


@pytest.fixture(scope="module")
def figmod():
    spec = importlib.util.spec_from_file_location(
        "generate_publication_figures", ROOT / "scripts" / "generate_publication_figures.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.apply_publication_style(STYLE)
    return module


def test_auc_of_perfect_and_chance_classifiers(figmod):
    assert figmod._auc([0.0, 0.0, 1.0], [0.0, 1.0, 1.0]) == pytest.approx(1.0)
    assert figmod._auc([0.0, 1.0], [0.0, 1.0]) == pytest.approx(0.5)


def test_make_figure4_writes_outputs(figmod, tmp_path: Path):
    dataset = pd.DataFrame(
        {
            "pdb_id": ["1ZY7", "5ICN", "1MAI", "1BQ3"],
            "sasa": [0.5, 12.0, 99.0, 355.0],
            "classification": ["Cryptic", "Semi-cryptic", "Surface", "Surface"],
            "resolution": [1.7, 2.8, 2.4, 2.7],
            "organism": ["NA", "NA", "NA", "NA"],
        }
    )
    controls = pd.DataFrame(
        {
            "protein": ["ADAR2", "PLCd1_PH"],
            "control_type": ["positive", "negative"],
            "score": [0.77, 0.24],
        }
    )
    dataset_csv = tmp_path / "dataset.csv"
    controls_csv = tmp_path / "controls.csv"
    dataset.to_csv(dataset_csv, index=False)
    controls.to_csv(controls_csv, index=False)

    config = {
        "figure4": {
            "dataset_csv": str(dataset_csv),
            "control_benchmark_csv": str(controls_csv),
        }
    }
    figmod.make_figure4(config, tmp_path)

    assert (tmp_path / "Figure4_Validation_Benchmark.png").exists()
    assert (tmp_path / "Figure4_Validation_Benchmark.pdf").exists()


def test_make_figure4_without_controls(figmod, tmp_path: Path):
    dataset = pd.DataFrame(
        {
            "pdb_id": ["1ZY7", "1MAI"],
            "sasa": [0.5, 99.0],
            "classification": ["Cryptic", "Surface"],
            "resolution": [1.7, 2.4],
        }
    )
    dataset_csv = tmp_path / "dataset.csv"
    dataset.to_csv(dataset_csv, index=False)

    config = {"figure4": {"dataset_csv": str(dataset_csv)}}
    figmod.make_figure4(config, tmp_path)

    assert (tmp_path / "Figure4_Validation_Benchmark.png").exists()
