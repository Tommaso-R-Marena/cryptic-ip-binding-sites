"""End-to-end test of the proteome-screen CLI with the structure analysis stubbed.

The analyzer itself needs fpocket and is exercised elsewhere; this checks the
plumbing that turns per-structure results into hit rates - catalogue, shards,
the all-pockets record, failure accounting and aggregation.
"""

from __future__ import annotations

import json

import pandas as pd
import pytest

from scripts import proteome_screen
from tests.test_proteome_catalog import _model_text, _write


class _StubAnalyzer:
    """Two pockets per structure; P0 has a strict hit, P2 raises."""

    def __init__(self, pdb_path, **_):
        self.pdb_path = str(pdb_path)

    def detect_pockets(self):
        if "P00002" in self.pdb_path:
            raise RuntimeError("fpocket produced no output")

    def score_all_pockets(self):
        hit = "P00000" in self.pdb_path
        return pd.DataFrame(
            [
                {
                    "pocket_id": 1,
                    "composite_score": 0.82 if hit else 0.55,
                    "sasa": 4.0,
                    "basic_residues": 6,
                    "volume": 1500.0,
                    "plddt_mean": 92.0,
                    "burial_depth": 6.0,
                    "hull_depth": 14.0,
                    "center": (1.0, 2.0, 3.0),
                    "pocket_residue_numbers": [10, 11, 12],
                },
                {
                    "pocket_id": 2,
                    "composite_score": 0.40,
                    "sasa": 30.0,
                    "basic_residues": 1,
                    "volume": 200.0,
                    "plddt_mean": 55.0,
                    "burial_depth": 1.0,
                    "hull_depth": 2.0,
                    "center": (0.0, 0.0, 0.0),
                    "pocket_residue_numbers": [1],
                },
            ]
        )


def test_catalog_screen_aggregate(tmp_path, monkeypatch):
    import cryptic_ip.analysis as analysis

    monkeypatch.setattr(analysis, "ProteinAnalyzer", _StubAnalyzer, raising=False)
    models = tmp_path / "models"
    models.mkdir()
    for i in range(4):
        _write(models, f"AF-P0000{i}-F1-model_v4.pdb", _model_text(300 + i))
    catalog_dir = tmp_path / "catalog"
    shards_dir = tmp_path / "shards"
    summary_dir = tmp_path / "summary"

    assert proteome_screen.main(
        ["catalog", "--organism", "yeast", "--structures-dir", str(models), "--output-dir", str(catalog_dir)]
    ) == 0
    qc = json.loads((catalog_dir / "yeast_qc.json").read_text())
    assert qc["screenable_f1_models"] == 4

    for index in range(2):
        assert proteome_screen.main(
            [
                "screen", "--organism", "yeast",
                "--catalog", str(catalog_dir / "yeast_catalog.csv"),
                "--shard-index", str(index), "--shard-count", "2",
                "--workers", "1", "--output-dir", str(shards_dir),
            ]
        ) == 0

    pockets = pd.concat(pd.read_csv(p) for p in shards_dir.glob("*_pockets_part*.csv.gz"))
    # Every pocket is kept, including the one that fails every gate.
    assert len(pockets) == 6
    assert {"center_x", "pocket_residues", "organism_key"} <= set(pockets.columns)
    # The pocket's hull depth is its descriptor, not recomputed by the worker.
    assert sorted(pockets["hull_depth"].unique()) == [2.0, 14.0]
    assert pockets["protein_max_hull_depth"].notna().all()
    status = pd.concat(pd.read_csv(p) for p in shards_dir.glob("*_status.csv"))
    assert len(status) == 4
    assert status.loc[status["uniprot_id"] != "P00002", "fpocket_seconds"].notna().all()

    # The stub's scores are not derived from its descriptors, so they are kept.
    assert proteome_screen.main(
        ["aggregate", "--shards-dir", str(shards_dir), "--catalog-dir", str(catalog_dir),
         "--output-dir", str(summary_dir), "--keep-screen-scores"]
    ) == 0
    summary = json.loads((summary_dir / "screen_summary.json").read_text())
    assert summary["scoring"] == {"rescored": False}
    rates = summary["plan"]["hit_rates"][0]
    # The failed structure is reported, and excluded from the denominator.
    assert summary["structures"]["failed"] == 1
    assert rates["screened"] == 3
    assert rates["hits"] == 1
    assert rates["hit_rate"] == pytest.approx(1 / 3)
    candidates = pd.read_csv(summary_dir / "plan" / "candidates.csv")
    assert candidates["uniprot_id"].tolist() == ["P00000"]
    failed = pd.read_csv(summary_dir / "failed_structures.csv")
    assert failed["uniprot_id"].tolist() == ["P00002"]
    assert "fpocket" in failed["error"].iloc[0]
    digest = (summary_dir / "DIGEST.md").read_text()
    assert digest.startswith("## Proteome screen")
    assert summary["structures"]["unscreened"] == {}
    assert "Incomplete" not in digest


def test_models_a_shard_never_reached_are_counted(tmp_path, monkeypatch):
    """A shard that timed out leaves models neither screened nor failed."""
    import cryptic_ip.analysis as analysis

    monkeypatch.setattr(analysis, "ProteinAnalyzer", _StubAnalyzer, raising=False)
    catalog = _catalog(tmp_path, n=4)
    shards, summary_dir = tmp_path / "shards", tmp_path / "summary"
    # Only shard 0 of 2 runs.
    assert proteome_screen.main(
        ["screen", "--organism", "yeast", "--catalog", str(catalog), "--shard-index", "0",
         "--shard-count", "2", "--workers", "1", "--output-dir", str(shards)]
    ) == 0
    assert proteome_screen.main(
        ["aggregate", "--shards-dir", str(shards), "--catalog-dir", str(catalog.parent),
         "--output-dir", str(summary_dir), "--keep-screen-scores"]
    ) == 0
    summary = json.loads((summary_dir / "screen_summary.json").read_text())
    assert summary["structures"]["unscreened"] == {"yeast": 2}
    screened = pd.read_csv(next(shards.glob("*_status.csv")))["uniprot_id"]
    unscreened = pd.read_csv(summary_dir / "unscreened_structures.csv")["uniprot_id"]
    assert set(screened).isdisjoint(unscreened) and len(unscreened) == 2
    assert "**Incomplete:**" in (summary_dir / "DIGEST.md").read_text()


def _catalog(tmp_path, n=3):
    models = tmp_path / "models"
    models.mkdir()
    for i in range(n):
        _write(models, f"AF-P0000{i}-F1-model_v4.pdb", _model_text(300 + i))
    catalog_dir = tmp_path / "catalog"
    assert proteome_screen.main(
        ["catalog", "--organism", "yeast", "--structures-dir", str(models), "--output-dir", str(catalog_dir)]
    ) == 0
    return catalog_dir / "yeast_catalog.csv"


def test_slow_models_are_named_in_the_log(tmp_path, monkeypatch, capsys):
    import cryptic_ip.analysis as analysis

    monkeypatch.setattr(analysis, "ProteinAnalyzer", _StubAnalyzer, raising=False)
    monkeypatch.setattr(proteome_screen, "SLOW_S", 0.0)
    catalog = _catalog(tmp_path)
    assert proteome_screen.main(
        ["screen", "--organism", "yeast", "--catalog", str(catalog), "--workers", "1",
         "--output-dir", str(tmp_path / "shards")]
    ) == 0
    log = capsys.readouterr().out
    assert "slow model P00000" in log and "(fpocket " in log
    assert "slow model P00002" in log and "fpocket produced no output" in log


def test_parallel_screen_records_every_model(tmp_path, capsys):
    """The process-pool path: every model returns a status, failures included.

    Without fpocket (or with these toy models) each structure fails, which is
    what is checked: a failed worker is a recorded error, never a lost model.
    """
    catalog = _catalog(tmp_path)
    shards = tmp_path / "shards"
    assert proteome_screen.main(
        ["screen", "--organism", "yeast", "--catalog", str(catalog), "--workers", "2",
         "--output-dir", str(shards)]
    ) == 0
    status = pd.read_csv(next(shards.glob("*_status.csv")))
    assert sorted(status["uniprot_id"]) == ["P00000", "P00001", "P00002"]
    assert "3/3" in capsys.readouterr().out


def test_aggregate_rescores_pockets_with_the_current_scorer(tmp_path, monkeypatch):
    """Stored scores are replaced by the current scorer's, and kept alongside."""
    import cryptic_ip.analysis as analysis
    from cryptic_ip.analysis.scorer import PocketScorer

    monkeypatch.setattr(analysis, "ProteinAnalyzer", _StubAnalyzer, raising=False)
    catalog = _catalog(tmp_path, n=3)
    shards, summary_dir = tmp_path / "shards", tmp_path / "summary"
    assert proteome_screen.main(
        ["screen", "--organism", "yeast", "--catalog", str(catalog), "--workers", "1",
         "--output-dir", str(shards)]
    ) == 0
    stored = pd.concat(pd.read_csv(p) for p in shards.glob("*_pockets_part*.csv.gz"))
    expected = PocketScorer().score_frame(stored)
    assert not (expected == stored["composite_score"].to_numpy()).all()  # stale by construction

    assert proteome_screen.main(
        ["aggregate", "--shards-dir", str(shards), "--catalog-dir", str(catalog.parent),
         "--output-dir", str(summary_dir)]
    ) == 0
    summary = json.loads((summary_dir / "screen_summary.json").read_text())
    assert summary["scoring"]["rescored"] is True
    assert summary["scoring"]["parameters"]["depth_measure"] == PocketScorer().parameters.depth_measure
    assert "recomputed at aggregation" in (summary_dir / "DIGEST.md").read_text()
    proteins = pd.read_csv(summary_dir / "plan" / "proteins.csv.gz")
    best = stored.assign(rescored=expected).groupby("uniprot_id")["rescored"].max()
    # The protein table is built from the rescored pockets.
    for uniprot_id, value in best.items():
        row = proteins[proteins["uniprot_id"] == uniprot_id]
        assert not row.empty
        assert row["best_score"].iloc[0] == pytest.approx(value)
