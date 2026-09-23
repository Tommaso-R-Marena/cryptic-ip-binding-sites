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

    def run_pipeline(self, include_electrostatics=False):
        if "P00002" in self.pdb_path:
            raise RuntimeError("fpocket produced no output")
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

    assert proteome_screen.main(
        ["aggregate", "--shards-dir", str(shards_dir), "--catalog-dir", str(catalog_dir),
         "--output-dir", str(summary_dir)]
    ) == 0
    summary = json.loads((summary_dir / "screen_summary.json").read_text())
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
    assert (summary_dir / "DIGEST.md").read_text().startswith("## Proteome screen")
