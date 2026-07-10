"""Tests for pocket detection via fpocket integration."""

import shutil

import pytest
from pathlib import Path

from cryptic_ip.analysis import ProteinAnalyzer


@pytest.fixture
def adar2_path(adar2_structure_path):
    if shutil.which("fpocket") is None:
        pytest.skip("fpocket not installed")
    return adar2_structure_path


class TestPocketDetection:
    def test_detect_pockets_adar2(self, adar2_path):
        analyzer = ProteinAnalyzer(str(adar2_path), skip_electrostatics=True)
        pockets = analyzer.detect_pockets()
        assert len(pockets) > 0
        assert "pocket_id" in pockets.columns
        assert "volume" in pockets.columns

    def test_score_pockets(self, adar2_path):
        analyzer = ProteinAnalyzer(str(adar2_path), skip_electrostatics=True)
        scored = analyzer.run_pipeline(include_electrostatics=False)
        assert not scored.empty
        assert scored["composite_score"].max() > 0
