"""Tests for pocket detection via fpocket integration."""

import pytest
from pathlib import Path

from cryptic_ip.analysis import ProteinAnalyzer


@pytest.fixture
def adar2_path():
    path = Path("data/validation/1ZY7.pdb")
    if not path.exists():
        pytest.skip("ADAR2 structure not available")
    return path


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
