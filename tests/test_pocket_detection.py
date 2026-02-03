"""Tests for pocket detection module."""

import pytest
from pathlib import Path
from cryptic_ip.analysis.pocket_detection import PocketDetector


class TestPocketDetector:
    """Test suite for PocketDetector class."""

    def test_initialization(self):
        """Test PocketDetector initialization."""
        detector = PocketDetector()
        assert detector is not None
        assert detector.fpocket_path is not None

    def test_detect_pockets_adar2(self, adar2_pdb, temp_output_dir):
        """Test pocket detection on ADAR2 structure."""
        if not adar2_pdb.exists():
            pytest.skip("ADAR2 structure not available")

        detector = PocketDetector()
        pockets = detector.detect_pockets(
            str(adar2_pdb),
            output_dir=str(temp_output_dir)
        )

        assert pockets is not None
        assert len(pockets) > 0
        # ADAR2 should have multiple pockets detected
        assert len(pockets) >= 3

    def test_pocket_volume_calculation(self, adar2_pdb, temp_output_dir):
        """Test that pocket volumes are calculated correctly."""
        if not adar2_pdb.exists():
            pytest.skip("ADAR2 structure not available")

        detector = PocketDetector()
        pockets = detector.detect_pockets(
            str(adar2_pdb),
            output_dir=str(temp_output_dir)
        )

        for pocket in pockets:
            assert 'volume' in pocket
            assert pocket['volume'] > 0
            # Reasonable volume range for protein pockets
            assert pocket['volume'] < 5000

    def test_pocket_depth_calculation(self, adar2_pdb, temp_output_dir):
        """Test pocket depth metrics."""
        if not adar2_pdb.exists():
            pytest.skip("ADAR2 structure not available")

        detector = PocketDetector()
        pockets = detector.detect_pockets(
            str(adar2_pdb),
            output_dir=str(temp_output_dir)
        )

        for pocket in pockets:
            assert 'depth' in pocket
            assert pocket['depth'] >= 0

    def test_invalid_pdb_handling(self, temp_output_dir):
        """Test handling of invalid PDB file."""
        invalid_pdb = temp_output_dir / "invalid.pdb"
        invalid_pdb.write_text("INVALID PDB CONTENT")

        detector = PocketDetector()
        with pytest.raises(Exception):
            detector.detect_pockets(
                str(invalid_pdb),
                output_dir=str(temp_output_dir)
            )

    def test_basic_residue_identification(self, adar2_pdb, temp_output_dir):
        """Test identification of basic residues in pockets."""
        if not adar2_pdb.exists():
            pytest.skip("ADAR2 structure not available")

        detector = PocketDetector()
        pockets = detector.detect_pockets(
            str(adar2_pdb),
            output_dir=str(temp_output_dir)
        )

        # At least one pocket should have basic residues
        found_basic_residues = False
        for pocket in pockets:
            if 'basic_residues' in pocket and pocket['basic_residues'] > 0:
                found_basic_residues = True
                break

        assert found_basic_residues, "No pockets with basic residues detected"
