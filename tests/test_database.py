"""Tests for database management module."""

import pytest
from pathlib import Path
from cryptic_ip.database.manager import ProteomeDatabase


class TestProteomeDatabase:
    """Test suite for ProteomeDatabase class."""

    def test_initialization(self, temp_output_dir):
        """Test database initialization."""
        db = ProteomeDatabase(db_path=str(temp_output_dir / "test.db"))
        assert db is not None

    def test_structure_download_url_generation(self):
        """Test AlphaFold URL generation."""
        db = ProteomeDatabase()
        url = db.get_alphafold_url("P78563")  # ADAR2 UniProt ID
        assert "alphafold" in url.lower()
        assert "P78563" in url

    def test_proteome_metadata(self):
        """Test proteome metadata retrieval."""
        db = ProteomeDatabase()

        yeast_meta = db.get_proteome_metadata("yeast")
        assert yeast_meta is not None
        assert yeast_meta['uniprot_id'] == "UP000002311"
        assert yeast_meta['organism'] == "Saccharomyces cerevisiae"

        human_meta = db.get_proteome_metadata("human")
        assert human_meta is not None
        assert human_meta['uniprot_id'] == "UP000005640"

    def test_structure_registration(self, temp_output_dir, mock_proteome_structures):
        """Test registering structures in database."""
        db = ProteomeDatabase(db_path=str(temp_output_dir / "test.db"))

        for i, structure in enumerate(mock_proteome_structures):
            db.register_structure(
                uniprot_id=f"P0000{i}",
                organism="test",
                file_path=str(structure),
                avg_plddt=85.0
            )

        # Verify structures were registered
        structures = db.get_structures_by_organism("test")
        assert len(structures) == len(mock_proteome_structures)

    def test_quality_metrics_storage(self, temp_output_dir):
        """Test storing quality metrics."""
        db = ProteomeDatabase(db_path=str(temp_output_dir / "test.db"))

        db.register_structure(
            uniprot_id="TEST001",
            organism="test",
            file_path="/fake/path.pdb",
            avg_plddt=87.5,
            protein_length=342
        )

        structure = db.get_structure("TEST001")
        assert structure is not None
        assert structure['avg_plddt'] == 87.5
        assert structure['protein_length'] == 342
