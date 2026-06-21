"""Tests for database API clients."""

import pytest
import requests
from pathlib import Path

from cryptic_ip.database.alphafold_client import AlphaFoldClient
from cryptic_ip.database.pdb_client import PDBClient
from cryptic_ip.database.uniprot_client import UniProtClient

pytestmark = pytest.mark.requires_network


class TestAlphaFoldClient:
    """Test AlphaFold API client."""

    def test_init(self, tmp_path):
        client = AlphaFoldClient(cache_dir=tmp_path)
        assert client.cache_dir == tmp_path
        assert client.cache_dir.exists()

    def test_fetch_structure_adar2(self, tmp_path):
        client = AlphaFoldClient(cache_dir=tmp_path)
        structure_path = client.fetch_structure("P78563")
        assert structure_path.exists()
        assert structure_path.stat().st_size > 10000
        assert structure_path.name.startswith("AF-P78563-F1-model_v")

    def test_get_metadata(self):
        client = AlphaFoldClient()
        metadata = client.get_metadata("P78563")
        assert metadata["uniprot_id"] == "P78563"
        assert metadata["organism"].lower() == "homo sapiens"
        assert metadata["sequence_length"] > 0
        assert metadata["global_metric_value"] > 0

    def test_invalid_uniprot_id(self, tmp_path):
        client = AlphaFoldClient(cache_dir=tmp_path)
        with pytest.raises((ValueError, requests.HTTPError)):
            client.fetch_structure("INVALID123")


class TestPDBClient:
    """Test PDB API client."""

    def test_fetch_structure_1zy7(self, tmp_path):
        client = PDBClient(cache_dir=tmp_path)
        structure_path = client.fetch_structure("1ZY7")
        assert structure_path.exists()
        assert structure_path.stat().st_size > 10000
        assert structure_path.name == "1ZY7.pdb"

    def test_get_entry_info(self):
        client = PDBClient()
        info = client.get_entry_info("1ZY7")
        assert info["pdb_id"] == "1ZY7"
        assert info["resolution"] is not None
        assert "X-RAY" in info["method"].upper()
        assert "ADAR" in info["title"].upper() or "DEAMINASE" in info["title"].upper()

    def test_search_by_uniprot(self):
        client = PDBClient()
        pdb_ids = client.search_by_uniprot("P78563")
        if len(pdb_ids) == 0:
            pytest.skip("RCSB search returned no PDB IDs for P78563")
        assert "1ZY7" in [pid.upper() for pid in pdb_ids]


class TestUniProtClient:
    """Test UniProt API client."""

    def test_get_protein_info(self):
        client = UniProtClient()
        info = client.get_protein_info("P78563")
        assert info["accession"] == "P78563"
        assert info["gene_name"].upper() == "ADARB1"
        assert "editase" in info["protein_name"].lower() or "deaminase" in info["protein_name"].lower()
        assert info["organism"] == "Homo sapiens"
        assert info["sequence_length"] > 700

    def test_get_sequence(self):
        client = UniProtClient()
        sequence = client.get_sequence("P78563")
        assert sequence.startswith(">")
        assert "P78563" in sequence
        assert len(sequence) > 700

    def test_search_by_gene(self):
        client = UniProtClient()
        results = client.search_by_gene("ADARB1", organism="Homo sapiens")
        assert len(results) > 0
        assert "P78563" in results


def test_complete_data_retrieval_workflow(tmp_path):
    """Test complete workflow: UniProt -> AlphaFold + PDB."""
    uniprot_client = UniProtClient()
    protein_info = uniprot_client.get_protein_info("P78563")
    assert protein_info["accession"] == "P78563"

    af_client = AlphaFoldClient(cache_dir=tmp_path / "alphafold")
    af_structure = af_client.fetch_structure("P78563")
    assert af_structure.exists()

    pdb_client = PDBClient(cache_dir=tmp_path / "pdb")
    pdb_structure = pdb_client.fetch_structure("1ZY7")
    assert pdb_structure.exists()
