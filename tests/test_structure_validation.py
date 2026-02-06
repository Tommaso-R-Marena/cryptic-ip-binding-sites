"""Tests for structure validation."""

import pytest
from pathlib import Path
import tempfile
from Bio import PDB
import requests


class TestStructureDownload:
    """Test downloading real structures."""
    
    def test_download_adar2_alphafold(self):
        """Download ADAR2 from AlphaFold."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            url = 'https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb'
            output = tmpdir / 'adar2.pdb'
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            output.write_bytes(response.content)
            
            assert output.exists()
            assert output.stat().st_size > 10000
    
    def test_download_adar2_crystal(self):
        """Download ADAR2 crystal from PDB."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            url = 'https://files.rcsb.org/download/1ZY7.pdb'
            output = tmpdir / '1zy7.pdb'
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            output.write_bytes(response.content)
            
            assert output.exists()
            assert output.stat().st_size > 10000


class TestStructureParsing:
    """Test parsing real structures with BioPython."""
    
    @pytest.fixture
    def adar2_structure(self):
        """Download and parse ADAR2 AlphaFold structure."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Download
            url = 'https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb'
            pdb_file = tmpdir / 'adar2.pdb'
            
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            pdb_file.write_bytes(response.content)
            
            # Parse
            parser = PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('ADAR2', str(pdb_file))
            
            yield structure
    
    def test_structure_has_residues(self, adar2_structure):
        """Check structure has amino acid residues."""
        model = adar2_structure[0]
        residues = [r for r in model.get_residues() if PDB.is_aa(r)]
        
        assert len(residues) > 700  # ADAR2 is ~700+ residues
        assert len(residues) < 800  # Sanity check
    
    def test_structure_has_atoms(self, adar2_structure):
        """Check structure has atoms."""
        model = adar2_structure[0]
        atoms = list(model.get_atoms())
        
        assert len(atoms) > 5000  # Should have thousands of atoms
    
    def test_plddt_scores_exist(self, adar2_structure):
        """Check pLDDT scores in B-factor column."""
        model = adar2_structure[0]
        atoms = list(model.get_atoms())
        
        b_factors = [atom.bfactor for atom in atoms]
        
        assert len(b_factors) > 0
        assert all(0 <= bf <= 100 for bf in b_factors)  # pLDDT range
        
        avg_plddt = sum(b_factors) / len(b_factors)
        assert avg_plddt > 50  # ADAR2 should have decent confidence
    
    def test_find_basic_residues(self, adar2_structure):
        """Find basic residues (Arg, Lys, His)."""
        model = adar2_structure[0]
        residues = [r for r in model.get_residues() if PDB.is_aa(r)]
        
        basic_residues = [r for r in residues if r.resname in ['ARG', 'LYS', 'HIS']]
        
        assert len(basic_residues) > 50  # Should have many basic residues
        
        # Check known IP6-binding residues exist
        known_ip6_residues = [376, 519, 522, 651, 672]
        residue_numbers = [r.id[1] for r in residues]
        
        for res_num in known_ip6_residues:
            assert res_num in residue_numbers, f"Missing known IP6 residue {res_num}"


def test_api_connectivity():
    """Test that all required APIs are accessible."""
    
    # Test AlphaFold API
    response = requests.get(
        'https://alphafold.ebi.ac.uk/api/prediction/P78563',
        timeout=10
    )
    assert response.status_code == 200, "AlphaFold API not accessible"
    
    # Test PDB API
    response = requests.get(
        'https://data.rcsb.org/rest/v1/core/entry/1ZY7',
        timeout=10
    )
    assert response.status_code == 200, "PDB API not accessible"
    
    # Test UniProt API
    response = requests.get(
        'https://rest.uniprot.org/uniprotkb/P78563',
        timeout=10
    )
    assert response.status_code == 200, "UniProt API not accessible"
    
    print("\n✓ All APIs are accessible!")
