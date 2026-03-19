import pytest
import json
from pathlib import Path
from pipeline.prepare_structure import clean_structure_and_extract_plddt

def test_plddt_extraction(tmp_path):
    # Create a mock PDB file
    mock_pdb = tmp_path / "mock.pdb"
    with open(mock_pdb, "w") as f:
        # ATOM      1  N   MET A   1      27.340  24.430   2.314  1.00 95.00           N  
        f.write("ATOM      1  N   MET A   1      27.340  24.430   2.314  1.00 95.00           N  \n")
        # ATOM      2  CA  MET A   1      28.240  25.430   2.314  1.00 90.00           C  
        f.write("ATOM      2  CA  MET A   1      28.240  25.430   2.314  1.00 90.50           C  \n")
        # HETATM 1000  P1  IHP B 200      10.000  10.000  10.000  1.00 20.00           P
        f.write("HETATM 1000  P1  IHP B 200      10.000  10.000  10.000  1.00 20.00           P  \n")
        f.write("TER\n")
        
    cleaned, plddt = clean_structure_and_extract_plddt(mock_pdb, tmp_path)
    
    assert cleaned.exists()
    assert plddt.exists()
    
    with open(plddt, "r") as f:
        data = json.load(f)
        
    assert "A_1" in data # We parsed "1" id, chain "A"
    assert data["A_1"] == 95.00 # First atom of residue is kept in our parser, which is fine

    with open(cleaned, "r") as f:
        content = f.read()
        assert "MET" in content
        assert "IHP" in content # Kept because it's in the IP list
