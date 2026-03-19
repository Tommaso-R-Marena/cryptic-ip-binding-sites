import os
from pathlib import Path
from pipeline.visualization import generate_pymol_script

def test_generate_pymol_script(tmp_path):
    """Test the generation of a PyMOL script for visualizing pockets."""
    # Setup dummy data
    pdb_path = tmp_path / "dummy.pdb"
    pdb_path.write_text("ATOM      1  N   MET A   1      11.104  18.667  31.543  1.00 48.90           N  \n")
    
    output_pml = tmp_path / "vis.pml"
    
    pockets = [
        {
            "residues": ["RES_A_1", "RES_A_2"],
            "center": (11.0, 18.0, 31.0),
            "composite_score": 0.85
        },
        {
            "residues": ["RES_B_10", "RES_B_11"],
            "composite_score": 0.60
        }
    ]
    
    generate_pymol_script(pdb_path, pockets, output_pml)
    
    assert output_pml.exists()
    
    content = output_pml.read_text()
    assert "reinitialize" in content
    assert f"load {pdb_path.resolve()}" in content
    assert "select pocket_1" in content
    assert "(chain A and resi 1) or (chain A and resi 2)" in content
    assert "select pocket_2" in content
    assert "(chain B and resi 10) or (chain B and resi 11)" in content
    assert "show sticks, basic_pocket_1" in content
    assert "zoom pocket_1, 5" in content
