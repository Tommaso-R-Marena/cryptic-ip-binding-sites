import pytest
from pipeline.score_pockets import compute_composite_score, calculate_basic_residue_count
from pipeline.utils import load_config

def test_composite_score():
    config = load_config()
    
    # Mock a perfect pocket scoring 1.0
    perfect_pocket = {
        "depth": 30.0,
        "mean_sasa": 0.0,
        "electrostatic_potential": 15.0,
        "basic_residue_count": 8
    }
    
    score = compute_composite_score(perfect_pocket, config)
    assert score == 1.0
    
    # Mock a terrible pocket scoring 0.0
    terrible_pocket = {
        "depth": 0.0,
        "mean_sasa": 50.0,
        "electrostatic_potential": 0.0,
        "basic_residue_count": 0
    }
    
    score = compute_composite_score(terrible_pocket, config)
    assert score == 0.0
    
    # Mock an intermediate pocket
    mid_pocket = {
        "depth": 15.0, # 0.5 normalized
        "mean_sasa": 25.0, # 0.5 normalized
        "electrostatic_potential": 7.5, # 0.5 normalized
        "basic_residue_count": 4 # 0.5 normalized
    }
    
    score = compute_composite_score(mid_pocket, config)
    assert 0.49 <= score <= 0.51

def test_calculate_basic_residue_count(tmp_path):
    mock_pdb = tmp_path / "mock.pdb"
    with open(mock_pdb, "w") as f:
        # Distance to center (0,0,0) varies
        f.write("ATOM      1  CA  ARG A   1       1.000   1.000   1.000  1.00 50.00           C  \n") # d=sqrt(3) ~ 1.7 (count)
        f.write("ATOM      2  CA  LYS A   2       2.000   2.000   2.000  1.00 50.00           C  \n") # d=sqrt(12) ~ 3.4 (count)
        f.write("ATOM      3  CA  HIS B   1      10.000  10.000  10.000  1.00 50.00           C  \n") # d=sqrt(300) > 5.0 (no count)
        f.write("ATOM      4  CA  ALA A   3       1.000   1.000   1.000  1.00 50.00           C  \n") # not basic (no count)
        f.write("TER\n")
        
    count = calculate_basic_residue_count((0.0, 0.0, 0.0), mock_pdb, distance_cutoff=5.0)
    assert count == 2

