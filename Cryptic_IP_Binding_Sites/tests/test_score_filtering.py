import pytest
import json
from pathlib import Path
from unittest.mock import patch
from pipeline.score_pockets import score_and_filter_pockets

@patch("pipeline.score_pockets.load_config")
@patch("pipeline.score_pockets.calculate_basic_residue_count")
def test_score_and_filter_pockets(mock_calc_basic, mock_load_config, tmp_path):
    # Mocking dependencies
    mock_load_config.return_value = {
        "pipeline": {
            "scoring": {
                "weights": {"depth": 0.25, "sasa": 0.25, "electrostatic": 0.25, "basic_residues": 0.25},
                "min_basic_residues": 1
            },
            "freesasa": {"max_mean_sasa_A2": 45.0},
            "apbs": {"min_electrostatic_kTe": 1.0},
            "plddt": {"min_plddt_for_pocket_residues": 70.0}
        }
    }
    
    mock_calc_basic.return_value = 2 # Always return 2 basic residues
    
    plddt_path = tmp_path / "mock_plddt.json"
    with open(plddt_path, "w") as f:
        json.dump({"A_1": 80.0, "A_2": 90.0, "B_1": 50.0}, f)
        
    cleaned_pdb_path = tmp_path / "mock_cleaned.pdb"
    cleaned_pdb_path.touch()
    
    sasa_data = {"RES_A_1": 10.0, "RES_A_2": 20.0, "RES_B_1": 80.0}
    
    pockets = [
        {
            "id": 1,
            "residues": ["RES_A_1", "RES_A_2"],
            "center": (0, 0, 0),
            "depth": 20.0,
            "electrostatic_potential": 5.0
        },
        {
            "id": 2, # Will fail filters (SASA > 45, pLDDT < 70)
            "residues": ["RES_B_1"],
            "center": (10, 10, 10),
            "depth": 5.0,
            "electrostatic_potential": 0.0
        }
    ]
    
    scored = score_and_filter_pockets(pockets, sasa_data, plddt_path, cleaned_pdb_path)
    
    assert len(scored) == 2
    
    # Pocket 1 should pass all filters and have a high score
    p1 = next(p for p in scored if p["id"] == 1)
    assert p1["pass_all_filters"] is True
    assert p1["mean_sasa"] == 15.0
    assert p1["mean_plddt"] == 85.0
    
    # Pocket 2 should fail multiple filters
    p2 = next(p for p in scored if p["id"] == 2)
    assert p2["pass_all_filters"] is False
    assert p2["mean_sasa"] == 80.0
    assert p2["mean_plddt"] == 50.0
    
    # Sorted by composite score (p1 should be first)
    assert scored[0]["id"] == 1
