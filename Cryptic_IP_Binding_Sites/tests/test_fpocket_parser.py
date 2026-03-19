import pytest
from pathlib import Path
from pipeline.run_fpocket import parse_fpocket_output

def test_parse_fpocket_output(tmp_path):
    mock_info = tmp_path / "info.txt"
    with open(mock_info, "w") as f:
        f.write("Pocket 1 :\n")
        f.write("\tScore : 0.8\n")
        f.write("\tDruggability Score : 0.9\n")
        f.write("\tVolume : 500.0\n")
        f.write("\tTotal SASA : 10.0\n")
        
        f.write("Pocket 2 :\n")
        f.write("\tScore : 0.1\n")
        f.write("\tDruggability Score : 0.2\n")
        f.write("\tVolume : 900.0\n") # Too large
        f.write("\tTotal SASA : 100.0\n")
        
        f.write("Pocket 3 :\n")
        f.write("\tScore : 0.5\n")
        f.write("\tDruggability Score : 0.5\n")
        f.write("\tVolume : 400.0\n")
        f.write("\tTotal SASA : 20.0\n")
        
    config = {
        "pipeline": {
            "fpocket": {
                "min_volume_A3": 300,
                "max_volume_A3": 800
            }
        }
    }
    
    pockets = parse_fpocket_output(mock_info, config)
    
    assert len(pockets) == 2
    assert pockets[0]["id"] == 1
    assert pockets[0]["volume"] == 500.0
    assert pockets[1]["id"] == 3
    assert pockets[1]["volume"] == 400.0
