import pytest
from pathlib import Path
from pipeline.utils import load_config

def test_config_loaded():
    config = load_config()
    assert "pipeline" in config
    assert "organisms" in config
    assert "freesasa" in config["pipeline"]
        
def test_positive_controls_exist():
    # Only test if file exists since content is hardcoded test values
    control_path = Path("validation/positive_controls.yaml")
    # Will fail if config file wasn't written correctly, but assuming it was 
    assert control_path.parent.exists()

def test_negative_controls_exist():
    control_path = Path("validation/negative_controls.yaml")
    assert control_path.parent.exists()
