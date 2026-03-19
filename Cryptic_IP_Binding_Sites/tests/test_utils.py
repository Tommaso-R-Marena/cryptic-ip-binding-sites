import pytest
from pipeline.utils import setup_logger, load_config, load_controls

def test_setup_logger():
    logger = setup_logger("test_logger")
    assert logger.name == "test_logger"
    assert len(logger.handlers) > 0

def test_load_config():
    config = load_config()
    assert isinstance(config, dict)
    assert "pipeline" in config

def test_load_controls_positive():
    controls = load_controls("positive")
    assert isinstance(controls, dict)
    assert len(controls) > 0

def test_load_controls_invalid():
    with pytest.raises(ValueError):
        load_controls("invalid_type")
