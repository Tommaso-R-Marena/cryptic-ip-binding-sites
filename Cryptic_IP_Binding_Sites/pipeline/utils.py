import logging
import yaml
from pathlib import Path
from typing import Dict, Any, List
from pydantic import BaseModel, Field, ValidationError


class OrganismConfig(BaseModel):
    name: str = Field(..., description="Target organism name")
    target_count: int = Field(..., gt=0)


class ScoringWeights(BaseModel):
    depth: float = Field(..., ge=0)
    sasa: float = Field(..., ge=0)
    electrostatic: float = Field(..., ge=0)
    basic_residues: float = Field(..., ge=0)


class PipelineSettings(BaseModel):
    fpocket: Dict[str, Any]
    freesasa: Dict[str, Any]
    apbs: Dict[str, Any]
    plddt: Dict[str, float]
    scoring: Dict[str, Any]


class FullConfigSchema(BaseModel):
    pipeline: PipelineSettings
    organisms: List[OrganismConfig]


def setup_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Sets up a logger with standard formatting."""
    logger = logging.getLogger(name)
    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(level)
    return logger


def load_config() -> Dict[str, Any]:
    """Loads and cryptographically validates the main configuration from config/config.yaml."""
    config_path = Path(__file__).parent.parent / "config" / "config.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found at {config_path}")
    with open(config_path, "r") as f:
        raw_config = yaml.safe_load(f)

    try:
        # Strictly validate the structure, types, and mathematical boundaries using Pydantic
        FullConfigSchema(**raw_config)
    except ValidationError as e:
        raise ValueError(f"CRITICAL: Configuration Schema Validation Failed!\n{e}")

    return raw_config


def load_controls(control_type: str) -> Dict[str, Any]:
    """Loads validation controls (positive or negative)."""
    if control_type not in ["positive", "negative"]:
        raise ValueError("control_type must be 'positive' or 'negative'")
    control_path = (
        Path(__file__).parent.parent / "validation" / f"{control_type}_controls.yaml"
    )
    if not control_path.exists():
        raise FileNotFoundError(f"Control file not found at {control_path}")
    with open(control_path, "r") as f:
        return yaml.safe_load(f)["controls"]
