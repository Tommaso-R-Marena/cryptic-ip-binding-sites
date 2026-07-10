"""Helpers for end-to-end pipeline integration tests."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
VALIDATION_DIR = ROOT / "data" / "validation"
DATASET_CSV = VALIDATION_DIR / "ip_binding_validation_dataset.csv"
BUNDLED_FEATURES = ROOT / "results" / "publication" / "ml_training" / "validation_pocket_features.csv"


def seed_ml_structures(raw_dir: Path) -> int:
    """Copy any cached validation PDBs referenced by the dataset into raw/structures."""
    structures_dir = raw_dir / "structures"
    structures_dir.mkdir(parents=True, exist_ok=True)
    if not DATASET_CSV.exists():
        return 0

    dataset = pd.read_csv(DATASET_CSV)
    copied = 0
    for pdb_id in dataset["pdb_id"].astype(str).str.upper().unique():
        src = VALIDATION_DIR / f"{pdb_id}.pdb"
        if src.exists():
            shutil.copy2(src, structures_dir / f"{pdb_id}.pdb")
            copied += 1
    return copied


def seed_yeast_structures(structures_dir: Path, *, n_proteins: int = 2) -> list[Path]:
    """Seed a tiny yeast pilot set from local validation caches (no network)."""
    structures_dir.mkdir(parents=True, exist_ok=True)
    candidates = [
        VALIDATION_DIR / "1MAI.pdb",
        VALIDATION_DIR / "AF-P78563-F1-model_v6.pdb",
        VALIDATION_DIR / "1ZY7.pdb",
    ]
    sources = [path for path in candidates if path.exists()]
    if not sources:
        raise FileNotFoundError("No local structures available to seed yeast pilot")

    uniprot_ids = ["P39968", "P25567", "P38972", "P32324", "P25367"]
    written: list[Path] = []
    for index in range(n_proteins):
        src = sources[index % len(sources)]
        uniprot_id = uniprot_ids[index % len(uniprot_ids)]
        dest = structures_dir / f"AF-{uniprot_id}-F1-model_v4.pdb"
        shutil.copy2(src, dest)
        written.append(dest)
    return written


def resolve_ml_features_csv(raw_dir: Path) -> Path | None:
    """Return bundled pocket features when local structure coverage is too sparse for training."""
    copied = seed_ml_structures(raw_dir)
    if copied >= 10:
        return None
    if BUNDLED_FEATURES.exists():
        return BUNDLED_FEATURES
    return None


def ensure_tier1_structures() -> None:
    """Ensure tier-1 control structures exist under data/validation."""
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)
    required = ("1ZY7", "1MAI", "1BWN", "5HDT", "5ICN")
    missing = [pdb_id for pdb_id in required if not (VALIDATION_DIR / f"{pdb_id}.pdb").exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing tier-1 structures in {VALIDATION_DIR}: {', '.join(missing)}"
        )
