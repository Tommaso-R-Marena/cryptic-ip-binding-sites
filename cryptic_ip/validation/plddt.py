"""AlphaFold pLDDT confidence extraction and pocket-level filtering."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set

import numpy as np
from Bio.PDB import PDBParser

from .structure_context import load_structure


def extract_residue_plddt(path: Path) -> Dict[int, float]:
    """Parse per-residue pLDDT from AlphaFold B-factors (CA atoms)."""
    structure = load_structure(path)
    plddt: Dict[int, float] = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" in residue:
                    plddt[int(residue.id[1])] = float(residue["CA"].get_bfactor())
    return plddt


def mean_plddt(path: Path) -> Optional[float]:
    """Return mean CA pLDDT for a structure file."""
    values = list(extract_residue_plddt(path).values())
    if not values:
        return None
    return float(np.mean(values))


def pocket_plddt_confidence(
    path: Path,
    pocket_residues: Iterable[int],
    *,
    min_plddt: float = 70.0,
) -> Dict[str, float]:
    """Summarize pLDDT for pocket-lining residues."""
    residue_plddt = extract_residue_plddt(path)
    pocket_set: Set[int] = {int(r) for r in pocket_residues}
    values = [residue_plddt[r] for r in pocket_set if r in residue_plddt]
    if not values:
        return {
            "plddt_mean": float("nan"),
            "plddt_min": float("nan"),
            "plddt_fraction_above_cutoff": 0.0,
            "passes_confidence_gate": False,
        }
    arr = np.asarray(values, dtype=float)
    fraction = float(np.mean(arr >= min_plddt))
    return {
        "plddt_mean": float(np.mean(arr)),
        "plddt_min": float(np.min(arr)),
        "plddt_fraction_above_cutoff": fraction,
        "passes_confidence_gate": fraction >= 0.80 and float(np.min(arr)) >= min_plddt - 10.0,
    }
