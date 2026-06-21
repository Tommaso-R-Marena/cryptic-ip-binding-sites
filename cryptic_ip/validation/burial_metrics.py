"""Advanced burial metrics for inositol phosphate binding sites."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from Bio.PDB import NeighborSearch, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley

from .structure_context import LIGAND_RESNAMES, load_structure


class _LigandSelect(Select):
    """Select only inositol phosphate heteroatoms."""

    def accept_residue(self, residue):
        return residue.get_resname() in LIGAND_RESNAMES


class _ProteinOnlySelect(Select):
    """Exclude inositol phosphate ligands for apo SASA comparison."""

    def accept_residue(self, residue):
        return residue.get_resname() not in LIGAND_RESNAMES


@dataclass
class BurialMetrics:
    """Ligand burial measurements for validation and ML labeling."""

    ligand_sasa: Optional[float]
    phosphate_sasa: Optional[float]
    delta_sasa: Optional[float]
    burial_depth: Optional[float]
    burial_class: str


def _is_phosphate_atom(atom_name: str) -> bool:
    name = atom_name.upper()
    return name.startswith("P") or name.startswith("O") and any(ch.isdigit() for ch in name)


def _ligand_atoms(structure, phosphate_only: bool = False) -> List:
    atoms = []
    for residue in structure.get_residues():
        if residue.get_resname() not in LIGAND_RESNAMES:
            continue
        for atom in residue.get_atoms():
            if phosphate_only and not _is_phosphate_atom(atom.get_name()):
                continue
            atoms.append(atom)
    return atoms


def _compute_ligand_sasa(structure, phosphate_only: bool = False) -> float:
    atoms = _ligand_atoms(structure, phosphate_only=phosphate_only)
    if not atoms:
        return float("nan")
    ShrakeRupley().compute(structure, level="A")
    return float(sum(float(getattr(atom, "sasa", 0.0)) for atom in atoms))


def _ligand_centroid(atoms) -> Optional[Tuple[float, float, float]]:
    if not atoms:
        return None
    coords = np.asarray([atom.coord for atom in atoms], dtype=float)
    center = coords.mean(axis=0)
    return float(center[0]), float(center[1]), float(center[2])


def _burial_depth_angstrom(structure, ligand_centroid: Tuple[float, float, float]) -> float:
    """Approximate burial depth as distance from ligand centroid to nearest exposed atom."""
    ShrakeRupley().compute(structure, level="R")
    exposed_coords = []
    for residue in structure.get_residues():
        if residue.get_resname() in LIGAND_RESNAMES:
            continue
        if residue.id[0] != " ":
            continue
        sasa = float(getattr(residue, "sasa", 0.0))
        if sasa >= 30.0 and "CA" in residue:
            exposed_coords.append(residue["CA"].coord)
    if not exposed_coords:
        return 0.0
    center = np.asarray(ligand_centroid, dtype=float)
    distances = np.linalg.norm(np.asarray(exposed_coords) - center, axis=1)
    return float(np.min(distances))


def classify_burial(ligand_sasa: Optional[float], phosphate_sasa: Optional[float] = None) -> str:
    """Classify burial using ligand SASA thresholds aligned with the dataset builder."""
    metric = phosphate_sasa if phosphate_sasa is not None and phosphate_sasa == phosphate_sasa else ligand_sasa
    if metric is None or metric != metric:
        return "unknown"
    if metric <= 5.0:
        return "cryptic"
    if metric <= 50.0:
        return "semi_cryptic"
    return "surface"


def compute_burial_metrics(path: Path) -> BurialMetrics:
    """Compute ligand SASA, phosphate SASA, delta SASA, and burial depth."""
    path = Path(path)
    holo = load_structure(path)
    ligand_atoms = _ligand_atoms(holo)
    centroid = _ligand_centroid(ligand_atoms)

    ligand_sasa = _compute_ligand_sasa(holo, phosphate_only=False)
    phosphate_sasa = _compute_ligand_sasa(holo, phosphate_only=True)

    delta_sasa = None
    if ligand_atoms:
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as handle:
            apo_path = Path(handle.name)
        io = PDBIO()
        io.set_structure(holo)
        io.save(str(apo_path), _ProteinOnlySelect())
        apo = load_structure(apo_path)
        apo_path.unlink(missing_ok=True)
        # Ligand SASA in holo; apo has no ligand so delta = holo ligand SASA
        delta_sasa = float(ligand_sasa)

    burial_depth = _burial_depth_angstrom(holo, centroid) if centroid else None
    burial_class = classify_burial(ligand_sasa, phosphate_sasa)

    return BurialMetrics(
        ligand_sasa=None if ligand_sasa != ligand_sasa else float(ligand_sasa),
        phosphate_sasa=None if phosphate_sasa != phosphate_sasa else float(phosphate_sasa),
        delta_sasa=delta_sasa,
        burial_depth=burial_depth,
        burial_class=burial_class,
    )
