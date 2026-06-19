"""Structure parsing and ligand-context helpers for validation workflows."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Set, Tuple

from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.SASA import ShrakeRupley

LIGAND_RESNAMES = {
    "IP3",
    "IP4",
    "IP5",
    "IP6",
    "IHP",
    "I3P",
    "4IP",
    "6A0",
    "INS",
}
BASIC_RESNAMES = {"ARG", "LYS", "HIS"}


def load_structure(path: Path):
    """Load a PDB or mmCIF structure file."""
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def parse_ligand_resnames(path: Path, target_ligands: Optional[Set[str]] = None) -> Set[str]:
    """Return ligand residue names present in a structure file."""
    allowed = target_ligands or LIGAND_RESNAMES
    found: Set[str] = set()
    structure = load_structure(path)
    for residue in structure.get_residues():
        resname = residue.get_resname().upper()
        if resname in allowed:
            found.add(resname)
    return found


def ligand_context(path: Path) -> Tuple[Set[int], Optional[float], Optional[Tuple[float, float, float]]]:
    """Return coordinating basic residues, minimum ligand SASA, and ligand centroid."""
    structure = load_structure(path)
    ShrakeRupley().compute(structure, level="R")

    ligand_atoms = []
    ligand_sasa_values = []
    for residue in structure.get_residues():
        if residue.get_resname() not in LIGAND_RESNAMES:
            continue
        residue_sasa = float(getattr(residue, "sasa", 0.0))
        ligand_sasa_values.append(residue_sasa)
        ligand_atoms.extend(list(residue.get_atoms()))

    if not ligand_atoms:
        return set(), None, None

    from Bio.PDB import NeighborSearch

    coords = [atom.coord for atom in ligand_atoms]
    xs, ys, zs = zip(*coords)
    centroid = (float(sum(xs) / len(xs)), float(sum(ys) / len(ys)), float(sum(zs) / len(zs)))

    neighbor_search = NeighborSearch(list(structure.get_atoms()))
    coordinating_residues: Set[int] = set()
    for atom in ligand_atoms:
        for neighbor in neighbor_search.search(atom.coord, 5.0):
            residue = neighbor.get_parent()
            if residue.get_resname() in BASIC_RESNAMES:
                coordinating_residues.add(int(residue.id[1]))

    return coordinating_residues, float(min(ligand_sasa_values)), centroid


def ligand_total_sasa(path: Path, ligand_id: str) -> float:
    """Compute total SASA for all copies of a ligand residue name."""
    structure = load_structure(path)
    ShrakeRupley().compute(structure, level="R")

    total = 0.0
    found = False
    for residue in structure.get_residues():
        if residue.get_resname() == ligand_id:
            found = True
            total += float(getattr(residue, "sasa", 0.0))

    if not found:
        raise ValueError(f"Ligand {ligand_id} not found in {path.name}")

    return total
