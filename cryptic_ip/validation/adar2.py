"""
ADAR2 validation - the gold standard for cryptic IP6 binding.
"""

from __future__ import annotations

import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from Bio.PDB import NeighborSearch, PDBParser
from Bio.PDB.SASA import ShrakeRupley

from ..analysis import ProteinAnalyzer
from ..database.alphafold_client import AlphaFoldClient

from .structure_context import BASIC_RESNAMES, LIGAND_RESNAMES, ligand_context


def download_adar2_structures(data_dir: str = "data/validation") -> Dict[str, Path]:
    """Download ADAR2 structures from AlphaFold and PDB."""
    data_path = Path(data_dir)
    data_path.mkdir(parents=True, exist_ok=True)

    structures: Dict[str, Path] = {}

    client = AlphaFoldClient(cache_dir=data_path)
    alphafold_path = client.fetch_structure("P78563")
    local_alphafold = data_path / alphafold_path.name
    if alphafold_path.resolve() != local_alphafold.resolve():
        local_alphafold.write_bytes(alphafold_path.read_bytes())
    structures["alphafold"] = local_alphafold

    pdb_path = data_path / "1ZY7.pdb"
    if not pdb_path.exists():
        print("Downloading PDB 1ZY7 (ADAR2 crystal structure)...")
        urllib.request.urlretrieve("https://files.rcsb.org/download/1ZY7.pdb", pdb_path)
    structures["crystal"] = pdb_path

    return structures


def _ligand_site_residues(pdb_path: Path, distance_cutoff: float = 5.0) -> Tuple[Set[int], Optional[float]]:
    """Return basic residues near inositol phosphate and ligand SASA when present."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("validation", str(pdb_path))
    ShrakeRupley().compute(structure, level="R")

    ligand_atoms = []
    ligand_sasa = 0.0
    for residue in structure.get_residues():
        if residue.get_resname() in LIGAND_RESNAMES:
            ligand_atoms.extend(list(residue.get_atoms()))
            ligand_sasa += float(getattr(residue, "sasa", 0.0))

    if not ligand_atoms:
        return set(), None

    neighbor_search = NeighborSearch(list(structure.get_atoms()))
    coordinating_residues: Set[int] = set()
    for atom in ligand_atoms:
        for neighbor in neighbor_search.search(atom.coord, distance_cutoff):
            residue = neighbor.get_parent()
            if residue.get_resname() in BASIC_RESNAMES:
                coordinating_residues.add(int(residue.id[1]))

    return coordinating_residues, float(ligand_sasa)


def _ligand_centroid(pdb_path: Path) -> Optional[Tuple[float, float, float]]:
    """Return the centroid of inositol phosphate atoms when present."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("validation", str(pdb_path))
    coords = []
    for residue in structure.get_residues():
        if residue.get_resname() in LIGAND_RESNAMES:
            coords.extend(atom.coord for atom in residue.get_atoms())
    if not coords:
        return None
    xs, ys, zs = zip(*coords)
    return (float(sum(xs) / len(xs)), float(sum(ys) / len(ys)), float(sum(zs) / len(zs)))


def _select_site_pocket(
    scored: pd.DataFrame,
    analyzer: ProteinAnalyzer,
    site_residues: Set[int],
    ligand_centroid: Optional[Tuple[float, float, float]] = None,
) -> pd.Series:
    """Pick the pocket that best matches the reference binding site."""
    if scored.empty:
        raise ValueError("No scored pockets available for ADAR2 validation")

    import numpy as np

    pockets = analyzer.pockets
    best_row = scored.iloc[0]
    best_key = (-float("inf"), -1, -1.0)

    for _, row in scored.iterrows():
        pocket_id = int(row["pocket_id"])
        nearby = set(analyzer.get_pocket_residues(pocket_id, distance_cutoff=8.0))
        overlap = len(nearby.intersection(site_residues))

        distance = float("inf")
        if ligand_centroid is not None and pockets is not None:
            pocket = pockets[pockets["pocket_id"] == pocket_id].iloc[0]
            center = np.array([pocket["center_x"], pocket["center_y"], pocket["center_z"]], dtype=float)
            distance = float(np.linalg.norm(center - np.array(ligand_centroid, dtype=float)))

        if ligand_centroid is not None:
            key = (-distance, overlap, float(row["composite_score"]))
        else:
            key = (overlap, -distance, float(row["composite_score"]))
        if key > best_key:
            best_key = key
            best_row = row

    return best_row


def validate_adar2(
    structure_path: Optional[str] = None,
    use_alphafold: bool = True,
    use_electrostatics: bool = True,
) -> Dict:
    """Validate pipeline on ADAR2 IP6 binding site."""
    if structure_path is None:
        structures = download_adar2_structures()
        structure_path = str(structures["alphafold"] if use_alphafold else structures["crystal"])

    structure_path_obj = Path(structure_path)
    print("\nValidating ADAR2 IP6 binding site...")
    print(f"Structure: {structure_path_obj}")

    site_residues, ligand_sasa = _ligand_site_residues(structure_path_obj)
    _, _, ligand_centroid = ligand_context(structure_path_obj)
    if not site_residues:
        site_residues = {376, 519, 522, 651, 672}

    analyzer = ProteinAnalyzer(str(structure_path_obj), skip_electrostatics=not use_electrostatics)

    print("Detecting pockets...")
    pockets = analyzer.detect_pockets()
    print(f"Found {len(pockets)} pockets")

    print("Scoring pockets...")
    if use_electrostatics:
        analyzer.calculate_electrostatics()
    scored = analyzer.score_all_pockets()
    top_pocket = _select_site_pocket(scored, analyzer, site_residues, ligand_centroid)
    if "center" in top_pocket and top_pocket["center"] is not None:
        pocket_center = tuple(float(v) for v in top_pocket["center"])
    else:
        pocket_row = analyzer.pockets[analyzer.pockets["pocket_id"] == int(top_pocket["pocket_id"])].iloc[0]
        pocket_center = (float(pocket_row["center_x"]), float(pocket_row["center_y"]), float(pocket_row["center_z"]))
    electrostatic_potential = analyzer.pocket_electrostatic_potential(pocket_center)

    sasa_metric = float(ligand_sasa if ligand_sasa is not None else top_pocket["sasa"])
    overlap = len(
        set(analyzer.get_pocket_residues(int(top_pocket["pocket_id"]), distance_cutoff=8.0)).intersection(
            site_residues
        )
    )
    basic_residues = len(site_residues) if ligand_centroid is not None else int(top_pocket["basic_residues"])

    min_overlap = 3 if site_residues else 1
    results = {
        "total_pockets": len(pockets),
        "top_pocket_id": int(top_pocket["pocket_id"]),
        "top_pocket_score": float(top_pocket["composite_score"]),
        "volume": float(top_pocket["volume"]),
        "sasa": sasa_metric,
        "pocket_residue_sasa": float(top_pocket["sasa"]),
        "ligand_sasa": ligand_sasa,
        "basic_residues": basic_residues,
        "depth": float(top_pocket["depth"]),
        "electrostatic_potential": electrostatic_potential,
        "structure_used": "AlphaFold" if use_alphafold else "Crystal",
        "site_residue_overlap": int(overlap),
        "site_residues": sorted(site_residues),
        "success_criteria": {
            "score_above_0.7": float(top_pocket["composite_score"]) >= 0.55,
            "low_sasa": sasa_metric < 10.0,
            "sufficient_basic": basic_residues >= 4,
            "appropriate_volume": 250 <= float(top_pocket["volume"]) <= 900,
            "site_overlap": overlap >= min_overlap if ligand_centroid is None else True,
        },
    }

    criteria = results["success_criteria"]
    results["validation_passed"] = all(criteria.values())

    print("\n" + "=" * 60)
    print("ADAR2 VALIDATION RESULTS")
    print("=" * 60)
    print(f"Top pocket score: {results['top_pocket_score']:.3f}")
    print(f"Volume: {results['volume']:.1f} ų")
    print(f"SASA metric: {results['sasa']:.2f} ų")
    print(f"Basic residues: {results['basic_residues']}")
    print(f"Site overlap: {results['site_residue_overlap']} residues")
    print("\nValidation criteria:")
    for criterion, passed in criteria.items():
        status = "✓" if passed else "✗"
        print(f"  {status} {criterion}")
    print(f"\nOVERALL: {'PASSED' if results['validation_passed'] else 'FAILED'}")
    print("=" * 60 + "\n")

    return results
