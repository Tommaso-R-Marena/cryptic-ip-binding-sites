"""Identify inositol phosphate ligands from coordinates rather than from a list.

Rationale
---------
:mod:`cryptic_ip.database.ip_ligands` already refuses to trust a hand-written
list of chemical component identifiers when building the ground-truth dataset:
it discovers components at run time and derives the inositol phosphate series by
counting phosphorus atoms in the reported formula, because an identifier is a
label and a formula is evidence.

The structure-parsing half of the pipeline had no equivalent. It matched residue
names against a hard-coded set, which carries the same two defects the database
module was written to eliminate:

1. **Incompleteness.** A binding site whose ligand identifier is absent from the
   set is invisible - not scored as a negative, simply never seen. Since that set
   also defines which pockets become positive training labels, anything missing
   from it biases the model rather than merely shrinking the dataset.
2. **Wrong species.** The set included ``INS`` (*myo*-inositol), which carries no
   phosphate at all. :mod:`~cryptic_ip.database.ip_ligands` records ``INS``
   explicitly as an unphosphorylated *negative reference*, so the two halves of
   the pipeline disagreed about what an inositol phosphate is. Because burial is
   measured on the *most buried* matching copy, a structure containing both a
   free inositol and a genuine inositol phosphate could have its burial measured
   on the wrong molecule entirely.

This module applies the database module's principle to coordinates. A residue is
an inositol phosphate when its *atoms* say so:

* six carbons joined in a ring, at cyclohexane bonding distances;
* an oxygen on essentially every ring carbon (inositol is a
  hexahydroxycyclohexane, so its ring is fully substituted);
* at least one phosphorus bonded through one of those oxygens.

The phosphorus count then gives the series (InsP1 ... InsP8) by the same rule
:func:`~cryptic_ip.database.ip_ligands.ip_series_label` uses on formulae. No
network access, no identifier vocabulary, and a deoxy or otherwise modified
analogue is recognised on its structure rather than on whether anyone typed its
code into a set.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

from ..database.ip_ligands import ip_series_label
from .structure_arrays import ResidueKey, StructureArrays

LOGGER = logging.getLogger(__name__)

#: Carbon-carbon single bond is ~1.54 A; the cutoff allows for coordinate error
#: and low-resolution refinement without reaching the ~2.5 A 1,3-distance across
#: a ring vertex, which would create spurious ring chords.
CC_BOND_MAX = 1.75

#: Carbon-oxygen single bond is ~1.43 A.
CO_BOND_MAX = 1.65

#: Phosphorus-oxygen bonds run ~1.5 A (P=O) to ~1.6 A (P-O-C).
PO_BOND_MAX = 1.90

#: The inositol ring is six carbons.
INOSITOL_RING_SIZE = 6

#: Ring carbons that must bear an oxygen. Inositol is fully substituted, but
#: deoxy analogues and partially modelled ligands occur, so one bare position is
#: tolerated; two or more means the ring is not an inositol.
MIN_SUBSTITUTED_RING_CARBONS = 5

#: Upper bound on residue carbon count, matching
#: :data:`cryptic_ip.database.ip_ligands.MAX_INOSITOL_CARBONS`. Small substituted
#: derivatives are admitted; lipid-linked species are not.
MAX_INOSITOL_CARBONS = 12


@dataclass(frozen=True)
class InositolResidue:
    """An inositol-cored residue identified from its coordinates.

    Attributes:
        residue_key: Chain-aware key ``(model, chain, resseq, icode)``.
        comp_id: Residue name as deposited. Recorded for provenance only - it
            plays no part in the identification.
        atom_indices: Indices into the structure arrays for this residue's atoms.
        n_phosphorus: Phosphorus atoms bonded through a ring-oxygen.
        n_ring_carbons: Carbons in the detected ring (always six when detected).
        n_substituted_ring_carbons: Ring carbons bearing an oxygen.
        series: Series label such as ``"InsP6"``.
        is_phosphorylated: ``True`` when at least one phosphate is attached.
    """

    residue_key: ResidueKey
    comp_id: str
    atom_indices: np.ndarray
    n_phosphorus: int
    n_ring_carbons: int
    n_substituted_ring_carbons: int
    series: str
    is_phosphorylated: bool


def _find_six_carbon_ring(carbon_coords: np.ndarray) -> Optional[List[int]]:
    """Return indices of six carbons forming a ring, or ``None``.

    Args:
        carbon_coords: Coordinates of the residue's carbon atoms.

    Returns:
        Local indices of a six-membered carbon cycle, or ``None`` when the
        carbon skeleton contains no such ring.
    """
    n_carbons = int(carbon_coords.shape[0])
    if n_carbons < INOSITOL_RING_SIZE:
        return None

    tree = cKDTree(carbon_coords)
    neighbours: List[List[int]] = [[] for _ in range(n_carbons)]
    for left, right in tree.query_pairs(CC_BOND_MAX):
        neighbours[left].append(right)
        neighbours[right].append(left)

    # A cyclohexane carbon has two ring neighbours, so any carbon with fewer than
    # two bonded carbons cannot lie on the ring.
    if all(len(adjacent) < 2 for adjacent in neighbours):
        return None

    # Depth-first search for a cycle of exactly six distinct carbons. Ligand
    # residues hold a few tens of atoms, so the bounded search is cheap.
    def walk(start: int, current: int, path: List[int]) -> Optional[List[int]]:
        if len(path) == INOSITOL_RING_SIZE:
            return list(path) if start in neighbours[current] else None
        for neighbour in neighbours[current]:
            if neighbour in path:
                continue
            path.append(neighbour)
            found = walk(start, neighbour, path)
            if found is not None:
                return found
            path.pop()
        return None

    for start in range(n_carbons):
        if len(neighbours[start]) < 2:
            continue
        ring = walk(start, start, [start])
        if ring is not None:
            return ring
    return None


def _classify_residue(
    coords: np.ndarray,
    elements: np.ndarray,
    atom_indices: np.ndarray,
    residue_key: ResidueKey,
    comp_id: str,
) -> Optional[InositolResidue]:
    """Classify one residue's atoms as an inositol core, or reject it.

    Args:
        coords: Coordinates of this residue's atoms.
        elements: Element symbols of this residue's atoms.
        atom_indices: Indices of these atoms in the parent structure arrays.
        residue_key: Chain-aware residue key.
        comp_id: Residue name, recorded for provenance.

    Returns:
        The classified residue, or ``None`` when it is not inositol-cored.
    """
    is_carbon = elements == "C"
    n_carbons = int(is_carbon.sum())
    if n_carbons < INOSITOL_RING_SIZE or n_carbons > MAX_INOSITOL_CARBONS:
        return None

    carbon_positions = np.flatnonzero(is_carbon)
    ring_local = _find_six_carbon_ring(coords[carbon_positions])
    if ring_local is None:
        return None
    ring_positions = carbon_positions[np.asarray(ring_local, dtype=int)]

    oxygen_positions = np.flatnonzero(elements == "O")
    if oxygen_positions.size == 0:
        return None

    # Inositol's ring is fully hydroxylated, so demand an oxygen on essentially
    # every ring carbon. This is what separates an inositol from an unrelated
    # six-carbon ring such as a phenyl or a sugar-free cyclohexane.
    oxygen_tree = cKDTree(coords[oxygen_positions])
    substituent_oxygens: set = set()
    n_substituted = 0
    for ring_position in ring_positions:
        attached = oxygen_tree.query_ball_point(coords[ring_position], CO_BOND_MAX)
        if attached:
            n_substituted += 1
            substituent_oxygens.update(int(oxygen_positions[i]) for i in attached)
    if n_substituted < MIN_SUBSTITUTED_RING_CARBONS:
        return None

    # Count only phosphorus bonded through one of those oxygens, so a phosphate
    # merely sitting nearby in the model is not credited to this ligand.
    phosphorus_positions = np.flatnonzero(elements == "P")
    n_phosphorus = 0
    if phosphorus_positions.size and substituent_oxygens:
        substituent_list = sorted(substituent_oxygens)
        substituent_tree = cKDTree(coords[substituent_list])
        for phosphorus_position in phosphorus_positions:
            if substituent_tree.query_ball_point(coords[phosphorus_position], PO_BOND_MAX):
                n_phosphorus += 1

    return InositolResidue(
        residue_key=residue_key,
        comp_id=comp_id,
        atom_indices=atom_indices,
        n_phosphorus=n_phosphorus,
        n_ring_carbons=INOSITOL_RING_SIZE,
        n_substituted_ring_carbons=n_substituted,
        series=ip_series_label(n_phosphorus),
        is_phosphorylated=n_phosphorus >= 1,
    )


def detect_inositol_residues(
    arrays: StructureArrays,
    *,
    require_phosphate: bool = True,
    include_polymer: bool = False,
) -> List[InositolResidue]:
    """Find inositol-cored ligands in a structure from their coordinates.

    Args:
        arrays: Parsed structure arrays.
        require_phosphate: Return only phosphorylated inositols. With this off,
            unphosphorylated inositols are returned too, carrying
            ``series="InsP0"`` - useful for auditing what a structure holds
            without admitting free inositol as an inositol phosphate.
        include_polymer: Also search polymer residues. Off by default; inositol
            phosphates are non-polymer ligands.

    Returns:
        Matching residues, ordered by descending phosphorus count so that the
        most heavily phosphorylated species - the one a screen is looking for -
        comes first when a structure holds several.
    """
    candidate_mask = arrays.is_hetero if not include_polymer else ~arrays.is_solvent
    if not np.any(candidate_mask):
        return []

    hydrogen_mask = arrays.elements == "H"
    candidate_mask = candidate_mask & ~hydrogen_mask

    found: List[InositolResidue] = []
    for residue_position in np.unique(arrays.residue_index[candidate_mask]):
        atom_indices = np.flatnonzero(
            candidate_mask & (arrays.residue_index == residue_position)
        )
        if atom_indices.size < INOSITOL_RING_SIZE:
            continue
        classified = _classify_residue(
            arrays.coords[atom_indices],
            arrays.elements[atom_indices],
            atom_indices,
            arrays.residue_keys[int(residue_position)],
            str(arrays.resnames[atom_indices[0]]),
        )
        if classified is None:
            continue
        if require_phosphate and not classified.is_phosphorylated:
            LOGGER.debug(
                "%s is an unphosphorylated inositol (%s); not an IP ligand",
                classified.comp_id,
                classified.series,
            )
            continue
        found.append(classified)

    found.sort(key=lambda residue: (-residue.n_phosphorus, residue.residue_key))
    return found


def summarise_detection(residues: Sequence[InositolResidue]) -> Dict[str, int]:
    """Count detected residues by series label, for provenance manifests.

    Args:
        residues: Detected inositol residues.

    Returns:
        Mapping from series label to the number of copies found.
    """
    counts: Dict[str, int] = {}
    for residue in residues:
        counts[residue.series] = counts.get(residue.series, 0) + 1
    return counts


def detected_comp_ids(residues: Sequence[InositolResidue]) -> Tuple[str, ...]:
    """Return the distinct component identifiers among detected residues.

    Identification does not depend on these; they are reported so a run can
    record which chemical components its structural test actually matched.

    Args:
        residues: Detected inositol residues.

    Returns:
        Sorted distinct component identifiers.
    """
    return tuple(sorted({residue.comp_id for residue in residues}))
