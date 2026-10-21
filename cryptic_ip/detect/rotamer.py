"""Side-chain rotation, because that is how most cryptic sites actually open.

Why this exists
---------------
`docs/DETECTOR.md` states the sharpest limit of the conformer ensemble: an
elastic network moves atoms rigidly with their residue, so it opens cavities by
moving domains and cannot open a site that is occluded by a single side chain.
That limit matters because side-chain occlusion is a large share of real
cryptic sites - a tryptophan or a leucine lying across a cleft is the
textbook case, and no amount of low-frequency backbone motion will move it.

So this module rotates side chains about their own torsions. It is the second
tier of the ensemble: backbone modes open what domain motion opens, rotamers
open what a chi flip opens, and the two compose.

The method
----------
For each residue whose side chain lines the pocket, rotate every atom beyond
the rotation axis by a set of fixed chi1 offsets about the CA-CB bond, and for
the long side chains also about CB-CG. The offsets are the rotamer library's
own spacing - staggered by 120 degrees for sp3 torsions, plus the +/-30 degree
neighbours that cover the width of a rotamer well - rather than a continuous
scan, because the point is to span the accessible states cheaply, not to
minimise an energy.

A rotation that puts the side chain inside another atom is discarded. The
clash test is a hard sphere overlap against everything that did not move,
which is crude and deliberately so: a soft-potential screen would need a force
field and a minimisation, and this has to run hundreds of times per protein.

What this is not
----------------
It is not a rotamer library and it does not know rotamer probabilities. A
conformer it produces may be strained in ways the clash test cannot see, and
its probability is unmodelled - so the ensemble it feeds reports *whether* a
pocket can open, never how often. It also cannot do a coupled flip, where two
side chains have to move together; each residue is rotated independently.

Those limits are why `docs/DETECTOR_PLAN.md` asks only whether crypticity
ranks known cryptic binders above matched controls, and not what the opening
free energy is.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field as dataclasses_field
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np

LOGGER = logging.getLogger(__name__)

#: Side-chain torsions, as the atom pair defining the axis and the atoms that
#: move with it. Only the heavy atoms a PDB file carries are listed; hydrogens
#: move with whatever they are attached to and are not needed for a clash test
#: at this resolution.
CHI1_AXIS = ("CA", "CB")
CHI2_AXIS = ("CB", "CG")

#: Residues with a rotatable chi1, and the atoms distal to the CA-CB axis.
#: Glycine has no side chain, alanine's is a single methyl with no distinct
#: rotamer, and proline's is a closed ring.
CHI1_MOVING: Dict[str, Tuple[str, ...]] = {
    "ARG": ("CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "ASN": ("CG", "OD1", "ND2"),
    "ASP": ("CG", "OD1", "OD2"),
    "CYS": ("SG",),
    "GLN": ("CG", "CD", "OE1", "NE2"),
    "GLU": ("CG", "CD", "OE1", "OE2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
    "ILE": ("CG1", "CG2", "CD1"),
    "LEU": ("CG", "CD1", "CD2"),
    "LYS": ("CG", "CD", "CE", "NZ"),
    "MET": ("CG", "SD", "CE"),
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "SER": ("OG",),
    "THR": ("OG1", "CG2"),
    "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "VAL": ("CG1", "CG2"),
}

#: Residues long enough for a second torsion to matter, and what moves with it.
CHI2_MOVING: Dict[str, Tuple[str, ...]] = {
    "ARG": ("CD", "NE", "CZ", "NH1", "NH2"),
    "GLN": ("CD", "OE1", "NE2"),
    "GLU": ("CD", "OE1", "OE2"),
    "LYS": ("CD", "CE", "NZ"),
    "MET": ("SD", "CE"),
    "PHE": ("CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "TYR": ("CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "HIS": ("ND1", "CD2", "CE1", "NE2"),
}

#: Torsion offsets in degrees. The staggered sp3 minima are 120 degrees apart;
#: the +/-40 degree neighbours cover the width of a rotamer well without
#: pretending to be a continuous scan.
DEFAULT_OFFSETS: Tuple[float, ...] = (-120.0, -80.0, -40.0, 40.0, 80.0, 120.0)

#: How close a moved atom may come to an atom that did not move, in angstrom.
#: 2.9 is a hard-sphere contact between two carbons, slightly permissive so a
#: conformer is not discarded for a tenth of an angstrom.
CLASH_DISTANCE = 2.9

#: A clash against the residue's own backbone within this many bonds is
#: geometric rather than real, so neighbours in sequence are exempt.
SEQUENCE_EXEMPT = 1

#: Physical bond length range for the torsion axis, in angstrom. A CA-CB bond
#: is 1.53 A and a CB-CG 1.52; anything outside this window is not a covalent
#: bond and no torsion about it means anything.
#:
#: This exists because `cryptic_ip.testing.synthetic` places atoms on a lattice
#: and labels them by priority - its own header says NOT A REAL STRUCTURE - so
#: its "CA-CB bonds" are 3.3, 11.3 and 14.7 A long. Rotating about those axes
#: produced contacts of 0.17 A, two atoms superimposed, and every conformer
#: was silently rejected as a clash. Refusing the input says what is wrong
#: instead.
BOND_LENGTH_RANGE = (1.2, 2.1)

#: How much closer than the input a rotation may bring a side chain before it
#: counts as a new clash, in angstrom. Only applies when the input is already
#: strained; a clean input is held to CLASH_DISTANCE outright.
CLASH_TOLERANCE = 0.2


class RotamerError(ValueError):
    """A residue or structure from which no rotation can be built."""


@dataclass(frozen=True)
class SideChain:
    """One rotatable side chain, as indices into the structure's atom arrays."""

    residue: int
    resname: str
    axis_from: int
    axis_to: int
    moving: Tuple[int, ...]
    torsion: str

    @property
    def n_moving(self) -> int:
        return len(self.moving)


def rotatable_side_chains(
    arrays,
    residues: Optional[Sequence[int]] = None,
    *,
    include_chi2: bool = True,
) -> List[SideChain]:
    """Every rotatable torsion of the named residues, in structure order.

    ``residues`` are residue indices as `StructureArrays.residue_index` numbers
    them; ``None`` means every polymer residue, which is usually far more than
    a pocket needs.
    """
    polymer = np.asarray(arrays.is_polymer, dtype=bool)
    names = np.asarray(arrays.atom_names, dtype=object)
    resnames = np.asarray(arrays.resnames, dtype=object)
    residue_index = np.asarray(arrays.residue_index, dtype=int)
    wanted: Optional[Set[int]] = None if residues is None else {int(r) for r in residues}

    out: List[SideChain] = []
    unphysical: List[Tuple[int, str, str, float]] = []
    for residue in sorted(set(residue_index[polymer].tolist())):
        if wanted is not None and residue not in wanted:
            continue
        atoms = np.nonzero(polymer & (residue_index == residue))[0]
        if atoms.size == 0:
            continue
        resname = str(resnames[atoms[0]]).upper()
        by_name = {str(names[i]).upper(): int(i) for i in atoms}
        for torsion, axis, table in (
            ("chi1", CHI1_AXIS, CHI1_MOVING),
            ("chi2", CHI2_AXIS, CHI2_MOVING),
        ):
            if torsion == "chi2" and not include_chi2:
                continue
            moving_names = table.get(resname)
            if not moving_names:
                continue
            if axis[0] not in by_name or axis[1] not in by_name:
                continue
            moving = tuple(by_name[name] for name in moving_names if name in by_name)
            if not moving:
                continue
            length = float(np.linalg.norm(
                arrays.coords[by_name[axis[1]]] - arrays.coords[by_name[axis[0]]]
            ))
            if not BOND_LENGTH_RANGE[0] <= length <= BOND_LENGTH_RANGE[1]:
                unphysical.append((int(residue), resname, torsion, length))
                continue
            out.append(SideChain(
                residue=int(residue), resname=resname,
                axis_from=by_name[axis[0]], axis_to=by_name[axis[1]],
                moving=moving, torsion=torsion,
            ))
    if unphysical and not out:
        worst = max(unphysical, key=lambda item: item[3])
        raise RotamerError(
            f"no residue has a physical torsion axis: {len(unphysical)} candidate(s) rejected, "
            f"the worst being {worst[1]} {worst[0]} {worst[2]} with an axis {worst[3]:.2f} A long "
            f"(a covalent bond is {BOND_LENGTH_RANGE[0]}-{BOND_LENGTH_RANGE[1]} A). "
            "This structure has no covalent geometry - rotating about such an axis flings atoms "
            "through the protein and every conformer is rejected as a clash, which looks like a "
            "packed pocket rather than like bad input."
        )
    if unphysical:
        LOGGER.warning(
            "%d torsion axis/axes are not covalent bonds and were skipped; worst %.2f A",
            len(unphysical), max(item[3] for item in unphysical),
        )
    return out


def rotate_about_axis(
    coords: np.ndarray,
    origin: np.ndarray,
    axis: np.ndarray,
    indices: Sequence[int],
    degrees: float,
) -> np.ndarray:
    """``coords`` with ``indices`` rotated about the ``origin``-``axis`` line."""
    coords = np.asarray(coords, dtype=float)
    axis = np.asarray(axis, dtype=float)
    norm = float(np.linalg.norm(axis))
    if norm == 0.0:
        raise RotamerError("a torsion axis of zero length: two atoms share a position")
    unit = axis / norm
    angle = np.deg2rad(float(degrees))
    cos, sin = np.cos(angle), np.sin(angle)
    cross = np.array([
        [0.0, -unit[2], unit[1]],
        [unit[2], 0.0, -unit[0]],
        [-unit[1], unit[0], 0.0],
    ])
    matrix = cos * np.eye(3) + sin * cross + (1.0 - cos) * np.outer(unit, unit)
    out = coords.copy()
    picked = np.asarray(list(indices), dtype=int)
    if picked.size == 0:
        return out
    out[picked] = (coords[picked] - origin) @ matrix.T + origin
    return out


def closest_contact(
    coords: np.ndarray,
    moved: Sequence[int],
    residue_index: np.ndarray,
    residue_of_moved: int,
) -> float:
    """The shortest distance from a moved atom to an atom that did not move.

    The residue's own atoms and its sequence neighbours are exempt: a side
    chain is always close to its own backbone, and that is geometry rather than
    a clash. ``inf`` when there is nothing to collide with.
    """
    picked = np.asarray(list(moved), dtype=int)
    if picked.size == 0:
        return float("inf")
    residue_index = np.asarray(residue_index, dtype=int)
    exempt = np.abs(residue_index - residue_of_moved) <= SEQUENCE_EXEMPT
    static = np.ones(len(coords), dtype=bool)
    static[picked] = False
    static &= ~exempt
    if not np.any(static):
        return float("inf")
    delta = coords[picked][:, None, :] - coords[static][None, :, :]
    return float(np.sqrt(np.min(np.sum(delta * delta, axis=2))))


def clashes(
    coords: np.ndarray,
    moved: Sequence[int],
    residue_index: np.ndarray,
    residue_of_moved: int,
    *,
    distance: float = CLASH_DISTANCE,
    reference: Optional[float] = None,
    tolerance: float = CLASH_TOLERANCE,
) -> bool:
    """Whether a rotation introduces a clash the input did not already have.

    An absolute threshold is the wrong test and measuring it that way produced
    an empty ensemble with no explanation. A predicted model routinely carries
    marginal contacts - a side chain 1.9 A from a wall - and against a fixed
    2.9 A cutoff *every* rotation of that side chain is a clash, including the
    one that opens the pocket. So the test is relative: a rotation is rejected
    if it brings the side chain closer than the input already had it, by more
    than ``tolerance``, or closer than ``distance`` when the input was clean.

    ``reference`` is the input's own closest contact for this side chain. With
    ``None`` the test falls back to the absolute threshold, which is correct
    only for a structure known to be unstrained.
    """
    contact = closest_contact(coords, moved, residue_index, residue_of_moved)
    if reference is None:
        return contact < distance
    allowed = min(distance, reference - tolerance) if reference < distance else distance
    return contact < allowed


@dataclass
class RotamerEnsemble:
    """The conformers one side chain's rotation produced, and what was rejected."""

    conformers: List[np.ndarray]
    accepted: List[Tuple[int, str, float]]
    rejected: List[Tuple[int, str, float]]
    strained_inputs: List[Tuple[int, str, float]] = dataclasses_field(default_factory=list)

    @property
    def n_accepted(self) -> int:
        return len(self.conformers)

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_accepted": self.n_accepted,
            "n_rejected": len(self.rejected),
            "accepted": [[r, t, float(d)] for r, t, d in self.accepted],
            "rejected": [[r, t, float(d)] for r, t, d in self.rejected],
            # Side chains the input already had in contact. Published because
            # their rotations were screened against a looser bar, and a reader
            # has to be able to see which ones.
            "strained_inputs": [[r, t, float(d)] for r, t, d in self.strained_inputs],
        }


def rotamer_conformers(
    arrays,
    side_chains: Sequence[SideChain],
    *,
    offsets: Sequence[float] = DEFAULT_OFFSETS,
    clash_distance: float = CLASH_DISTANCE,
    max_conformers: Optional[int] = None,
) -> RotamerEnsemble:
    """One conformer per accepted (side chain, offset) pair.

    Each side chain is rotated independently from the input structure, not
    cumulatively. A coupled flip, where two side chains have to move together
    to open a site, is therefore out of reach - stated here because it is the
    obvious next limitation and it is not addressed.
    """
    coords = np.asarray(arrays.coords, dtype=float)
    residue_index = np.asarray(arrays.residue_index, dtype=int)
    conformers: List[np.ndarray] = []
    accepted: List[Tuple[int, str, float]] = []
    rejected: List[Tuple[int, str, float]] = []

    strained: List[Tuple[int, str, float]] = []
    for chain in side_chains:
        origin = coords[chain.axis_from]
        axis = coords[chain.axis_to] - origin
        # What the input already had, so a strained start does not reject
        # every rotation including the one that opens the pocket.
        reference = closest_contact(coords, chain.moving, residue_index, chain.residue)
        if reference < clash_distance:
            strained.append((chain.residue, chain.torsion, float(reference)))
        for offset in offsets:
            if max_conformers is not None and len(conformers) >= max_conformers:
                return RotamerEnsemble(conformers, accepted, rejected, strained)
            try:
                moved = rotate_about_axis(coords, origin, axis, chain.moving, offset)
            except RotamerError:
                rejected.append((chain.residue, chain.torsion, offset))
                continue
            if clashes(
                moved, chain.moving, residue_index, chain.residue,
                distance=clash_distance, reference=reference,
            ):
                rejected.append((chain.residue, chain.torsion, offset))
                continue
            conformers.append(moved)
            accepted.append((chain.residue, chain.torsion, float(offset)))
    return RotamerEnsemble(conformers, accepted, rejected, strained)


def lining_residues(arrays, centre: Sequence[float], radius: float = 8.0) -> List[int]:
    """Polymer residues with any atom within ``radius`` of a point.

    Deliberately generous and deliberately not requiring contact with a
    detected cavity surface: a side chain that *occludes* a cryptic site is
    not lining a cavity, because there is no cavity yet.
    """
    polymer = np.asarray(arrays.is_polymer, dtype=bool)
    residue_index = np.asarray(arrays.residue_index, dtype=int)
    centre = np.asarray(centre, dtype=float).reshape(3)
    distance = np.linalg.norm(np.asarray(arrays.coords, dtype=float) - centre[None, :], axis=1)
    near = polymer & (distance <= float(radius))
    return sorted(set(residue_index[near].tolist()))
