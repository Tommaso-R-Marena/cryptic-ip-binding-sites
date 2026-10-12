"""Synthetic structures with real covalent geometry, for benchmarking.

Why this exists
---------------
`cryptic_ip.testing.synthetic` builds structures by placing atoms on a lattice
and labelling them by priority. Its own header says NOT A REAL STRUCTURE, and
it is right: in the structures it writes, a CA-CB "bond" can be 3.3, 11.3 or
14.7 angstrom long. That is fine for what it was built for - exercising
descriptor plumbing and burial measurement, which only need positions.

It is not fine for anything that depends on covalent geometry. Rotating a side
chain about a 14.7 angstrom "bond" flings its atoms through the protein; every
conformer is then rejected as a clash, and the result looks exactly like a
tightly packed pocket rather than like invalid input. Benchmarking the rotamer
ensemble or the elastic network on those structures would measure nothing and
would look like a measurement.

So this module builds backbones with ideal geometry - real bond lengths, real
angles, real phi/psi - and attaches side chains with real chi torsions. The
cavities are carved by leaving residues out rather than by deleting atoms from
a lattice, so what remains is a protein-like shell whose torsions mean
something.

What these are and are not
--------------------------
They are geometrically valid and biologically arbitrary. The sequences are
chosen to put a particular residue in a particular place, the folds are
idealised helices and sheets rather than real topologies, and no energy has
been minimised. They can validate that a geometric method does what it claims
on geometry it can trust. They cannot validate that the method finds real
binding sites, and `docs/DETECTOR_PLAN.md` does not ask them to.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: Ideal backbone geometry, in angstrom and degrees. These are the standard
#: values used by every structure-building program.
N_CA = 1.458
CA_C = 1.525
C_N = 1.329
CA_CB = 1.530
CB_CG = 1.520
ANGLE_N_CA_C = 111.0
ANGLE_CA_C_N = 116.2
ANGLE_C_N_CA = 121.7
ANGLE_N_CA_CB = 110.5
ANGLE_CA_CB_CG = 114.0
OMEGA = 180.0

#: Alpha helix and beta strand, as (phi, psi) in degrees.
HELIX = (-57.0, -47.0)
STRAND = (-139.0, 135.0)

#: Van der Waals radii, matching cryptic_ip.detect.cavity.
ELEMENT_OF = {"N": "N", "CA": "C", "C": "C", "O": "O", "CB": "C", "CG": "C",
              "CD": "C", "CD1": "C", "CD2": "C", "CE": "C", "NE": "N", "CZ": "C",
              "NH1": "N", "NH2": "N", "NZ": "N", "OD1": "O", "OD2": "O",
              "OE1": "O", "OE2": "O", "SG": "S", "SD": "S", "OG": "O", "OG1": "O",
              "CG1": "C", "CG2": "C", "ND1": "N", "NE1": "N", "NE2": "N",
              "CE1": "C", "CE2": "C", "CE3": "C", "CZ2": "C", "CZ3": "C",
              "CH2": "C", "OH": "O"}


class FixtureError(ValueError):
    """A fixture specification that cannot be built."""


def place_atom(a: np.ndarray, b: np.ndarray, c: np.ndarray,
               bond: float, angle: float, torsion: float) -> np.ndarray:
    """The position of atom D given A-B-C and the D internal coordinates.

    Standard NeRF placement: ``bond`` is C-D, ``angle`` is B-C-D in degrees,
    ``torsion`` is A-B-C-D in degrees. This is what makes the geometry real
    rather than approximate, and it is the one piece of arithmetic every
    structure builder shares.
    """
    angle_rad = np.deg2rad(angle)
    torsion_rad = np.deg2rad(torsion)
    # Local frame on C, with x along B->C.
    bc = c - b
    bc /= np.linalg.norm(bc)
    ab = b - a
    n = np.cross(ab, bc)
    norm = np.linalg.norm(n)
    if norm < 1e-9:
        raise FixtureError("three collinear atoms give no torsion reference frame")
    n /= norm
    m = np.cross(n, bc)
    d2 = np.array([
        -bond * np.cos(angle_rad),
        bond * np.cos(torsion_rad) * np.sin(angle_rad),
        bond * np.sin(torsion_rad) * np.sin(angle_rad),
    ])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


@dataclass
class Residue:
    """One built residue: atom names, coordinates and elements."""

    index: int
    resname: str
    names: List[str]
    coords: List[np.ndarray]

    def atom(self, name: str) -> np.ndarray:
        return self.coords[self.names.index(name)]

    def has(self, name: str) -> bool:
        return name in self.names


#: Standard deviation of phi/psi in real secondary structure, in degrees.
#: Helices in deposited structures scatter by roughly this much about their
#: ideal torsions, so it is the scale at which a fixture stops being idealised
#: and starts resembling a protein.
REALISTIC_JITTER = 8.0

#: Torsions of a single-residue bulge: a residue that interrupts a helix and
#: shifts everything after it out of register. Taken from the alpha-L region,
#: which is where a classic helix bulge sits.
BULGE = (57.0, 40.0)


def build_backbone(
    n_residues: int,
    *,
    phi_psi: Sequence[Tuple[float, float]] = (HELIX,),
    jitter: float = 0.0,
    bulges: Sequence[int] = (),
    rng: Optional[np.random.Generator] = None,
) -> List[Residue]:
    """A backbone of ``n_residues`` with ideal bond lengths and angles.

    ``phi_psi`` is cycled, so a single pair gives a uniform helix or strand and
    a list gives a mixed fold.

    ``jitter`` adds a Gaussian perturbation, in degrees, to each residue's phi
    and psi. This is the knob that turns an idealised fixture into something
    protein-like: it leaves every bond length and bond angle **exactly** ideal,
    because only torsions move, while the helix axis wanders and the face
    periodicity stops being exactly 3.6. `REALISTIC_JITTER` is the scale real
    structures show.

    ``bulges`` lists residue indices given `BULGE` torsions instead, which
    interrupts the helix and shifts the register of everything after it - the
    structural irregularity most likely to break a periodicity argument.
    """
    if n_residues < 3:
        raise FixtureError(f"need at least 3 residues, got {n_residues}")
    if jitter < 0.0:
        raise FixtureError(f"jitter must not be negative, got {jitter}")
    bulge_at = {int(b) for b in bulges}
    if any(b < 1 or b >= n_residues for b in bulge_at):
        raise FixtureError(f"bulges must lie in 1..{n_residues - 1}, got {sorted(bulge_at)}")
    generator = rng if rng is not None else np.random.default_rng(0)
    residues: List[Residue] = []
    # Seed the first residue in a fixed frame.
    n0 = np.array([0.0, 0.0, 0.0])
    ca0 = np.array([N_CA, 0.0, 0.0])
    c0 = ca0 + CA_C * np.array([np.cos(np.deg2rad(180 - ANGLE_N_CA_C)),
                                np.sin(np.deg2rad(180 - ANGLE_N_CA_C)), 0.0])
    residues.append(Residue(0, "GLY", ["N", "CA", "C"], [n0, ca0, c0]))

    for index in range(1, n_residues):
        if index in bulge_at:
            phi, psi = BULGE
        else:
            phi, psi = phi_psi[(index - 1) % len(phi_psi)]
        if jitter > 0.0:
            phi += float(generator.normal(0.0, jitter))
            psi += float(generator.normal(0.0, jitter))
        previous = residues[-1]
        n_prev, ca_prev, c_prev = previous.atom("N"), previous.atom("CA"), previous.atom("C")
        # The next N is placed by psi of the previous residue.
        n_next = place_atom(n_prev, ca_prev, c_prev, C_N, ANGLE_CA_C_N, psi)
        ca_next = place_atom(ca_prev, c_prev, n_next, N_CA, ANGLE_C_N_CA, OMEGA)
        c_next = place_atom(c_prev, n_next, ca_next, CA_C, ANGLE_N_CA_C, phi)
        residues.append(Residue(index, "GLY", ["N", "CA", "C"], [n_next, ca_next, c_next]))
    return residues


def add_carbonyl(residues: Sequence[Residue]) -> None:
    """Place each residue's carbonyl oxygen, in place."""
    for index, residue in enumerate(residues):
        following = residues[index + 1] if index + 1 < len(residues) else None
        reference = following.atom("N") if following is not None else residue.atom("N")
        oxygen = place_atom(residue.atom("N"), residue.atom("CA"), residue.atom("C"),
                            1.231, 120.5, 180.0 if following is not None else 0.0)
        del reference
        residue.names.append("O")
        residue.coords.append(oxygen)


#: Side-chain build plans: each atom placed from three already-placed atoms,
#: with (bond, angle, torsion) in angstrom and degrees. The torsions are the
#: chi values of the most common rotamer, so the built side chain is in a
#: plausible state before anything rotates it.
SIDE_CHAIN_PLANS: Dict[str, List[Tuple[str, Tuple[str, str, str], Tuple[float, float, float]]]] = {
    "ALA": [("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5))],
    "LEU": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (CB_CG, ANGLE_CA_CB_CG, -60.0)),
        ("CD1", ("CA", "CB", "CG"), (1.521, 110.7, 174.0)),
        ("CD2", ("CA", "CB", "CG"), (1.521, 110.7, -65.0)),
    ],
    "ARG": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (CB_CG, ANGLE_CA_CB_CG, -60.0)),
        ("CD", ("CA", "CB", "CG"), (1.520, 111.3, -177.0)),
        ("NE", ("CB", "CG", "CD"), (1.460, 112.0, -179.0)),
        ("CZ", ("CG", "CD", "NE"), (1.329, 124.2, -178.0)),
        ("NH1", ("CD", "NE", "CZ"), (1.326, 120.0, 0.0)),
        ("NH2", ("CD", "NE", "CZ"), (1.326, 120.0, 180.0)),
    ],
    "LYS": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (CB_CG, ANGLE_CA_CB_CG, -60.0)),
        ("CD", ("CA", "CB", "CG"), (1.520, 111.3, -177.0)),
        ("CE", ("CB", "CG", "CD"), (1.520, 111.3, -179.0)),
        ("NZ", ("CG", "CD", "CE"), (1.489, 111.9, -179.0)),
    ],
    "TRP": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (1.498, 113.6, -60.0)),
        ("CD1", ("CA", "CB", "CG"), (1.365, 126.9, 95.0)),
        ("CD2", ("CA", "CB", "CG"), (1.433, 126.7, -85.0)),
        ("NE1", ("CB", "CG", "CD1"), (1.375, 110.2, 180.0)),
        ("CE2", ("CB", "CG", "CD2"), (1.409, 107.2, 180.0)),
        ("CE3", ("CB", "CG", "CD2"), (1.398, 133.9, 0.0)),
        ("CZ2", ("CG", "CD2", "CE2"), (1.394, 122.4, 180.0)),
        ("CZ3", ("CG", "CD2", "CE3"), (1.382, 118.6, 180.0)),
        ("CH2", ("CD2", "CE2", "CZ2"), (1.368, 117.5, 0.0)),
    ],
    "ASP": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (1.516, 112.6, -60.0)),
        ("OD1", ("CA", "CB", "CG"), (1.249, 118.4, -20.0)),
        ("OD2", ("CA", "CB", "CG"), (1.249, 118.4, 160.0)),
    ],
    "GLU": [
        ("CB", ("C", "N", "CA"), (CA_CB, ANGLE_N_CA_CB, 122.5)),
        ("CG", ("N", "CA", "CB"), (1.520, 113.4, -60.0)),
        ("CD", ("CA", "CB", "CG"), (1.516, 112.6, -177.0)),
        ("OE1", ("CB", "CG", "CD"), (1.249, 118.4, -20.0)),
        ("OE2", ("CB", "CG", "CD"), (1.249, 118.4, 160.0)),
    ],
    "GLY": [],
}


def add_side_chain(residue: Residue, resname: str) -> None:
    """Build ``resname``'s side chain onto an existing backbone residue."""
    plan = SIDE_CHAIN_PLANS.get(resname.upper())
    if plan is None:
        raise FixtureError(f"no build plan for {resname}; known: {sorted(SIDE_CHAIN_PLANS)}")
    residue.resname = resname.upper()
    for name, (a_name, b_name, c_name), (bond, angle, torsion) in plan:
        for required in (a_name, b_name, c_name):
            if not residue.has(required):
                raise FixtureError(f"{resname} {name} needs {required}, which is not built yet")
        position = place_atom(
            residue.atom(a_name), residue.atom(b_name), residue.atom(c_name),
            bond, angle, torsion,
        )
        residue.names.append(name)
        residue.coords.append(position)


@dataclass
class Fixture:
    """A built structure, in the shape the detector's functions expect."""

    coords: np.ndarray
    atom_names: np.ndarray
    resnames: np.ndarray
    residue_index: np.ndarray
    chain_ids: np.ndarray
    resseqs: np.ndarray
    elements: np.ndarray
    is_polymer: np.ndarray
    is_solvent: np.ndarray
    bfactors: np.ndarray
    radii: np.ndarray
    occupancies: np.ndarray
    icodes: np.ndarray
    model_ids: np.ndarray
    residue_keys: List[Tuple[str, int, str]]
    source_path: Optional[str] = None

    @property
    def n_atoms(self) -> int:
        return int(len(self.coords))

    @property
    def is_hetero(self) -> np.ndarray:
        return ~self.is_polymer & ~self.is_solvent


def assemble(residues: Sequence[Residue], *, plddt: float = 90.0) -> Fixture:
    """Turn built residues into the parallel arrays the detector reads."""
    coords: List[np.ndarray] = []
    names: List[str] = []
    resnames: List[str] = []
    residue_index: List[int] = []
    resseqs: List[int] = []
    for residue in residues:
        for name, position in zip(residue.names, residue.coords):
            coords.append(position)
            names.append(name)
            resnames.append(residue.resname)
            residue_index.append(residue.index)
            resseqs.append(residue.index + 1)
    n = len(coords)
    elements = [ELEMENT_OF.get(name, "C") for name in names]
    from cryptic_ip.detect.cavity import VDW_RADII, DEFAULT_RADIUS
    return Fixture(
        coords=np.asarray(coords, dtype=float),
        atom_names=np.asarray(names, dtype=object),
        resnames=np.asarray(resnames, dtype=object),
        residue_index=np.asarray(residue_index, dtype=int),
        chain_ids=np.asarray(["A"] * n, dtype=object),
        resseqs=np.asarray(resseqs, dtype=int),
        elements=np.asarray(elements, dtype=object),
        is_polymer=np.ones(n, dtype=bool),
        is_solvent=np.zeros(n, dtype=bool),
        bfactors=np.full(n, float(plddt)),
        radii=np.asarray([VDW_RADII.get(e, DEFAULT_RADIUS) for e in elements], dtype=float),
        occupancies=np.ones(n, dtype=float),
        icodes=np.asarray([""] * n, dtype=object),
        model_ids=np.zeros(n, dtype=int),
        residue_keys=[("A", residue.index + 1, "") for residue in residues],
    )


def helix_barrel(
    n_helices: int = 6,
    residues_per_helix: int = 18,
    barrel_radius: float = 9.0,
    *,
    sequence: str = "ALA",
    lining: Optional[Dict[int, str]] = None,
    plddt: float = 90.0,
    jitter: float = 0.0,
    bulges: Optional[Dict[int, Sequence[int]]] = None,
    seed: int = 0,
) -> Fixture:
    """A ring of ideal helices around a central channel.

    This is the workhorse fixture: a barrel has a real cavity down its axis,
    every residue has valid covalent geometry, and a residue placed on the
    inner face can be made to occlude the channel. ``lining`` maps a helix
    index to the residue name placed at its mid-point, facing inward.
    """
    if n_helices < 3:
        raise FixtureError(f"a barrel needs at least 3 helices, got {n_helices}")
    lining = dict(lining or {})
    bulges = dict(bulges or {})
    generator = np.random.default_rng(seed)
    all_residues: List[Residue] = []
    offset = 0
    for helix in range(n_helices):
        built = build_backbone(residues_per_helix, phi_psi=(HELIX,), jitter=jitter,
                               bulges=bulges.get(helix, ()), rng=generator)
        add_carbonyl(built)
        for position, residue in enumerate(built):
            name = lining.get(helix) if position == residues_per_helix // 2 else None
            add_side_chain(residue, name or sequence)
        coords = np.asarray([c for r in built for c in r.coords])
        centre = coords.mean(axis=0)
        # Lay the helix axis along z, then place it on the barrel.
        axis = _principal_axis(coords - centre)
        rotation = _align(axis, np.array([0.0, 0.0, 1.0]))
        angle = 2.0 * np.pi * helix / n_helices
        shift = np.array([barrel_radius * np.cos(angle), barrel_radius * np.sin(angle), 0.0])
        for residue in built:
            residue.index += offset
            residue.coords = [rotation @ (c - centre) + shift for c in residue.coords]
        offset += residues_per_helix
        all_residues.extend(built)
    return assemble(all_residues, plddt=plddt)


def strand_sheet(
    n_strands: int = 5,
    residues_per_strand: int = 12,
    strand_spacing: float = 4.8,
    *,
    sequence: str = "ALA",
    lining: Optional[Dict[int, str]] = None,
    antiparallel: bool = True,
    plddt: float = 90.0,
) -> Fixture:
    """A beta sheet: strands laid side by side, one face pointing the same way.

    The second fold family in this module, and it exists so that a method which
    claims to read architecture can be asked to tell two folds apart. A sheet's
    structure is different from a barrel's in exactly the way that matters
    here: an extended strand alternates faces every residue, so the residues on
    one face of the sheet are two apart along the chain rather than three or
    four.

    ``strand_spacing`` is 4.8 angstrom, the inter-strand distance of a real
    sheet. ``antiparallel`` flips alternate strands, which is what most sheets
    do and which changes the sign of the chain direction without changing the
    face periodicity.
    """
    if n_strands < 2:
        raise FixtureError(f"a sheet needs at least 2 strands, got {n_strands}")
    lining = dict(lining or {})
    all_residues: List[Residue] = []
    offset = 0
    for strand in range(n_strands):
        built = build_backbone(residues_per_strand, phi_psi=(STRAND,))
        add_carbonyl(built)
        for position, residue in enumerate(built):
            name = lining.get(strand) if position == residues_per_strand // 2 else None
            add_side_chain(residue, name or sequence)
        coords = np.asarray([c for r in built for c in r.coords])
        centre = coords.mean(axis=0)
        axis = _principal_axis(coords - centre)
        rotation = _align(axis, np.array([0.0, 0.0, 1.0]))
        if antiparallel and strand % 2 == 1:
            # Flip about x, which reverses the chain direction along z while
            # leaving the sheet plane and the face periodicity alone.
            rotation = np.diag([1.0, -1.0, -1.0]) @ rotation
        shift = np.array([strand * strand_spacing, 0.0, 0.0])
        for residue in built:
            residue.index += offset
            residue.coords = [rotation @ (c - centre) + shift for c in residue.coords]
        offset += residues_per_strand
        all_residues.extend(built)
    return assemble(all_residues, plddt=plddt)


def _principal_axis(centred: np.ndarray) -> np.ndarray:
    """The direction of greatest variance, which for a helix is its axis."""
    _, _, vectors = np.linalg.svd(centred, full_matrices=False)
    return vectors[0]


def _align(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """A rotation matrix taking ``source`` to ``target``."""
    a = source / np.linalg.norm(source)
    b = target / np.linalg.norm(target)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def bond_length_report(fixture: Fixture) -> Dict[str, float]:
    """Measured CA-CB and CB-CG lengths, so a fixture can prove it is valid.

    A benchmark that depends on covalent geometry has to be able to show its
    inputs have some, which is the check the lattice fixtures fail.
    """
    out: Dict[str, List[float]] = {"CA_CB": [], "CB_CG": [], "N_CA": [], "CA_C": []}
    for residue in sorted(set(fixture.residue_index.tolist())):
        mask = fixture.residue_index == residue
        names = {str(n): fixture.coords[i] for i, n in
                 zip(np.nonzero(mask)[0], fixture.atom_names[mask])}
        for key, (a, b) in (("CA_CB", ("CA", "CB")), ("CB_CG", ("CB", "CG")),
                            ("N_CA", ("N", "CA")), ("CA_C", ("CA", "C"))):
            if a in names and b in names:
                out[key].append(float(np.linalg.norm(names[b] - names[a])))
    return {
        f"{key}_{stat}": float(value)
        for key, values in out.items() if values
        for stat, value in (("mean", np.mean(values)), ("max_deviation",
                            np.max(np.abs(np.asarray(values) - np.mean(values)))))
    }
