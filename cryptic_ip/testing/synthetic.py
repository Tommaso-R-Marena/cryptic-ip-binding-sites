"""Deterministic synthetic protein-ligand structures with known ground truth.

Purpose
-------
Validating a structural pipeline requires inputs whose correct answers are known
in advance. Real PDB entries do not provide that: the "true" burial of a ligand
is itself a measurement. The generators here construct structures where the
answer follows from the construction:

* A ligand placed at the centre of a closed protein shell **must** have relative
  SASA of essentially zero, enclosure of essentially one, and burial depth close
  to the shell radius.
* The same ligand placed outside the shell **must** have relative SASA close to
  one, enclosure near zero, and burial depth near zero.

Any implementation that fails those statements is wrong, independently of
biology. That makes these structures a usable regression suite for the geometry,
labelling and training code, and it makes the full end-to-end pipeline runnable
in environments with no access to the RCSB.

The generated coordinates are physically plausible - correct bond lengths,
realistic packing density, valid PDB records that external tools such as fpocket
parse - but they are *not* real proteins and must never be used to support a
biological claim.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.spatial import cKDTree

#: Body-centred cubic lattice spacing giving the heavy-atom number density of a
#: folded protein interior (2 atoms per a^3 with a = 3.83 Å is ~0.036 atoms/Å³).
PROTEIN_ATOM_SPACING = 3.83

#: Residue composition used for filler positions, in rough proportion to average
#: amino acid frequencies so that hydrophobicity and charge features take
#: realistic values.
FILLER_RESIDUES: Tuple[str, ...] = (
    "ALA", "LEU", "GLY", "VAL", "SER", "THR", "ILE", "GLU", "ASP", "ASN",
    "PHE", "PRO", "GLN", "TYR", "MET", "TRP", "CYS", "ALA", "LEU", "VAL",
)

#: Side-chain atom names emitted per residue type (beyond the backbone).
SIDE_CHAIN_ATOMS: Dict[str, Tuple[Tuple[str, str], ...]] = {
    "ARG": (("CB", "C"), ("CG", "C"), ("CD", "C"), ("NE", "N"), ("CZ", "C"), ("NH1", "N"), ("NH2", "N")),
    "LYS": (("CB", "C"), ("CG", "C"), ("CD", "C"), ("CE", "C"), ("NZ", "N")),
    "HIS": (("CB", "C"), ("CG", "C"), ("ND1", "N"), ("CD2", "C"), ("NE2", "N")),
    "SER": (("CB", "C"), ("OG", "O")),
    "THR": (("CB", "C"), ("OG1", "O"), ("CG2", "C")),
    "TYR": (("CB", "C"), ("CG", "C"), ("CD1", "C"), ("CE1", "C"), ("CZ", "C"), ("OH", "O")),
    "ASP": (("CB", "C"), ("CG", "C"), ("OD1", "O"), ("OD2", "O")),
    "GLU": (("CB", "C"), ("CG", "C"), ("CD", "C"), ("OE1", "O"), ("OE2", "O")),
    "ASN": (("CB", "C"), ("CG", "C"), ("OD1", "O"), ("ND2", "N")),
    "GLN": (("CB", "C"), ("CG", "C"), ("CD", "C"), ("OE1", "O"), ("NE2", "N")),
    "LEU": (("CB", "C"), ("CG", "C"), ("CD1", "C"), ("CD2", "C")),
    "ILE": (("CB", "C"), ("CG1", "C"), ("CG2", "C"), ("CD1", "C")),
    "VAL": (("CB", "C"), ("CG1", "C"), ("CG2", "C")),
    "PHE": (("CB", "C"), ("CG", "C"), ("CD1", "C"), ("CE1", "C"), ("CZ", "C")),
    "TRP": (("CB", "C"), ("CG", "C"), ("CD1", "C"), ("NE1", "N"), ("CZ2", "C")),
    "MET": (("CB", "C"), ("CG", "C"), ("SD", "S"), ("CE", "C")),
    "CYS": (("CB", "C"), ("SG", "S")),
    "PRO": (("CB", "C"), ("CG", "C"), ("CD", "C")),
    "ALA": (("CB", "C"),),
    "GLY": (),
}


@dataclass
class SyntheticStructureSpec:
    """Specification of one synthetic protein-ligand structure.

    Attributes:
        name: Identifier used for the file stem and PDB header.
        protein_radius: Radius of the protein sphere (Å).
        site_offset: Distance of the ligand site from the protein centre (Å).
            ``0`` places the site in the core (buried); a value at or beyond
            ``protein_radius`` places it outside (surface).
        cavity_clearance: Extra clearance added to the 3.2 Å van der Waals
            contact distance when carving the cavity around the ligand (Å).
            Clearance loosens the fit; more than a few tenths of an Ångström
            opens gaps wide enough for the 1.4 Å solvent probe, which would make
            a nominally buried ligand register as partly exposed.
        n_basic_lining: Number of cavity-lining residues converted to Arg/Lys.
        ligand_comp_id: Ligand residue name to emit.
        plddt_mean: Mean B-factor written for protein atoms, standing in for
            AlphaFold pLDDT.
        plddt_spread: Standard deviation of the B-factor distribution.
        seed: Random seed; identical seeds give byte-identical files.
        expected_burial_class: Ground-truth burial class implied by the geometry.
        include_ligand: Write the ligand. Set ``False`` for a decoy: the cavity is
            still carved, so a pocket is detected, but no ligand occupies it and
            every pocket of the structure is a true negative.
        acidic_lining: Line the cavity with Asp/Glu instead of Arg/Lys. Decoys need
            a genuinely different chemistry, otherwise the benchmark asks the model
            to separate two identical distributions.
    """

    name: str
    protein_radius: float = 22.0
    site_offset: float = 0.0
    cavity_clearance: float = 0.0
    n_basic_lining: int = 6
    ligand_comp_id: str = "IHP"
    plddt_mean: float = 88.0
    plddt_spread: float = 4.0
    seed: int = 0
    expected_burial_class: str = "cryptic"
    include_ligand: bool = True
    acidic_lining: bool = False

    def to_dict(self) -> Dict[str, object]:
        """Return a JSON-serialisable representation."""
        return {
            "name": self.name,
            "protein_radius": self.protein_radius,
            "site_offset": self.site_offset,
            "cavity_clearance": self.cavity_clearance,
            "n_basic_lining": self.n_basic_lining,
            "ligand_comp_id": self.ligand_comp_id,
            "include_ligand": self.include_ligand,
            "acidic_lining": self.acidic_lining,
            "seed": self.seed,
            "expected_burial_class": self.expected_burial_class,
        }


@dataclass
class _Atom:
    """One PDB atom record under construction."""

    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    coord: np.ndarray
    element: str
    bfactor: float
    hetatm: bool = False
    occupancy: float = 1.0


def inositol_hexakisphosphate_coords(
    centre: Sequence[float], *, seed: int = 0
) -> List[Tuple[str, str, np.ndarray]]:
    """Build heavy-atom coordinates for an inositol hexakisphosphate-like ligand.

    The geometry reproduces the connectivity and bond lengths of InsP6 - a
    puckered six-carbon ring, one ester oxygen per carbon, one phosphorus per
    ester oxygen and three terminal oxygens per phosphorus - giving 36 heavy
    atoms, matching the C6O24P6 composition of the real component. Correct bond
    lengths matter because van der Waals radii and SASA depend on them.

    Args:
        centre: Ring centre coordinate.
        seed: Seed for the small ring-orientation jitter.

    Returns:
        ``(atom_name, element, coordinate)`` triples.
    """
    rng = np.random.default_rng(seed)
    centre = np.asarray(centre, dtype=float)

    ring_radius = 1.45  # C-C 1.52 Å around a six-membered ring
    ester_length = 1.43  # C-O
    p_length = 1.60  # O-P
    terminal_length = 1.50  # P=O / P-O

    # Random but deterministic ring orientation.
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis)
    reference = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(axis, reference)) > 0.9:
        reference = np.array([0.0, 1.0, 0.0])
    u = np.cross(axis, reference)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)

    atoms: List[Tuple[str, str, np.ndarray]] = []
    for i in range(6):
        angle = 2.0 * np.pi * i / 6.0
        radial = np.cos(angle) * u + np.sin(angle) * v
        # Alternating pucker approximates the chair conformation.
        pucker = 0.25 * (1.0 if i % 2 == 0 else -1.0)
        c_pos = centre + ring_radius * radial + pucker * axis
        atoms.append((f"C{i + 1}", "C", c_pos))

        o_pos = c_pos + ester_length * radial
        atoms.append((f"O{i + 1}", "O", o_pos))

        p_pos = o_pos + p_length * radial
        atoms.append((f"P{i + 1}", "P", p_pos))

        # Three terminal oxygens arranged tetrahedrally about the phosphorus.
        basis_a = np.cross(radial, axis)
        norm_a = np.linalg.norm(basis_a)
        basis_a = basis_a / norm_a if norm_a > 1e-9 else u
        basis_b = np.cross(radial, basis_a)
        for k in range(3):
            phi = 2.0 * np.pi * k / 3.0
            direction = (
                0.33 * radial + 0.94 * (np.cos(phi) * basis_a + np.sin(phi) * basis_b)
            )
            direction /= np.linalg.norm(direction)
            atoms.append((f"O{i + 1}{chr(ord('A') + k)}", "O", p_pos + terminal_length * direction))
    return atoms


def _lattice_points_in_sphere(radius: float, spacing: float) -> np.ndarray:
    """Return body-centred cubic lattice points inside a sphere.

    A BCC lattice at the C-alpha packing distance reproduces the atom density of
    a folded protein core far better than uniform random sampling, which leaves
    voids that fpocket would report as spurious pockets.

    Args:
        radius: Sphere radius (Å).
        spacing: Lattice spacing (Å).

    Returns:
        Coordinates of shape ``(n, 3)``, ordered by distance from the origin.
    """
    n = int(np.ceil(radius / spacing)) + 1
    grid = np.arange(-n, n + 1) * spacing
    xs, ys, zs = np.meshgrid(grid, grid, grid, indexing="ij")
    primary = np.column_stack((xs.ravel(), ys.ravel(), zs.ravel()))
    offset = primary + spacing / 2.0
    points = np.vstack((primary, offset))
    distances = np.linalg.norm(points, axis=1)
    keep = distances <= radius
    points = points[keep]
    order = np.argsort(distances[keep])
    return points[order]


def _lining_shell_positions(
    ligand_coords: np.ndarray,
    *,
    bulk_positions: np.ndarray,
    protein_centre: np.ndarray,
    protein_radius: float,
    contact_distance: float,
    shell_thickness: float = 1.6,
    min_separation: float = 2.9,
    sample_spacing: float = 0.8,
) -> np.ndarray:
    """Pack atoms at van der Waals contact against the ligand surface.

    Candidate positions are drawn from a fine grid, restricted to the thin shell
    just outside contact distance from the ligand, and then thinned greedily so
    no two retained positions are closer than ``min_separation``. Only positions
    inside the protein sphere are kept, so a ligand docked on the exterior gains
    a lining on the protein side only and stays solvent exposed on the other.

    Args:
        ligand_coords: Ligand heavy-atom coordinates.
        bulk_positions: Already-placed protein atoms to avoid clashing with.
        protein_centre: Centre of the protein sphere.
        protein_radius: Radius of the protein sphere.
        contact_distance: Minimum allowed ligand-protein atom distance.
        shell_thickness: Thickness of the candidate shell beyond contact (Å).
        min_separation: Minimum distance between retained lining atoms (Å).
        sample_spacing: Candidate grid spacing (Å).

    Returns:
        Retained lining positions, shape ``(n, 3)``.
    """
    ligand_coords = np.asarray(ligand_coords, dtype=float)
    if ligand_coords.size == 0:
        return np.empty((0, 3), dtype=float)

    lower = ligand_coords.min(axis=0) - (contact_distance + shell_thickness + sample_spacing)
    upper = ligand_coords.max(axis=0) + (contact_distance + shell_thickness + sample_spacing)
    axes = [np.arange(lo, hi + sample_spacing, sample_spacing) for lo, hi in zip(lower, upper)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)

    ligand_tree = cKDTree(ligand_coords)
    distance, _ = ligand_tree.query(grid, k=1)
    in_shell = (distance > contact_distance) & (distance <= contact_distance + shell_thickness)
    grid = grid[in_shell]
    if grid.size == 0:
        return np.empty((0, 3), dtype=float)

    inside = np.linalg.norm(grid - np.asarray(protein_centre, dtype=float), axis=1) <= protein_radius
    grid = grid[inside]
    if grid.size == 0:
        return np.empty((0, 3), dtype=float)

    if bulk_positions is not None and len(bulk_positions):
        bulk_distance, _ = cKDTree(np.asarray(bulk_positions, dtype=float)).query(grid, k=1)
        grid = grid[bulk_distance > min_separation]
        if grid.size == 0:
            return np.empty((0, 3), dtype=float)

    # Greedy thinning, innermost first, so the layer touching the ligand is
    # populated before the outer part of the shell.
    order = np.argsort(ligand_tree.query(grid, k=1)[0])
    grid = grid[order]
    kept: List[np.ndarray] = []
    for candidate in grid:
        if not kept or np.min(np.linalg.norm(np.asarray(kept) - candidate, axis=1)) > min_separation:
            kept.append(candidate)
    return np.asarray(kept, dtype=float) if kept else np.empty((0, 3), dtype=float)


def _residue_atom_plan(resname: str) -> List[Tuple[str, str]]:
    """Return ``(atom_name, element)`` records for one residue, backbone first."""
    backbone = [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]
    return backbone + [(name, element) for name, element in SIDE_CHAIN_ATOMS.get(resname, ())]


def _site_priority(resname: str, atom_name: str) -> int:
    """Ordering key placing phosphate-coordinating atoms nearest the ligand.

    Basic side-chain nitrogens are laid down closest to the site so that a
    synthetic "cryptic" pocket really is lined by coordinating groups, which is
    what the coordination features are supposed to detect.
    """
    from ..analysis.structure_arrays import BASIC_NITROGEN_ATOMS, BASIC_RESIDUES

    if resname in BASIC_RESIDUES and atom_name in BASIC_NITROGEN_ATOMS:
        return 0
    if atom_name in {"OG", "OG1", "OH"}:
        return 1
    if atom_name in {"CB", "CG", "CD", "CE", "CZ"}:
        return 2
    return 3


def build_synthetic_atoms(spec: SyntheticStructureSpec) -> Tuple[List[_Atom], np.ndarray]:
    """Build the atom records for one synthetic structure.

    Atoms are laid on a body-centred cubic lattice at 3.83 Å, which reproduces
    the heavy-atom number density of a folded protein interior (about
    0.035 atoms per Å³). Placing *atoms* rather than residue centres with
    projected side chains matters: a sparse shell leaves gaps wider than the
    1.4 Å solvent probe, so a nominally buried ligand would register as partly
    exposed and the benchmark would encode the wrong ground truth.

    Args:
        spec: Structure specification.

    Returns:
        ``(atoms, site_centre)`` where ``site_centre`` is the ligand centroid.
    """
    rng = np.random.default_rng(spec.seed)
    site_centre = np.array([float(spec.site_offset), 0.0, 0.0])

    # Build the ligand first: the cavity is sized from the ligand's own extent so
    # the protein packs snugly against it.
    ligand_atoms = inositol_hexakisphosphate_coords(site_centre, seed=spec.seed)
    # The cavity is carved against the ligand's own surface, not against its
    # bounding sphere. A spherical carve leaves voids between the phosphate
    # groups and around the ring plane that are wider than the solvent probe, so
    # a fully enclosed ligand would still report substantial accessible surface.
    # Molding the protein to the ligand shape - as a real binding site does -
    # gives the snug packing the ground truth assumes.
    contact_distance = 3.2 + float(spec.cavity_clearance)
    ligand_coords = np.asarray([coord for _, _, coord in ligand_atoms], dtype=float)

    positions = _lattice_points_in_sphere(spec.protein_radius, PROTEIN_ATOM_SPACING)
    positions = positions + rng.normal(scale=0.12, size=positions.shape)
    ligand_tree = cKDTree(ligand_coords)
    nearest_ligand_distance, _ = ligand_tree.query(positions, k=1)
    positions = positions[nearest_ligand_distance > contact_distance]
    if positions.shape[0] == 0:
        raise ValueError(f"Specification {spec.name!r} leaves no protein atoms")

    # The bulk lattice is coarse relative to the solvent probe, so it alone
    # leaves crevices around the ligand through which the probe reaches the
    # ligand surface. An explicit lining shell, packed at van der Waals contact
    # against the ligand, closes them - which is exactly what the first
    # coordination shell of a real binding site does.
    lining = _lining_shell_positions(
        ligand_coords,
        bulk_positions=positions,
        protein_centre=np.zeros(3),
        protein_radius=spec.protein_radius,
        contact_distance=contact_distance,
    )
    if lining.shape[0]:
        positions = np.vstack((lining, positions))

    # Order by distance from the site so the cavity lining is filled first, which
    # makes the lining residues the basic ones.
    order = np.argsort(np.linalg.norm(positions - site_centre, axis=1))
    positions = positions[order]

    atoms: List[_Atom] = []
    serial = 1
    resseq = 1
    cursor = 0
    n_positions = positions.shape[0]
    while cursor < n_positions:
        is_lining = resseq <= max(0, spec.n_basic_lining)
        if is_lining and spec.acidic_lining:
            resname = "ASP" if resseq % 2 == 1 else "GLU"
        elif is_lining:
            resname = "ARG" if resseq % 2 == 1 else "LYS"
        else:
            resname = FILLER_RESIDUES[resseq % len(FILLER_RESIDUES)]

        plan = _residue_atom_plan(resname)
        chunk = positions[cursor : cursor + len(plan)]
        cursor += len(plan)
        if chunk.shape[0] == 0:
            break
        plan = plan[: chunk.shape[0]]

        # Within a residue, place coordinating atoms closest to the site.
        chunk_order = np.argsort(np.linalg.norm(chunk - site_centre, axis=1))
        plan_order = sorted(range(len(plan)), key=lambda i: _site_priority(resname, plan[i][0]))
        bfactor = float(np.clip(rng.normal(spec.plddt_mean, spec.plddt_spread), 20.0, 100.0))
        for slot, plan_index in enumerate(plan_order):
            name, element = plan[plan_index]
            coord = chunk[chunk_order[slot]]
            atoms.append(_Atom(serial, name, resname, "A", resseq, coord, element, bfactor))
            serial += 1
        resseq += 1

    if not spec.include_ligand:
        return atoms, site_centre

    ligand_resseq = resseq + 1
    for name, element, coord in ligand_atoms:
        atoms.append(
            _Atom(
                serial,
                name,
                spec.ligand_comp_id,
                "A",
                ligand_resseq,
                coord,
                element,
                float(np.clip(rng.normal(spec.plddt_mean, spec.plddt_spread), 20.0, 100.0)),
                hetatm=True,
            )
        )
        serial += 1

    return atoms, site_centre


def write_synthetic_structure(spec: SyntheticStructureSpec, out_dir: Path) -> Path:
    """Write a synthetic structure to a PDB file.

    Args:
        spec: Structure specification.
        out_dir: Destination directory, created when missing.

    Returns:
        Path to the written ``.pdb`` file.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    atoms, _ = build_synthetic_atoms(spec)

    lines: List[str] = [
        f"HEADER    SYNTHETIC BENCHMARK STRUCTURE           {spec.name}",
        "REMARK   1 GENERATED BY cryptic_ip.testing.synthetic -- NOT A REAL STRUCTURE.",
        "REMARK   1 FOR PIPELINE VERIFICATION ONLY; NOT VALID FOR BIOLOGICAL CLAIMS.",
        f"REMARK   2 EXPECTED BURIAL CLASS: {spec.expected_burial_class.upper()}",
    ]
    for atom in atoms:
        record = "HETATM" if atom.hetatm else "ATOM  "
        # Atom names occupy columns 13-16; names of fewer than four characters
        # are conventionally offset by one so the element aligns in column 14.
        name = atom.name if len(atom.name) >= 4 else f" {atom.name:<3s}"
        lines.append(
            f"{record}{atom.serial:5d} {name:<4s}{'':1s}{atom.resname:>3s} {atom.chain}"
            f"{atom.resseq:4d}{'':1s}   "
            f"{atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}"
            f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}          {atom.element:>2s}"
        )
    lines.append("END")

    out_path = out_dir / f"{spec.name}.pdb"
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def default_benchmark_specs(
    *,
    n_buried: int = 8,
    n_surface: int = 8,
    n_decoy: int = 8,
    seed: int = 20240601,
) -> List[SyntheticStructureSpec]:
    """Return a balanced set of buried and surface-site specifications.

    Protein radius, cavity radius and lining composition are varied across the
    set so that a classifier trained on it cannot succeed on a single trivial
    cue, and so feature distributions have non-degenerate spread.

    Args:
        n_buried: Number of buried (cryptic) structures.
        n_surface: Number of surface structures.
        n_decoy: Number of ligand-free decoy structures. Decoys supply the
            negative pockets: without them every detected pocket holds a ligand
            and there is no negative class to train against.
        seed: Base seed; each structure gets a distinct derived seed.

    Returns:
        The specifications, buried first, then surface, then decoys.
    """
    rng = np.random.default_rng(seed)
    specs: List[SyntheticStructureSpec] = []

    for i in range(n_buried):
        radius = float(rng.uniform(20.0, 26.0))
        specs.append(
            SyntheticStructureSpec(
                name=f"SYN_BURIED_{i:02d}",
                protein_radius=radius,
                site_offset=0.0,
                cavity_clearance=float(rng.uniform(0.0, 0.25)),
                n_basic_lining=int(rng.integers(5, 9)),
                plddt_mean=float(rng.uniform(82.0, 94.0)),
                seed=seed + i,
                expected_burial_class="cryptic",
            )
        )

    for i in range(n_surface):
        radius = float(rng.uniform(20.0, 26.0))
        specs.append(
            SyntheticStructureSpec(
                name=f"SYN_SURFACE_{i:02d}",
                protein_radius=radius,
                # Site centre sits beyond the protein sphere, so the ligand is
                # docked against the exterior rather than enclosed.
                site_offset=radius + 2.0,
                cavity_clearance=float(rng.uniform(0.0, 0.25)),
                n_basic_lining=int(rng.integers(2, 6)),
                plddt_mean=float(rng.uniform(70.0, 88.0)),
                seed=seed + 100 + i,
                expected_burial_class="surface",
            )
        )

    for i in range(n_decoy):
        radius = float(rng.uniform(20.0, 26.0))
        specs.append(
            SyntheticStructureSpec(
                name=f"SYN_DECOY_{i:02d}",
                protein_radius=radius,
                site_offset=float(rng.uniform(0.0, radius)),
                cavity_clearance=float(rng.uniform(0.0, 0.25)),
                n_basic_lining=int(rng.integers(0, 3)),
                plddt_mean=float(rng.uniform(70.0, 92.0)),
                seed=seed + 200 + i,
                expected_burial_class="none",
                include_ligand=False,
                acidic_lining=True,
            )
        )
    return specs


@dataclass
class SyntheticBenchmark:
    """A written synthetic benchmark set.

    Attributes:
        structure_dir: Directory containing the written PDB files.
        paths: Written structure paths, in specification order.
        specs: The specifications used.
    """

    structure_dir: Path
    paths: List[Path]
    specs: List[SyntheticStructureSpec] = field(default_factory=list)

    def ground_truth(self) -> Dict[str, str]:
        """Map structure name to its expected burial class."""
        return {spec.name: spec.expected_burial_class for spec in self.specs}


def build_synthetic_benchmark(
    out_dir: Path,
    *,
    n_buried: int = 8,
    n_surface: int = 8,
    n_decoy: int = 8,
    seed: int = 20240601,
    specs: Optional[Sequence[SyntheticStructureSpec]] = None,
) -> SyntheticBenchmark:
    """Write a synthetic benchmark set to disk.

    Args:
        out_dir: Destination directory for structure files.
        n_buried: Number of buried structures (ignored when ``specs`` is given).
        n_surface: Number of surface structures (ignored when ``specs`` is given).
        n_decoy: Number of ligand-free decoys (ignored when ``specs`` is given).
        seed: Base random seed.
        specs: Explicit specifications, overriding the defaults.

    Returns:
        A :class:`SyntheticBenchmark` describing what was written.
    """
    chosen = (
        list(specs)
        if specs is not None
        else default_benchmark_specs(
            n_buried=n_buried, n_surface=n_surface, n_decoy=n_decoy, seed=seed
        )
    )
    out_dir = Path(out_dir)
    paths = [write_synthetic_structure(spec, out_dir) for spec in chosen]
    return SyntheticBenchmark(structure_dir=out_dir, paths=paths, specs=chosen)
