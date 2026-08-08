"""Burial metrics for inositol phosphate ligands, measured per ligand copy.

What changed and why it matters
-------------------------------
The previous implementation summed solvent accessible surface area over *every
copy* of a ligand residue name in an entry and classified the entry from that
sum. Two consequences followed, and both were visible in the shipped dataset:

* A structure with six InsP6 copies reported roughly six times the SASA of a
  structure with one, so copy number - a crystallisation artefact - determined
  the burial class. In the previously shipped ground-truth table this pushed
  almost every entry into ``Surface``, leaving 5 positive pockets out of 12 190
  and a classifier at chance performance (AUROC 0.50).
* Absolute SASA is not comparable across ligands of different size: InsP3 has
  far less surface than InsP6, so the same threshold means different things.

This module fixes both by measuring **each ligand copy separately** and
normalising:

``relative_sasa = SASA(ligand copy in complex) / SASA(same copy in isolation)``

That ratio is the fraction of the ligand's own surface still reachable by
solvent. It is dimensionless, independent of ligand size and copy count, and
directly expresses the quantity of interest - the ADAR2 InsP6 site sits near
0.0, a PH-domain surface site near 0.5-0.8.

Four complementary measurements are reported per copy, because no single number
distinguishes a buried cofactor site from a deep surface groove:

1. ``relative_sasa`` - how much ligand surface solvent still reaches.
2. ``relative_phosphate_sasa`` - the same for phosphate groups only, which is
   what a polyanion-binding site must sequester.
3. ``burial_depth`` - how far below the molecular surface the copy sits.
4. ``enclosure`` - what fraction of directions out of the site are blocked.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from ..analysis.geometry import (
    DEFAULT_N_POINTS,
    burial_depth as _burial_depth,
    count_within,
    enclosure_fraction,
    point_cloud_shape,
    shrake_rupley_sasa,
)
from ..analysis.structure_arrays import (
    BASIC_NITROGEN_ATOMS,
    BASIC_RESIDUES,
    HYDROXYL_RESIDUES,
    ResidueKey,
    StructureArrays,
    load_structure_arrays,
    phosphate_group_indices,
)
from .structure_context import LIGAND_RESNAMES

LOGGER = logging.getLogger(__name__)

#: Burial classification thresholds on ``relative_sasa`` (fraction of the
#: ligand's own surface still solvent-accessible).
#:
#: The boundaries are anchored on the paradigm cases rather than chosen for
#: convenience. Macbeth et al. (*Science* 309:1534-1539, 2005) describe the
#: ADAR2 InsP6 as completely encapsulated, with only a narrow window to the
#: exterior - a few per cent of its surface. Surface signalling sites such as the
#: PLC-delta-1 PH domain InsP3 complex leave roughly half the ligand solvent
#: exposed. A ``semi_cryptic`` band between the two captures interface sites
#: (e.g. inositol phosphate at a subunit interface), which are mechanistically
#: intermediate and are analysed separately rather than being forced into one of
#: the extremes.
CRYPTIC_RELATIVE_SASA_MAX = 0.05
SEMI_CRYPTIC_RELATIVE_SASA_MAX = 0.25

#: Legacy absolute-SASA thresholds (Å²) retained for the historical
#: :func:`classify_burial` signature.
CRYPTIC_ABSOLUTE_SASA_MAX = 5.0
SEMI_CRYPTIC_ABSOLUTE_SASA_MAX = 50.0

#: A ligand copy making fewer than this many protein heavy-atom contacts within
#: :data:`CONTACT_CUTOFF` is unlikely to occupy a genuine binding site and is
#: flagged as a probable crystallisation additive rather than scored.
MIN_CONTACTS_FOR_SITE = 8
CONTACT_CUTOFF = 4.5

#: Distance within which a basic side-chain nitrogen is treated as coordinating.
COORDINATION_CUTOFF = 4.0

#: Per-atom SASA above which an atom counts as part of the molecular surface for
#: burial-depth measurement (Å²). See :func:`cryptic_ip.analysis.geometry.burial_depth`.
EXPOSED_ATOM_SASA = 5.0

BurialClass = str


@dataclass
class LigandInstanceBurial:
    """Burial measurements for a single ligand copy.

    Attributes:
        comp_id: Ligand chemical component identifier.
        chain_id: Author chain identifier of the copy.
        resseq: Author residue sequence number of the copy.
        icode: Insertion code (``""`` when absent).
        model_id: Model identifier.
        n_atoms: Heavy-atom count of the copy.
        n_phosphorus: Phosphorus atom count.
        sasa_complex: SASA of the copy inside the complex (Å²).
        sasa_isolated: SASA of the same coordinates in isolation (Å²).
        relative_sasa: ``sasa_complex / sasa_isolated`` in ``[0, 1]``.
        phosphate_sasa_complex: SASA of phosphate-group atoms in the complex (Å²).
        phosphate_sasa_isolated: The same in isolation (Å²).
        relative_phosphate_sasa: Ratio of the two.
        burial_depth: Distance from the copy centroid to the nearest exposed
            atom (Å).
        enclosure: Fraction of directions out of the site that are blocked.
        n_protein_contacts: Protein heavy atoms within :data:`CONTACT_CUTOFF`.
        n_contact_chains: Distinct chains contacting the copy.
        n_basic_residues: Basic residues with a side-chain nitrogen within
            :data:`COORDINATION_CUTOFF` of a phosphate oxygen.
        n_basic_nitrogens: Individual coordinating side-chain nitrogens.
        n_hydroxyl_residues: Ser/Thr/Tyr residues in contact.
        mean_occupancy: Mean atom occupancy of the copy.
        mean_bfactor: Mean atom B-factor of the copy.
        radius_of_gyration: Radius of gyration of the copy (Å).
        burial_class: One of ``cryptic``, ``semi_cryptic``, ``surface``,
            ``crystal_artifact`` or ``unknown``.
        is_probable_artifact: ``True`` when contacts are too few for a site.
    """

    comp_id: str
    chain_id: str
    resseq: int
    icode: str
    model_id: int
    n_atoms: int
    n_phosphorus: int
    sasa_complex: float
    sasa_isolated: float
    relative_sasa: float
    phosphate_sasa_complex: float
    phosphate_sasa_isolated: float
    relative_phosphate_sasa: float
    burial_depth: float
    enclosure: float
    n_protein_contacts: int
    n_contact_chains: int
    n_basic_residues: int
    n_basic_nitrogens: int
    n_hydroxyl_residues: int
    mean_occupancy: float
    mean_bfactor: float
    radius_of_gyration: float
    burial_class: BurialClass
    is_probable_artifact: bool

    @property
    def instance_id(self) -> str:
        """Stable identifier of the copy, e.g. ``"IHP_A_501"``."""
        suffix = f"{self.resseq}{self.icode}".strip()
        return f"{self.comp_id}_{self.chain_id}_{suffix}"

    @property
    def centroid_key(self) -> ResidueKey:
        """The residue key of this copy."""
        return (self.model_id, self.chain_id, self.resseq, self.icode)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON/CSV-serialisable representation."""
        payload = asdict(self)
        payload["instance_id"] = self.instance_id
        return payload


@dataclass
class BurialMetrics:
    """Structure-level burial summary.

    The legacy field names are preserved so existing callers keep working, but
    the values are now taken from the **most buried ligand copy** rather than
    from a sum over copies. The most buried copy is the biologically meaningful
    one: an entry demonstrates a cryptic site if *any* copy is cryptic, whereas
    an average is dominated by surface-bound copies of the same additive.

    Attributes:
        ligand_sasa: Absolute SASA of the most buried copy (Å²).
        phosphate_sasa: Absolute phosphate SASA of that copy (Å²).
        delta_sasa: Ligand surface buried on complex formation (Å²), i.e.
            ``sasa_isolated - sasa_complex``.
        burial_depth: Burial depth of that copy (Å).
        burial_class: Burial class of that copy.
        relative_sasa: Relative SASA of that copy.
        relative_phosphate_sasa: Relative phosphate SASA of that copy.
        enclosure: Enclosure of that copy.
        n_basic_residues: Coordinating basic residues of that copy.
        n_instances: Total ligand copies measured.
        n_cryptic_instances: Copies classified ``cryptic``.
        instances: All per-copy measurements.
    """

    ligand_sasa: Optional[float]
    phosphate_sasa: Optional[float]
    delta_sasa: Optional[float]
    burial_depth: Optional[float]
    burial_class: BurialClass
    relative_sasa: Optional[float] = None
    relative_phosphate_sasa: Optional[float] = None
    enclosure: Optional[float] = None
    n_basic_residues: int = 0
    n_instances: int = 0
    n_cryptic_instances: int = 0
    instances: List[LigandInstanceBurial] = field(default_factory=list)

    @property
    def most_buried(self) -> Optional[LigandInstanceBurial]:
        """The copy with the lowest relative SASA, if any."""
        return self.instances[0] if self.instances else None

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable representation."""
        return {
            "ligand_sasa": self.ligand_sasa,
            "phosphate_sasa": self.phosphate_sasa,
            "delta_sasa": self.delta_sasa,
            "burial_depth": self.burial_depth,
            "burial_class": self.burial_class,
            "relative_sasa": self.relative_sasa,
            "relative_phosphate_sasa": self.relative_phosphate_sasa,
            "enclosure": self.enclosure,
            "n_basic_residues": self.n_basic_residues,
            "n_instances": self.n_instances,
            "n_cryptic_instances": self.n_cryptic_instances,
            "instances": [inst.to_dict() for inst in self.instances],
        }


def classify_burial_relative(
    relative_sasa: Optional[float],
    *,
    relative_phosphate_sasa: Optional[float] = None,
    is_probable_artifact: bool = False,
) -> BurialClass:
    """Classify burial from relative (size-normalised) SASA.

    When both are available the **larger** of the whole-ligand and phosphate
    ratios is used. That is the conservative choice: a site only counts as
    cryptic when neither the inositol ring nor its phosphates remain exposed, so
    a ligand buried to the ring but with phosphates projecting into solvent is
    not mistaken for a sequestered cofactor.

    Args:
        relative_sasa: Fraction of ligand surface accessible in the complex.
        relative_phosphate_sasa: The same for phosphate groups.
        is_probable_artifact: When ``True``, return ``"crystal_artifact"``.

    Returns:
        A burial class label.
    """
    if is_probable_artifact:
        return "crystal_artifact"

    values = [
        value
        for value in (relative_sasa, relative_phosphate_sasa)
        if value is not None and np.isfinite(value)
    ]
    if not values:
        return "unknown"
    metric = float(max(values))
    if metric <= CRYPTIC_RELATIVE_SASA_MAX:
        return "cryptic"
    if metric <= SEMI_CRYPTIC_RELATIVE_SASA_MAX:
        return "semi_cryptic"
    return "surface"


def classify_burial(
    ligand_sasa: Optional[float], phosphate_sasa: Optional[float] = None
) -> BurialClass:
    """Classify burial from absolute SASA (legacy interface).

    Retained for backward compatibility with callers and stored datasets that
    only have absolute SASA. New code should prefer
    :func:`classify_burial_relative`, which is size- and copy-number invariant.

    Args:
        ligand_sasa: Absolute ligand SASA in Å².
        phosphate_sasa: Absolute phosphate SASA in Å².

    Returns:
        A burial class label.
    """
    metric = (
        phosphate_sasa
        if phosphate_sasa is not None and np.isfinite(phosphate_sasa)
        else ligand_sasa
    )
    if metric is None or not np.isfinite(metric):
        return "unknown"
    if metric <= CRYPTIC_ABSOLUTE_SASA_MAX:
        return "cryptic"
    if metric <= SEMI_CRYPTIC_ABSOLUTE_SASA_MAX:
        return "semi_cryptic"
    return "surface"


def find_ligand_instances(
    arrays: StructureArrays, comp_ids: Optional[Iterable[str]] = None
) -> List[Tuple[ResidueKey, str, np.ndarray]]:
    """Locate every copy of the requested ligand components.

    Args:
        arrays: Parsed structure arrays.
        comp_ids: Component identifiers to match. Defaults to
            :data:`cryptic_ip.validation.structure_context.LIGAND_RESNAMES`.

    Returns:
        ``(residue_key, comp_id, atom_indices)`` per copy, ordered by key.
    """
    wanted = {str(c).upper() for c in (comp_ids or LIGAND_RESNAMES)}
    out: List[Tuple[ResidueKey, str, np.ndarray]] = []
    for slot, key in enumerate(arrays.residue_keys):
        atom_indices = np.flatnonzero(arrays.residue_index == slot)
        if atom_indices.size == 0:
            continue
        resname = str(arrays.resnames[atom_indices[0]])
        if resname in wanted:
            out.append((key, resname, atom_indices))
    return out


def compute_instance_burial(
    arrays: StructureArrays,
    ligand_key: ResidueKey,
    comp_id: str,
    ligand_atoms: np.ndarray,
    *,
    other_ligand_atoms: Optional[np.ndarray] = None,
    n_points: int = DEFAULT_N_POINTS,
) -> LigandInstanceBurial:
    """Measure burial for one ligand copy.

    The measurement context is protein plus non-solvent heteroatoms, **excluding
    other copies of the inositol phosphate ligands**. Excluding sibling copies
    matters because in crystals they frequently pack against one another; if a
    neighbouring copy is allowed to occlude, a surface-bound ligand can appear
    buried by an artefact of packing rather than by the protein.

    Args:
        arrays: Structure arrays for the whole entry.
        ligand_key: Residue key of the copy.
        comp_id: Component identifier of the copy.
        ligand_atoms: Atom indices of the copy.
        other_ligand_atoms: Atom indices of sibling ligand copies to exclude.
        n_points: SASA sample points per atom.

    Returns:
        The per-copy burial record.
    """
    n_atoms_total = arrays.n_atoms
    exclude = np.zeros(n_atoms_total, dtype=bool)
    if other_ligand_atoms is not None and len(other_ligand_atoms):
        exclude[np.asarray(other_ligand_atoms, dtype=int)] = True
    exclude[ligand_atoms] = False

    context_mask = ~arrays.is_solvent & ~exclude
    context_indices = np.flatnonzero(context_mask)
    # Position of the ligand atoms inside the context sub-array.
    position_of = {int(atom): pos for pos, atom in enumerate(context_indices)}
    ligand_positions = np.asarray(
        [position_of[int(atom)] for atom in ligand_atoms if int(atom) in position_of], dtype=int
    )

    context_coords = arrays.coords[context_indices]
    context_radii = arrays.radii[context_indices]
    ligand_coords = arrays.coords[ligand_atoms]
    ligand_radii = arrays.radii[ligand_atoms]

    sasa_complex_atoms = shrake_rupley_sasa(
        context_coords, context_radii, subset=ligand_positions, n_points=n_points
    )
    sasa_isolated_atoms = shrake_rupley_sasa(ligand_coords, ligand_radii, n_points=n_points)
    sasa_complex = float(np.sum(sasa_complex_atoms))
    sasa_isolated = float(np.sum(sasa_isolated_atoms))
    relative_sasa = sasa_complex / sasa_isolated if sasa_isolated > 0 else float("nan")

    # Phosphate subset, resolved by element and P-O bonding geometry.
    ligand_mask = np.zeros(n_atoms_total, dtype=bool)
    ligand_mask[ligand_atoms] = True
    phosphate_atoms = phosphate_group_indices(arrays, ligand_mask)
    if phosphate_atoms.size:
        phosphate_positions = np.asarray(
            [position_of[int(a)] for a in phosphate_atoms if int(a) in position_of], dtype=int
        )
        local = {int(atom): pos for pos, atom in enumerate(ligand_atoms)}
        phosphate_local = np.asarray(
            [local[int(a)] for a in phosphate_atoms if int(a) in local], dtype=int
        )
        phosphate_complex = float(
            np.sum(
                shrake_rupley_sasa(
                    context_coords, context_radii, subset=phosphate_positions, n_points=n_points
                )
            )
        )
        phosphate_isolated = float(np.sum(sasa_isolated_atoms[phosphate_local]))
        relative_phosphate = (
            phosphate_complex / phosphate_isolated if phosphate_isolated > 0 else float("nan")
        )
    else:
        phosphate_complex = float("nan")
        phosphate_isolated = float("nan")
        relative_phosphate = float("nan")

    centroid = ligand_coords.mean(axis=0)

    # Surface definition uses the holo complex so the ligand-filled cavity is
    # not mistaken for solvent-accessible space (see geometry.burial_depth).
    protein_positions = np.flatnonzero(arrays.is_polymer[context_indices])
    protein_sasa = shrake_rupley_sasa(
        context_coords, context_radii, subset=protein_positions, n_points=n_points
    )
    depth = _burial_depth(
        centroid, context_coords[protein_positions], protein_sasa, min_atom_sasa=EXPOSED_ATOM_SASA
    )
    enclosure = enclosure_fraction(
        centroid, context_coords[protein_positions], context_radii[protein_positions]
    )

    protein_atom_indices = context_indices[protein_positions]
    n_contacts, contact_rows = count_within(
        ligand_coords, arrays.coords[protein_atom_indices], CONTACT_CUTOFF
    )
    contact_atom_indices = protein_atom_indices[contact_rows] if contact_rows.size else contact_rows
    contact_chains = (
        set(arrays.chain_ids[contact_atom_indices].tolist()) if contact_atom_indices.size else set()
    )

    n_basic_residues, n_basic_nitrogens, n_hydroxyl = _coordination_counts(
        arrays, ligand_atoms, phosphate_atoms, protein_atom_indices
    )

    is_artifact = n_contacts < MIN_CONTACTS_FOR_SITE
    burial_class = classify_burial_relative(
        relative_sasa,
        relative_phosphate_sasa=relative_phosphate,
        is_probable_artifact=is_artifact,
    )

    return LigandInstanceBurial(
        comp_id=comp_id,
        chain_id=str(ligand_key[1]),
        resseq=int(ligand_key[2]),
        icode=str(ligand_key[3]),
        model_id=int(ligand_key[0]),
        n_atoms=int(ligand_atoms.size),
        n_phosphorus=int(np.count_nonzero(arrays.elements[ligand_atoms] == "P")),
        sasa_complex=sasa_complex,
        sasa_isolated=sasa_isolated,
        relative_sasa=float(relative_sasa),
        phosphate_sasa_complex=phosphate_complex,
        phosphate_sasa_isolated=phosphate_isolated,
        relative_phosphate_sasa=float(relative_phosphate),
        burial_depth=float(depth),
        enclosure=float(enclosure),
        n_protein_contacts=int(n_contacts),
        n_contact_chains=len(contact_chains),
        n_basic_residues=int(n_basic_residues),
        n_basic_nitrogens=int(n_basic_nitrogens),
        n_hydroxyl_residues=int(n_hydroxyl),
        mean_occupancy=float(np.mean(arrays.occupancies[ligand_atoms])),
        mean_bfactor=float(np.mean(arrays.bfactors[ligand_atoms])),
        radius_of_gyration=point_cloud_shape(ligand_coords).radius_of_gyration,
        burial_class=burial_class,
        is_probable_artifact=bool(is_artifact),
    )


def _coordination_counts(
    arrays: StructureArrays,
    ligand_atoms: np.ndarray,
    phosphate_atoms: np.ndarray,
    protein_atom_indices: np.ndarray,
) -> Tuple[int, int, int]:
    """Count coordinating basic residues, basic nitrogens and hydroxyl residues.

    Coordination is measured from **phosphate oxygens** where available rather
    than from the whole ligand, since it is the phosphate groups that accept the
    hydrogen bonds and salt bridges. Residues are identified by chain-aware keys,
    so the same sequence number in different chains counts once per chain
    instead of collapsing into a single residue.

    Args:
        arrays: Structure arrays.
        ligand_atoms: Atom indices of the ligand copy.
        phosphate_atoms: Atom indices of its phosphate groups.
        protein_atom_indices: Candidate protein atom indices.

    Returns:
        ``(n_basic_residues, n_basic_nitrogens, n_hydroxyl_residues)``.
    """
    from scipy.spatial import cKDTree

    if protein_atom_indices.size == 0:
        return 0, 0, 0

    anchor_atoms = phosphate_atoms if phosphate_atoms.size else ligand_atoms
    anchor_coords = arrays.coords[anchor_atoms]
    tree = cKDTree(anchor_coords)

    basic_residue_keys = set()
    hydroxyl_residue_keys = set()
    n_basic_nitrogens = 0

    for atom_index in protein_atom_indices:
        resname = str(arrays.resnames[atom_index])
        if resname not in BASIC_RESIDUES and resname not in HYDROXYL_RESIDUES:
            continue
        neighbours = tree.query_ball_point(arrays.coords[atom_index], COORDINATION_CUTOFF)
        if not neighbours:
            continue
        key = (
            int(arrays.model_ids[atom_index]),
            str(arrays.chain_ids[atom_index]),
            int(arrays.resseqs[atom_index]),
            str(arrays.icodes[atom_index]),
        )
        atom_name = str(arrays.atom_names[atom_index]).upper()
        if resname in BASIC_RESIDUES:
            if atom_name in BASIC_NITROGEN_ATOMS:
                basic_residue_keys.add(key)
                n_basic_nitrogens += 1
        if resname in HYDROXYL_RESIDUES and atom_name in {"OG", "OG1", "OH"}:
            hydroxyl_residue_keys.add(key)

    return len(basic_residue_keys), n_basic_nitrogens, len(hydroxyl_residue_keys)


def compute_ligand_burial(
    path: Path,
    *,
    comp_ids: Optional[Sequence[str]] = None,
    n_points: int = DEFAULT_N_POINTS,
) -> List[LigandInstanceBurial]:
    """Measure burial for every inositol phosphate copy in a structure file.

    Args:
        path: Path to a PDB or mmCIF file.
        comp_ids: Component identifiers to measure. Defaults to the built-in
            inositol phosphate residue names.
        n_points: SASA sample points per atom.

    Returns:
        Per-copy records sorted by ascending relative SASA, i.e. most buried
        first.
    """
    arrays = load_structure_arrays(Path(path))
    instances = find_ligand_instances(arrays, comp_ids)
    if not instances:
        return []

    all_ligand_atoms = np.concatenate([atoms for _, _, atoms in instances])
    records: List[LigandInstanceBurial] = []
    for key, comp_id, atoms in instances:
        siblings = np.setdiff1d(all_ligand_atoms, atoms, assume_unique=False)
        try:
            records.append(
                compute_instance_burial(
                    arrays, key, comp_id, atoms, other_ligand_atoms=siblings, n_points=n_points
                )
            )
        except Exception as exc:  # noqa: BLE001 - one bad copy must not abort the entry
            LOGGER.warning("Burial measurement failed for %s in %s: %s", comp_id, path, exc)

    records.sort(
        key=lambda rec: (
            rec.relative_sasa if np.isfinite(rec.relative_sasa) else np.inf,
            rec.instance_id,
        )
    )
    return records


def compute_burial_metrics(
    path: Path,
    *,
    comp_ids: Optional[Sequence[str]] = None,
    n_points: int = DEFAULT_N_POINTS,
) -> BurialMetrics:
    """Summarise burial for a structure, reporting the most buried ligand copy.

    Args:
        path: Path to a PDB or mmCIF file.
        comp_ids: Component identifiers to measure.
        n_points: SASA sample points per atom.

    Returns:
        A :class:`BurialMetrics` summary. When the structure contains no matching
        ligand, all numeric fields are ``None`` and ``burial_class`` is
        ``"unknown"``.
    """
    instances = compute_ligand_burial(path, comp_ids=comp_ids, n_points=n_points)
    if not instances:
        return BurialMetrics(
            ligand_sasa=None,
            phosphate_sasa=None,
            delta_sasa=None,
            burial_depth=None,
            burial_class="unknown",
            instances=[],
        )

    best = instances[0]
    delta = (
        best.sasa_isolated - best.sasa_complex
        if np.isfinite(best.sasa_isolated) and np.isfinite(best.sasa_complex)
        else None
    )
    return BurialMetrics(
        ligand_sasa=best.sasa_complex,
        phosphate_sasa=(
            best.phosphate_sasa_complex if np.isfinite(best.phosphate_sasa_complex) else None
        ),
        delta_sasa=delta,
        burial_depth=best.burial_depth if np.isfinite(best.burial_depth) else None,
        burial_class=best.burial_class,
        relative_sasa=best.relative_sasa if np.isfinite(best.relative_sasa) else None,
        relative_phosphate_sasa=(
            best.relative_phosphate_sasa if np.isfinite(best.relative_phosphate_sasa) else None
        ),
        enclosure=best.enclosure if np.isfinite(best.enclosure) else None,
        n_basic_residues=best.n_basic_residues,
        n_instances=len(instances),
        n_cryptic_instances=sum(1 for inst in instances if inst.burial_class == "cryptic"),
        instances=instances,
    )
