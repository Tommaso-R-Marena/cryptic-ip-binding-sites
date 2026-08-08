"""Array view of a macromolecular structure for vectorised analysis.

Motivation
----------
Object-graph traversal of a parsed structure is convenient but makes two classes
of error easy, and both were present in earlier versions of this pipeline:

1. **Chain-blind residue keys.** Indexing residues by sequence number alone
   collides across chains. In a homodimer, residue 42 of chain A and residue 42
   of chain B become one entry, so per-residue SASA values overwrite each other
   and basic-residue counts are inflated by the number of chains. Every residue
   here is keyed by ``(model, chain, resseq, icode)``.
2. **Alternate-location double counting.** Structures with alternate
   conformations contribute several coordinate sets for the same atom. Summing
   over them inflates atom counts, volumes and SASA. This module keeps only the
   highest-occupancy altloc per atom name.

The resulting :class:`StructureArrays` exposes parallel NumPy arrays, so
downstream geometry is vectorised and unambiguous.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .geometry import element_radius

#: The 20 standard amino acids plus common modified residues that still form
#: part of the polymer backbone.
STANDARD_AMINO_ACIDS = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        "MSE", "SEC", "PYL", "HYP", "CSO", "CSD", "PTR", "SEP", "TPO", "MLY",
        "KCX", "ALY", "M3L", "HIC", "CME", "CAS", "ASX", "GLX", "UNK",
    }
)

#: Basic side chains capable of coordinating a phosphate group. Histidine is
#: included because it is protonated and phosphate-coordinating in several
#: inositol phosphate sites, but it is also tracked separately so its weaker,
#: pH-dependent contribution can be modelled apart from Arg/Lys.
BASIC_RESIDUES = frozenset({"ARG", "LYS", "HIS"})
STRONG_BASIC_RESIDUES = frozenset({"ARG", "LYS"})

#: Side-chain nitrogen atoms that donate hydrogen bonds to phosphate oxygens.
BASIC_NITROGEN_ATOMS = frozenset({"NZ", "NE", "NH1", "NH2", "ND1", "NE2"})

#: Hydroxyl-bearing residues; Ser/Thr/Tyr hydroxyls are frequent phosphate
#: hydrogen-bond donors alongside the basic side chains.
HYDROXYL_RESIDUES = frozenset({"SER", "THR", "TYR"})

#: Aromatic residues, which stack against the inositol ring.
AROMATIC_RESIDUES = frozenset({"PHE", "TYR", "TRP", "HIS"})

#: Acidic residues; their presence disfavours polyanion binding.
ACIDIC_RESIDUES = frozenset({"ASP", "GLU"})

#: Solvent and cryoprotectant residues excluded from every geometric context.
SOLVENT_RESNAMES = frozenset({"HOH", "DOD", "WAT", "D2O"})

#: Kyte-Doolittle hydropathy index (*J Mol Biol* 157:105-132, 1982).
KYTE_DOOLITTLE: Dict[str, float] = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}

#: Maximum residue SASA in the Gly-X-Gly reference state (Tien et al.,
#: *PLoS ONE* 8:e80635, 2013, theoretical values). Used to convert absolute
#: residue SASA into relative solvent accessibility, which is comparable across
#: residue types - a fully exposed tryptophan has far more surface than a fully
#: exposed glycine, so raw SASA is not.
MAX_RESIDUE_SASA: Dict[str, float] = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

#: Formal side-chain charge at pH 7.4 used for the Coulombic electrostatic
#: surrogate. Histidine carries a fractional charge from the Henderson-
#: Hasselbalch equation at pKa 6.0, which is more faithful than treating it as
#: either fully neutral or fully charged.
FORMAL_CHARGES: Dict[str, float] = {
    "ARG": 1.0,
    "LYS": 1.0,
    "HIS": 0.1,
    "ASP": -1.0,
    "GLU": -1.0,
}

ResidueKey = Tuple[int, str, int, str]


@dataclass
class StructureArrays:
    """Parallel arrays describing every retained atom of a structure.

    Attributes:
        coords: Atom coordinates, shape ``(n_atoms, 3)``.
        radii: Van der Waals radii in Å.
        elements: Upper-case element symbols.
        atom_names: Atom names as deposited.
        resnames: Upper-case residue names.
        chain_ids: Author chain identifiers.
        resseqs: Author residue sequence numbers.
        icodes: Insertion codes (``""`` when absent).
        model_ids: Model identifiers.
        occupancies: Atom occupancies.
        bfactors: Atom B-factors (pLDDT for AlphaFold models).
        is_polymer: ``True`` for standard polymer residues.
        is_solvent: ``True`` for water and related solvent residues.
        residue_index: Index into :attr:`residue_keys` for each atom.
        residue_keys: Ordered unique residue keys.
        source_path: Path the structure was read from.
    """

    coords: np.ndarray
    radii: np.ndarray
    elements: np.ndarray
    atom_names: np.ndarray
    resnames: np.ndarray
    chain_ids: np.ndarray
    resseqs: np.ndarray
    icodes: np.ndarray
    model_ids: np.ndarray
    occupancies: np.ndarray
    bfactors: np.ndarray
    is_polymer: np.ndarray
    is_solvent: np.ndarray
    residue_index: np.ndarray
    residue_keys: List[ResidueKey]
    source_path: Optional[str] = None

    @property
    def n_atoms(self) -> int:
        """Number of retained atoms."""
        return int(self.coords.shape[0])

    @property
    def is_hetero(self) -> np.ndarray:
        """``True`` for non-polymer, non-solvent atoms (ligands, ions)."""
        return ~self.is_polymer & ~self.is_solvent

    def mask_resnames(self, resnames: Iterable[str]) -> np.ndarray:
        """Return a boolean atom mask selecting the given residue names.

        Args:
            resnames: Residue names to select (case-insensitive).

        Returns:
            Boolean array of length :attr:`n_atoms`.
        """
        wanted = {str(name).upper() for name in resnames}
        if not wanted:
            return np.zeros(self.n_atoms, dtype=bool)
        return np.isin(self.resnames, list(wanted))

    def mask_residue(self, key: ResidueKey) -> np.ndarray:
        """Return a boolean atom mask selecting one residue by key.

        Args:
            key: ``(model, chain, resseq, icode)``.

        Returns:
            Boolean array of length :attr:`n_atoms`.
        """
        try:
            index = self.residue_keys.index(key)
        except ValueError:
            return np.zeros(self.n_atoms, dtype=bool)
        return self.residue_index == index

    def residue_atom_indices(self) -> Dict[ResidueKey, np.ndarray]:
        """Map each residue key to the indices of its atoms."""
        out: Dict[ResidueKey, np.ndarray] = {}
        for index, key in enumerate(self.residue_keys):
            out[key] = np.flatnonzero(self.residue_index == index)
        return out

    def residue_names_by_key(self) -> Dict[ResidueKey, str]:
        """Map each residue key to its residue name."""
        out: Dict[ResidueKey, str] = {}
        for index, key in enumerate(self.residue_keys):
            atoms = np.flatnonzero(self.residue_index == index)
            if atoms.size:
                out[key] = str(self.resnames[atoms[0]])
        return out

    def subset(self, mask: np.ndarray) -> "StructureArrays":
        """Return a new :class:`StructureArrays` restricted to ``mask``.

        Residue keys are recomputed so the result remains self-consistent.

        Args:
            mask: Boolean atom mask.

        Returns:
            The restricted structure view.
        """
        mask = np.asarray(mask, dtype=bool)
        keys_in_order: List[ResidueKey] = []
        remap: Dict[int, int] = {}
        new_residue_index = np.empty(int(mask.sum()), dtype=int)
        for new_pos, old_pos in enumerate(np.flatnonzero(mask)):
            old_res = int(self.residue_index[old_pos])
            if old_res not in remap:
                remap[old_res] = len(keys_in_order)
                keys_in_order.append(self.residue_keys[old_res])
            new_residue_index[new_pos] = remap[old_res]

        return StructureArrays(
            coords=self.coords[mask],
            radii=self.radii[mask],
            elements=self.elements[mask],
            atom_names=self.atom_names[mask],
            resnames=self.resnames[mask],
            chain_ids=self.chain_ids[mask],
            resseqs=self.resseqs[mask],
            icodes=self.icodes[mask],
            model_ids=self.model_ids[mask],
            occupancies=self.occupancies[mask],
            bfactors=self.bfactors[mask],
            is_polymer=self.is_polymer[mask],
            is_solvent=self.is_solvent[mask],
            residue_index=new_residue_index,
            residue_keys=keys_in_order,
            source_path=self.source_path,
        )


def _guess_element(atom_name: str, declared: Optional[str]) -> str:
    """Infer an element symbol, falling back to the atom name.

    Args:
        atom_name: Atom name as deposited.
        declared: Element column value, if present.

    Returns:
        Upper-case element symbol (possibly ``""``).
    """
    if declared and declared.strip():
        return declared.strip().upper()
    name = (atom_name or "").strip()
    if not name:
        return ""
    # PDB atom names put the element in columns 13-14; a leading digit is a
    # branch index, not an element.
    if name[0].isdigit():
        name = name[1:]
    return name[:1].upper()


def load_structure_arrays(
    path: Path,
    *,
    model_index: int = 0,
    keep_solvent: bool = False,
    keep_hydrogens: bool = False,
) -> StructureArrays:
    """Parse a structure file into :class:`StructureArrays`.

    Args:
        path: Path to a PDB or mmCIF file.
        model_index: Which model to read. NMR ensembles and multi-model
            predictions contain several; mixing them would superimpose distinct
            conformers into one atom cloud, so exactly one is used.
        keep_solvent: Retain water molecules. Excluded by default because
            crystallographic waters are partially modelled and would occlude
            ligand surface inconsistently between entries.
        keep_hydrogens: Retain hydrogen atoms. Excluded by default so that
            structures with and without riding hydrogens give comparable SASA.

    Returns:
        The parsed array view.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        ValueError: If the file contains no usable atoms.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Structure file not found: {path}")

    structure = _parse_structure(path)
    models = list(structure)
    if not models:
        raise ValueError(f"No models found in {path}")
    model = models[min(model_index, len(models) - 1)]

    coords: List[Sequence[float]] = []
    radii: List[float] = []
    elements: List[str] = []
    atom_names: List[str] = []
    resnames: List[str] = []
    chain_ids: List[str] = []
    resseqs: List[int] = []
    icodes: List[str] = []
    model_ids: List[int] = []
    occupancies: List[float] = []
    bfactors: List[float] = []
    is_polymer: List[bool] = []
    is_solvent: List[bool] = []
    residue_index: List[int] = []
    residue_keys: List[ResidueKey] = []

    model_id = int(getattr(model, "id", 0))
    for chain in model:
        chain_id = str(getattr(chain, "id", "A"))
        for residue in chain:
            het_flag, resseq, icode = residue.id
            resname = str(residue.get_resname()).strip().upper()
            solvent = resname in SOLVENT_RESNAMES
            if solvent and not keep_solvent:
                continue
            polymer = het_flag == " " and resname in STANDARD_AMINO_ACIDS

            key: ResidueKey = (model_id, chain_id, int(resseq), str(icode).strip())
            selected = _select_altloc_atoms(residue)
            if not selected:
                continue
            residue_slot = len(residue_keys)
            residue_keys.append(key)

            for atom in selected:
                element = _guess_element(atom.get_name(), getattr(atom, "element", None))
                if element in {"H", "D"} and not keep_hydrogens:
                    continue
                coords.append(atom.get_coord())
                radii.append(element_radius(element))
                elements.append(element)
                atom_names.append(str(atom.get_name()).strip())
                resnames.append(resname)
                chain_ids.append(chain_id)
                resseqs.append(int(resseq))
                icodes.append(str(icode).strip())
                model_ids.append(model_id)
                occupancy = atom.get_occupancy()
                occupancies.append(1.0 if occupancy is None else float(occupancy))
                bfactor = atom.get_bfactor()
                bfactors.append(0.0 if bfactor is None else float(bfactor))
                is_polymer.append(polymer)
                is_solvent.append(solvent)
                residue_index.append(residue_slot)

    if not coords:
        raise ValueError(f"No usable atoms parsed from {path}")

    # A residue whose every atom was filtered out (e.g. hydrogen-only) leaves an
    # unreferenced key. Slots are assigned monotonically, so compacting the key
    # list and remapping keeps residue_index a valid, gap-free index array.
    used = sorted(set(residue_index))
    remap = {old: new for new, old in enumerate(used)}
    compact_keys = [residue_keys[old] for old in used]

    return StructureArrays(
        coords=np.asarray(coords, dtype=float),
        radii=np.asarray(radii, dtype=float),
        elements=np.asarray(elements, dtype=object),
        atom_names=np.asarray(atom_names, dtype=object),
        resnames=np.asarray(resnames, dtype=object),
        chain_ids=np.asarray(chain_ids, dtype=object),
        resseqs=np.asarray(resseqs, dtype=int),
        icodes=np.asarray(icodes, dtype=object),
        model_ids=np.asarray(model_ids, dtype=int),
        occupancies=np.asarray(occupancies, dtype=float),
        bfactors=np.asarray(bfactors, dtype=float),
        is_polymer=np.asarray(is_polymer, dtype=bool),
        is_solvent=np.asarray(is_solvent, dtype=bool),
        residue_index=np.asarray([remap[i] for i in residue_index], dtype=int),
        residue_keys=compact_keys,
        source_path=str(path),
    )


def _parse_structure(path: Path):
    """Parse a coordinate file with the parser matching its extension."""
    from Bio.PDB import MMCIFParser, PDBParser

    suffixes = {suffix.lower() for suffix in path.suffixes}
    if {".cif", ".mmcif"} & suffixes:
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def _select_altloc_atoms(residue) -> List:
    """Return one atom per atom name, keeping the highest-occupancy altloc.

    Args:
        residue: A parsed residue.

    Returns:
        Atoms with unique names.
    """
    best: Dict[str, object] = {}
    for atom in residue.get_unpacked_list() if hasattr(residue, "get_unpacked_list") else residue:
        name = str(atom.get_name()).strip()
        occupancy = atom.get_occupancy()
        occupancy = 1.0 if occupancy is None else float(occupancy)
        current = best.get(name)
        if current is None:
            best[name] = atom
            continue
        current_occ = current.get_occupancy()  # type: ignore[attr-defined]
        current_occ = 1.0 if current_occ is None else float(current_occ)
        if occupancy > current_occ:
            best[name] = atom
    return list(best.values())


def phosphate_group_indices(
    arrays: StructureArrays,
    atom_mask: np.ndarray,
    *,
    p_o_bond_cutoff: float = 1.9,
) -> np.ndarray:
    """Identify phosphate-group atoms within a masked atom set.

    Phosphate atoms are found by **element and bonding geometry**, not by atom
    name. Name-prefix matching is unreliable: an inositol phosphate's ring
    oxygens, ester oxygens and terminal phosphate oxygens all have names
    beginning with ``O``, and a naive ``startswith("P")`` test also matches
    ``PA``/``PB`` bridging atoms of pyrophosphates while missing unconventional
    names. Selecting phosphorus atoms by element and then taking the oxygens
    within a P-O bond distance identifies exactly the phosphate groups.

    Args:
        arrays: Structure arrays.
        atom_mask: Boolean mask restricting the search (e.g. one ligand copy).
        p_o_bond_cutoff: Maximum P-O distance treated as a bond (Å). P-O bonds
            are 1.5-1.6 Å; 1.9 Å tolerates refinement error without reaching
            non-bonded oxygens.

    Returns:
        Sorted indices of phosphorus and phosphate-oxygen atoms.
    """
    from scipy.spatial import cKDTree

    mask = np.asarray(atom_mask, dtype=bool)
    candidate_indices = np.flatnonzero(mask)
    if candidate_indices.size == 0:
        return np.empty(0, dtype=int)

    elements = arrays.elements[candidate_indices]
    p_local = np.flatnonzero(elements == "P")
    if p_local.size == 0:
        return np.empty(0, dtype=int)

    o_local = np.flatnonzero(elements == "O")
    selected = {int(candidate_indices[i]) for i in p_local}
    if o_local.size:
        tree = cKDTree(arrays.coords[candidate_indices[o_local]])
        for i in p_local:
            for hit in tree.query_ball_point(
                arrays.coords[candidate_indices[i]], float(p_o_bond_cutoff)
            ):
                selected.add(int(candidate_indices[o_local[hit]]))
    return np.asarray(sorted(selected), dtype=int)
