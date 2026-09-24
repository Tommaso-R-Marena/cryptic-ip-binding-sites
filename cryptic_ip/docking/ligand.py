"""Ligands for docking, built from the Chemical Component Dictionary (docs/REDOCKING_PLAN.md).

The docked ligand is never built from crystal coordinates: those carry the
answer. It is built from the component's CCD SMILES, embedded by RDKit (ETKDGv3,
fixed seed), rotated at random, and checked to start more than 2 Å from the
crystal pose. The crystal copy is used only as the RMSD reference, after its
bond orders are assigned from the same CCD template and its configuration is
checked against the template's.

Protonation (the plan's two states):

* ``primary`` - every terminal phosphate oxygen deprotonated, then
  ``floor(n_P / 2)`` protons put back, one per phosphorus, on the phosphorus
  atoms of lowest canonical rank, each on that phosphorus's single-bonded
  terminal oxygen of lowest canonical rank. InsP6 carries -9, InsP5 -8, InsP4 -6,
  InsP3 -5: close to the -9 to -10 measured for InsP6 near pH 7.4.
* ``deprotonated`` - every terminal phosphate oxygen deprotonated (InsP6 -12).

Vina's scoring has no electrostatic term and ignores partial charges; the
states differ only by which oxygens carry a donor hydrogen. They matter more
for AD4, whose scoring uses charges.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .rmsd import symmetric_rmsd, terminal_phosphate_oxygens

#: CCD descriptor preference: isomeric SMILES first.
SMILES_PREFERENCE: Tuple[Tuple[str, str], ...] = (
    ("SMILES_CANONICAL", "CACTVS"),
    ("SMILES_CANONICAL", "OpenEye OEToolkits"),
    ("SMILES", "CACTVS"),
    ("SMILES", "OpenEye OEToolkits"),
)
PROTONATION_STATES = ("primary", "deprotonated")
#: A start pose closer than this to the crystal pose is re-embedded with the next seed.
MIN_START_RMSD = 2.0
MAX_EMBED_ATTEMPTS = 20


class LigandError(ValueError):
    """The ligand cannot be built or checked; the copy is excluded with this reason."""


# ------------------------------------------------------------------ CCD
def ccd_smiles(descriptors: Sequence[Mapping[str, str]]) -> Tuple[str, str]:
    """Pick the preferred SMILES from CCD ``_pdbx_chem_comp_descriptor`` rows.

    Each row has ``type``, ``program`` and ``descriptor``. Returns
    ``(smiles, "type/program")``.
    """
    for kind, program in SMILES_PREFERENCE:
        for row in descriptors:
            if row.get("type") == kind and str(row.get("program", "")).startswith(program):
                return str(row["descriptor"]).strip().strip('"'), f"{kind}/{program}"
    raise LigandError("no SMILES in the CCD entry")


def parse_ccd_cif(text: str) -> Dict[str, object]:
    """Descriptors and formal charge from a CCD mmCIF (``files.rcsb.org/ligands/download/X.cif``)."""
    import gemmi

    doc = gemmi.cif.read_string(text)
    block = doc.sole_block()
    table = block.find("_pdbx_chem_comp_descriptor.", ["type", "program", "descriptor"])
    rows = [{"type": gemmi.cif.as_string(r[0]), "program": gemmi.cif.as_string(r[1]),
             "descriptor": gemmi.cif.as_string(r[2])} for r in table]
    comp_id = gemmi.cif.as_string(block.find_value("_chem_comp.id") or block.name)
    return {"comp_id": comp_id, "descriptors": rows}


# ------------------------------------------------------------ chemistry
def parent_molecule(smiles: str):
    """Heavy-atom molecule from SMILES, with stereocentres checked to be fully specified."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise LigandError(f"RDKit cannot parse the CCD SMILES {smiles!r}")
    centres = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    unassigned = [idx for idx, label in centres if label == "?"]
    if unassigned:
        raise LigandError(f"CCD SMILES leaves {len(unassigned)} stereocentre(s) unassigned")
    return mol


def protonate(parent, state: str):
    """The ligand in one of :data:`PROTONATION_STATES`, with explicit hydrogens (no coordinates)."""
    from rdkit import Chem

    if state not in PROTONATION_STATES:
        raise ValueError(f"unknown protonation state {state!r}")
    rw = Chem.RWMol(Chem.RemoveHs(parent))
    terminal = terminal_phosphate_oxygens(rw)
    single_bonded: Dict[int, List[int]] = {}
    for p_idx, oxygens in terminal.items():
        singles = [o for o in oxygens
                   if rw.GetBondBetweenAtoms(p_idx, o).GetBondType() == Chem.BondType.SINGLE]
        single_bonded[p_idx] = singles
        for o in singles:
            atom = rw.GetAtomWithIdx(o)
            atom.SetFormalCharge(-1)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    if state == "primary":
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=True))
        phosphorus = sorted((p for p in single_bonded if single_bonded[p]), key=lambda p: ranks[p])
        n_protons = len(terminal) // 2
        rw = Chem.RWMol(mol)
        for p_idx in phosphorus[:n_protons]:
            o = min(single_bonded[p_idx], key=lambda i: ranks[i])
            atom = rw.GetAtomWithIdx(o)
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(1)
        mol = rw.GetMol()
        Chem.SanitizeMol(mol)
    return Chem.AddHs(mol)


def net_charge(mol) -> int:
    from rdkit import Chem

    return int(Chem.GetFormalCharge(mol))


# ---------------------------------------------------------- crystal copy
def crystal_molecule(elements: Sequence[str], coords: np.ndarray, atom_names: Sequence[str], template):
    """The crystal copy as an RDKit molecule with CCD bond orders and 3-D stereo.

    Connectivity comes from interatomic distances; bond orders from the CCD
    template (``AssignBondOrdersFromTemplate``), so graph automorphisms and
    stereochemistry are those of the real compound.

    Returns ``(mol, complete)``: ``complete`` is False when the copy lacks heavy
    atoms of the template, in which case bond orders are not assigned and the
    configuration cannot be checked.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDetermineBonds

    rw = Chem.RWMol()
    conf = Chem.Conformer(len(elements))
    for i, (el, xyz) in enumerate(zip(elements, coords)):
        symbol = el.strip().capitalize()
        atom = Chem.Atom(symbol)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
        conf.SetAtomPosition(i, [float(v) for v in xyz])
    rw.AddConformer(conf, assignId=True)
    mol = rw.GetMol()
    rdDetermineBonds.DetermineConnectivity(mol, useHueckel=False)
    heavy_template = Chem.RemoveHs(template)
    complete = mol.GetNumAtoms() == heavy_template.GetNumAtoms()
    if not complete:
        return mol, False
    try:
        generic_template = Chem.Mol(heavy_template)
        for atom in generic_template.GetAtoms():
            atom.SetFormalCharge(0)
            atom.SetNoImplicit(False)
            atom.SetNumExplicitHs(0)
        Chem.SanitizeMol(generic_template)
        assigned = AllChem.AssignBondOrdersFromTemplate(generic_template, mol)
    except (ValueError, RuntimeError) as exc:
        raise LigandError(f"crystal copy does not match the CCD template: {exc}") from exc
    Chem.AssignStereochemistryFrom3D(assigned)
    return assigned, True


def _neutral_isomeric_smiles(mol) -> str:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    heavy = Chem.RemoveHs(mol)
    neutral = rdMolStandardize.Uncharger().uncharge(heavy)
    return Chem.MolToSmiles(neutral, isomericSmiles=True)


def configuration_matches(crystal, template) -> bool:
    """True when the crystal copy's configuration (myo, scyllo, ...) equals the template's."""
    return _neutral_isomeric_smiles(crystal) == _neutral_isomeric_smiles(template)


# -------------------------------------------------------------- embedding
@dataclass
class StartPose:
    mol: object  # RDKit molecule with explicit H and one conformer
    seed: int
    start_rmsd: float
    attempts: int


def _random_rotation(seed: int) -> np.ndarray:
    from scipy.spatial.transform import Rotation

    return Rotation.random(random_state=seed).as_matrix()


def start_pose(ligand, centre: Sequence[float], seed: int, crystal=None,
               min_start_rmsd: float = MIN_START_RMSD) -> StartPose:
    """ETKDGv3 embedding with a fixed seed, a random rigid rotation, centroid on ``centre``.

    Re-embeds with the next seed while the start pose is within ``min_start_rmsd``
    of the crystal pose, so the input never encodes the answer.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem

    for attempt in range(MAX_EMBED_ATTEMPTS):
        s = seed + attempt
        mol = Chem.Mol(ligand)
        params = AllChem.ETKDGv3()
        params.randomSeed = s
        if AllChem.EmbedMolecule(mol, params) != 0:
            params.useRandomCoords = True
            if AllChem.EmbedMolecule(mol, params) != 0:
                continue
        xyz = mol.GetConformer().GetPositions()
        xyz = (xyz - xyz.mean(axis=0)) @ _random_rotation(s).T + np.asarray(centre, dtype=float)
        conf = mol.GetConformer()
        for i, p in enumerate(xyz):
            conf.SetAtomPosition(i, p.tolist())
        rmsd = symmetric_rmsd(crystal, mol) if crystal is not None else math.inf
        if rmsd > min_start_rmsd:
            return StartPose(mol, s, float(rmsd), attempt + 1)
    raise LigandError(f"no start pose more than {min_start_rmsd} Å from the crystal pose "
                      f"after {MAX_EMBED_ATTEMPTS} seeds")


def pose_on_crystal(ligand, crystal):
    """The docked ligand's atoms placed on the crystal coordinates (hydrogens rebuilt).

    Used for the scoring-versus-sampling control: the crystal pose, scored and
    locally minimised in the prepared receptor.
    """
    from rdkit import Chem

    from .rmsd import _element_graph, atom_mappings, heavy_coordinates

    heavy = Chem.RemoveHs(ligand)
    lig_graph, _ = _element_graph(heavy)
    cry_graph, cry_xyz = heavy_coordinates(crystal)
    mappings = atom_mappings(cry_graph, lig_graph)
    if not mappings:
        raise LigandError("crystal copy does not map onto the ligand")
    core_map, cry_terminal, lig_terminal = mappings[0]
    full = dict(core_map)
    for p_cry, oxygens in cry_terminal.items():
        partners = lig_terminal[core_map[p_cry]]
        if len(partners) != len(oxygens):
            raise LigandError("terminal oxygen counts differ")
        full.update(dict(zip(oxygens, partners)))
    if len(full) != heavy.GetNumAtoms():
        raise LigandError("crystal copy is incomplete")
    conf = Chem.Conformer(heavy.GetNumAtoms())
    for c_idx, l_idx in full.items():
        conf.SetAtomPosition(l_idx, cry_xyz[c_idx].tolist())
    placed = Chem.Mol(heavy)
    placed.RemoveAllConformers()
    placed.AddConformer(conf, assignId=True)
    return Chem.AddHs(placed, addCoords=True)


# ------------------------------------------------------------------ PDBQT
def to_pdbqt(mol) -> str:
    """Meeko ligand preparation: rigid rings, rotatable hydroxyls, Gasteiger charges."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    preparation = MoleculePreparation(rigid_macrocycles=True)
    setups = preparation.prepare(mol)
    text, ok, error = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise LigandError(f"Meeko could not write the ligand: {error}")
    return text


def poses_from_pdbqt(text: str) -> List[object]:
    """RDKit molecules (with hydrogens) for every pose in a Vina output PDBQT."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem

    pdbqt = PDBQTMolecule(text, is_dlg=False, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pdbqt)
    mol = mols[0]
    out = []
    for conf in mol.GetConformers():
        single = Chem.Mol(mol)
        single.RemoveAllConformers()
        single.AddConformer(Chem.Conformer(conf), assignId=True)
        out.append(single)
    return out
