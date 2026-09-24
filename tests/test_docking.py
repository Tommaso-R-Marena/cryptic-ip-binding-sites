"""The docking core (docs/REDOCKING_PLAN.md): RMSD, ligand states, configuration, receptor typing."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

Chem = pytest.importorskip("rdkit.Chem")
AllChem = pytest.importorskip("rdkit.Chem.AllChem")

from cryptic_ip.docking import rmsd as R  # noqa: E402
from cryptic_ip.docking import ligand as L  # noqa: E402

MYO_IP6 = ("O[P](=O)(O)O[C@H]1[C@H](O[P](=O)(O)O)[C@@H](O[P](=O)(O)O)[C@H](O[P](=O)(O)O)"
           "[C@@H](O[P](=O)(O)O)[C@@H]1O[P](=O)(O)O")
MYO_IP3 = "O[C@H]1[C@@H](O)[C@H](OP(=O)(O)O)[C@@H](OP(=O)(O)O)[C@H](O)[C@H]1OP(=O)(O)O"
ATP = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O"


def embedded(smiles: str, seed: int = 3):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return Chem.RemoveHs(mol)


def translated(mol, shift):
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, (p.x + shift[0], p.y + shift[1], p.z + shift[2]))
    return out


def brute_force_rmsd(ref, probe) -> float:
    """Minimum over every element-graph automorphism, enumerated by RDKit (small molecules only)."""
    g_ref, x_ref = R.heavy_coordinates(ref)
    g_probe, x_probe = R.heavy_coordinates(probe)
    best = np.inf
    for match in g_probe.GetSubstructMatches(g_ref, uniquify=False, useChirality=False, maxMatches=10_000_000):
        best = min(best, float(np.sqrt(np.mean(np.sum((x_ref - x_probe[list(match)]) ** 2, axis=1)))))
    return best


# ------------------------------------------------------------------ RMSD
def test_relabelled_identical_pose_is_zero_while_naive_rmsd_is_not():
    ref = embedded(MYO_IP6)
    order = list(np.random.default_rng(0).permutation(ref.GetNumAtoms()))
    relabelled = Chem.RenumberAtoms(ref, [int(i) for i in order])
    assert R.symmetric_rmsd(ref, relabelled) == pytest.approx(0.0, abs=1e-6)
    assert R.naive_rmsd(ref, relabelled) > 1.0


def test_swapping_terminal_phosphate_oxygens_is_zero():
    ref = embedded(MYO_IP6)
    graph, _ = R.heavy_coordinates(ref)
    p_idx, oxygens = next(iter(R.terminal_phosphate_oxygens(graph).items()))
    a, b = oxygens[:2]
    swapped = Chem.Mol(ref)
    conf = swapped.GetConformer()
    pa, pb = conf.GetAtomPosition(a), conf.GetAtomPosition(b)
    conf.SetAtomPosition(a, pb)
    conf.SetAtomPosition(b, pa)
    assert R.naive_rmsd(ref, swapped) > 0.5
    assert R.symmetric_rmsd(ref, swapped) == pytest.approx(0.0, abs=1e-6)


def test_rigid_translation_gives_its_length():
    ref = embedded(MYO_IP6)
    assert R.symmetric_rmsd(ref, translated(ref, (3.0, 0.0, 0.0))) == pytest.approx(3.0, abs=1e-6)
    assert R.symmetric_rmsd(ref, translated(ref, (1.0, 2.0, 2.0))) == pytest.approx(3.0, abs=1e-6)


def test_no_superposition_is_done():
    """A rotated pose about its own centroid is not aligned back: RMSD stays large."""
    from scipy.spatial.transform import Rotation

    ref = embedded(MYO_IP6)
    xyz = ref.GetConformer().GetPositions()
    c = xyz.mean(axis=0)
    rot = Rotation.from_euler("z", 90, degrees=True).as_matrix()
    moved = Chem.Mol(ref)
    for i, p in enumerate((xyz - c) @ rot.T + c):
        moved.GetConformer().SetAtomPosition(i, p.tolist())
    assert R.symmetric_rmsd(ref, moved) > 1.0


@pytest.mark.parametrize("smiles", [MYO_IP3, "CC(C)OP(=O)(O)O", "OC1C(O)C(OP(=O)(O)O)C(O)C(O)C1O"])
def test_matches_brute_force_on_a_different_conformer(smiles):
    ref = embedded(smiles, seed=3)
    probe = translated(embedded(smiles, seed=11), (0.7, -0.4, 0.2))
    assert R.symmetric_rmsd(ref, probe) == pytest.approx(brute_force_rmsd(ref, probe), abs=1e-6)


def test_ip6_different_conformer_agrees_with_explicit_enumeration():
    """IP6 has 12 core automorphisms x 6^6 oxygen permutations; check against a direct enumeration."""
    ref = embedded(MYO_IP6, seed=3)
    probe = embedded(MYO_IP6, seed=21)
    g_ref, x_ref = R.heavy_coordinates(ref)
    g_probe, x_probe = R.heavy_coordinates(probe)
    best = np.inf
    for core_map, ref_t, probe_t in R.atom_mappings(g_ref, g_probe):
        total = sum(np.sum((x_ref[i] - x_probe[j]) ** 2) for i, j in core_map.items())
        for p, oxygens in ref_t.items():
            cands = probe_t[core_map[p]]
            total += min(sum(np.sum((x_ref[o] - x_probe[c]) ** 2) for o, c in zip(oxygens, perm))
                         for perm in itertools.permutations(cands))
        best = min(best, total)
    expected = float(np.sqrt(best / ref.GetNumAtoms()))
    assert R.symmetric_rmsd(ref, probe) == pytest.approx(expected, abs=1e-6)
    assert len(R.atom_mappings(g_ref, g_probe)) == 12


def test_incomplete_reference_is_matched_as_a_substructure():
    probe = embedded(MYO_IP6)
    graph, _ = R.heavy_coordinates(probe)
    oxygen = next(iter(R.terminal_phosphate_oxygens(graph).values()))[0]
    rw = Chem.RWMol(probe)
    rw.RemoveAtom(oxygen)
    partial = rw.GetMol()
    assert R.symmetric_rmsd(partial, probe) == pytest.approx(0.0, abs=1e-6)


def test_phosphorus_only_rmsd():
    ref = embedded(MYO_IP6)
    assert R.symmetric_rmsd(ref, translated(ref, (0, 0, 2.0)), phosphorus_only=True) == pytest.approx(2.0)


def test_different_compounds_are_refused():
    with pytest.raises(ValueError):
        R.symmetric_rmsd(embedded(MYO_IP6), embedded(ATP))


# ------------------------------------------------------------ protonation
@pytest.mark.parametrize("smiles,primary,deprotonated", [(MYO_IP6, -9, -12), (MYO_IP3, -5, -6), (ATP, -3, -4)])
def test_protonation_states(smiles, primary, deprotonated):
    parent = L.parent_molecule(smiles)
    assert L.net_charge(L.protonate(parent, "primary")) == primary
    assert L.net_charge(L.protonate(parent, "deprotonated")) == deprotonated


def test_protonation_is_deterministic():
    parent = L.parent_molecule(MYO_IP6)
    a = Chem.MolToSmiles(L.protonate(parent, "primary"))
    b = Chem.MolToSmiles(L.protonate(L.parent_molecule(MYO_IP6), "primary"))
    assert a == b


def test_unassigned_stereocentre_is_refused():
    with pytest.raises(L.LigandError):
        L.parent_molecule("OC1C(O)C(OP(=O)(O)O)C(O)C(O)C1O")


# ---------------------------------------------------------- configuration
def _crystal_from(smiles_for_coords: str, template_smiles: str, seed: int = 5):
    mol = embedded(smiles_for_coords, seed)
    elements = [a.GetSymbol() for a in mol.GetAtoms()]
    names = [f"{e}{i}" for i, e in enumerate(elements)]
    template = L.parent_molecule(template_smiles)
    return L.crystal_molecule(elements, mol.GetConformer().GetPositions(), names, template), template


def test_configuration_matches_itself():
    (crystal, complete), template = _crystal_from(MYO_IP6, MYO_IP6)
    assert complete
    assert L.configuration_matches(crystal, template)


def test_epimer_is_detected():
    epimer = MYO_IP6.replace("[C@H]1", "[C@@H]1", 1)  # invert one ring carbon
    (crystal, complete), template = _crystal_from(epimer, MYO_IP6)
    assert complete
    assert not L.configuration_matches(crystal, template)


def test_incomplete_copy_is_flagged():
    mol = embedded(MYO_IP6)
    graph, _ = R.heavy_coordinates(mol)
    drop = next(iter(R.terminal_phosphate_oxygens(graph).values()))[0]
    keep = [i for i in range(mol.GetNumAtoms()) if i != drop]
    xyz = mol.GetConformer().GetPositions()[keep]
    elements = [mol.GetAtomWithIdx(i).GetSymbol() for i in keep]
    _, complete = L.crystal_molecule(elements, xyz, [str(i) for i in keep], L.parent_molecule(MYO_IP6))
    assert not complete


def test_ccd_smiles_preference():
    rows = [{"type": "SMILES", "program": "CACTVS", "descriptor": "C"},
            {"type": "SMILES_CANONICAL", "program": "OpenEye OEToolkits", "descriptor": "CC"},
            {"type": "SMILES_CANONICAL", "program": "CACTVS", "descriptor": "CCC"}]
    assert L.ccd_smiles(rows) == ("CCC", "SMILES_CANONICAL/CACTVS")
    with pytest.raises(L.LigandError):
        L.ccd_smiles([{"type": "InChI", "program": "x", "descriptor": "y"}])


def test_start_pose_is_far_from_the_crystal_and_reproducible():
    parent = L.parent_molecule(MYO_IP6)
    ligand = L.protonate(parent, "primary")
    crystal = embedded(MYO_IP6, seed=9)
    centre = crystal.GetConformer().GetPositions().mean(axis=0)
    a = L.start_pose(ligand, centre, 1, crystal=crystal)
    b = L.start_pose(ligand, centre, 1, crystal=crystal)
    assert a.start_rmsd > L.MIN_START_RMSD
    assert a.start_rmsd == pytest.approx(b.start_rmsd)
    xyz = Chem.RemoveHs(a.mol).GetConformer().GetPositions()
    assert np.allclose(xyz.mean(axis=0), centre, atol=1.5)


def test_start_pose_reembeds_when_too_close():
    parent = L.parent_molecule(MYO_IP6)
    ligand = L.protonate(parent, "primary")
    crystal = embedded(MYO_IP6, seed=9)
    centre = crystal.GetConformer().GetPositions().mean(axis=0)
    with pytest.raises(L.LigandError):
        L.start_pose(ligand, centre, 1, crystal=crystal, min_start_rmsd=1e6)


def test_pose_on_crystal_places_every_heavy_atom():
    parent = L.parent_molecule(MYO_IP6)
    ligand = L.protonate(parent, "primary")
    crystal = embedded(MYO_IP6, seed=4)
    placed = L.pose_on_crystal(ligand, crystal)
    assert R.symmetric_rmsd(crystal, placed) == pytest.approx(0.0, abs=1e-6)
    assert placed.GetNumAtoms() > crystal.GetNumAtoms()  # hydrogens rebuilt


# ------------------------------------------------------------- engine
def test_box_size_rule():
    from cryptic_ip.docking.engine import box_size, greedy_clusters

    assert box_size(np.array([[0, 0, 0], [1, 1, 1.0]])) == [22.0, 22.0, 22.0]
    assert box_size(np.array([[0, 0, 0], [10, 1, 1.0]]))[0] == pytest.approx(26.0)
    m = np.array([[0, 1, 5], [1, 0, 5], [5, 5, 0.0]])
    assert greedy_clusters(m, 2.0) == [[0, 1], [2]]


def test_pdbqt_types_reads_the_type_column():
    from cryptic_ip.docking.engine import pdbqt_types
    from cryptic_ip.docking.receptor import Atom, pdbqt_line

    line = pdbqt_line(1, Atom("ATOM", "NZ", "LYS", "A", 12, "", np.array([1.0, 2.0, -3.5]), "N", 0.33), "N")
    assert pdbqt_types(line) == ["N"]
    assert float(line[30:38]) == 1.0 and float(line[46:54]) == -3.5
    assert float(line[70:76]) == pytest.approx(0.33)


# ----------------------------------------------------------- receptor
def _atom(name, resname, element, xyz, charge=0.0, resseq=1):
    from cryptic_ip.docking.receptor import Atom

    return Atom("ATOM", name, resname, "A", resseq, "", np.array(xyz, dtype=float), element, charge)


def test_autodock_typing():
    from cryptic_ip.docking.receptor import autodock_types

    atoms = [
        _atom("CG", "HID", "C", [0, 0, 0]), _atom("ND1", "HID", "N", [1.3, 0, 0]),
        _atom("HD1", "HID", "H", [2.3, 0, 0], 0.3), _atom("NE2", "HID", "N", [0, 1.3, 0]),
        _atom("CB", "LYS", "C", [5, 5, 5], -0.1, 2), _atom("HB2", "LYS", "H", [5.9, 5, 5], 0.05, 2),
        _atom("CZ", "PHE", "C", [9, 9, 9], 0, 3), _atom("SD", "MET", "S", [12, 0, 0], 0, 4),
        _atom("O", "GLY", "O", [15, 0, 0], -0.5, 5),
    ]
    typed = {(a.name, a.resname): t for a, t in autodock_types(atoms)}
    assert typed[("CG", "HID")] == "A"
    assert typed[("ND1", "HID")] == "N"      # carries HD1: donor, not acceptor
    assert typed[("NE2", "HID")] == "NA"     # no hydrogen: acceptor
    assert typed[("HD1", "HID")] == "HD"
    assert ("HB2", "LYS") not in typed       # non-polar hydrogen merged
    assert typed[("CZ", "PHE")] == "A" and typed[("SD", "MET")] == "SA" and typed[("O", "GLY")] == "OA"
    charges = {a.name: a.charge for a, _ in autodock_types(atoms)}
    assert charges["CB"] == pytest.approx(-0.05)  # merged hydrogen charge


def test_nucleotide_typing():
    from cryptic_ip.docking.receptor import autodock_types

    atoms = [_atom("N9", "DA", "N", [0, 0, 0]), _atom("N7", "DA", "N", [3, 0, 0]),
             _atom("N1", "DC", "N", [6, 0, 0]), _atom("N3", "DC", "N", [9, 0, 0]),
             _atom("C8", "DA", "C", [12, 0, 0]), _atom("C1'", "DA", "C", [15, 0, 0])]
    typed = {(a.name, a.resname): t for a, t in autodock_types(atoms)}
    assert typed[("N9", "DA")] == "N" and typed[("N7", "DA")] == "NA"
    assert typed[("N1", "DC")] == "N" and typed[("N3", "DC")] == "NA"
    assert typed[("C8", "DA")] == "A" and typed[("C1'", "DA")] == "C"


def test_unknown_element_is_refused():
    from cryptic_ip.docking.receptor import ReceptorError, autodock_types

    with pytest.raises(ReceptorError):
        autodock_types([_atom("X", "UNK", "XE", [0, 0, 0])])


def test_strip_keeps_every_polymer_chain_and_nothing_else(tmp_path):
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.receptor import nearby_metals, write_polymer_pdb

    lines = [
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00  0.00           C",
        "ATOM      3  N   GLY B   1       5.000   0.000   0.000  1.00  0.00           N",
        "ATOM      4  CA  GLY B   1       6.450   0.000   0.000  1.00  0.00           C",
        "HETATM    5  N   MSE B   2       8.000   0.000   0.000  1.00  0.00           N",
        "HETATM    6 SE   MSE B   2       9.000   0.000   0.000  1.00  0.00          SE",
        "HETATM    7 MG    MG C   1      20.000   0.000   0.000  1.00  0.00          MG",
        "HETATM    8  O   HOH C   2      30.000   0.000   0.000  1.00  0.00           O",
        "HETATM    9  C1  GOL C   3      40.000   0.000   0.000  1.00  0.00           C",
        "END",
    ]
    src = tmp_path / "x.pdb"
    src.write_text("\n".join(lines) + "\n")
    arrays = load_structure_arrays(src)
    report = write_polymer_pdb(arrays, tmp_path / "poly.pdb")
    text = (tmp_path / "poly.pdb").read_text()
    assert " A " in text and " B " in text
    assert "MG" not in text and "GOL" not in text and "HOH" not in text
    assert "MET" in text and " SD " in text and report.selenomethionines == 2
    metals = nearby_metals(arrays, np.array([[21.0, 0.0, 0.0]]))
    assert [m.element for m in metals] == ["MG"]
    assert nearby_metals(arrays, np.array([[25.0, 0.0, 0.0]])) == []


# --------------------------------------------------------------- stats
def test_group_equal_estimand_weights_groups_not_copies():
    from cryptic_ip.docking.stats import mean_estimates, per_group_mean

    y = np.array([1, 1, 1, 1, 0.0])
    groups = np.array(["adar", "adar", "adar", "adar", "other"])
    codes = np.array([0, 0, 0, 0, 1])
    assert per_group_mean(y, codes, np.ones(5)) == pytest.approx(0.5)
    est = mean_estimates(y, groups, n_bootstrap=50)
    assert est["per_copy"]["point"] == pytest.approx(0.8)
    assert est["per_group"]["point"] == pytest.approx(0.5)
    assert est["evidence"] is False  # 2 groups


def test_bootstrap_resamples_groups_not_copies():
    from cryptic_ip.docking.stats import mean_estimates

    # One group of 100 successes and one group of 100 failures: resampling copies
    # would give a tight interval around 0.5; resampling groups cannot.
    y = np.r_[np.ones(100), np.zeros(100)]
    groups = np.r_[["a"] * 100, ["b"] * 100]
    est = mean_estimates(y, groups, n_bootstrap=400)
    assert est["per_copy"]["low"] == 0.0 and est["per_copy"]["high"] == 1.0


def test_difference_and_auc():
    from cryptic_ip.docking.stats import auc_estimate, difference_estimates

    y = np.array([1, 1, 0, 0.0])
    d = difference_estimates(y, [True, True, False, False], [False, False, True, True], ["a", "b", "c", "d"],
                             n_bootstrap=50)
    assert d["per_group"]["point"] == pytest.approx(1.0)
    a = auc_estimate([1, 1, 0, 0], [3, 4, 1, 2], ["a", "b", "a", "b"], n_bootstrap=50)
    assert a["roc_auc"]["point"] == pytest.approx(1.0)


def test_phosphodiester_phosphorus_does_not_break_the_configuration_check():
    (crystal, complete), template = _crystal_from(ATP, ATP)
    assert complete
    assert L.configuration_matches(crystal, template)


def test_multicharacter_chains_are_kept_apart(tmp_path):
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.receptor import write_polymer_pdb

    cif = """data_x
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . GLY A 1 1 ? 0.0 0.0 0.0 1.0 0.0 1 GLY AA N 1
ATOM 2 C CA . GLY A 1 1 ? 1.4 0.0 0.0 1.0 0.0 1 GLY AA CA 1
ATOM 3 N N . GLY B 1 1 ? 5.0 0.0 0.0 1.0 0.0 1 GLY AB N 1
ATOM 4 C CA . GLY B 1 1 ? 6.4 0.0 0.0 1.0 0.0 1 GLY AB CA 1
"""
    src = tmp_path / "x.cif"
    src.write_text(cif)
    write_polymer_pdb(load_structure_arrays(src), tmp_path / "out.pdb")
    chains = {line[21] for line in (tmp_path / "out.pdb").read_text().splitlines() if line.startswith("ATOM")}
    assert len(chains) == 2
