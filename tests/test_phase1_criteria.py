"""Tests for the Phase 1 criteria report (parts that do not need fpocket)."""

from __future__ import annotations

import numpy as np
import pytest
from Bio.PDB import PDBIO, PDBParser

from cryptic_ip.analysis.inositol_detection import detect_inositol_residues
from cryptic_ip.analysis.structure_arrays import load_structure_arrays
from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure
from cryptic_ip.validation.phase1_criteria import (
    CONTROLS,
    Control,
    binding_region_superposition,
    contact_residues,
    critical_test,
    judge,
    kabsch,
    select_ligand_copy,
    write_apo_chain,
    write_holo_site,
)


@pytest.fixture(scope="module")
def holo(tmp_path_factory):
    root = tmp_path_factory.mktemp("phase1")
    spec = default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=0)[0]
    return write_synthetic_structure(spec, root)


def _random_rotation(seed: int) -> np.ndarray:
    q, r = np.linalg.qr(np.random.default_rng(seed).normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] = -q[:, 0]
    return q


def test_kabsch_recovers_a_rigid_motion():
    rng = np.random.default_rng(1)
    mobile = rng.normal(size=(40, 3)) * 10
    rotation = _random_rotation(2)
    target = mobile @ rotation + np.array([5.0, -3.0, 12.0])
    fitted, mc, tc, rmsd = kabsch(mobile, target)
    assert rmsd == pytest.approx(0.0, abs=1e-6)
    assert np.allclose((mobile - mc) @ fitted + tc, target, atol=1e-6)


def test_ligand_copy_and_chain(holo):
    ligand = select_ligand_copy(holo)
    assert ligand is not None
    assert ligand["chain"] == "A"
    assert ligand["n_contacts"] > 0
    assert ligand["coords"].shape[1] == 3


def test_apo_chain_has_no_ligand_and_holo_site_keeps_only_the_ligand(holo, tmp_path):
    ligand = select_ligand_copy(holo)
    apo = load_structure_arrays(write_apo_chain(holo, "A", tmp_path / "apo.pdb"))
    assert detect_inositol_residues(apo) == []
    assert apo.is_hetero.sum() == 0
    site = load_structure_arrays(write_holo_site(holo, "A", ligand["residue_key"], tmp_path / "site.pdb"))
    assert len(detect_inositol_residues(site)) == 1
    assert not np.isin(site.resnames, ["HOH", "WAT"]).any()


def test_contact_residues_are_the_basic_lining(holo):
    ligand = select_ligand_copy(holo)
    contacts = contact_residues(holo, ligand["coords"], "A", 4.0)
    names = {name for _, _, name in contacts}
    assert names & {"ARG", "LYS"}


def test_superposition_pairs_by_sequence_not_number(holo, tmp_path):
    """A renumbered, moved and N-terminally truncated copy superposes exactly."""
    ligand = select_ligand_copy(holo)
    structure = PDBParser(QUIET=True).get_structure("m", str(holo))
    rotation, shift = _random_rotation(7), np.array([30.0, -12.0, 4.0])
    chain = structure[0]["A"]
    for residue in list(chain):
        if residue.id[0] != " ":
            chain.detach_child(residue.id)
        elif residue.id[1] <= 3:
            chain.detach_child(residue.id)  # truncate the N-terminus
    for residue in list(chain):
        residue.id = (" ", residue.id[1] + 500, " ")  # renumber
    for atom in structure.get_atoms():
        atom.set_coord(atom.coord @ rotation + shift)
    io = PDBIO()
    io.set_structure(structure)
    model_path = tmp_path / "model.pdb"
    io.save(str(model_path))

    result = binding_region_superposition(holo, "A", model_path, ligand["coords"])
    assert result["ok"]
    assert result["rmsd"] == pytest.approx(0.0, abs=1e-2)
    assert result["sequence_identity"] == pytest.approx(1.0)
    expected = ligand["coords"] @ rotation + shift
    assert np.allclose(result["ligand_in_model_frame"], expected, atol=1e-2)


def _site(rank, n, score, potential=None):
    return {"ok": True, "site_rank": rank, "n_pockets": n, "site_rank_fraction": rank / n,
            "site_score": score, "apbs_potential_kT": potential}


def test_judge_positive_and_unmeasured():
    adar2 = CONTROLS[0]
    crystal = {
        "site": _site(2, 20, 0.81, potential=None),
        "residues": {"basic_coordinating_sasa_holo_mean": 1.2, "n_basic_within_5A": 7},
    }
    criteria = judge(adar2, crystal, None)
    assert criteria["pocket_rank"]["passes"] is True
    assert criteria["site_sasa_holo"]["passes"] is True
    assert criteria["basic_residues"]["passes"] is True
    assert criteria["composite_score"]["passes"] is True
    # No APBS value: unmeasured, never reported as met.
    assert criteria["apbs_potential"]["passes"] is None
    assert "alphafold_vs_crystal_rmsd" not in criteria


def test_judge_negative_rank():
    negative = Control("PH", "XXXX", "negative")
    assert judge(negative, {"site": _site(15, 20, 0.3)}, None)["pocket_rank_bottom_half"]["passes"] is True
    assert judge(negative, {"site": _site(2, 20, 0.6)}, None)["pocket_rank_bottom_half"]["passes"] is False


def test_critical_test():
    def entry(name, role, site):
        return {"control": {"name": name, "role": role}, "crystal": {"site": site}}

    results = [
        entry("ADAR2", "positive", _site(1, 20, 0.85)),
        entry("Pds5B", "positive", _site(4, 30, 0.55)),
        entry("PLCd1_PH", "negative", _site(12, 15, 0.45)),
        entry("Btk_PH_IP4", "negative", _site(9, 10, 0.30)),
    ]
    test = critical_test(results)
    assert test["adar2_site_in_top_3"] is True
    assert test["plc_and_btk_sites_in_bottom_half"] is True
    # Lowest positive (Pds5B 0.55) is above the highest negative (PLCd1 0.45).
    assert test["no_overlap_all_positives_vs_negatives"] is True
    results[1]["crystal"]["site"]["site_score"] = 0.40
    assert critical_test(results)["no_overlap_all_positives_vs_negatives"] is False
    assert critical_test(results)["no_overlap_adar2_vs_negatives"] is True


def test_hull_depth_ignores_internal_cavities():
    """A point at the centre of a hollow ball is deep; one near its skin is not."""
    from cryptic_ip.validation.phase1_criteria import hull_depth

    rng = np.random.default_rng(3)
    directions = rng.normal(size=(4000, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    radii = rng.uniform(12.0, 20.0, size=(4000, 1))  # a shell with an empty core
    shell = directions * radii
    # The hull's flat facets sit just inside the 20 A outer radius.
    assert 18.0 < hull_depth(shell, (0.0, 0.0, 0.0)) <= 20.0
    assert 0.0 < hull_depth(shell, (0.0, 0.0, 18.0)) < 2.0
    assert hull_depth(shell, (0.0, 0.0, 25.0)) < 0
