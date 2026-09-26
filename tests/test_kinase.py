"""scripts/kinase.py (docs/KINASE_PLAN.md) on synthetic structures and pose lists."""

from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def kinase():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("kinase", ROOT / "scripts" / "kinase.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["kinase"] = module
    spec.loader.exec_module(module)
    return module


def _arrays(resnames, elements, coords, polymer=None, chains=None, resseqs=None):
    n = len(resnames)
    return types.SimpleNamespace(
        resnames=np.array(resnames), elements=np.array(elements),
        coords=np.asarray(coords, dtype=float),
        is_polymer=np.array(polymer if polymer is not None else [False] * n),
        chain_ids=np.array(chains or ["A"] * n),
        resseqs=np.array(resseqs or list(range(1, n + 1))),
        icodes=np.array([""] * n),
        atom_names=np.array([f"X{i}" for i in range(n)]))


# ------------------------------------------------------- the cofactor census
def test_a_nearby_nucleotide_is_found(kinase):
    arrays = _arrays(["ATP"] * 3, ["P", "O", "N"], [[0, 0, 0], [1, 0, 0], [2, 0, 0]],
                     resseqs=[10, 10, 10])
    found = kinase.cofactor_atoms(arrays, np.array([[3.0, 0.0, 0.0]]), cutoff=6.0)
    assert len(found) == 1 and found[0]["resname"] == "ATP" and found[0]["n_atoms"] == 3
    assert found[0]["min_distance"] == pytest.approx(1.0)


def test_a_distant_nucleotide_does_not_qualify(kinase):
    arrays = _arrays(["ATP"], ["P"], [[0, 0, 0]], resseqs=[10])
    assert kinase.cofactor_atoms(arrays, np.array([[50.0, 0.0, 0.0]]), cutoff=6.0) == []


def test_only_nucleotide_residues_count(kinase):
    arrays = _arrays(["HOH", "GOL", "IHP"], ["O", "C", "P"], [[0, 0, 0]] * 3, resseqs=[1, 2, 3])
    assert kinase.cofactor_atoms(arrays, np.array([[1.0, 0.0, 0.0]])) == []


def test_a_polymer_atom_is_never_a_cofactor(kinase):
    """A residue named ATP inside the chain is a modelling artefact, not a cosubstrate."""
    arrays = _arrays(["ATP"], ["P"], [[0, 0, 0]], polymer=[True], resseqs=[10])
    assert kinase.cofactor_atoms(arrays, np.array([[1.0, 0.0, 0.0]])) == []


def test_hydrogens_are_not_counted(kinase):
    arrays = _arrays(["ATP", "ATP"], ["P", "H"], [[0, 0, 0], [0.5, 0, 0]], resseqs=[10, 10])
    assert kinase.cofactor_atoms(arrays, np.array([[1.0, 0.0, 0.0]]))[0]["n_atoms"] == 1


def test_two_cofactor_copies_are_separate_residues(kinase):
    arrays = _arrays(["ATP", "ADP"], ["P", "P"], [[0, 0, 0], [1, 0, 0]], resseqs=[10, 11])
    found = kinase.cofactor_atoms(arrays, np.array([[2.0, 0.0, 0.0]]))
    assert sorted(c["resname"] for c in found) == ["ADP", "ATP"]


# ------------------------------------------------------------ the IPK stratum
def test_the_family_is_recognised_from_any_accession_in_the_field(kinase):
    assert kinase.ipk_label("N9UNA8") == "EhIP6KA"
    assert kinase.ipk_label("P12345;O43314") == "PPIP5K2"
    assert kinase.ipk_label("P12345") is None and kinase.ipk_label("") is None


def test_ip6k_itself_is_in_the_stratum(kinase):
    assert {"Q92551", "Q9UHH9", "Q96PC2"} <= set(kinase.IPK_ACCESSIONS)


# --------------------------------------------------------- the holo receptor
def _apo(tmp_path):
    path = tmp_path / "apo.pdbqt"
    path.write_text("ATOM      1  N   ALA A   1      "
                    "  0.000   0.000   0.000  1.00  0.00    -0.100 N \n")
    return path


def test_the_cofactor_is_appended_and_the_apo_lines_are_untouched(kinase, tmp_path):
    apo = _apo(tmp_path)
    before = apo.read_text().splitlines()
    out = tmp_path / "holo.pdbqt"
    residue = {"resname": "ATP", "chain": "A", "resseq": 10, "icode": "",
               "atoms": [{"name": "PA", "element": "P", "xyz": [1.0, 2.0, 3.0]},
                         {"name": "O1A", "element": "O", "xyz": [2.0, 2.0, 3.0]}]}
    info = kinase.holo_receptor(apo, [residue], out)
    lines = out.read_text().splitlines()
    assert lines[:1] == before                      # the apo receptor is copied verbatim
    assert info["appended"] == 2 and len(lines) == 3
    # the type sits in the last field, left-justified in two columns, as the protein atoms are
    assert lines[1].split()[-1] == "P" and lines[2].split()[-1] == "OA"


def test_an_untypable_element_is_skipped_and_named(kinase, tmp_path):
    out = tmp_path / "holo.pdbqt"
    residue = {"resname": "ATP", "chain": "A", "resseq": 1, "icode": "",
               "atoms": [{"name": "U1", "element": "U", "xyz": [0.0, 0.0, 0.0]}]}
    info = kinase.holo_receptor(_apo(tmp_path), [residue], out)
    assert info["appended"] == 0 and info["untyped"] == ["U"]


def test_no_cofactor_leaves_the_receptor_identical(kinase, tmp_path):
    """The pin: with nothing to append, the holo path must reproduce the apo receptor exactly."""
    apo = _apo(tmp_path)
    out = tmp_path / "holo.pdbqt"
    kinase.holo_receptor(apo, [], out)
    assert out.read_text() == apo.read_text()


def test_charges_default_to_zero_without_a_ccd_template(kinase, tmp_path):
    assert kinase.cofactor_charges("ATP", [], None) == {}
    assert kinase.cofactor_charges("ATP", [], tmp_path) == {}       # directory has no ATP.cif


# ------------------------------------------------------------------ decisions
def _census(keys):
    return pd.DataFrame({"copy_key": keys, "pdb_id": [k.split(":")[0] for k in keys],
                         "homology_group_strict": [f"G{i}" for i in range(len(keys))],
                         "burial_class": "surface", "species": "InsP6",
                         "uniprot_ids": "N9UNA8", "primary_set": "True"})


def _records(keys, rmsd_top, seeds=(1, 2, 3)):
    """Pose lists where the first pose is at ``rmsd_top`` and a near-native pose is present."""
    out = []
    for k in keys:
        runs = []
        for s in seeds:
            runs.append({"seed": s, "vina": [-9.0, -8.0, -7.0], "rmsd": [rmsd_top, 1.0, 9.0],
                         "eel": [0.0, -5.0, 0.0], "seconds": 1.0})
        out.append({"copy_key": k, "runs": runs})
    return out


def test_the_cofactor_rescuing_every_copy_is_detected(kinase):
    keys = [f"P{i}:A:1" for i in range(12)]
    r = kinase.build(_census(keys), {k: {"qualifies": True} for k in keys},
                     _records(keys, 8.0), _records(keys, 1.0), None, n_bootstrap=300)
    assert r["K1"]["decision"] == "better with the cofactor"
    assert r["K1"]["estimate"]["point"] == pytest.approx(1.0)


def test_no_change_is_reported_as_no_change(kinase):
    keys = [f"P{i}:A:1" for i in range(12)]
    r = kinase.build(_census(keys), {k: {"qualifies": True} for k in keys},
                     _records(keys, 8.0), _records(keys, 8.0), None, n_bootstrap=300)
    assert r["K1"]["decision"] == "no material change"
    assert r["K1"]["estimate"]["point"] == pytest.approx(0.0)


def test_the_ceiling_is_reported_separately_from_the_top_pose(kinase):
    """A cofactor that changes ranking but not sampling must not be called a sampling win."""
    keys = [f"P{i}:A:1" for i in range(12)]
    r = kinase.build(_census(keys), {k: {"qualifies": True} for k in keys},
                     _records(keys, 8.0), _records(keys, 1.0), None, n_bootstrap=300)
    assert r["K2"]["estimate"]["point"] == pytest.approx(0.0)      # the ceiling was already 1
    assert r["K1"]["estimate"]["point"] > 0.9


def test_no_shared_copy_is_not_evaluable(kinase):
    keys = [f"P{i}:A:1" for i in range(6)]
    r = kinase.build(_census(keys), {}, _records(keys, 8.0), [], None, n_bootstrap=100)
    assert r["K1"]["decision"].startswith("not evaluable")


def test_few_groups_are_not_evaluable(kinase):
    keys = ["P0:A:1", "P1:A:1"]
    r = kinase.build(_census(keys), {k: {"qualifies": True} for k in keys},
                     _records(keys, 8.0), _records(keys, 1.0), None, n_bootstrap=100)
    assert r["K1"]["decision"].startswith("not evaluable")


def test_k4_says_so_when_no_cofactor_carried_charges(kinase):
    keys = [f"P{i}:A:1" for i in range(12)]
    r = kinase.build(_census(keys), {k: {"qualifies": True} for k in keys},
                     _records(keys, 8.0), _records(keys, 8.0), None, n_bootstrap=200)
    assert r["K4_electrostatics"]["copies_with_cofactor_charges"] == 0
    assert "not evidence about electrostatics" in r["K4_electrostatics"]["note"]


# ------------------------------------------------------- descriptive re-analysis
def test_the_family_table_counts_the_ipk_superfamily(kinase):
    copies = pd.DataFrame({"copy_key": ["A:A:1", "B:A:1", "C:A:1"],
                           "uniprot_ids": ["N9UNA8", "O43314", "P99999"],
                           "success": [0.0, 0.0, 1.0], "burial_class": "surface",
                           "failure_kind": ["scoring", "scoring", ""], "species": "InsP6",
                           "homology_group_strict": ["G1", "G2", "G3"]})
    out = kinase.family_table(pd.DataFrame(), copies)
    assert out["copies"] == 2 and out["successes"] == 0
    assert sorted(out["proteins"]) == ["EhIP6KA", "PPIP5K2"]
    assert out["failure_kinds"] == {"scoring": 2}


def test_the_pyrophosphates_are_pooled_into_one_stratum(kinase):
    copies = pd.DataFrame({"species": ["InsP7"] * 4 + ["InsP8"] * 2 + ["InsP6"] * 6,
                           "success": [0.0] * 6 + [1.0] * 6,
                           "homology_group_strict": [f"G{i}" for i in range(12)]})
    out = kinase.species_table(pd.DataFrame(), copies, n_bootstrap=200)
    assert out["PP-IP (InsP7 + InsP8)"]["copies"] == 6
    assert out["PP-IP (InsP7 + InsP8)"]["estimate"]["per_group"]["point"] == pytest.approx(0.0)
    assert out["InsP6"]["estimate"]["per_group"]["point"] == pytest.approx(1.0)


def test_report_writes_its_outputs(kinase, tmp_path):
    keys = [f"P{i}:A:1" for i in range(12)]
    _census(keys).to_csv(tmp_path / "census.csv", index=False)
    cof = {k: {"qualifies": True, "cofactors": [{"resname": "ATP"}]} for k in keys}
    (tmp_path / "cof.json").write_text(json.dumps(cof))
    for name, recs in (("apo", _records(keys, 8.0)), ("holo", _records(keys, 1.0))):
        d = tmp_path / name
        d.mkdir()
        (d / f"{name}.jsonl").write_text("\n".join(json.dumps(r) for r in recs))
    assert kinase.main(["report", "--census", str(tmp_path / "census.csv"),
                        "--cofactors", str(tmp_path / "cof.json"),
                        "--apo-dir", str(tmp_path / "apo"), "--holo-dir", str(tmp_path / "holo"),
                        "--out-dir", str(tmp_path / "out"), "--n-bootstrap", "100"]) == 0
    text = (tmp_path / "out" / "KINASE.md").read_text()
    assert "**K1 (primary):** better with the cofactor" in text
    assert "ATP" in text
