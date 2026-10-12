"""scripts/template_fit.py (docs/TEMPLATE_PLAN.md) on synthetic geometry and records."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def tf():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("template_fit", ROOT / "scripts" / "template_fit.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["template_fit"] = module
    spec.loader.exec_module(module)
    return module


def _rotation(seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q if np.linalg.det(q) > 0 else -q


# ------------------------------------------------------------------ geometry
def test_kabsch_recovers_an_exact_rigid_transform(tf):
    rng = np.random.default_rng(1)
    points = rng.normal(size=(4, 3)) * 5
    rot, shift = _rotation(2), np.array([3.0, -1.0, 7.0])
    moved = points @ rot + shift
    got_rot, got_shift, rmsd = tf.kabsch_batch(points[None], moved[None])
    assert rmsd[0] == pytest.approx(0.0, abs=1e-8)
    assert np.allclose(points @ got_rot[0] + got_shift[0], moved, atol=1e-8)


def test_kabsch_never_returns_a_reflection(tf):
    """A mirrored point set must fit badly, not perfectly: proteins are chiral."""
    points = np.array([[0.0, 0, 0], [3, 0, 0], [0, 4, 0], [0, 0, 5]])
    mirrored = points * np.array([1.0, 1.0, -1.0])
    _, _, rmsd = tf.kabsch_batch(points[None], mirrored[None])
    assert rmsd[0] > 1.0
    rot, _, _ = tf.kabsch_batch(points[None], mirrored[None])
    assert np.linalg.det(rot[0]) > 0


def test_kabsch_is_batched_over_independent_correspondences(tf):
    rng = np.random.default_rng(3)
    mobile = rng.normal(size=(7, 4, 3))
    target = rng.normal(size=(7, 4, 3))
    _, _, batched = tf.kabsch_batch(mobile, target)
    for i in range(7):
        _, _, one = tf.kabsch_batch(mobile[i][None], target[i][None])
        assert batched[i] == pytest.approx(one[0])


# ------------------------------------------------------------------ matching
def test_correspondences_find_the_true_map_of_a_moved_constellation(tf):
    template = np.array([[0.0, 0, 0], [6, 0, 0], [0, 7, 0], [0, 0, 8]])
    query = template @ _rotation(4) + np.array([10.0, 2, -3])
    maps = tf.correspondences(template, query, 4)
    assert ((0, 1, 2, 3), (0, 1, 2, 3)) in maps


def test_correspondences_prune_distance_incompatible_maps(tf):
    """Pruning is exact for a rigid map, so a scrambled query must yield nothing."""
    template = np.array([[0.0, 0, 0], [6, 0, 0], [0, 7, 0], [0, 0, 8]])
    query = np.array([[0.0, 0, 0], [20, 0, 0], [0, 30, 0], [0, 0, 40]])
    assert tf.correspondences(template, query, 4) == []


def test_correspondences_need_enough_points(tf):
    pts = np.zeros((3, 3))
    assert tf.correspondences(pts, pts, 4) == []


# ------------------------------------------------------------------ placement
def _tree(points):
    from scipy.spatial import cKDTree

    return cKDTree(np.asarray(points, dtype=float))


def test_best_fit_places_a_template_on_a_clear_pocket(tf):
    template = np.array([[0.0, 0, 0], [6, 0, 0], [0, 7, 0], [0, 0, 8]])
    ligand = np.array([[2.0, 2, 2], [2.5, 2, 2]])
    rot, shift = _rotation(5), np.array([30.0, 30, 30])
    query = template @ rot + shift
    far = _tree([[-100.0, -100, -100]])
    fit = tf.best_fit(template, ligand, query, far)
    assert fit is not None
    assert fit["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert fit["n_anchors"] == 4


def test_best_fit_rejects_a_placement_that_clashes(tf):
    template = np.array([[0.0, 0, 0], [6, 0, 0], [0, 7, 0], [0, 0, 8]])
    ligand = np.array([[2.0, 2, 2]])
    query = template.copy()
    # A protein atom sitting exactly where the ligand would land.
    assert tf.best_fit(template, ligand, query, _tree([[2.0, 2, 2]])) is None
    assert tf.best_fit(template, ligand, query, _tree([[-50.0, 0, 0]])) is not None


def test_best_fit_returns_the_best_admissible_placement_not_the_first_tried(tf):
    """A near-perfect placement that clashes must be passed over for the next best one."""
    template = np.array([[0.0, 0, 0], [6, 0, 0], [0, 7, 0], [0, 0, 8]])
    ligand = np.array([[1.0, 1.0, 1.0]])
    # Two copies of the site: one exact, one slightly distorted and 60 A away.
    exact = template.copy()
    distorted = template * 1.05 + np.array([60.0, 0, 0])
    query = np.vstack([exact, distorted])
    blocked = tf.best_fit(template, ligand, query, _tree([[1.0, 1.0, 1.0]]))
    assert blocked is not None
    assert blocked["rmsd"] > 0.0           # the exact fit was rejected for clashing
    assert set(blocked["query_anchors"]) == {4, 5, 6, 7}
    clear = tf.best_fit(template, ligand, query, _tree([[-99.0, 0, 0]]))
    assert clear["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert set(clear["query_anchors"]) == {0, 1, 2, 3}


def test_clash_free_uses_the_plan_s_cutoff(tf):
    ligand = np.array([[0.0, 0, 0]])
    assert not tf.clash_free(ligand, _tree([[0.0, 0, tf.CLASH_DISTANCE - 0.1]]))
    assert tf.clash_free(ligand, _tree([[0.0, 0, tf.CLASH_DISTANCE + 0.1]]))


# ------------------------------------------------------------------ anchors
def _arrays(resnames, atom_names, coords, chains=None, resseqs=None, elements=None, polymer=True):
    n = len(resnames)
    resseqs = list(range(1, n + 1)) if resseqs is None else resseqs
    keys = [f"{c}:{r}" for c, r in zip(chains or ["A"] * n, resseqs)]
    order = list(dict.fromkeys(keys))
    return types.SimpleNamespace(
        coords=np.asarray(coords, dtype=float),
        resnames=np.array(resnames),
        atom_names=np.array(atom_names),
        chain_ids=np.array(chains or ["A"] * n),
        resseqs=np.array(resseqs),
        elements=np.array(elements if elements is not None else ["N"] * n),
        is_polymer=np.full(n, polymer),
        residue_index=np.array([order.index(k) for k in keys]),
    )


def test_residue_points_average_a_residue_s_own_nitrogens(tf):
    arrays = _arrays(["ARG", "ARG", "ARG", "LYS", "ALA"],
                     ["NE", "NH1", "NH2", "NZ", "CA"],
                     [[0, 0, 0], [2, 0, 0], [4, 0, 0], [10, 0, 0], [20, 0, 0]],
                     resseqs=[1, 1, 1, 2, 3])
    points, labels = tf.residue_points(arrays, ())
    assert labels == ["A:ARG1", "A:LYS2"]
    assert points[0] == pytest.approx([2.0, 0, 0])


def test_residue_points_honour_a_requested_residue_set(tf):
    arrays = _arrays(["LYS", "LYS"], ["NZ", "NZ"], [[0, 0, 0], [9, 0, 0]], resseqs=[1, 2])
    points, labels = tf.residue_points(arrays, [("A", 2)])
    assert labels == ["A:LYS2"] and points[0][0] == pytest.approx(9.0)


def test_nearest_keeps_only_the_closest_points(tf):
    points = np.array([[0.0, 0, 0], [1, 0, 0], [50, 0, 0]])
    kept, labels = tf.nearest(points, ["a", "b", "c"], np.zeros(3), 2)
    assert labels == ["a", "b"] and len(kept) == 2


def test_template_anchors_use_phosphate_oxygens_not_the_whole_ligand(tf):
    """A basic residue packed against the ring, far from every phosphate, is not an anchor."""
    coords = [[0, 0, 0], [1.5, 0, 0],            # P and its oxygen
              [10, 0, 0],                        # a ring carbon, far from the phosphate
              [3.0, 0, 0],                       # LYS near the phosphate oxygen
              [11.0, 0, 0]]                      # LYS near the ring carbon only
    arrays = _arrays(["IHP", "IHP", "IHP", "LYS", "LYS"],
                     ["P1", "O11", "C1", "NZ", "NZ"], coords,
                     resseqs=[100, 100, 100, 1, 2],
                     elements=["P", "O", "C", "N", "N"])
    points, labels, ligand = tf.template_anchors(arrays, np.array([0, 1, 2]))
    assert labels == ["A:LYS1"]
    assert len(ligand) == 3


# ------------------------------------------------------------------ leakage
def test_eligible_templates_drop_the_query_s_own_protein_and_cluster(tf):
    templates = [{"copy_key": "a", "uniprot_ids": "P11111"},
                 {"copy_key": "b", "uniprot_ids": "P22222,P33333"},
                 {"copy_key": "c", "uniprot_ids": "Q99999"}]
    kept, removed = tf.eligible_templates(templates, "P11111", "P22222")
    assert [t["copy_key"] for t in kept] == ["c"]
    assert removed == {"same_accession": 1, "same_cluster": 1}


def test_eligible_templates_keep_everything_without_a_cluster(tf):
    templates = [{"copy_key": "a", "uniprot_ids": "P11111"}]
    kept, removed = tf.eligible_templates(templates, "Q00000", float("nan"))
    assert len(kept) == 1 and removed == {"same_accession": 0, "same_cluster": 0}


# ------------------------------------------------------------------ decisions
def _records(binder_scores, binder_control_scores, cand_scores=(), control_scores=()):
    out = []
    for i, s in enumerate(binder_scores):
        out.append({"arm": "annotated", "uniprot_id": f"B{i}", "cluster": f"cb{i}", "score": s})
    for i, s in enumerate(binder_control_scores):
        out.append({"arm": "annotated_controls", "uniprot_id": f"C{i}", "cluster": f"cc{i}",
                    "score": s, "matched_to": f"B{i}"})
    for i, s in enumerate(cand_scores):
        out.append({"arm": "candidates", "uniprot_id": f"K{i}", "cluster": f"ck{i}", "score": s})
    for i, s in enumerate(control_scores):
        out.append({"arm": "controls", "uniprot_id": f"M{i}", "cluster": f"cm{i}", "score": s,
                    "matched_to": f"K{i}"})
    return out


def test_guard_passes_when_binders_separate_from_their_controls(tf):
    records = _records([-0.5] * 8, [-3.5] * 8)
    result = tf.guard(records, n_bootstrap=200)
    assert result["verdict"] == "pass"
    assert result["auc"]["roc_auc"]["point"] == pytest.approx(1.0)


def test_guard_fails_when_the_score_cannot_tell_them_apart(tf):
    records = _records([-2.0] * 8, [-2.0] * 8)
    assert tf.guard(records, n_bootstrap=200)["verdict"] == "fail"


def test_guard_is_not_evaluable_below_the_project_s_group_floor(tf):
    records = _records([-0.5] * 2, [-3.5] * 2)
    result = tf.guard(records, n_bootstrap=200)
    assert result["verdict"] == "not evaluable"
    assert str(tf.MIN_GROUPS) in result["reason"]


def test_a_failed_guard_stops_the_study_before_any_candidate_is_ranked(tf):
    records = _records([-2.0] * 8, [-2.0] * 8, cand_scores=[-0.1] * 8, control_scores=[-3.0] * 8)
    result = tf.build(records, n_bootstrap=200)
    assert result["N1"]["verdict"] == "fail"
    assert result["N2"]["decision"] == "not run"
    assert result["N3"]["members"] == []
    assert "holm" not in result


def test_a_passed_guard_lets_n2_and_n3_run(tf):
    records = _records([-0.5] * 8, [-3.5] * 8, cand_scores=[-0.5] * 8, control_scores=[-3.0] * 8)
    result = tf.build(records, n_bootstrap=200)
    assert result["N1"]["verdict"] == "pass"
    assert result["N2"]["decision"] == "better in candidates"
    assert len(result["N3"]["members"]) == 8
    assert set(result["holm"]) == {"N1", "N2"}


def test_n2_calls_a_tiny_difference_no_difference(tf):
    scores = [-2.0, -2.05, -1.95, -2.02, -1.98, -2.01, -1.99, -2.0]
    records = _records([-0.5] * 8, [-3.5] * 8, cand_scores=scores, control_scores=scores)
    result = tf.build(records, n_bootstrap=400)
    assert result["N2"]["decision"] == "no difference"


def test_shortlist_uses_the_control_arm_s_percentile(tf):
    records = _records([], [], cand_scores=[-0.2, -3.9], control_scores=[-3.0] * 20)
    picked = tf.shortlist(records)
    assert [m["uniprot_id"] for m in picked["members"]] == ["K0"]
    assert picked["threshold"] == pytest.approx(-3.0)


def test_markdown_reports_a_failed_guard_without_a_candidate_table(tf):
    text = tf.markdown(tf.build(_records([-2.0] * 8, [-2.0] * 8), n_bootstrap=200))
    assert "Verdict: **fail**" in text
    assert "| accession |" not in text


def test_markdown_never_calls_a_shortlisted_protein_a_predicted_binder(tf):
    records = _records([-0.5] * 8, [-3.5] * 8, cand_scores=[-0.5] * 8, control_scores=[-3.0] * 8)
    text = tf.markdown(tf.build(records, n_bootstrap=200))
    assert "template-compatible" in text
    assert "predicted binder" not in text.replace("not predicted binders", "")


# ------------------------------------------------------------------ plumbing
def test_shard_partitions_without_loss_or_overlap(tf):
    items = list(range(23))
    parts = [tf.shard(items, i, 4) for i in range(4)]
    assert sorted(x for p in parts for x in p) == items


def test_parse_positions_reads_the_screen_s_residue_lists(tf):
    assert tf.parse_positions("12,15, 19") == [12, 15, 19]
    assert tf.parse_positions(float("nan")) == []
    assert tf.parse_positions(None) == []


def test_json_round_trips_through_gzip(tf, tmp_path):
    path = tmp_path / "x.json.gz"
    tf.write_json(path, {"a": np.float64(1.5), "b": np.array([1, 2])})
    assert tf.read_json(path) == {"a": 1.5, "b": [1, 2]}


# ------------------------------------------------------------------ end to end
def _pdb(residues, ligand=None) -> str:
    """A minimal PDB: ``residues`` as (resname, resseq, {atom: xyz}), plus one IHP copy."""
    lines, serial = [], 1
    for resname, resseq, atoms in residues:
        for name, (x, y, z) in atoms.items():
            lines.append(f"ATOM  {serial:5d} {name:<4s} {resname:>3s} A{resseq:4d}    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00          {name[0]:>2s}")
            serial += 1
    if ligand:
        for name, (x, y, z) in ligand.items():
            lines.append(f"HETATM{serial:5d} {name:<4s} IHP A 900    "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          {name[0]:>2s}")
            serial += 1
    return "\n".join(lines) + "\nEND\n"


def _backbone(resseq, centre):
    """Backbone atoms placed well away from the anchor, so they cannot cause a clash."""
    x, y, z = centre
    return {"N": (x, y, z + 12.0), "CA": (x, y, z + 13.0), "C": (x, y, z + 14.0), "O": (x, y, z + 15.0)}


def _site(anchors):
    residues = []
    for i, xyz in enumerate(anchors, start=1):
        atoms = dict(_backbone(i, xyz))
        atoms["NZ"] = tuple(xyz)
        residues.append(("LYS", i, atoms))
    return residues


ANCHORS = np.array([[0.0, 0, 0], [7.0, 0, 0], [0.0, 8.0, 0], [0.0, 0, 9.0]])


def _ligand(anchors=ANCHORS):
    """One phosphate per anchor: an oxygen 3 A inside it, its phosphorus 1.4 A further in.

    The distances matter - an anchor counts only within 4 A of a *phosphate* oxygen, and
    an oxygen counts as a phosphate oxygen only within 1.8 A of a phosphorus.
    """
    centre = anchors.mean(axis=0)
    atoms = {}
    for i, a in enumerate(anchors, start=1):
        direction = (centre - a) / np.linalg.norm(centre - a)
        oxygen = a + 3.0 * direction
        atoms[f"O{i}1"] = tuple(oxygen)
        atoms[f"P{i}"] = tuple(oxygen + 1.4 * direction)
    return atoms


def test_templates_and_fit_agree_on_a_synthetic_site(tf, tmp_path):
    """The whole path through the real PDB parser: extract a template, then re-find it."""
    import pandas as pd

    structures = tmp_path / "structures"
    structures.mkdir()
    (structures / "1TST.pdb").write_text(_pdb(_site(ANCHORS), _ligand()))
    census = pd.DataFrame([{"pdb_id": "1TST", "copy_key": "1TST:A:900", "comp_id": "IHP", "chain": "A",
                            "resseq": 900, "icode": "", "status": "eligible", "complete": True,
                            "symmetry_contact": False, "homology_group_strict": "G:1TST",
                            "uniprot_ids": "P00001", "burial_class": "cryptic"}])
    census_path = tmp_path / "census.csv"
    census.to_csv(census_path, index=False)
    out = tmp_path / "templates.json"
    assert tf.main(["templates", "--census", str(census_path), "--structures", str(structures),
                    "--out", str(out)]) == 0
    library = tf.read_json(out)
    assert library["n_kept"] == 1, library["dropped"]
    template = library["templates"][0]
    assert len(template["anchors"]) == 4 and len(template["ligand"]) == 8

    # The same site, rigidly moved, is a query: the transplant must land back on it.
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    moved = ANCHORS @ _rotation(9) + np.array([40.0, -20.0, 5.0])
    query_path = tmp_path / "AF-Q00001-F1-model_v4.pdb"
    query_path.write_text(_pdb(_site(moved)))
    arrays = load_structure_arrays(query_path)
    fit = tf.score_query(library["templates"], arrays, [1, 2, 3, 4])
    assert fit["censored"] is False
    assert fit["rmsd"] == pytest.approx(0.0, abs=1e-3)
    assert fit["template"] == "1TST:A:900"


def test_a_pocket_whose_basic_residues_are_scrambled_is_censored(tf, tmp_path):
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    template = {"copy_key": "1TST:A:900", "uniprot_ids": "P00001",
                "anchors": ANCHORS, "ligand": np.array([[2.0, 2.0, 2.0]])}
    scrambled = np.array([[0.0, 0, 0], [20.0, 0, 0], [0.0, 25.0, 0], [0.0, 0, 30.0]])
    path = tmp_path / "AF-Q00002-F1-model_v4.pdb"
    path.write_text(_pdb(_site(scrambled)))
    arrays = load_structure_arrays(path)
    fit = tf.score_query([template], arrays, [1, 2, 3, 4])
    assert fit["censored"] is True
    assert fit["score"] == pytest.approx(-tf.RMSD_CEILING)
