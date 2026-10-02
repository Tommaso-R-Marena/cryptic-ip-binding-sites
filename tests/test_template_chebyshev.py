"""scripts/template_chebyshev.py (docs/TEMPLATE2_PLAN.md): study O, the corrected transplant."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def tc():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(
        "template_chebyshev", ROOT / "scripts" / "template_chebyshev.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["template_chebyshev"] = module
    spec.loader.exec_module(module)
    return module


def _tree(points):
    from scipy.spatial import cKDTree

    return cKDTree(np.asarray(points, dtype=float))


def _rotation(seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    return q if np.linalg.det(q) > 0 else -q


ANCHORS = np.array([[0.0, 0, 0], [7.0, 0, 0], [0.0, 8.0, 0], [0.0, 0, 9.0]])
FAR = [[-200.0, -200.0, -200.0]]


# ------------------------------------------------------- the criterion itself
def test_the_pruning_tolerance_is_twice_eps(tc):
    """The lemma is 2*eps; anything smaller is the unsound setting study N shipped."""
    assert tc.PRUNE_TOLERANCE == pytest.approx(2.0 * tc.EPS)


def test_an_exact_transplant_is_found(tc):
    query = ANCHORS @ _rotation(3) + np.array([25.0, -10.0, 4.0])
    fit = tc.best_chebyshev(ANCHORS, np.array([[2.0, 2.0, 2.0]]), query, _tree(FAR))
    assert fit is not None
    assert fit["deviation"] == pytest.approx(0.0, abs=1e-6)


def test_a_placement_within_eps_survives_pruning_even_at_the_2eps_limit(tc):
    """The lossless property: distances may disagree by up to 2*eps and still be kept.

    Study N pruned at 1.5 A while accepting on RMSD <= 4.0, which discarded admissible
    correspondences. Here the tolerance is tied to the criterion, so a correspondence whose
    anchors are each within eps is never thrown away - even when that pushes a pairwise
    distance almost the full 2*eps out.
    """
    eps = tc.EPS
    query = ANCHORS.copy()
    # Move two anchors by eps in opposite directions along x: each residual is <= eps,
    # while the distance between them changes by nearly 2*eps.
    query[0] = query[0] + np.array([-eps, 0.0, 0.0])
    query[1] = query[1] + np.array([+eps, 0.0, 0.0])
    gap = abs(np.linalg.norm(ANCHORS[0] - ANCHORS[1]) - np.linalg.norm(query[0] - query[1]))
    assert gap == pytest.approx(2 * eps)                       # at the lemma's bound
    assert gap > tf_tolerance()                                # study N would have pruned it
    maps = tc.tf.correspondences(ANCHORS, query, tc.K_ANCHORS, tolerance=tc.PRUNE_TOLERANCE)
    assert ((0, 1, 2, 3), (0, 1, 2, 3)) in maps                # kept here
    fit = tc.best_chebyshev(ANCHORS, np.array([[3.0, 3.0, 3.0]]), query, _tree(FAR))
    assert fit is not None and fit["deviation"] <= eps + 1e-9


def tf_tolerance() -> float:
    spec = importlib.util.spec_from_file_location("template_fit", ROOT / "scripts" / "template_fit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.DISTANCE_TOLERANCE


def test_one_badly_placed_anchor_is_rejected_however_good_the_mean(tc):
    """Study N accepted this on RMSD; Chebyshev is the criterion that refuses it."""
    query = ANCHORS.copy()
    query[3] = query[3] + np.array([0.0, 0.0, 5.0])            # one anchor 5 A out
    _, _, rmsd = tc.tf.kabsch_batch(ANCHORS[None], query[None])
    assert rmsd[0] < 4.0                                       # inside study N's ceiling
    assert tc.best_chebyshev(ANCHORS, np.array([[2.0, 2.0, 2.0]]), query, _tree(FAR)) is None


def test_a_clashing_placement_is_rejected(tc):
    ligand = np.array([[2.0, 2.0, 2.0]])
    assert tc.best_chebyshev(ANCHORS, ligand, ANCHORS.copy(), _tree([[2.0, 2.0, 2.0]])) is None
    assert tc.best_chebyshev(ANCHORS, ligand, ANCHORS.copy(), _tree(FAR)) is not None


# ------------------------------------------------------- control selection
def _pool():
    rows = [{"uniprot_id": "Q1", "arm": "candidates", "matched_to": "Q1", "cluster": "c1"},
            {"uniprot_id": "Q2", "arm": "candidates", "matched_to": "Q2", "cluster": "c2"},
            {"uniprot_id": "B1", "arm": "annotated", "matched_to": "B1", "cluster": "b1"}]
    for acc, target in (("C1", "Q1"), ("C2", "Q1"), ("C3", "Q1"), ("D1", "Q2"), ("E1", "B1")):
        rows.append({"uniprot_id": acc, "arm": "controls" if target.startswith("Q")
                     else "annotated_controls", "matched_to": target, "cluster": acc})
    return pd.DataFrame(rows)


def test_select_controls_keeps_the_closest_basic_count(tc):
    counts = {"Q1": 9, "Q2": 4, "B1": 7, "C1": 3, "C2": 8, "C3": 20, "D1": 5, "E1": 7}
    chosen = tc.select_controls(_pool(), counts)
    picked = chosen[chosen["arm"] == "controls"].set_index("matched_to")["uniprot_id"].to_dict()
    assert picked["Q1"] == "C2"           # 8 is nearest to 9, not 3 and not 20
    assert picked["Q2"] == "D1"
    assert len(chosen[chosen["arm"] == "controls"]) == 2    # one control per candidate


def test_select_controls_breaks_ties_on_accession_not_row_order(tc):
    counts = {"Q1": 10, "Q2": 4, "B1": 7, "C1": 8, "C2": 12, "C3": 8, "D1": 4, "E1": 7}
    pool = _pool()
    first = tc.select_controls(pool, counts)
    shuffled = tc.select_controls(pool.sample(frac=1, random_state=5), counts)
    a = first[first["arm"] == "controls"].set_index("matched_to")["uniprot_id"].to_dict()
    b = shuffled[shuffled["arm"] == "controls"].set_index("matched_to")["uniprot_id"].to_dict()
    assert a == b
    assert a["Q1"] == "C1"                # C1 and C3 both at gap 2; C1 sorts first


def test_select_controls_skips_a_pool_with_no_measured_count(tc):
    counts = {"Q1": 9, "Q2": 4, "B1": 7, "D1": 5, "E1": 7}   # C1..C3 unmeasured
    chosen = tc.select_controls(_pool(), counts)
    picked = chosen[chosen["arm"] == "controls"]["matched_to"].tolist()
    assert "Q1" not in picked and "Q2" in picked


# ------------------------------------------------------- decisions
def _records(binders, binder_controls, cands=(), controls=(), buried=True):
    out = []
    for i, s in enumerate(binders):
        out.append({"arm": "annotated", "uniprot_id": f"B{i}", "cluster": f"cb{i}", "score": s,
                    "censored": False, "buried": True, "cryptic": False, "n_basic": 8})
    for i, s in enumerate(binder_controls):
        out.append({"arm": "annotated_controls", "uniprot_id": f"C{i}", "cluster": f"cc{i}",
                    "score": s, "matched_to": f"B{i}", "censored": False, "buried": True,
                    "cryptic": False, "n_basic": 8})
    for i, s in enumerate(cands):
        out.append({"arm": "candidates", "uniprot_id": f"K{i}", "cluster": f"ck{i}", "score": s,
                    "deviation": -s, "censored": False, "buried": buried, "cryptic": False,
                    "n_basic": 7, "relative_sasa": 0.2, "template": "1ABC:A:1"})
    for i, s in enumerate(controls):
        out.append({"arm": "controls", "uniprot_id": f"M{i}", "cluster": f"cm{i}", "score": s,
                    "deviation": -s, "matched_to": f"K{i}", "censored": False, "buried": False,
                    "cryptic": False, "n_basic": 7, "relative_sasa": 0.6, "template": "1ABC:A:1"})
    return out


def test_a_failed_guard_stops_the_study_before_any_candidate_is_ranked(tc):
    r = tc.build(_records([-1.0] * 8, [-1.0] * 8, cands=[-0.1] * 8, controls=[-2.0] * 8),
                 n_bootstrap=200)
    assert r["O1"]["verdict"] == "fail"
    assert r["O2"]["decision"] == "not run"
    assert r["O3"]["members"] == []
    assert "holm" not in r


def test_a_passed_guard_lets_o2_and_o3_run(tc):
    r = tc.build(_records([-0.3] * 8, [-2.4] * 8, cands=[-0.3] * 8, controls=[-2.0] * 8),
                 n_bootstrap=200)
    assert r["O1"]["verdict"] == "pass"
    assert r["O2"]["decision"] == "better in candidates"
    assert len(r["O3"]["members"]) == 8
    assert set(r["holm"]) == {"O1", "O2"}


def test_o3_uses_absolute_criteria_not_a_control_percentile(tc):
    """N3 collapsed because its cut was the control arm's 95th percentile. O3 cannot."""
    r = tc.build(_records([-0.3] * 8, [-2.4] * 8, cands=[-0.3] * 8, controls=[-0.3] * 8),
                 n_bootstrap=200)
    o3 = r["O3"]
    # Controls score identically to candidates here, yet only the buried arm passes, and
    # the count of passing controls is reported rather than used as a threshold.
    assert o3["n_candidates_passing"] == 8
    assert o3["n_controls_passing"] == 0
    assert o3["buried_max_relative_sasa"] == tc.BURIED_MAX_REL_SASA


def test_an_unburied_candidate_never_reaches_the_shortlist(tc):
    r = tc.build(_records([-0.3] * 8, [-2.4] * 8, cands=[-0.1] * 8, controls=[-2.0] * 8,
                          buried=False), n_bootstrap=200)
    assert r["O1"]["verdict"] == "pass"
    assert r["O3"]["members"] == []
    assert r["O3"]["n_candidates_passing"] == 0


def test_the_report_states_the_basic_count_balance(tc):
    r = tc.build(_records([-0.3] * 8, [-2.4] * 8, cands=[-0.3] * 8, controls=[-2.0] * 8),
                 n_bootstrap=200)
    assert r["balance"]["candidates"]["query_mean_basic"] == pytest.approx(7.0)
    text = tc.markdown(r)
    assert "Basic-residue balance after matching" in text
    assert "template-compatible" in text and "predicted binder" not in text.replace(
        "not predicted binders", "")


def test_markdown_reports_a_failed_guard_without_a_candidate_table(tc):
    text = tc.markdown(tc.build(_records([-1.0] * 8, [-1.0] * 8), n_bootstrap=200))
    assert "Verdict: **fail**" in text and "| accession |" not in text


def test_guard_is_not_evaluable_below_the_group_floor(tc):
    r = tc.guard(_records([-0.3] * 2, [-2.4] * 2), n_bootstrap=200)
    assert r["verdict"] == "not evaluable"


# ------------------------------------------------------- structures
def _pdb(residues, ligand=None) -> str:
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


def _site(anchors):
    residues = []
    for i, xyz in enumerate(anchors, start=1):
        x, y, z = xyz
        atoms = {"N": (x, y, z + 12.0), "CA": (x, y, z + 13.0), "C": (x, y, z + 14.0),
                 "O": (x, y, z + 15.0), "NZ": tuple(xyz)}
        residues.append(("LYS", i, atoms))
    return residues


def test_basic_count_counts_the_pockets_k_r_h(tc, tmp_path):
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    path = tmp_path / "AF-Q0-F1-model_v4.pdb"
    path.write_text(_pdb(_site(ANCHORS)))
    arrays = load_structure_arrays(path)
    assert tc.basic_count(arrays, [1, 2, 3, 4]) == 4
    assert tc.basic_count(arrays, [1, 2]) == 2


def test_ligand_identity_matches_the_coordinate_order(tc, tmp_path):
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import find_ligand_instances

    ligand = {"P1": (2.0, 2.0, 2.0), "O11": (3.2, 2.0, 2.0)}
    path = tmp_path / "1TST.pdb"
    path.write_text(_pdb(_site(ANCHORS), ligand))
    arrays = load_structure_arrays(path)
    atoms = next(a for _k, _c, a in find_ligand_instances(arrays, comp_ids=("IHP",)))
    names, elements = tc.ligand_identity(arrays, atoms)
    _points, _labels, xyz = tc.tf.template_anchors(arrays, atoms)
    assert len(names) == len(elements) == len(xyz)
    assert names == ["P1", "O11"] and elements == ["P", "O"]


def test_placed_relative_sasa_separates_a_buried_from_an_exposed_ligand(tc, tmp_path):
    """A ligand inside a shell of protein must read far more buried than one in free space."""
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    rng = np.random.default_rng(0)
    shell = []
    directions = rng.normal(size=(60, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    for i, d in enumerate(directions * 6.0, start=1):
        x, y, z = d
        shell.append(("ALA", i, {"N": (x, y, z), "CA": (x * 1.05, y * 1.05, z * 1.05),
                                 "C": (x * 1.1, y * 1.1, z * 1.1), "O": (x * 1.15, y * 1.15, z * 1.15)}))
    path = tmp_path / "shell.pdb"
    path.write_text(_pdb(shell))
    arrays = load_structure_arrays(path)

    names, elements = ["P1", "O11", "O12"], ["P", "O", "O"]
    inside = np.array([[0.0, 0, 0], [1.4, 0, 0], [-1.4, 0, 0]])
    outside = inside + np.array([60.0, 0, 0])
    buried = tc.placed_relative_sasa(arrays, inside, names, elements)
    exposed = tc.placed_relative_sasa(arrays, outside, names, elements)
    assert buried is not None and exposed is not None
    assert buried < exposed
    assert exposed > tc.BURIED_MAX_REL_SASA


def test_the_report_survives_a_run_with_nothing_scored(tc, tmp_path):
    """Study J's report crashed on an empty frame; this pins the sibling against it."""
    import json

    for payload in ([], [{"uniprot_id": "X", "arm": "candidates", "score": None,
                          "error": "no AlphaFold model"}]):
        fits = tmp_path / "fits.json"
        fits.write_text(json.dumps(payload))
        out, md = tmp_path / "o.json", tmp_path / "O.md"
        assert tc.main(["report", "--fits", str(fits), "--out", str(out), "--markdown", str(md)]) == 0
        result = json.loads(out.read_text())
        assert result["O1"]["verdict"] == "not evaluable"
        assert result["O2"]["decision"] == "not run"
        assert md.read_text().strip()
