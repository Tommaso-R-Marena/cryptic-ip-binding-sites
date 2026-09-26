"""cryptic_ip/rescoring and scripts/rerank*.py (docs/RERANK_PLAN.md) on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.rescoring import crossfit as cf
from cryptic_ip.rescoring.electrostatics import (COULOMB, KAPPA, ChargedAtoms, ReceptorField, read_pdbqt,
                                                 read_pdbqt_models)

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _line(serial, xyz, q, atype="OA"):
    return (f"ATOM  {serial:5d}  O   LIG A   1    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00"
            f"    {q:6.3f} {atype:<2s}")


# ------------------------------------------------------------ electrostatics
def test_screened_coulomb_of_two_charges():
    field = ReceptorField(ChargedAtoms(np.array([[0.0, 0.0, 0.0]]), np.array([1.0])))
    e = field.energy(ChargedAtoms(np.array([[3.0, 0.0, 0.0]]), np.array([-1.0])))
    assert e == pytest.approx(-COULOMB * np.exp(-KAPPA * 3.0) / (4.0 * 9.0))
    assert field.energy(ChargedAtoms(np.array([[13.0, 0.0, 0.0]]), np.array([-1.0]))) == 0.0  # beyond 12 Å
    close = field.energy(ChargedAtoms(np.array([[0.5, 0.0, 0.0]]), np.array([-1.0])))
    assert close == pytest.approx(-COULOMB * np.exp(-KAPPA * 1.5) / (4.0 * 1.5 ** 2))  # floored at 1.5 Å


def test_pdbqt_parsing_single_and_models():
    text = "\n".join([_line(1, (0, 0, 0), -0.5), _line(2, (1, 2, 3), 0.25)])
    atoms = read_pdbqt(text)
    assert atoms.xyz.shape == (2, 3) and atoms.net == pytest.approx(-0.25)
    models = read_pdbqt_models("MODEL 1\n" + text + "\nENDMDL\nMODEL 2\n" + _line(1, (5, 5, 5), -1.0) + "\nENDMDL\n")
    assert len(models) == 2 and models[1].xyz[0, 0] == pytest.approx(5.0) and models[1].net == pytest.approx(-1.0)


# ---------------------------------------------------------------- crossfit
def _runs(n_groups=12, per_group=2, informative=True, seed=0):
    """Vina ranks a wrong pose first; the correct pose (rmsd 1) is second. E_el marks it when informative."""
    rng = np.random.default_rng(seed)
    runs = []
    for g in range(n_groups):
        for c in range(per_group):
            for s in (1, 2, 3):
                n = 10
                vina = np.sort(rng.normal(-6, 1, n))
                rmsd = rng.uniform(4, 10, n)
                rmsd[1] = 1.0
                eel = rng.normal(0, 1, n)
                if informative:
                    eel[1] = -200.0
                runs.append(cf.Run(f"G{g}C{c}", f"G{g}", s, vina, rmsd, eel))
    return runs


def test_top_index_and_success_matrix():
    run = cf.Run("a", "g", 1, np.array([-7.0, -6.0]), np.array([5.0, 1.0]), np.array([0.0, -100.0]))
    assert cf.top_index(run, 0.0) == 0 and cf.top_index(run, 0.05) == 1
    succ = cf.success_matrix([run], ["a"], grid=(0.0, 0.05))
    assert succ.tolist() == [[0.0, 1.0]]


def test_folds_keep_groups_together():
    codes = np.array([0, 0, 1, 2, 2, 2, 3, 4, 5, 6])
    folds = cf.fold_of_groups(codes)
    for g in np.unique(codes):
        assert len(set(folds[codes == g])) == 1
    assert set(folds) <= set(range(cf.N_FOLDS))


def test_informative_term_improves_and_permutation_removes_it():
    res = cf.evaluate(_runs(informative=True), n_bootstrap=200, n_permutations=40)
    assert res["F1"]["decision"] == "improves"
    assert res["F1"]["group"]["vina"]["point"] == pytest.approx(0.0)
    assert res["F1"]["group"]["reranked"]["point"] == pytest.approx(1.0)
    assert res["F3"]["mean"] < 0.3 and res["F3"]["fraction_at_least_observed"] == 0.0
    shares = res["F2"]["per_run_shares"]
    assert shares["vina"]["scoring_failure"] == pytest.approx(1.0) and shares["reranked"]["ok"] == pytest.approx(1.0)


def test_uninformative_term_is_not_an_improvement():
    """Vina is always wrong here, so noise alone lifts success: the permutation gate must catch it."""
    res = cf.evaluate(_runs(informative=False, seed=3), n_bootstrap=200, n_permutations=40)
    assert res["F1"]["decision"] in ("gain not specific to electrostatics", "no detectable difference")
    assert res["F3"]["p_value"] >= 0.05


def test_decision_rule_as_amended():
    ci = {"low": 0.01, "high": 0.2}
    assert cf.decide(ci, 0.01, 10) == "improves"
    assert cf.decide(ci, 0.2, 10) == "gain not specific to electrostatics"
    assert cf.decide({"low": -0.2, "high": -0.01}, 0.9, 10) == "worsens"
    assert cf.decide({"low": -0.1, "high": 0.1}, 0.01, 10) == "no detectable difference"
    assert cf.decide(ci, 0.01, 4).startswith("not evaluable")


def test_few_groups_are_not_evaluable():
    res = cf.evaluate(_runs(n_groups=3), n_bootstrap=50, n_permutations=10)
    assert res["F1"]["decision"].startswith("not evaluable")


def test_crossfit_chooses_zero_when_nothing_helps():
    succ = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 0.0], [1.0, 1.0]])
    codes = np.array([0, 1, 2, 3])
    oof, chosen = cf.crossfit(succ, codes, np.array([0, 1, 2, 3]), n_folds=4)
    assert chosen == [0, 0, 0, 0] and oof.tolist() == succ[:, 0].tolist()


# ---------------------------------------------------------------- end to end
@pytest.mark.skipif(shutil.which("pdb2pqr") is None and shutil.which("pdb2pqr30") is None,
                    reason="pdb2pqr not installed")
def test_dock_and_report_end_to_end(tmp_path):
    pytest.importorskip("vina")
    pytest.importorskip("meeko")
    pytest.importorskip("gemmi")
    sys.path.insert(0, str(ROOT / "scripts"))
    helpers = _load("test_redocking_script_helpers", ROOT / "tests" / "test_redocking_script.py")
    redock = _load("redocking", ROOT / "scripts" / "redocking.py")
    rerank = _load("rerank", ROOT / "scripts" / "rerank.py")
    report = _load("rerank_report", ROOT / "scripts" / "rerank_report.py")
    helpers._synthetic_entry(tmp_path)
    census = tmp_path / "census.csv"
    assert redock.main(["census", "--table", str(tmp_path / "table.csv.gz"), "--entries",
                        str(tmp_path / "entries.csv"), "--structures-dir", str(tmp_path), "--ccd-dir",
                        str(tmp_path / "ccd"), "--output", str(census)]) == 0
    assert rerank.main(["dock", "--census", str(census), "--structures-dir", str(tmp_path), "--ccd-dir",
                        str(tmp_path / "ccd"), "--out-dir", str(tmp_path / "arms"), "--work-dir",
                        str(tmp_path / "work"), "--exhaustiveness", "1"]) == 0
    record = json.loads((tmp_path / "arms" / "rerank_0.jsonl").read_text().splitlines()[0])
    assert "error" not in record, record
    seeds = [r for r in record["runs"] if r["seed"] >= 1]
    assert len(seeds) == 3
    for r in seeds:
        assert len(r["vina"]) == len(r["rmsd"]) == len(r["eel"]) >= 1
        assert r["ligand_net_charge"] == pytest.approx(-9.0, abs=0.05)  # IP6, primary state
    # a -9 ligand next to a Lys/Arg-rich peptide: attractive on balance
    assert np.median([e for r in seeds for e in r["eel"]]) < 0
    out = tmp_path / "rep"
    assert report.main(["--census", str(census), "--arms-dir", str(tmp_path / "arms"), "--out-dir", str(out),
                        "--n-bootstrap", "20", "--n-permutations", "10"]) == 0
    res = json.loads((out / "rerank.json").read_text())
    assert res["copies"] == 1 and res["F1"]["decision"].startswith("not evaluable")
    assert (out / "report.html").read_text().startswith("<!DOCTYPE html>")
    assert pd.read_csv(census)["selected"].iloc[0]


# ------------------------------------------------------- more electrostatics
def test_receptor_field_matches_brute_force():
    rng = np.random.default_rng(4)
    rec = ChargedAtoms(rng.uniform(-15, 15, (400, 3)), rng.normal(0, 0.4, 400))
    rec.charge[::7] = 0.0
    lig = ChargedAtoms(rng.uniform(-4, 4, (30, 3)), rng.normal(-0.3, 0.3, 30))
    r = np.linalg.norm(lig.xyz[:, None, :] - rec.xyz[None, :, :], axis=2)
    within = r <= 12.0
    rf = np.maximum(r, 1.5)
    brute = COULOMB * np.sum(np.where(within, lig.charge[:, None] * rec.charge[None, :] * np.exp(-KAPPA * rf)
                                      / (4.0 * rf ** 2), 0.0))
    assert ReceptorField(rec).energy(lig) == pytest.approx(brute, rel=1e-9)


def test_energy_is_additive_and_sign_symmetric():
    rng = np.random.default_rng(5)
    field = ReceptorField(ChargedAtoms(rng.uniform(-8, 8, (100, 3)), rng.normal(0, 0.5, 100)))
    a = ChargedAtoms(rng.uniform(-3, 3, (5, 3)), rng.normal(0, 1, 5))
    b = ChargedAtoms(rng.uniform(-3, 3, (7, 3)), rng.normal(0, 1, 7))
    both = ChargedAtoms(np.vstack([a.xyz, b.xyz]), np.r_[a.charge, b.charge])
    assert field.energy(both) == pytest.approx(field.energy(a) + field.energy(b))
    assert field.energy(ChargedAtoms(a.xyz, -a.charge)) == pytest.approx(-field.energy(a))


def test_uncharged_receptor_or_ligand_gives_zero():
    field = ReceptorField(ChargedAtoms(np.zeros((3, 3)), np.zeros(3)))
    assert field.energy(ChargedAtoms(np.ones((2, 3)), np.array([-1.0, -1.0]))) == 0.0
    charged = ReceptorField(ChargedAtoms(np.zeros((1, 3)), np.array([1.0])))
    assert charged.energy(ChargedAtoms(np.ones((2, 3)), np.zeros(2))) == 0.0
    assert charged.energy(ChargedAtoms(np.zeros((0, 3)), np.zeros(0))) == 0.0


def test_models_without_model_records_read_as_one_pose():
    text = "\n".join([_line(1, (0, 0, 0), -0.5), _line(2, (1, 0, 0), 0.5)])
    models = read_pdbqt_models(text)
    assert len(models) == 1 and models[0].xyz.shape == (2, 3)
    assert read_pdbqt(text + "\nENDMDL\n" + _line(3, (9, 9, 9), 1.0)).xyz.shape == (2, 3)  # first model only


# ------------------------------------------------------------ more crossfit
def test_zero_weight_groups_do_not_train():
    succ = np.array([[0.0, 1.0], [0.0, 1.0], [1.0, 0.0], [1.0, 0.0]])
    codes = np.array([0, 1, 2, 3])
    folds = np.array([0, 0, 1, 1])
    _, chosen = cf.crossfit(succ, codes, folds, weights=np.array([1.0, 1.0, 0.0, 0.0]), n_folds=2)
    assert chosen[1] == 1  # fold 1 trains on groups 0 and 1 only, where w = grid[1] wins
    assert chosen[0] == 0  # fold 0's training groups carry no weight: fall back to Vina (w = 0)


def test_success_matrix_averages_seeds_and_skips_nothing():
    runs = [cf.Run("a", "g", s, np.array([-7.0, -6.0]), np.array([5.0, 1.0]), np.array([0.0, -100.0 * (s == 1)]))
            for s in (1, 2, 3)]
    succ = cf.success_matrix(runs, ["a"], grid=(0.0, 0.05))
    assert succ[0, 0] == 0.0 and succ[0, 1] == pytest.approx(1 / 3)


def test_evaluate_drops_runs_with_missing_rmsd():
    good = _runs(n_groups=6)
    bad = [cf.Run("X", "GX", 1, np.array([-7.0]), np.array([np.nan]), np.array([0.0]))]
    res = cf.evaluate(good + bad, n_bootstrap=20, n_permutations=10)
    assert res["copies"] == len({r.copy_key for r in good})


def test_stratum_marks_small_strata():
    runs = _runs(n_groups=3)
    s = cf.stratum(runs, {r.copy_key for r in runs}, {r.copy_key: 0.0 for r in runs})
    assert s["groups"] == 3 and not s["evidence"] and s["vina"] == pytest.approx(0.0)
    assert cf.stratum(runs, set(), {})["copies"] == 0


def test_permutation_control_is_reproducible():
    runs = _runs(n_groups=6)
    a = cf.evaluate(runs, n_bootstrap=20, n_permutations=10)["F3"]
    b = cf.evaluate(runs, n_bootstrap=20, n_permutations=10)["F3"]
    assert a["values"] == b["values"] and 1 / 11 <= a["p_value"] <= 1.0


# ---------------------------------------------------------------- report
@pytest.fixture()
def rerank_report():
    sys.path.insert(0, str(ROOT / "scripts"))
    return _load("rerank_report", ROOT / "scripts" / "rerank_report.py")


def _records_from_runs(runs, extra=()):
    by_copy = {}
    for r in runs:
        by_copy.setdefault(r.copy_key, []).append({"seed": r.seed, "vina": r.vina.tolist(), "rmsd": r.rmsd.tolist(),
                                                   "eel": r.eel.tolist()})
    recs = [{"copy_key": k, "pdb_id": k[:4], "runs": v + [{"seed": 0, "crystal_minimised_vina": -8.0,
                                                           "crystal_minimised_eel": -150.0}]}
            for k, v in by_copy.items()]
    return recs + list(extra)


def test_report_on_synthetic_records(rerank_report, tmp_path):
    runs = _runs(n_groups=10)
    extra = [{"copy_key": "ZZZZ:A:1", "pdb_id": "ZZZZ", "error": "not reached: shard time budget"},
             {"copy_key": "YYYY:A:1", "pdb_id": "YYYY", "error": "RuntimeError: boom"}]
    census = pd.DataFrame({"copy_key": sorted({r.copy_key for r in runs}) + ["ZZZZ:A:1", "YYYY:A:1"]})
    census["homology_group_strict"] = [k.split("C")[0] if k[0] == "G" else "GZ" for k in census["copy_key"]]
    census["burial_class"] = ["cryptic" if k.startswith("G1") else "surface" for k in census["copy_key"]]
    census["icode"] = ""
    census.to_csv(tmp_path / "census.csv", index=False)
    (tmp_path / "arms").mkdir()
    (tmp_path / "arms" / "rerank_0.jsonl").write_text(
        "\n".join(json.dumps(r) for r in _records_from_runs(runs, extra)) + "\n")
    first = {"arm": "vina_s1", "seed": 1, "top_rmsd": float(runs[0].rmsd[0])}
    primary = [{"copy_key": runs[0].copy_key, "runs": [first, {"arm": "crystal_control", "crystal_score": -8.0}]}]
    (tmp_path / "prim").mkdir()
    (tmp_path / "prim" / "primary_0.jsonl").write_text(json.dumps(primary[0]) + "\n")
    assert rerank_report.main(["--census", str(tmp_path / "census.csv"), "--arms-dir", str(tmp_path / "arms"),
                               "--primary-dir", str(tmp_path / "prim"), "--out-dir", str(tmp_path / "out"),
                               "--n-bootstrap", "50", "--n-permutations", "20"]) == 0
    res = json.loads((tmp_path / "out" / "rerank.json").read_text())
    assert res["F1"]["decision"] == "improves"
    assert res["accounting"] == {"records": len(census), "errors": 2, "not_reached": 1,
                                 "error_examples": ["not reached: shard time budget", "RuntimeError: boom"]}
    assert res["reproducibility_vs_redocking_primary"] == {"pairs": 1, "same_success": 1.0, "within_0.5A": 1.0}
    assert res["crystal_minimised"]["copies"] == 20 and res["F4"]["cryptic"]["groups"] == 1
    assert not res["F4"]["cryptic"]["evidence"] and res["F4"]["classic arrestins"]["copies"] == 0
    md = (tmp_path / "out" / "RERANK.md").read_text()
    assert "**F1:** improves." in md and "not evidence" in md
    assert (tmp_path / "out" / "report.html").read_text().startswith("<!DOCTYPE html>")
