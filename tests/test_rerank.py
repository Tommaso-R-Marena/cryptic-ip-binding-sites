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
