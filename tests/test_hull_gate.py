"""scripts/hull_gate.py (docs/HULL_GATE_PLAN.md) on synthetic screens."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def hg():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("hull_gate", ROOT / "scripts" / "hull_gate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _screen(n=600, binders_shallow=True, seed=0):
    rng = np.random.default_rng(seed)
    rows, prots = [], []
    for i in range(n):
        acc = f"P{i:05d}"
        binder = i < 30
        # binders have high composite scores; their hull depth decides whether the gate hides them
        score = 0.8 + 0.1 * rng.random() if binder else 0.3 + 0.5 * rng.random()
        hull = (6.0 + rng.random() if binders_shallow else 15.0) if binder else 3 + 14 * rng.random()
        rows.append({"uniprot_id": acc, "pocket_id": 1, "composite_score": score, "hull_depth": hull,
                     "sasa": 20.0, "basic_residues": 6, "volume": 800.0, "plddt_mean": 90.0})
        prots.append({"uniprot_id": acc, "seen": i % 97 == 0 and not binder, "annotated": binder,
                      "cluster": f"C{i // 2}", "organism_key": "human"})
    return pd.DataFrame(rows), pd.DataFrame(prots)


def test_arms_differ_only_in_the_hull_gate(hg):
    a, b, c = (hg.arm_criteria(h) for h in (10.0, 5.0, None))
    for x in (b, c):
        da, dx = dict(a.__dict__), dict(x.__dict__)
        da.pop("min_hull_depth")
        dx.pop("min_hull_depth")
        assert da == dx
    assert (a.min_hull_depth, b.min_hull_depth, c.min_hull_depth) == (10.0, 5.0, None)


def test_gate_that_hides_binders_is_removed(hg):
    pockets, proteins = _screen(binders_shallow=True)
    arms = {n: hg.protein_scores(pockets, hg.arm_criteria(h)) for n, h in hg.ARMS.items()}
    result = hg.evaluate(proteins, arms, n_bootstrap=200)
    assert result["arms"]["none"]["roc_auc"]["point"] > result["arms"]["gate10"]["roc_auc"]["point"]
    assert hg.decide(result)["decision"] == "remove"


def test_gate_that_helps_is_kept(hg):
    pockets, proteins = _screen(binders_shallow=False)
    # make shallow non-binders score as high as binders: only the gate separates them
    binders = proteins.loc[proteins["annotated"], "uniprot_id"]
    shallow = (pockets["hull_depth"] < 10) & ~pockets["uniprot_id"].isin(binders)
    pockets.loc[shallow, "composite_score"] = 0.95
    arms = {n: hg.protein_scores(pockets, hg.arm_criteria(h)) for n, h in hg.ARMS.items()}
    result = hg.evaluate(proteins, arms, n_bootstrap=200)
    assert result["delta_none"]["point"] < 0
    assert hg.decide(result)["decision"] == "keep"


def test_seen_proteins_are_not_evaluated(hg):
    pockets, proteins = _screen()
    arms = {n: hg.protein_scores(pockets, hg.arm_criteria(h)) for n, h in hg.ARMS.items()}
    result = hg.evaluate(proteins, arms, n_bootstrap=20)
    assert result["unseen_proteins"] == int((~proteins["seen"]).sum())


def test_protein_without_eligible_pocket_ranks_last(hg):
    pockets = pd.DataFrame([{"uniprot_id": "A", "pocket_id": 1, "composite_score": 0.9, "hull_depth": 4.0,
                             "sasa": 1.0, "basic_residues": 6, "volume": 800.0, "plddt_mean": 90.0}])
    gated = hg.protein_scores(pockets, hg.arm_criteria(10.0))
    assert np.isnan(gated.loc["A", "rank_score"]) and not gated.loc["A", "hit"]
    free = hg.protein_scores(pockets, hg.arm_criteria(None))
    assert free.loc["A", "rank_score"] == pytest.approx(0.9) and free.loc["A", "hit"]
