"""scripts/sampling.py and scripts/sampling_report.py (docs/SAMPLING_PLAN.md) on synthetic data."""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path=None):
    path = path or ROOT / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def sampling():
    sys.path.insert(0, str(ROOT / "scripts"))
    return _load("sampling")


@pytest.fixture(scope="module")
def report():
    sys.path.insert(0, str(ROOT / "scripts"))
    return _load("sampling_report")


# ------------------------------------------------------------ the subset
def test_subset_is_deterministic_stratified_and_outcome_blind(sampling):
    frame = pd.DataFrame({"copy_key": [f"{i:04d}:A:1" for i in range(250)],
                          "burial_class": ["surface"] * 186 + ["semi_cryptic"] * 51 + ["cryptic"] * 13})
    keys = sampling.subset_keys(frame, size=80)
    assert keys == sampling.subset_keys(frame.sample(frac=1, random_state=3), size=80)  # row order cannot matter
    assert len(keys) == len(set(keys)) and 70 <= len(keys) <= 90
    by_class = frame.set_index("copy_key")["burial_class"]
    counts = by_class.loc[keys].value_counts()
    assert counts["surface"] > counts["semi_cryptic"] > counts["cryptic"]  # keeps the class order
    assert all(counts[c] >= 1 for c in ("surface", "semi_cryptic", "cryptic"))


def test_subset_handles_a_class_smaller_than_its_share(sampling):
    frame = pd.DataFrame({"copy_key": [f"{i:03d}:A:1" for i in range(12)],
                          "burial_class": ["surface"] * 10 + ["cryptic"] * 2})
    keys = sampling.subset_keys(frame, size=100)  # asks for more than exist
    assert len(keys) == 12


def test_arm_definitions_match_the_plan(sampling):
    assert sampling.ARMS["e128"] == (128, (1, 2, 3), False)
    assert sampling.ARMS["e512"] == (512, (1,), True)


# --------------------------------------------------------- run metrics
def test_run_metrics_uses_the_frozen_weight(report):
    run = {"seed": 1, "vina": [-7.0, -6.0, -5.0], "rmsd": [5.0, 1.0, 9.0], "eel": [0.0, -100.0, 0.0],
           "seconds": 12.0}
    m = report.run_metrics(run)
    assert m["vina"] == 0.0 and m["reranked"] == 1.0 and m["ceiling"] == 1.0 and m["poses"] == 3.0
    # a weight of zero reproduces Vina exactly
    assert report.run_metrics(run, w=0.0)["reranked"] == 0.0
    assert report.FROZEN_W == 0.1


def test_run_metrics_rejects_unusable_lists(report):
    assert report.run_metrics({"seed": 1, "vina": [], "rmsd": [], "eel": []}) is None
    assert report.run_metrics({"seed": 1, "vina": [-7.0], "rmsd": [float("nan")], "eel": [0.0]}) is None


def test_ceiling_is_the_best_of_the_list(report):
    m = report.run_metrics({"seed": 1, "vina": [-7.0, -6.0], "rmsd": [9.0, 1.9], "eel": [0.0, 0.0]})
    assert m["ceiling"] == 1.0 and m["vina"] == 0.0  # sampled but not ranked first: a scoring failure


# ------------------------------------------------------------- decisions
def _records(keys, *, seeds=(1, 2, 3), ceiling=0.0, top=0.0, rng=None, seconds=100.0):
    """Pose lists engineered so that a copy's ceiling/top-pose outcome is exactly as asked."""
    rng = rng or np.random.default_rng(0)
    out = []
    for key in keys:
        runs = []
        for s in seeds:
            n = 8
            vina = np.sort(rng.normal(-6, 0.5, n))
            rmsd = rng.uniform(4, 9, n)
            eel = np.zeros(n)
            if rng.random() < top:
                rmsd[0] = 1.0            # the top-ranked pose is near native
            elif rng.random() < ceiling:
                rmsd[3] = 1.0            # a near-native pose exists but is not ranked first
            runs.append({"seed": s, "vina": vina.tolist(), "rmsd": rmsd.tolist(), "eel": eel.tolist(),
                         "seconds": seconds})
        out.append({"copy_key": key, "pdb_id": key[:4], "runs": runs})
    return out


def _census(n=60, n_groups=12):
    return pd.DataFrame({"copy_key": [f"{i:04d}:A:1" for i in range(n)],
                         "homology_group_strict": [f"G{i % n_groups}" for i in range(n)],
                         "burial_class": [("surface", "semi_cryptic", "cryptic")[i % 3] for i in range(n)],
                         "icode": ""})


def test_a_higher_ceiling_is_budget_limited(report):
    census = _census()
    keys = census["copy_key"].tolist()
    e32 = _records(keys, ceiling=0.0, top=0.0, rng=np.random.default_rng(1))
    e128 = _records(keys, ceiling=1.0, top=0.0, rng=np.random.default_rng(2))
    r = report.build(census, e32, {"e128": e128}, n_bootstrap=200)
    assert r["decisions"]["G1"]["decision"] == "budget-limited"
    assert r["decisions"]["G1"]["estimate"]["point"] > 0.5
    # the ceiling rose while the top pose did not: the report must not call that an improvement
    assert r["decisions"]["G2"]["decision"] in ("no gain", "inconclusive")


def test_an_unchanged_ceiling_is_search_saturated(report):
    census = _census()
    keys = census["copy_key"].tolist()
    e32 = _records(keys, ceiling=0.3, top=0.1, rng=np.random.default_rng(5))
    e128 = _records(keys, ceiling=0.3, top=0.1, rng=np.random.default_rng(5))  # identical draws
    r = report.build(census, e32, {"e128": e128}, n_bootstrap=200)
    assert r["decisions"]["G1"]["decision"] == "search-saturated"
    assert r["decisions"]["G1"]["estimate"]["point"] == pytest.approx(0.0)


def test_few_groups_are_not_evaluable(report):
    census = _census(n=9, n_groups=3)
    keys = census["copy_key"].tolist()
    r = report.build(census, _records(keys), {"e128": _records(keys)}, n_bootstrap=50)
    assert all(d["decision"].startswith("not evaluable") for d in r["decisions"].values())


def test_holm_covers_the_three_primary_questions(report):
    census = _census()
    keys = census["copy_key"].tolist()
    r = report.build(census, _records(keys, ceiling=0.2, top=0.1, rng=np.random.default_rng(7)),
                     {"e128": _records(keys, ceiling=0.9, top=0.6, rng=np.random.default_rng(8))},
                     n_bootstrap=200)
    assert set(r["decisions"]) == {"G1", "G2", "G3"}
    for d in r["decisions"].values():
        est = d.get("estimate")
        if d.get("holm_p") is not None and est:
            assert d["holm_p"] >= est["p_value"] - 1e-12


def test_only_copies_in_both_arms_are_compared(report):
    census = _census()
    keys = census["copy_key"].tolist()
    e32 = _records(keys)
    e128 = _records(keys[:30]) + [{"copy_key": keys[31], "pdb_id": "x", "error": "not reached: shard time budget"}]
    r = report.build(census, e32, {"e128": e128}, n_bootstrap=50)
    assert r["copies_compared"] == 30
    assert r["accounting"]["e128"]["not_reached"] == 1


def test_rescue_rate_and_ladder_and_outputs(report, tmp_path):
    census = _census()
    keys = census["copy_key"].tolist()
    e32 = _records(keys, ceiling=0.0, top=0.0, rng=np.random.default_rng(11))
    e128 = _records(keys, ceiling=1.0, top=0.0, rng=np.random.default_rng(12))
    e512 = _records(keys[:20], seeds=(1,), ceiling=1.0, top=0.5, rng=np.random.default_rng(13))
    (tmp_path / "e32").mkdir()
    (tmp_path / "arms").mkdir()
    (tmp_path / "e32" / "rerank_0.jsonl").write_text("\n".join(json.dumps(r) for r in e32) + "\n")
    (tmp_path / "arms" / "e128_0.jsonl").write_text("\n".join(json.dumps(r) for r in e128) + "\n")
    (tmp_path / "arms" / "e512_0.jsonl").write_text("\n".join(json.dumps(r) for r in e512) + "\n")
    census.to_csv(tmp_path / "census.csv", index=False)
    assert report.main(["--census", str(tmp_path / "census.csv"), "--e32-dir", str(tmp_path / "e32"),
                        "--arms-dir", str(tmp_path / "arms"), "--out-dir", str(tmp_path / "out"),
                        "--n-bootstrap", "50"]) == 0
    r = json.loads((tmp_path / "out" / "sampling.json").read_text())
    g4 = r["G4_rescue"]
    assert g4["e32_seed_runs_without_a_near_native_pose"] == 180  # every run missed at E32
    assert g4["e128"]["rescued"]["per_group"]["point"] == pytest.approx(1.0)
    assert "e512" in g4 and r["E512_ladder_seed1"]["e512"]["copies"] == 20
    assert set(r["strata"]) == {"surface", "semi_cryptic", "cryptic"}
    text = (tmp_path / "out" / "SAMPLING.md").read_text()
    assert "### Decisions (Holm across G1-G3)" in text and "| G1 |" in text
    assert (tmp_path / "out" / "report.html").read_text().startswith("<!DOCTYPE html>")


def test_no_shared_copies_is_reported_not_crashed(report):
    census = _census()
    r = report.build(census, _records(["9999:A:1"]), {}, n_bootstrap=20)
    assert r["decisions"]["G1"]["decision"].startswith("not evaluable")


# ------------------------------------------------------------ end to end
@pytest.mark.skipif(shutil.which("pdb2pqr") is None and shutil.which("pdb2pqr30") is None,
                    reason="pdb2pqr not installed")
def test_dock_matches_study_f_on_the_same_settings(tmp_path):
    """Study G's dock routine mirrors study F's; this pins them on one synthetic copy."""
    pytest.importorskip("vina")
    pytest.importorskip("meeko")
    pytest.importorskip("gemmi")
    sys.path.insert(0, str(ROOT / "scripts"))
    helpers = _load("helpers_for_sampling", ROOT / "tests" / "test_redocking_script.py")
    redock = _load("redocking", ROOT / "scripts" / "redocking.py")
    rerank = _load("rerank", ROOT / "scripts" / "rerank.py")
    sampling = _load("sampling", ROOT / "scripts" / "sampling.py")
    helpers._synthetic_entry(tmp_path)
    census = tmp_path / "census.csv"
    assert redock.main(["census", "--table", str(tmp_path / "table.csv.gz"), "--entries",
                        str(tmp_path / "entries.csv"), "--structures-dir", str(tmp_path), "--ccd-dir",
                        str(tmp_path / "ccd"), "--output", str(census)]) == 0
    row = pd.read_csv(census, dtype={"icode": str}).to_dict(orient="records")[0]
    row["icode"] = ""
    ctx = redock.Context(row, tmp_path, tmp_path / "ccd", tmp_path / "work")
    mine = sampling.dock_copy(ctx, exhaustiveness=1, seeds=(1,))
    rerank.EXHAUSTIVENESS[0] = 1
    theirs = rerank.dock_copy(ctx)
    seed1 = [r for r in theirs if r.get("seed") == 1][0]
    assert set(mine[0]) >= {"seed", "vina", "rmsd", "eel"}
    assert mine[0]["vina"] == seed1["vina"]      # same receptor, ligand, box, starting pose and seed
    assert mine[0]["rmsd"] == seed1["rmsd"]
    assert mine[0]["eel"] == pytest.approx(seed1["eel"])
