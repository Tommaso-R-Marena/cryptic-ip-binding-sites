"""scripts/redocking_report.py on a realistic synthetic benchmark: every arm, failures and re-dispatches.

The report job has one shot at the real arms, so its bookkeeping is exercised here
on 48 copies in 9 strict groups with every kind of record the dock arms write.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CLASSES = ("surface", "semi_cryptic", "cryptic")


@pytest.fixture(scope="module")
def rr():
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location("redocking_report", ROOT / "scripts" / "redocking_report.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["redocking_report"] = module
    spec.loader.exec_module(module)
    return module


def _census(n=48, n_groups=9):
    rows = []
    for i in range(n):
        rows.append({"copy_key": f"{i:04d}:A:{i}", "pdb_id": f"{i:04d}", "selected": True,
                     "primary_set": i % 11 != 0, "status": "eligible", "homology_group_strict": f"G{i % n_groups}",
                     "burial_class": CLASSES[i % 3], "resolution": 1.5 + (i % 5) * 0.5,
                     "experimental_method": "X-RAY DIFFRACTION" if i % 7 else "ELECTRON MICROSCOPY",
                     "species": ("InsP6", "InsP3", "InsP5", "IP7")[i % 4], "metal": i % 2 == 0,
                     "interface": i % 3 == 0, "icode": ""})
    rows.append({**rows[1], "copy_key": "9999:A:1", "pdb_id": "9999", "selected": False,
                 "status": "excluded: configuration differs from the CCD"})
    return pd.DataFrame(rows)


def _vina(rng, arm, seed, rmsd=None, score=None):
    r = float(rng.uniform(0.3, 8.0)) if rmsd is None else rmsd
    return {"arm": arm, "seed": seed, "top_rmsd": r, "best_rmsd": min(r, float(rng.uniform(0.3, 4.0))),
            "top_rmsd_p": 0.8 * r, "top_score": float(rng.normal(-7, 1)) if score is None else score,
            "spearman": float(rng.uniform(-1, 1)), "scores": [], "rmsd": []}


def _records(census, seed=0):
    rng = np.random.default_rng(seed)
    recs = {"primary": [], "secondary": [], "pockets": [], "alphafold": []}
    for i, row in census[census["selected"]].iterrows():
        key, pdb = row["copy_key"], row["pdb_id"]
        if i == 3:
            recs["primary"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "primary",
                                    "error": "not reached: shard time budget"})
            continue
        runs = [_vina(rng, f"vina_s{s}", s) for s in (1, 2, 3)]
        if i == 5:
            runs[1]["top_rmsd"] = None
        runs[0]["receptor"] = {"his": {"A:10": "HID"}}
        runs.append({"arm": "crystal_control", "crystal_score": -8.0,
                     "crystal_minimised_score": float(rng.normal(-8, 1)), "minimised_rmsd": 0.4})
        recs["primary"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "primary", "runs": runs})
        sec = [_vina(rng, "vinardo_s1", 1),
               {"arm": "ad4_s1", "error": "autogrid4 failed"} if i == 7 else _vina(rng, "ad4_s1", 1),
               _vina(rng, "vina_deprotonated_s1", 1)]
        if row["metal"]:
            sec += [_vina(rng, f"vina_metals_s{s}", s) for s in (1, 2, 3)]
        recs["secondary"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "secondary", "runs": sec})
        if row["primary_set"]:
            recs["pockets"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "pockets", "runs": [
                _vina(rng, "decoy_s1", 1),
                {"arm": "site_finding", "lands_in_true_site": bool(i % 2),
                 "any_positive_pocket": None if i == 9 else bool(i % 3)}]})
            if i % 4 == 0:
                recs["alphafold"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "alphafold",
                                          "error": "no AlphaFold model"})
            else:
                recs["alphafold"].append({"copy_key": key, "pdb_id": pdb, "arm_group": "alphafold", "runs": [
                    {**_vina(rng, "alphafold_s1", 1), "site_ca_rmsd": 1.0, "accession": "P1", "identity": 1.0}]})
    return recs


def _write(directory: Path, recs):
    directory.mkdir(parents=True, exist_ok=True)
    for group, rs in recs.items():
        (directory / f"{group}_0.jsonl").write_text("\n".join(json.dumps(r) for r in rs) + "\n")


@pytest.fixture(scope="module")
def run(rr, tmp_path_factory):
    tmp = tmp_path_factory.mktemp("redock_report")
    census = _census()
    census.to_csv(tmp / "census.csv", index=False)
    _write(tmp / "arms", _records(census))
    assert rr.main(["--census", str(tmp / "census.csv"), "--results-dir", str(tmp / "arms"),
                    "--out-dir", str(tmp / "out"), "--n-bootstrap", "100"]) == 0
    return {"census": census, "dir": tmp, "report": json.loads((tmp / "out" / "redocking.json").read_text()),
            "copies": pd.read_csv(tmp / "out" / "copies.csv")}


def test_accounting(run):
    acc = run["report"]["accounting"]
    census = run["census"]
    assert acc["copies_found"] == len(census) and acc["copies_selected"] == int(census["selected"].sum())
    assert acc["primary_set"] == int((census["selected"] & census["primary_set"]).sum())
    assert acc["incomplete_selected"] == int((census["selected"] & ~census["primary_set"]).sum())
    assert acc["arms"]["primary"]["failed_or_not_reached"] == {"not reached: shard time budget": 1}
    assert acc["arms"]["alphafold"]["failed_or_not_reached"]["no AlphaFold model"] > 0
    assert "9999:A:1" not in set(run["copies"]["copy_key"])  # unselected copies never enter the table


def test_per_copy_outcomes(run):
    copies = run["copies"].set_index("copy_key")
    assert copies.loc["0003:A:3", "primary_error"].startswith("not reached")
    assert np.isnan(copies.loc["0003:A:3", "success"])
    assert copies.loc["0005:A:5", "n_seeds"] == 2  # one seed without an RMSD counts as missing, not a failure
    row = copies.loc["0001:A:1"]
    seeds = [row[f"top_rmsd_s{s}"] <= 2.0 for s in (1, 2, 3)]
    assert row["success"] == pytest.approx(np.mean(seeds))
    assert copies.loc["0007:A:7", "ad4_note"] == "autogrid4 failed" and np.isnan(copies.loc["0007:A:7", "ad4_success"])
    assert np.isnan(copies.loc["0001:A:1", "metals_success"]) and np.isfinite(copies.loc["0002:A:2", "metals_success"])


def test_failure_decomposition_is_complete(run):
    copies = run["copies"]
    prim = copies[copies["primary_set"]]
    failed = prim[pd.to_numeric(prim["top_rmsd_s1"], errors="coerce") > 2.0]
    fd = run["report"]["controls"]["failure_decomposition"]
    assert fd["failed_seed1"] == len(failed) == fd.get("sampling", 0) + fd.get("scoring", 0)
    for _, r in failed.iterrows():
        expected = "sampling" if r["crystal_minimised_score"] < r["top_score_s1"] else "scoring"
        assert r["failure_kind"] == expected


def test_decisions_holm_and_evidence(run):
    dec = run["report"]["decisions"]
    assert set(dec) == {"R1", "R2", "R3", "R4"}
    assert dec["R2"]["decision"].startswith("not evaluable")  # burial class is constant within a group here
    for d in dec.values():
        est = d.get("estimate")
        if d.get("holm_p") is not None and est and "p_value" in est:
            assert d["holm_p"] >= est["p_value"] - 1e-12
    assert dec["R3"]["copies"] == int(run["copies"]["af_success"].notna().sum())


def test_strata_flag_small_groups(run):
    for levels in run["report"]["strata"].values():
        for e in levels.values():
            assert e["evidence"] == (e["groups"] >= 5)


def test_outputs_render(run):
    out = run["dir"] / "out"
    assert (out / "report.html").read_text().startswith("<!DOCTYPE html>")
    text = (out / "REPORT.md").read_text()
    assert "### Decisions" in text and "| R1 |" in text


def test_redispatch_replaces_not_reached_but_never_a_docked_record(rr, tmp_path):
    earlier, later = tmp_path / "earlier", tmp_path / "later"
    docked = {"copy_key": "k1", "arm_group": "primary", "runs": [{"arm": "vina_s1", "top_rmsd": 1.0}]}
    missing = {"copy_key": "k1", "arm_group": "primary", "error": "not reached: shard time budget"}
    other = {"copy_key": "k2", "arm_group": "primary", "error": "not reached: shard time budget"}
    _write(earlier, {"primary": [missing, other]})
    _write(later, {"primary": [docked]})
    recs = rr.read_records(later, earlier)
    by_key = {r["copy_key"]: r for r in recs["primary"]}
    assert "runs" in by_key["k1"] and "error" in by_key["k2"]
    _write(earlier, {"primary": [docked]})
    _write(later, {"primary": [missing]})
    assert "runs" in {r["copy_key"]: r for r in rr.read_records(later, earlier)["primary"]}["k1"]


def test_missing_arm_groups_do_not_break_the_report(rr, tmp_path):
    census = _census(n=12, n_groups=6)
    recs = _records(census)
    _write(tmp_path / "arms", {"primary": recs["primary"]})
    census.to_csv(tmp_path / "census.csv", index=False)
    assert rr.main(["--census", str(tmp_path / "census.csv"), "--results-dir", str(tmp_path / "arms"),
                    "--out-dir", str(tmp_path / "out"), "--n-bootstrap", "30"]) == 0
    rep = json.loads((tmp_path / "out" / "redocking.json").read_text())
    assert rep["decisions"]["R4"]["decision"] == "not evaluable: no decoy scores"
    assert rep["decisions"]["R3"]["decision"].startswith("not evaluable")
