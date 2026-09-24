"""The benchmark workflow's run matrix: complete, paired, and blind in smoke mode."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

WORKFLOW = Path(__file__).resolve().parents[1] / ".github" / "workflows" / "benchmark.yml"


def _plan(mode, tmp_path):
    workflow = yaml.safe_load(WORKFLOW.read_text())
    script = workflow["jobs"]["plan"]["steps"][0]["run"]
    out = tmp_path / "out"
    subprocess.run([sys.executable, "-c", script], check=True, env={**os.environ, "MODE": mode, "GITHUB_OUTPUT": str(out)})
    values = dict(line.split("=", 1) for line in out.read_text().splitlines())
    return json.loads(values["runs"])["include"], json.loads(values["shards"]), values


@pytest.mark.parametrize("mode", ["smoke", "full"])
def test_every_run_has_its_paired_arm(mode, tmp_path):
    runs, _, _ = _plan(mode, tmp_path)
    keys = {(r["task"], r["arm"], r["grouping"], r["repeat"], r["permute"]) for r in runs if not (r["permute"] and mode == "full")}
    for task, arm, grouping, repeat, permute in keys:
        other = "no_hull_depth" if arm == "full" else "full"
        assert (task, other, grouping, repeat, permute) in keys


def test_full_mode_matches_the_plan(tmp_path):
    runs, shards, values = _plan("full", tmp_path)
    real = [r for r in runs if not r["permute"]]
    for task in ("ip_site", "cryptic_ip_site", "burial"):
        seq = {r["repeat"] for r in real if r["task"] == task and r["grouping"] == "sequence"}
        strict = {r["repeat"] for r in real if r["task"] == task and r["grouping"] == "strict"}
        assert seq == {0, 1, 2} and strict == {0}
        # One permutation control per task, on the full arm.
        perms = [r for r in runs if r["permute"] and r["task"] == task]
        assert len(perms) == 1 and perms[0]["arm"] == "full"
    # The locked model scores the holdout once per task and arm.
    assert all(r["repeat"] == 0 and r["grouping"] == "sequence" and not r["permute"] for r in runs if r["holdout"])
    assert values["max_entries"] == "0" and values["draws"] == "10"
    assert len(runs) <= 256  # GitHub's matrix limit


def test_smoke_mode_is_blind(tmp_path):
    runs, _, values = _plan("smoke", tmp_path)
    assert all(r["permute"] for r in runs)
    assert values["max_entries"] != "0"


def test_push_runs_the_full_analysis_and_smoke_is_on_demand():
    workflow = yaml.safe_load(WORKFLOW.read_text())
    assert "'full'" in workflow["env"]["MODE"]
    assert "mode" in workflow[True]["workflow_dispatch"]["inputs"]
