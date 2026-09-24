"""Structure of the proteome screen and aggregation workflows.

The screen takes hours of runner time, so its wiring is checked here rather
than discovered broken at the end of a run.
"""

from __future__ import annotations

from pathlib import Path

import yaml

WORKFLOWS = Path(__file__).resolve().parents[1] / ".github" / "workflows"


def _load(name):
    data = yaml.safe_load((WORKFLOWS / name).read_text())
    # PyYAML reads the bare key `on` as True.
    data["on"] = data.pop(True, data.get("on"))
    return data


def test_screen_aggregates_its_own_run_through_the_reusable_workflow():
    screen = _load("proteome-screen.yml")
    job = screen["jobs"]["aggregate"]
    assert job["uses"] == "./.github/workflows/proteome-aggregate.yml"
    assert job["with"]["run_id"] == "${{ github.run_id }}"
    assert job["permissions"]["actions"] == "read"
    assert "cancelled()" in job["if"]


def test_aggregation_reads_artifacts_of_the_requested_run():
    aggregate = _load("proteome-aggregate.yml")
    assert {"workflow_call", "workflow_dispatch"} <= set(aggregate["on"])
    downloads = [
        step for step in aggregate["jobs"]["aggregate"]["steps"]
        if str(step.get("uses", "")).startswith("actions/download-artifact")
    ]
    assert {step["with"]["pattern"] for step in downloads} == {"catalog-*", "shard-*"}
    for step in downloads:
        assert step["with"]["run-id"] == "${{ inputs.run_id }}"
        assert step["with"]["github-token"] == "${{ github.token }}"


def test_only_changes_to_what_is_measured_trigger_a_rescreen():
    """Scoring and hit statistics are applied at aggregation, not in the screen."""
    paths = set(_load("proteome-screen.yml")["on"]["push"]["paths"])
    assert "scripts/proteome_screen.py" in paths
    assert "cryptic_ip/analysis/proteome_stats.py" not in paths
    assert "cryptic_ip/analysis/scorer.py" not in paths


def test_screen_step_has_a_timeout_inside_the_job_timeout():
    job = _load("proteome-screen.yml")["jobs"]["screen"]
    step = next(s for s in job["steps"] if s.get("name") == "Screen shard")
    assert step["timeout-minutes"] < job["timeout-minutes"]
    upload = next(s for s in job["steps"] if s.get("name") == "Upload shard")
    assert upload["if"] == "always()"
