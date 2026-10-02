"""Every workflow parses, and every job id is a string.

YAML reads bare ``null``, ``true``, ``on`` and the like as non-strings; a job
named ``null`` makes GitHub reject the whole file at push time.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

WORKFLOWS = sorted((Path(__file__).resolve().parents[1] / ".github" / "workflows").glob("*.yml"))


@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_workflow_jobs_are_named_by_strings(path):
    data = yaml.safe_load(path.read_text())
    jobs = data["jobs"]
    assert jobs and all(isinstance(name, str) for name in jobs), list(jobs)
    for name, job in jobs.items():
        needs = job.get("needs", [])
        for dependency in [needs] if isinstance(needs, str) else needs:
            assert dependency in jobs, f"{path.name}: {name} needs unknown job {dependency!r}"


@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_no_step_input_parses_to_null(path):
    """``path: null`` is YAML null, not a directory called null: the input is silently dropped."""
    data = yaml.safe_load(path.read_text())
    for name, job in data["jobs"].items():
        for step in job.get("steps", []):
            for key, value in (step.get("with") or {}).items():
                assert value is not None, f"{path.name}: job {name!r} step input {key!r} is YAML null"


#: Workflows that run analyses (hours of compute) rather than tests.
ANALYSES = {
    "benchmark.yml", "burial-survey.yml", "explore.yml", "permutation-null.yml", "phase1-criteria.yml",
    "proteome-screen.yml", "train-real-data.yml", "yeast-rescreen.yml", "transfer.yml",
    "transfer-secondary.yml", "learned-screen.yml", "redocking.yml", "arrestin.yml", "specificity.yml",
    "hull-gate.yml", "redocking-audit.yml", "rerank.yml", "sampling.yml", "triage.yml",
    "orthogonal.yml", "coevolution.yml", "kinase.yml", "cofold-probe.yml",
    "template-fit.yml", "template2.yml", "flexible.yml",
}

#: Workflows written under the conventions below from the start; older ones are
#: not edited, because touching their files would re-run their analyses.
STRICT = {"redocking.yml", "arrestin.yml", "specificity.yml", "hull-gate.yml", "redocking-audit.yml",
          "rerank.yml", "sampling.yml", "triage.yml", "orthogonal.yml", "coevolution.yml", "kinase.yml",
          "cofold-probe.yml", "template-fit.yml", "template2.yml", "flexible.yml"}


@pytest.mark.parametrize("path", [p for p in WORKFLOWS if p.name in ANALYSES], ids=lambda p: p.name)
def test_analysis_workflows_never_auto_run_on_main(path):
    data = yaml.safe_load(path.read_text())
    on = data.get(True, data.get("on"))
    push = on.get("push") if isinstance(on, dict) else None
    assert push is None or "main" in push.get("branches-ignore", []), path.name
    assert "workflow_dispatch" in on, f"{path.name} must stay runnable on demand"


@pytest.mark.parametrize("path", [p for p in WORKFLOWS if p.name in STRICT], ids=lambda p: p.name)
def test_conda_jobs_run_every_step_in_the_conda_shell(path):
    """A step without ``bash -l {0}`` in a setup-miniconda job runs outside the environment."""
    data = yaml.safe_load(path.read_text())
    for name, job in data["jobs"].items():
        steps = job.get("steps", [])
        if not any("setup-miniconda" in str(step.get("uses", "")) for step in steps):
            continue
        for step in steps:
            if "run" in step:
                assert step.get("shell") == "bash -l {0}", f"{path.name}: {name}: {step.get('name')}"


@pytest.mark.parametrize("path", [p for p in WORKFLOWS if p.name in STRICT], ids=lambda p: p.name)
def test_analysis_workflows_trigger_only_on_their_own_files(path):
    data = yaml.safe_load(path.read_text())
    on = data.get(True, data.get("on"))
    paths = on["push"]["paths"]
    assert f".github/workflows/{path.name}" in paths
    assert all(not p.startswith("tests/") and p not in ("requirements.txt", "setup.py") for p in paths)


@pytest.mark.parametrize("path", [p for p in WORKFLOWS if p.name in STRICT], ids=lambda p: p.name)
def test_results_are_printed_between_markers(path):
    text = path.read_text()
    assert "extract_log_block.py emit" in text, f"{path.name} never prints its results for extraction"


@pytest.mark.parametrize("path", [p for p in WORKFLOWS if p.name in STRICT], ids=lambda p: p.name)
def test_matrices_stay_small(path):
    data = yaml.safe_load(path.read_text())
    for name, job in data["jobs"].items():
        matrix = (job.get("strategy") or {}).get("matrix") or {}
        size = 1
        for key, values in matrix.items():
            if isinstance(values, list):
                size *= len(values)
        assert size <= 30, f"{path.name}: {name} has {size} matrix entries"


#: GitHub generates at most this many jobs from one matrix. Past it the job is dropped
#: without failing the run, which is how study J's 272-shard matrix produced a green run
#: with no docking in it (run 36883632819).
MATRIX_JOB_LIMIT = 256


def _flexible():
    path = Path(__file__).resolve().parents[1] / ".github" / "workflows" / "flexible.yml"
    return path, yaml.safe_load(path.read_text())


def test_the_shard_halves_tile_every_shard_exactly_once():
    """The split must cover range(DOCK_SHARDS) with no gap and no overlap.

    A gap silently drops copies from the study; an overlap docks them twice.
    """
    _path, data = _flexible()
    n = int(data["env"]["DOCK_SHARDS"])
    half = (n + 1) // 2
    first, second = list(range(half)), list(range(half, n))
    assert first + second == list(range(n))
    assert not set(first) & set(second)


def test_neither_shard_matrix_can_reach_the_matrix_job_limit():
    _path, data = _flexible()
    n = int(data["env"]["DOCK_SHARDS"])
    half = (n + 1) // 2
    assert half <= MATRIX_JOB_LIMIT, f"first half is {half} jobs"
    assert n - half <= MATRIX_JOB_LIMIT, f"second half is {n - half} jobs"


def test_both_shard_matrices_are_wired_to_their_own_half():
    _path, data = _flexible()
    jobs = data["jobs"]
    for job, output in (("dock_a", "dock_a"), ("dock_b", "dock_b")):
        matrix = jobs[job]["strategy"]["matrix"]["shard"]
        assert f"needs.plan.outputs.{output}" in matrix, f"{job} reads {matrix}"
    assert set(jobs["report"]["needs"]) == {"dock_a", "dock_b"}


def test_the_report_fails_when_no_shard_produced_a_file():
    """A dropped matrix must not look like a run whose copies all failed."""
    _path, data = _flexible()
    merge = [s for s in data["jobs"]["report"]["steps"] if s.get("name") == "Merge the shards"]
    assert merge, "the report no longer has a merge step"
    assert "no shard produced an output file" in merge[0]["run"]
