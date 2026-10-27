"""The manifest layer: stages, resources, dependency order, and the plan gate."""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

from cryptic_ip.hpc.manifest import (
    ManifestError, Resources, Stage, StudyManifest, load_manifest, total_tasks,
)
from cryptic_ip.hpc.site import Partition, load_site

STUDIES = Path("hpc/studies")
WORKFLOW = Path(".github/workflows/flexible.yml")


def _stage(**overrides):
    raw = {
        "name": "work",
        "command": ["python", "scripts/flexible.py", "dock"],
        "resources": {"partition": "shared", "cpus": 2, "mem_gb": 8, "walltime": "01:00:00"},
    }
    raw.update(overrides)
    return raw


def _manifest(**overrides):
    raw = {"study": "X", "title": "t", "plan": ["README.md"], "stages": [_stage()]}
    raw.update(overrides)
    return raw


# ------------------------------------------------------------------- resources

def test_resources_need_a_partition_role():
    with pytest.raises(ManifestError, match="partition"):
        Resources.from_dict({"cpus": 1})


@pytest.mark.parametrize("field,value,word", [("cpus", 0, "at least 1"), ("mem_gb", 0, "positive")])
def test_resources_reject_explicit_zeros(field, value, word):
    with pytest.raises(ManifestError, match=word):
        Resources.from_dict({"partition": "shared", field: value})


def test_resources_reject_negative_gpus():
    with pytest.raises(ManifestError, match="negative"):
        Resources.from_dict({"partition": "gpu", "gpus": -1})


def test_a_stage_over_the_walltime_cap_is_reported_not_silently_clipped():
    partition = Partition.from_dict("cpu", {"name": "defq", "max_walltime": "02:00:00"})
    problems = Resources.from_dict({"partition": "cpu", "walltime": "05:55:00"}).problems_against(partition)
    assert len(problems) == 1
    assert "05:55:00" in problems[0] and "02:00:00" in problems[0]


def test_asking_a_cpu_partition_for_a_gpu_is_caught_before_submission():
    partition = Partition.from_dict("cpu", {"name": "defq", "max_walltime": "72:00:00"})
    problems = Resources.from_dict({"partition": "cpu", "gpus": 1}).problems_against(partition)
    assert any("has none" in problem for problem in problems)


def test_more_gpus_than_a_node_holds_is_caught():
    partition = Partition.from_dict(
        "gpu", {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4, "gpu_vram_gb": 40},
    )
    problems = Resources.from_dict({"partition": "gpu", "gpus": 8}).problems_against(partition)
    assert any("per node" in problem for problem in problems)


def test_cpu_and_memory_over_a_node_are_both_reported():
    partition = Partition.from_dict(
        "cpu",
        {"name": "defq", "max_walltime": "72:00:00", "max_cpus_per_node": 48, "max_mem_gb_per_node": 192},
    )
    problems = Resources.from_dict(
        {"partition": "cpu", "cpus": 96, "mem_gb": 500, "walltime": "01:00:00"},
    ).problems_against(partition)
    assert len(problems) == 2


# ----------------------------------------------------------------------- stage

def test_a_command_must_be_tokens_not_one_string():
    # One string would be passed to the shell, and quoting would decide what ran.
    with pytest.raises(ManifestError, match="list of tokens"):
        Stage.from_dict(_stage(command="python scripts/flexible.py dock"))


def test_an_empty_command_is_refused():
    with pytest.raises(ManifestError, match="command is empty"):
        Stage.from_dict(_stage(command=[]))


def test_stage_names_are_constrained_so_they_can_be_job_names():
    with pytest.raises(ManifestError, match="lowercase"):
        Stage.from_dict(_stage(name="Dock Stage"))


def test_an_unknown_placeholder_is_a_typo_and_is_refused():
    with pytest.raises(ManifestError, match="unknown placeholder"):
        Stage.from_dict(_stage(command=["python", "x.py", "--out", "{outdir}/x.json"]))


def test_shard_placeholder_without_shards_is_refused():
    with pytest.raises(ManifestError, match=re.escape("uses {shard} but declares no shards")):
        Stage.from_dict(_stage(command=["python", "x.py", "--shard", "{shard}"]))


def test_a_sharded_stage_whose_output_does_not_vary_would_overwrite_itself():
    with pytest.raises(ManifestError, match="write over the others"):
        Stage.from_dict(_stage(shards=4, output="out.json"))


@pytest.mark.parametrize("bad", ["/abs/out.json", "../out-{shard}.json", "${STAGE_DIR}/out-{shard}.json"])
def test_outputs_stay_inside_the_stage_directory(bad):
    with pytest.raises(ManifestError):
        Stage.from_dict(_stage(shards=2, output=bad))


def test_an_output_may_only_vary_by_shard():
    with pytest.raises(ManifestError, match=re.escape("only use {shard}")):
        Stage.from_dict(_stage(shards=2, output="{study}-{shard}.json"))


def test_shards_and_throttle_must_be_positive():
    with pytest.raises(ManifestError, match="shards must be at least 1"):
        Stage.from_dict(_stage(shards=0))
    with pytest.raises(ManifestError, match="throttle must be at least 1"):
        Stage.from_dict(_stage(shards=4, output="o-{shard}", throttle=0))


def test_after_must_be_a_kind_slurm_understands():
    with pytest.raises(ManifestError, match="after must be one of"):
        Stage.from_dict(_stage(depends_on=["other"], after="whenever"))


def test_after_without_a_dependency_means_nothing_and_is_refused():
    with pytest.raises(ManifestError, match="means nothing without depends_on"):
        Stage.from_dict(_stage(after="any"))


def test_output_for_resolves_the_shard():
    stage = Stage.from_dict(_stage(shards=3, output="flexible-{shard}.json"))
    assert stage.output_for(2) == "flexible-2.json"
    assert stage.task_count == 3 and stage.is_array


def test_a_single_task_stage_counts_as_one():
    stage = Stage.from_dict(_stage())
    assert not stage.is_array and stage.task_count == 1


# -------------------------------------------------------------------- manifest

def test_duplicate_stage_names_are_refused():
    with pytest.raises(ManifestError, match="duplicate stage"):
        StudyManifest.from_dict(_manifest(stages=[_stage(), _stage()]))


def test_a_dependency_on_an_unknown_stage_is_refused():
    with pytest.raises(ManifestError, match="unknown stage"):
        StudyManifest.from_dict(_manifest(stages=[_stage(depends_on=["ghost"])]))


def test_a_stage_cannot_depend_on_itself():
    with pytest.raises(ManifestError, match="depends on itself"):
        StudyManifest.from_dict(_manifest(stages=[_stage(depends_on=["work"])]))


def test_a_dependency_cycle_is_refused_rather_than_deadlocking_the_queue():
    stages = [
        _stage(name="a", depends_on=["b"]),
        _stage(name="b", depends_on=["c"]),
        _stage(name="c", depends_on=["a"]),
    ]
    with pytest.raises(ManifestError, match="cycle"):
        StudyManifest.from_dict(_manifest(stages=stages))


def test_order_puts_every_dependency_first():
    stages = [
        _stage(name="report", depends_on=["dock"]),
        _stage(name="dock", depends_on=["inputs"]),
        _stage(name="inputs"),
    ]
    manifest = StudyManifest.from_dict(_manifest(stages=stages))
    assert [s.name for s in manifest.order()] == ["inputs", "dock", "report"]


def test_a_manifest_with_no_plan_is_refused():
    # An empty plan list is exactly how an un-pre-registered study would get run.
    with pytest.raises(ManifestError, match="pre-registered plan"):
        StudyManifest.from_dict(_manifest(plan=[]))


def test_missing_plan_docs_are_named(tmp_path):
    manifest = StudyManifest.from_dict(_manifest(plan=["docs/THERE.md", "docs/NOT_THERE.md"]))
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "THERE.md").write_text("x")
    assert manifest.missing_plan_docs(tmp_path) == ["docs/NOT_THERE.md"]


def test_total_tasks_counts_array_tasks_not_stages():
    stages = [_stage(name="a"), _stage(name="b", shards=10, output="o-{shard}")]
    assert total_tasks(StudyManifest.from_dict(_manifest(stages=stages))) == 11


def test_problems_against_a_site_name_the_stage():
    site = load_site("hpc/sites/rockfish.yml")
    stages = [_stage(resources={"partition": "nosuch", "cpus": 1, "mem_gb": 1, "walltime": "01:00:00"})]
    problems = StudyManifest.from_dict(_manifest(stages=stages)).problems_against(site)
    assert problems and "work" in problems[0]


# ------------------------------------------------------- the shipped manifests

def test_every_shipped_manifest_loads_and_fits_the_site():
    site = load_site("hpc/sites/rockfish.yml")
    paths = sorted(STUDIES.glob("*.yml"))
    assert paths, "no study manifests found"
    for path in paths:
        manifest = load_manifest(path)
        assert manifest.problems_against(site) == [], f"{path.name} does not fit the site"


def test_the_gpu_studies_are_gated_on_plans_that_do_not_exist_yet():
    # Study M (co-folding) and study P (MD) are engineered but not
    # pre-registered. docs/GPU_REQUIREMENTS.md commits to pre-registering study
    # M before any candidate is scored, so the gate has to be holding them.
    for name in ("cofold.yml", "md.yml"):
        manifest = load_manifest(STUDIES / name)
        assert manifest.missing_plan_docs(Path(".")), f"{name} should still be plan-gated"


def test_study_j_is_runnable_because_its_plan_and_amendments_are_committed():
    manifest = load_manifest(STUDIES / "flexible.yml")
    assert manifest.missing_plan_docs(Path(".")) == []
    assert len(manifest.plan) == 4  # the plan plus three amendments


def test_the_study_j_port_shards_exactly_as_the_workflow_it_came_from():
    # The port has to be the same study, not a similar one. The workflow's
    # DOCK_SHARDS is the reference; one copy per shard is amendment 3.
    workflow = yaml.safe_load(WORKFLOW.read_text())
    expected = int(workflow["env"]["DOCK_SHARDS"])
    dock = load_manifest(STUDIES / "flexible.yml").stage("dock")
    assert dock.task_count == expected == 272


def test_the_study_j_dock_walltime_leaves_room_for_its_own_shard_budget():
    # The script records its own censoring when it runs out of budget. If the
    # scheduler kills the task first, that record is never written - which is
    # the difference between a reported timeout and a vanished shard.
    dock = load_manifest(STUDIES / "flexible.yml").stage("dock")
    budget = int(dock.command[dock.command.index("--shard-budget") + 1])
    assert dock.resources.walltime > budget
