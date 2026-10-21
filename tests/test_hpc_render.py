"""The sbatch renderer, including running the scripts it writes.

The tests that matter most here execute the generated script and check what it
leaves behind. Study J lost two docking copies to reclaimed runners and no
record was written for either, so they were absent from the census rather than
present among the failures. A rendered task has to record its own death.
"""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import time
from pathlib import Path

import pytest

from cryptic_ip.hpc.manifest import StudyManifest
from cryptic_ip.hpc.render import expand_token, render_stage, render_study, write_scripts
from cryptic_ip.hpc.shards import chunk_shards
from cryptic_ip.hpc.site import SiteConfig

SITE = SiteConfig.from_dict({
    "name": "test",
    "scratch_root": "/scratch",
    "account": "tm1",
    "module_loads": ["anaconda3"],
    "env_activate": ["conda activate cryptic-ip"],
    "max_array_size": 1001,
    "max_array_throttle": 40,
    "partitions": {
        "shared": {"name": "shared", "max_walltime": "72:00:00", "max_cpus_per_node": 48,
                   "max_mem_gb_per_node": 192},
        "gpu": {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4,
                "gpu_model": "A100", "gpu_vram_gb": 40},
    },
    "verified": {"date": "2026-10-10", "by": "test"},
})


def _manifest(stage_overrides=None, **kwargs):
    stage = {
        "name": "work",
        "command": ["python", "x.py", "--out", "{stage_dir}/out.json"],
        "resources": {"partition": "shared", "cpus": 2, "mem_gb": 8, "walltime": "05:55:00"},
    }
    stage.update(stage_overrides or {})
    raw = {"study": "J", "title": "t", "plan": ["README.md"], "stages": [stage]}
    raw.update(kwargs)
    return StudyManifest.from_dict(raw)


def _render(stage_overrides=None, shards=None, **kwargs):
    manifest = _manifest(stage_overrides, **kwargs)
    stage = manifest.stage("work")
    chunk = None
    if stage.is_array:
        chunk, = chunk_shards(shards or range(stage.task_count), max_array_size=1001, throttle=40)
    return render_stage(manifest, stage, SITE, work_dir="/wd", repo_root="/repo", chunk=chunk)


# ---------------------------------------------------------- token expansion

def test_a_placeholder_becomes_a_shell_variable():
    assert expand_token("{stage_dir}/out.json", {"stage_dir": "${STAGE_DIR}"}) == '"${STAGE_DIR}/out.json"'


def test_a_literal_value_is_substituted_and_escaped():
    assert expand_token("--study={study}", {"study": "J"}) == '"--study=J"'


def test_a_brace_that_is_not_a_placeholder_stays_literal():
    # A JSON literal in a command line contains braces. A hand-rolled scanner
    # read `{"shard": %s}` as a placeholder named `"shard": %s` and crashed at
    # submit time, although the manifest had loaded cleanly.
    out = expand_token('printf \'[{"shard": %s}]\' {shard}', {"shard": "${SHARD}"})
    assert out == '"printf \'[{\\"shard\\": %s}]\' ${SHARD}"'


@pytest.mark.parametrize("hostile", ['a"b', "a$b", "a`b", "a\\b"])
def test_shell_metacharacters_in_a_literal_are_escaped(hostile):
    rendered = expand_token(hostile, {})
    # Round-trip it through bash: what comes out must equal what went in.
    result = subprocess.run(["bash", "-c", f"printf '%s' {rendered}"], capture_output=True, text=True)
    assert result.stdout == hostile


def test_an_unknown_placeholder_at_render_time_is_an_error():
    with pytest.raises(KeyError):
        expand_token("{nope}", {"shard": "${SHARD}"})


# ------------------------------------------------------------------ directives

def test_the_directives_carry_the_resources_the_stage_asked_for():
    script = _render()
    assert "#SBATCH --partition=shared" in script
    assert "#SBATCH --time=05:55:00" in script
    assert "#SBATCH --cpus-per-task=2" in script
    assert "#SBATCH --mem=8G" in script
    assert "#SBATCH --account=tm1" in script


def test_no_gres_line_unless_the_stage_wants_a_gpu():
    assert "--gres" not in _render()


def test_a_gpu_stage_gets_a_gres_line_naming_the_partitions_gres():
    script = _render({"resources": {"partition": "gpu", "cpus": 4, "mem_gb": 32,
                                    "walltime": "02:00:00", "gpus": 1}})
    assert "#SBATCH --gres=gpu:1" in script


def test_an_array_stage_logs_per_task_and_a_single_stage_per_job():
    array = _render({"shards": 4, "output": "out-{shard}.json"})
    # The throttle is clamped to the task count, so a 4-task array reads %4.
    assert "#SBATCH --array=0-3%4" in array
    assert "work-%A_%a.out" in array
    assert "work-%j.out" in _render()


def test_the_task_is_requeueable_and_warned_before_the_walltime_kill():
    # The USR1 warning is what lets the trap write a record while the task is
    # still alive, instead of the task simply disappearing.
    script = _render()
    assert "#SBATCH --requeue" in script
    assert "#SBATCH --signal=B:USR1@120" in script


# --------------------------------------------------------------- shard mapping

def test_a_contiguous_chunk_uses_arithmetic_not_a_lookup_table():
    script = _render({"shards": 272, "output": "o-{shard}.json"})
    assert "SHARD=$(( SLURM_ARRAY_TASK_ID + 0 ))" in script
    assert "SHARD_MAP" not in script


def test_a_sparse_chunk_carries_an_explicit_map():
    script = _render({"shards": 272, "output": "o-{shard}.json"}, shards=[117, 245])
    assert "SHARD_MAP=(117 245)" in script
    assert "#SBATCH --array=0-1%2" in script


def test_an_array_script_run_without_a_task_id_says_so():
    script = _render({"shards": 4, "output": "o-{shard}.json"})
    assert "SLURM_ARRAY_TASK_ID" in script
    assert "submit it with sbatch" in script


# ------------------------------------------------------------- the environment

def test_modules_and_activation_appear_before_the_command():
    script = _render()
    lines = script.splitlines()
    assert lines.index("module load anaconda3") < lines.index('cd "$REPO_ROOT"')
    assert "conda activate cryptic-ip" in script


def test_a_stage_that_pins_its_thread_count_is_not_overridden():
    # Study J pins OMP_NUM_THREADS=1. An unconditional default written after it
    # silently raised it to the core count, which changes what Vina does.
    script = _render({"env": {"OMP_NUM_THREADS": "1"}})
    exports = [line for line in script.splitlines() if "OMP_NUM_THREADS" in line]
    assert exports == ['export OMP_NUM_THREADS="1"']


def test_without_a_pin_the_thread_count_follows_the_allocation():
    script = _render()
    assert 'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"' in script


# ------------------------------------------------- what the scripts actually do

def _runnable(tmp_path: Path, body: str, shards=(0, 1, 2)) -> Path:
    """Render a script whose command we can run here, with modules stubbed out."""
    manifest = _manifest({
        "name": "work", "shards": len(shards), "output": "out-{shard}.txt",
        "command": ["bash", "-c", body],
        "resources": {"partition": "shared", "cpus": 1, "mem_gb": 1, "walltime": "00:05:00"},
    })
    chunk, = chunk_shards(shards, max_array_size=1001, throttle=2)
    text = render_stage(
        manifest, manifest.stage("work"), SITE,
        work_dir=str(tmp_path / "wd"), repo_root=str(tmp_path), chunk=chunk,
    )
    text = re.sub(r"^(module |conda activate|source ).*$", "true", text, flags=re.MULTILINE)
    path = tmp_path / "work.sbatch"
    path.write_text(text)
    return path


def test_the_rendered_script_is_valid_bash(tmp_path):
    subprocess.run(["bash", "-n", str(_runnable(tmp_path, "true"))], check=True)


def _status(tmp_path: Path, shard: int):
    return json.loads((tmp_path / "wd" / "work" / "_status" / f"shard-{shard}.json").read_text())


def test_a_task_that_succeeds_records_ok_and_its_shard(tmp_path):
    script = _runnable(tmp_path, 'echo hi > "{stage_dir}/out-{shard}.txt"')
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "2"}
    subprocess.run(["bash", str(script)], check=True, env=env)
    record = _status(tmp_path, 2)
    assert record["status"] == "ok" and record["shard"] == 2 and record["exit_code"] == 0
    assert (tmp_path / "wd" / "work" / "out-2.txt").is_file()


def test_a_task_that_fails_records_the_exit_code(tmp_path):
    script = _runnable(tmp_path, "exit 7")
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "1"}
    result = subprocess.run(["bash", str(script)], env=env)
    assert result.returncode == 7
    record = _status(tmp_path, 1)
    assert record["status"] == "failed" and record["exit_code"] == 7


def test_a_task_the_scheduler_kills_records_that_it_was_killed(tmp_path):
    # This is the study J case: two copies were lost when their runners were
    # reclaimed, nothing was written, and they vanished from the census instead
    # of appearing among the failures.
    script = _runnable(tmp_path, "sleep 30")
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "0"}
    process = subprocess.Popen(["bash", str(script)], env=env)
    deadline = time.time() + 10
    status_file = tmp_path / "wd" / "work" / "_status" / "shard-0.json"
    while time.time() < deadline and not status_file.is_file():
        time.sleep(0.05)
    process.send_signal(signal.SIGTERM)
    process.wait(timeout=10)
    record = _status(tmp_path, 0)
    assert record["status"] == "killed", "a killed task must not be indistinguishable from one that never ran"
    assert record["exit_code"] == 143


def test_a_running_task_is_visible_as_running_before_it_finishes(tmp_path):
    script = _runnable(tmp_path, "sleep 5")
    env = {**os.environ, "SLURM_ARRAY_TASK_ID": "0"}
    process = subprocess.Popen(["bash", str(script)], env=env)
    status_file = tmp_path / "wd" / "work" / "_status" / "shard-0.json"
    deadline = time.time() + 10
    while time.time() < deadline and not status_file.is_file():
        time.sleep(0.05)
    assert json.loads(status_file.read_text())["status"] == "running"
    process.kill()
    process.wait(timeout=10)


# ----------------------------------------------------------------- whole study

def test_render_study_returns_one_script_per_chunk_in_dependency_order():
    manifest = StudyManifest.from_dict({
        "study": "J", "title": "t", "plan": ["README.md"],
        "stages": [
            {"name": "report", "depends_on": ["dock"], "after": "any",
             "command": ["echo", "done"],
             "resources": {"partition": "shared", "walltime": "01:00:00"}},
            {"name": "dock", "shards": 5, "output": "o-{shard}.json",
             "command": ["echo", "{shard}"],
             "resources": {"partition": "shared", "walltime": "01:00:00"}},
        ],
    })
    small = SiteConfig.from_dict({
        "name": "s", "scratch_root": "/scratch", "max_array_size": 3,
        "partitions": {"shared": {"name": "shared", "max_walltime": "72:00:00"}},
        "verified": {"date": "2026-10-10", "by": "t"},
    })
    rendered = render_study(manifest, small, work_dir="/wd", repo_root="/repo")
    # max_array_size 3 means at most 2 tasks a chunk, so 5 shards need 3 chunks.
    assert [item.script_name for item in rendered] == [
        "j-dock.sbatch", "j-dock-c1.sbatch", "j-dock-c2.sbatch", "j-report.sbatch",
    ]
    assert [len(item.shards) for item in rendered] == [2, 2, 1, 1]
    assert rendered[-1].after == "any" and rendered[-1].depends_on == ("dock",)


def test_only_shards_restricts_a_resume_to_the_shards_that_did_not_land():
    manifest = _manifest({"shards": 272, "output": "o-{shard}.json"})
    rendered = render_study(
        manifest, SITE, work_dir="/wd", repo_root="/repo", only_shards={"work": [117, 245]},
    )
    assert len(rendered) == 1
    assert rendered[0].shards == (117, 245)


def test_a_stage_with_no_shards_to_rerun_is_skipped_entirely():
    manifest = _manifest({"shards": 4, "output": "o-{shard}.json"})
    assert render_study(manifest, SITE, work_dir="/wd", repo_root="/repo", only_shards={"work": []}) == []


def test_written_scripts_are_executable(tmp_path):
    manifest = _manifest()
    rendered = render_study(manifest, SITE, work_dir="/wd", repo_root="/repo")
    paths = write_scripts(rendered, tmp_path / "scripts")
    assert paths and all(os.access(path, os.X_OK) for path in paths)


def test_rendering_a_sharded_stage_without_a_chunk_is_refused():
    manifest = _manifest({"shards": 4, "output": "o-{shard}.json"})
    with pytest.raises(ValueError, match="one chunk at a time"):
        render_stage(manifest, manifest.stage("work"), SITE, work_dir="/wd", repo_root="/repo")
