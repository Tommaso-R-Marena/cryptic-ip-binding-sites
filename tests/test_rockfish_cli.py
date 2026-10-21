"""The rockfish CLI, including the three refusals that are the reason it exists."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("rockfish_cli", ROOT / "scripts" / "rockfish.py")
rockfish = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(rockfish)

SITE = "hpc/sites/rockfish.yml"
FLEXIBLE = "hpc/studies/flexible.yml"
COFOLD = "hpc/studies/cofold.yml"


def _site_file(tmp_path: Path, **extra) -> Path:
    body = {
        "name": "local",
        "scratch_root": str(tmp_path / "scratch"),
        "max_array_size": 1001,
        "max_array_throttle": 4,
        "partitions": {
            "shared": {"name": "shared", "max_walltime": "72:00:00", "max_cpus_per_node": 48,
                       "max_mem_gb_per_node": 187},
            "gpu": {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4,
                    "gpu_model": "A100", "gpu_vram_gb": 40},
        },
    }
    body.update(extra)
    import yaml
    path = tmp_path / "site.yml"
    path.write_text(yaml.safe_dump(body))
    return path


def _study_file(tmp_path: Path, **extra) -> Path:
    body = {
        "study": "T", "title": "toy", "plan": ["README.md"],
        "stages": [
            {"name": "work", "shards": 3, "output": "out-{shard}.json",
             "command": ["bash", "-c", 'printf "[]" > {stage_dir}/out-{shard}.json'],
             "resources": {"partition": "shared", "cpus": 1, "mem_gb": 1, "walltime": "00:05:00"}},
            {"name": "report", "depends_on": ["work"], "after": "any", "output": "r.json",
             "command": ["echo", "ok"],
             "resources": {"partition": "shared", "cpus": 1, "mem_gb": 1, "walltime": "00:05:00"}},
        ],
    }
    body.update(extra)
    import yaml
    path = tmp_path / "toy.yml"
    path.write_text(yaml.safe_dump(body))
    return path


# ------------------------------------------------------------------- validate

def test_validate_accepts_the_shipped_runnable_study():
    assert rockfish.main(["--site", SITE, "validate", "--study", FLEXIBLE]) == 0


def test_validate_reports_the_plan_gated_study_as_a_problem(capsys):
    code = rockfish.main(["--site", SITE, "validate", "--study", COFOLD])
    out = capsys.readouterr().out
    assert code == 1
    assert "missing pre-registered plan: docs/COFOLD_PLAN.md" in out


def test_validate_lists_the_stages_with_their_resources(capsys):
    rockfish.main(["--site", SITE, "validate", "--study", FLEXIBLE])
    out = capsys.readouterr().out
    assert "dock: 272 shards" in out
    assert "05:55:00" in out


def test_validate_warns_that_the_shipped_site_is_unverified(capsys):
    rockfish.main(["--site", SITE, "validate", "--study", FLEXIBLE])
    assert "never been checked against the live cluster" in capsys.readouterr().out


# --------------------------------------------------------------------- refusals

def test_submit_refuses_a_study_whose_plan_is_not_in_the_repository(capsys, tmp_path):
    site = _site_file(tmp_path, verified={"date": "2026-10-10", "by": "test"})
    code = rockfish.main(["--site", str(site), "submit", "--study", COFOLD, "--dry-run"])
    out = capsys.readouterr().out
    assert code == 2
    assert "pre-registered plan is not in the repository" in out
    assert "docs/COFOLD_PLAN.md" in out


def test_submit_refuses_an_unverified_site_unless_told_otherwise(capsys, tmp_path):
    site = _site_file(tmp_path)  # no verified block
    study = _study_file(tmp_path)
    code = rockfish.main(["--site", str(site), "submit", "--study", str(study), "--dry-run"])
    out = capsys.readouterr().out
    assert code == 2
    assert "never been checked against the live cluster" in out


def test_submit_proceeds_on_an_unverified_site_when_that_is_stated(capsys, tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    code = rockfish.main([
        "--site", str(site), "submit", "--study", str(study),
        "--dry-run", "--allow-unverified-site",
    ])
    assert code == 0
    assert "would run: sbatch" in capsys.readouterr().out


def test_a_dry_run_submits_nothing_and_writes_no_submission_record(tmp_path):
    site = _site_file(tmp_path, verified={"date": "2026-10-10", "by": "test"})
    study = _study_file(tmp_path)
    rockfish.main(["--site", str(site), "submit", "--study", str(study), "--dry-run"])
    assert not list((tmp_path / "scratch").rglob("submission.json"))


# ------------------------------------------------------------------------ plan

def test_plan_writes_one_script_per_stage_and_chunk(tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    out_dir = tmp_path / "scripts"
    assert rockfish.main([
        "--site", str(site), "plan", "--study", str(study), "--out-dir", str(out_dir),
    ]) == 0
    names = sorted(path.name for path in out_dir.glob("*.sbatch"))
    assert names == ["t-report.sbatch", "t-work.sbatch"]


def test_plan_blocks_and_says_why_when_the_plan_doc_is_missing(capsys, tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path, plan=["docs/NOT_WRITTEN.md"])
    code = rockfish.main(["--site", str(site), "plan", "--study", str(study)])
    assert code == 1
    assert "BLOCKED: missing pre-registered plan docs/NOT_WRITTEN.md" in capsys.readouterr().out


def test_the_dependency_is_wired_with_the_stages_own_after_kind(capsys, tmp_path):
    site = _site_file(tmp_path, verified={"date": "2026-10-10", "by": "test"})
    study = _study_file(tmp_path)
    rockfish.main(["--site", str(site), "submit", "--study", str(study), "--dry-run"])
    out = capsys.readouterr().out
    # The reporting stage must run even when shards failed, which is afterany.
    assert "--dependency=afterany:" in out


# ----------------------------------------------------------------- check-inputs

def test_check_inputs_passes_when_everything_is_staged(tmp_path, capsys):
    (tmp_path / "census").mkdir()
    (tmp_path / "census" / "census.csv").write_text("id\n")
    (tmp_path / "structures").mkdir()
    (tmp_path / "structures" / "1ZY7.pdb").write_text("ATOM\n")
    code = rockfish.main([
        "check-inputs", "--work-dir", str(tmp_path),
        "--require", "census/census.csv", "--require-dir", "structures",
        "--out", str(tmp_path / "inputs.json"),
    ])
    assert code == 0
    assert json.loads((tmp_path / "inputs.json").read_text())["ok"] is True
    assert "inputs ok" in capsys.readouterr().out


def test_check_inputs_fails_before_an_array_queues(tmp_path, capsys):
    code = rockfish.main([
        "check-inputs", "--work-dir", str(tmp_path), "--require", "census/census.csv",
    ])
    assert code == 1
    assert "missing file" in capsys.readouterr().out


def test_an_empty_input_directory_counts_as_missing(tmp_path, capsys):
    # A directory that exists but holds nothing is the shape a failed staging
    # step leaves behind, and 272 tasks would queue against it.
    (tmp_path / "structures").mkdir()
    code = rockfish.main(["check-inputs", "--work-dir", str(tmp_path), "--require-dir", "structures"])
    assert code == 1
    assert "empty directory" in capsys.readouterr().out


# --------------------------------------------------------------------- collect

def _run_toy_work(tmp_path: Path, shards=(0, 1, 2)):
    """Produce outputs for some shards of the toy study, by hand."""
    stage_dir = tmp_path / "scratch" / "t" / "work"
    stage_dir.mkdir(parents=True, exist_ok=True)
    for shard in shards:
        (stage_dir / f"out-{shard}.json").write_text(json.dumps([{"shard": shard}]))
    return stage_dir


def test_collect_merges_and_writes_the_census(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path)
    work_dir = tmp_path / "scratch" / "t"
    code = rockfish.main([
        "--site", str(site), "collect", "--study", str(study), "--stage", "work",
        "--work-dir", str(work_dir), "--records", str(work_dir / "records.json"),
        "--census-out", str(work_dir / "census.json"),
    ])
    assert code == 0
    assert json.loads((work_dir / "records.json").read_text()) == [
        {"shard": 0}, {"shard": 1}, {"shard": 2},
    ]
    assert json.loads((work_dir / "census.json").read_text())["reconciles"] is True
    assert "3 of 3 shards produced output" in capsys.readouterr().out


def test_collect_refuses_and_explains_when_a_shard_vanished(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path, shards=(0, 1))
    work_dir = tmp_path / "scratch" / "t"
    code = rockfish.main([
        "--site", str(site), "collect", "--study", str(study), "--stage", "work",
        "--work-dir", str(work_dir), "--records", str(work_dir / "records.json"),
    ])
    out = capsys.readouterr().out
    assert code == 1
    assert "neither output nor a status record" in out
    assert not (work_dir / "records.json").exists()


def test_collect_runs_the_follow_on_command_only_after_reconciling(tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path)
    work_dir = tmp_path / "scratch" / "t"
    marker = tmp_path / "report_ran.txt"
    code = rockfish.main([
        "--site", str(site), "collect", "--study", str(study), "--stage", "work",
        "--work-dir", str(work_dir), "--records", str(work_dir / "records.json"),
        "--then", f"touch {marker}",
    ])
    assert code == 0 and marker.exists()


def test_the_follow_on_command_does_not_run_when_collection_is_refused(tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path, shards=(0,))
    work_dir = tmp_path / "scratch" / "t"
    marker = tmp_path / "report_ran.txt"
    code = rockfish.main([
        "--site", str(site), "collect", "--study", str(study), "--stage", "work",
        "--work-dir", str(work_dir), "--then", f"touch {marker}",
    ])
    assert code == 1
    assert not marker.exists(), "a report must not run on an unreconciled stage"


def test_collect_writes_a_provenance_stub(tmp_path):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path)
    work_dir = tmp_path / "scratch" / "t"
    rockfish.main([
        "--site", str(site), "collect", "--study", str(study), "--stage", "work",
        "--work-dir", str(work_dir), "--provenance", str(work_dir / "PROVENANCE.md"),
    ])
    text = (work_dir / "PROVENANCE.md").read_text()
    assert "Shard census" in text and "README.md" in text


# ---------------------------------------------------------------------- status

def test_status_reports_the_census_per_stage(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path, shards=(0, 2))
    rockfish.main(["--site", str(site), "status", "--study", str(study), "--verbose"])
    out = capsys.readouterr().out
    assert "produced=2" in out and "missing=1" in out
    assert "missing: 1" in out


def test_status_on_an_untouched_work_dir_says_so(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    rockfish.main(["--site", str(site), "status", "--study", str(study)])
    assert "nothing there yet" in capsys.readouterr().out


# ---------------------------------------------------------------------- resume

def test_resume_targets_only_what_the_cluster_lost(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path, shards=(0, 2))
    code = rockfish.main([
        "--site", str(site), "resume", "--study", str(study), "--stage", "work", "--dry-run",
    ])
    out = capsys.readouterr().out
    assert code == 0
    assert "1 shard(s) to re-run: 1" in out


def test_resume_can_also_pick_up_shards_that_failed_in_the_code(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path, shards=(0,))
    status_dir = tmp_path / "scratch" / "t" / "work" / "_status"
    status_dir.mkdir(parents=True, exist_ok=True)
    (status_dir / "shard-1.json").write_text(json.dumps({"status": "failed", "exit_code": 7, "shard": 1}))
    rockfish.main([
        "--site", str(site), "resume", "--study", str(study), "--stage", "work",
        "--include-failed", "--dry-run",
    ])
    assert "2 shard(s) to re-run: 1, 2" in capsys.readouterr().out


def test_resume_says_nothing_to_do_when_every_shard_landed(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    _run_toy_work(tmp_path)
    assert rockfish.main([
        "--site", str(site), "resume", "--study", str(study), "--stage", "work", "--dry-run",
    ]) == 0
    assert "nothing to resume" in capsys.readouterr().out


def test_resuming_an_unsharded_stage_points_at_submit_instead(tmp_path, capsys):
    site = _site_file(tmp_path)
    study = _study_file(tmp_path)
    code = rockfish.main([
        "--site", str(site), "resume", "--study", str(study), "--stage", "report", "--dry-run",
    ])
    assert code == 1
    assert "not sharded" in capsys.readouterr().out


# ------------------------------------------------------------------- preflight

def test_preflight_on_a_machine_with_no_scheduler_does_not_claim_success(capsys):
    code = rockfish.main(["--site", SITE, "preflight"])
    out = capsys.readouterr().out
    assert code == 0
    assert "NOT VERIFIED" in out
    assert "login node" in out


def test_preflight_can_write_its_report(tmp_path):
    out_json = tmp_path / "preflight.json"
    out_md = tmp_path / "preflight.md"
    rockfish.main([
        "--site", SITE, "preflight", "--out", str(out_json), "--markdown", str(out_md),
    ])
    payload = json.loads(out_json.read_text())
    assert payload["site"] == "rockfish"
    assert payload["verdict"] == "indeterminate"
    assert "Preflight: rockfish" in out_md.read_text()


def test_an_unreadable_site_file_is_an_error_not_a_traceback(capsys, tmp_path):
    bad = tmp_path / "bad.yml"
    bad.write_text("name: x\n")  # no partitions, no scratch_root
    code = rockfish.main(["--site", str(bad), "validate"])
    assert code == 1
    assert "site file" in capsys.readouterr().out


@pytest.mark.parametrize("command", ["validate", "plan", "submit", "status", "collect", "resume"])
def test_every_command_is_registered(command):
    with pytest.raises(SystemExit):
        rockfish.main([command, "--help"])
