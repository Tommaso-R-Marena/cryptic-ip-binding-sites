"""Collection as reconciliation: a shard that vanished must not shrink the data.

Study J merged 270 records from a 272-shard run. The two absent copies had lost
their runners, so nothing was written for them and they appeared nowhere in the
failure census - the only trace was the count. These tests pin that the same
situation now stops the merge instead of producing a quietly shorter dataset.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from cryptic_ip.hpc.collect import (
    FAILED, KILLED, MISSING, PRODUCED, RUNNING, CollectError, collect_stage, merge_records,
    provenance, read_status, take_census,
)
from cryptic_ip.hpc.manifest import StudyManifest
from cryptic_ip.hpc.render import STATUS_DIR


def _manifest(shards=6, output="out-{shard}.json"):
    return StudyManifest.from_dict({
        "study": "J", "title": "t", "plan": ["README.md"],
        "stages": [{
            "name": "dock", "shards": shards, "output": output,
            "command": ["echo", "{shard}"],
            "resources": {"partition": "shared", "walltime": "01:00:00"},
        }],
    })


def _write_output(stage_dir: Path, shard: int, payload=None):
    stage_dir.mkdir(parents=True, exist_ok=True)
    body = payload if payload is not None else [{"shard": shard}]
    (stage_dir / f"out-{shard}.json").write_text(json.dumps(body))


def _write_status(stage_dir: Path, shard: int, status: str, exit_code: int = 0):
    directory = stage_dir / STATUS_DIR
    directory.mkdir(parents=True, exist_ok=True)
    (directory / f"shard-{shard}.json").write_text(
        json.dumps({"study": "J", "stage": "dock", "shard": shard, "status": status, "exit_code": exit_code})
    )


def test_every_shard_lands_in_exactly_one_state(tmp_path):
    stage_dir = tmp_path / "dock"
    manifest = _manifest(shards=5)
    _write_output(stage_dir, 0)
    _write_status(stage_dir, 0, "ok")
    _write_status(stage_dir, 1, "failed", 7)
    _write_status(stage_dir, 2, "killed", 143)
    _write_status(stage_dir, 3, "running", -1)
    # shard 4 leaves nothing at all
    census = take_census(manifest, manifest.stage("dock"), stage_dir)
    assert census.states == {0: PRODUCED, 1: FAILED, 2: KILLED, 3: RUNNING, 4: MISSING}
    assert census.counts() == {PRODUCED: 1, FAILED: 1, KILLED: 1, RUNNING: 1, MISSING: 1}
    assert not census.reconciles


def test_a_shard_with_output_counts_as_produced_whatever_its_record_says(tmp_path):
    # The output is the evidence. A stale "running" record next to a finished
    # file should not remove a copy from the analysis.
    stage_dir = tmp_path / "dock"
    manifest = _manifest(shards=1)
    _write_output(stage_dir, 0)
    _write_status(stage_dir, 0, "running", -1)
    assert take_census(manifest, manifest.stage("dock"), stage_dir).states == {0: PRODUCED}


def test_a_task_claiming_success_with_no_output_is_a_failure_not_an_absence(tmp_path):
    stage_dir = tmp_path / "dock"
    manifest = _manifest(shards=1)
    _write_status(stage_dir, 0, "ok")
    assert take_census(manifest, manifest.stage("dock"), stage_dir).states == {0: FAILED}


def test_an_unreadable_status_file_is_a_failure_not_a_crash(tmp_path):
    stage_dir = tmp_path / "dock"
    (stage_dir / STATUS_DIR).mkdir(parents=True)
    (stage_dir / STATUS_DIR / "shard-0.json").write_text("{not json")
    manifest = _manifest(shards=1)
    assert take_census(manifest, manifest.stage("dock"), stage_dir).states == {0: FAILED}
    assert read_status(stage_dir, 0)["status"] == "unreadable"


def test_records_from_every_produced_shard_are_merged_in_shard_order(tmp_path):
    stage_dir = tmp_path / "dock"
    manifest = _manifest(shards=3)
    for shard in range(3):
        _write_output(stage_dir, shard, [{"copy": f"c{shard}"}])
    census = take_census(manifest, manifest.stage("dock"), stage_dir)
    assert merge_records(census) == [{"copy": "c0"}, {"copy": "c1"}, {"copy": "c2"}]
    assert census.records_merged == 3


def test_a_shard_writing_one_object_rather_than_a_list_is_kept(tmp_path):
    stage_dir = tmp_path / "dock"
    manifest = _manifest(shards=1)
    _write_output(stage_dir, 0, {"copy": "single"})
    census = take_census(manifest, manifest.stage("dock"), stage_dir)
    assert merge_records(census) == [{"copy": "single"}]


def test_unreadable_output_is_reported_with_the_shard_that_wrote_it(tmp_path):
    stage_dir = tmp_path / "dock"
    stage_dir.mkdir()
    (stage_dir / "out-0.json").write_text("{truncated")
    manifest = _manifest(shards=1)
    census = take_census(manifest, manifest.stage("dock"), stage_dir)
    with pytest.raises(CollectError, match="shard 0"):
        merge_records(census)


def test_collect_refuses_when_a_shard_left_nothing(tmp_path):
    manifest = _manifest(shards=3)
    stage_dir = tmp_path / "dock"
    _write_output(stage_dir, 0)
    _write_output(stage_dir, 1)
    # shard 2 vanished
    with pytest.raises(CollectError, match="neither output nor a status record"):
        collect_stage(manifest, "dock", tmp_path, records_out=tmp_path / "records.json")
    assert not (tmp_path / "records.json").exists(), "a refused collection must not leave a merged file"


def test_the_census_is_written_even_when_the_merge_is_refused(tmp_path):
    # The census is the diagnosis. Withholding it on failure would leave the
    # person with a refusal and no list of which shards to re-run.
    manifest = _manifest(shards=2)
    _write_output(tmp_path / "dock", 0)
    with pytest.raises(CollectError):
        collect_stage(manifest, "dock", tmp_path, census_out=tmp_path / "census.json")
    census = json.loads((tmp_path / "census.json").read_text())
    assert census["missing"] == [1] and census["reconciles"] is False


def test_allow_missing_makes_the_shortfall_a_stated_decision(tmp_path):
    manifest = _manifest(shards=3)
    stage_dir = tmp_path / "dock"
    _write_output(stage_dir, 0)
    _write_output(stage_dir, 1)
    census = collect_stage(
        manifest, "dock", tmp_path, records_out=tmp_path / "records.json", allow_missing=1,
    )
    assert census.missing == [2]
    assert json.loads((tmp_path / "records.json").read_text()) == [{"shard": 0}, {"shard": 1}]


def test_allow_missing_does_not_cover_more_than_it_says(tmp_path):
    manifest = _manifest(shards=4)
    _write_output(tmp_path / "dock", 0)
    with pytest.raises(CollectError, match="allowing 1"):
        collect_stage(manifest, "dock", tmp_path, allow_missing=1)


def test_a_stage_where_nothing_produced_output_is_a_structural_failure(tmp_path):
    # Study J's report once ran on zero shards because the matrix had been
    # dropped, and reported "not evaluable" rather than failing. A stage where
    # every task ran and every copy failed still writes files.
    manifest = _manifest(shards=3)
    stage_dir = tmp_path / "dock"
    for shard in range(3):
        _write_status(stage_dir, shard, "failed", 1)
    with pytest.raises(CollectError, match="structural failure"):
        collect_stage(manifest, "dock", tmp_path, allow_missing=3)


def test_a_fully_failed_but_fully_accounted_stage_still_reconciles(tmp_path):
    manifest = _manifest(shards=2)
    stage_dir = tmp_path / "dock"
    _write_output(stage_dir, 0)
    _write_status(stage_dir, 1, "failed", 1)
    census = collect_stage(manifest, "dock", tmp_path)
    assert census.reconciles and census.counts()[FAILED] == 1


def test_collecting_an_unsharded_stage_is_refused():
    manifest = StudyManifest.from_dict({
        "study": "J", "title": "t", "plan": ["README.md"],
        "stages": [{"name": "report", "command": ["echo", "x"], "output": "r.json",
                    "resources": {"partition": "shared", "walltime": "01:00:00"}}],
    })
    with pytest.raises(CollectError, match="not sharded"):
        collect_stage(manifest, "report", Path("."))


def test_collecting_a_stage_with_no_declared_output_is_refused(tmp_path):
    manifest = _manifest(shards=2, output=None)
    with pytest.raises(CollectError, match="declares no output"):
        collect_stage(manifest, "dock", tmp_path)


def test_the_summary_explains_a_missing_shard_rather_than_only_counting_it(tmp_path):
    manifest = _manifest(shards=2)
    _write_output(tmp_path / "dock", 0)
    census = take_census(manifest, manifest.stage("dock"), tmp_path / "dock")
    summary = census.summary()
    assert "1 missing" in summary
    assert "shrink the denominator" in summary


def test_provenance_states_the_census_and_the_plan(tmp_path):
    manifest = _manifest(shards=2)
    _write_output(tmp_path / "dock", 0)
    _write_status(tmp_path / "dock", 1, "killed", 143)
    census = collect_stage(manifest, "dock", tmp_path)
    text = provenance(manifest, census, site_name="rockfish", job_ids=["123456"], commit="abc1234")
    assert "rockfish" in text and "123456" in text and "abc1234" in text
    assert "README.md" in text
    assert "| killed | 1 |" in text
