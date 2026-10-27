"""The retry policy that assembles a rate-limited structure set.

The strict yeast re-screen needs all 499 of its models and the endpoint
rate-limits by address, so the set is built over several attempts. Until now
the thing that caused the next attempt was a person re-running the workflow,
and the set sat incomplete for days between nudges.

The rules that matter here are the stopping rules. A loop that always sleeps
its full budget wastes a runner; one that gives up on the first empty attempt
never assembles anything. The tests pin all three exits - complete, stalled,
budget exhausted - plus the one deficit that must never be waited on.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "fetch_until_complete", ROOT / "scripts" / "fetch_until_complete.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["fetch_until_complete"] = script
SPEC.loader.exec_module(script)


def _manifest(path: Path, *, ok: int, failed: int = 0, not_found: int = 0) -> None:
    rows = ([{"key": f"ok{i}", "ok": True, "not_found": False} for i in range(ok)]
            + [{"key": f"bad{i}", "ok": False, "not_found": False} for i in range(failed)]
            + [{"key": f"gone{i}", "ok": False, "not_found": True} for i in range(not_found)])
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rows))


class _Endpoint:
    """A fetcher that hands over `per_attempt` models each time it is called."""

    def __init__(self, manifest: Path, total: int, per_attempt):
        self.manifest, self.total = manifest, total
        self.per_attempt = list(per_attempt)
        self.ok, self.calls = 0, 0

    def __call__(self, argv):
        gain = self.per_attempt[self.calls] if self.calls < len(self.per_attempt) else 0
        self.calls += 1
        self.ok = min(self.total, self.ok + gain)
        _manifest(self.manifest, ok=self.ok, failed=self.total - self.ok)
        return 0 if self.ok == self.total else 1


# ------------------------------------------------------------ reading state

def test_an_absent_manifest_is_not_an_error(tmp_path):
    assert script.read_manifest(tmp_path / "none.json") == {
        "total": 0, "ok": 0, "not_found": 0, "failed": 0}


def test_an_unparseable_manifest_is_not_an_error(tmp_path):
    """A job killed mid-write leaves a truncated file; the next run must start."""
    path = tmp_path / "m.json"
    path.write_text("[{\"ok\": tru")
    assert script.read_manifest(path)["total"] == 0


def test_a_cached_model_counts_as_present(tmp_path):
    """The screen needs the file, not the fetch: `from_cache` is still present."""
    path = tmp_path / "m.json"
    path.write_text(json.dumps([{"ok": True, "from_cache": True, "not_found": False},
                                {"ok": True, "from_cache": False, "not_found": False}]))
    counts = script.read_manifest(path)
    assert counts == {"total": 2, "ok": 2, "not_found": 0, "failed": 0}


def test_absent_upstream_is_counted_apart_from_owed(tmp_path):
    path = tmp_path / "m.json"
    _manifest(path, ok=3, failed=2, not_found=1)
    assert script.read_manifest(path) == {"total": 6, "ok": 3, "not_found": 1, "failed": 2}


# ----------------------------------------------------------- stopping rules

def test_it_stops_as_soon_as_the_set_is_complete(tmp_path):
    manifest = tmp_path / "alphafold_manifest.json"
    endpoint = _Endpoint(manifest, total=10, per_attempt=[4, 6, 99])
    slept = []
    report = script.assemble([], manifest=manifest, attempts=5, sleep_s=30,
                             fetch=endpoint, sleeper=slept.append)
    assert report["reason"] == script.COMPLETE and report["complete"] is True
    assert endpoint.calls == 2, "it kept fetching after the set was complete"
    assert slept == [30.0], "it slept once, between the two attempts"


def test_it_does_not_sleep_after_the_attempt_that_completes_it(tmp_path):
    manifest = tmp_path / "alphafold_manifest.json"
    slept = []
    script.assemble([], manifest=manifest, attempts=3, sleep_s=60,
                    fetch=_Endpoint(manifest, total=5, per_attempt=[5]),
                    sleeper=slept.append)
    assert slept == [], "a completed set must not wait"


def test_two_barren_attempts_stop_the_loop(tmp_path):
    """The window is longer than this job can wait, so stop rather than burn it."""
    manifest = tmp_path / "alphafold_manifest.json"
    endpoint = _Endpoint(manifest, total=10, per_attempt=[3, 0, 0, 7])
    slept = []
    report = script.assemble([], manifest=manifest, attempts=9, sleep_s=10,
                             fetch=endpoint, sleeper=slept.append)
    assert report["reason"] == script.STALLED and report["complete"] is False
    assert endpoint.calls == 3
    assert report["ok"] == 3, "the partial set is kept and reported"


def test_one_barren_attempt_does_not_stop_the_loop(tmp_path):
    """A single empty pass is ordinary against a throttling endpoint."""
    manifest = tmp_path / "alphafold_manifest.json"
    endpoint = _Endpoint(manifest, total=10, per_attempt=[5, 0, 5])
    report = script.assemble([], manifest=manifest, attempts=9, sleep_s=1,
                             fetch=endpoint, sleeper=lambda s: None)
    assert report["reason"] == script.COMPLETE
    assert endpoint.calls == 3


def test_the_attempt_budget_is_respected(tmp_path):
    manifest = tmp_path / "alphafold_manifest.json"
    endpoint = _Endpoint(manifest, total=100, per_attempt=[1, 1, 1, 1, 1, 1])
    report = script.assemble([], manifest=manifest, attempts=3, sleep_s=1,
                             fetch=endpoint, sleeper=lambda s: None)
    assert report["reason"] == script.EXHAUSTED and report["complete"] is False
    assert endpoint.calls == 3 and report["attempts_used"] == 3


def test_a_deficit_of_absent_models_is_never_waited_on(tmp_path):
    """No amount of waiting gives the endpoint a model it does not have."""
    manifest = tmp_path / "alphafold_manifest.json"
    calls = []

    def fetch(argv):
        calls.append(1)
        _manifest(manifest, ok=8, failed=0, not_found=2)
        return 1

    slept = []
    report = script.assemble([], manifest=manifest, attempts=6, sleep_s=99,
                             fetch=fetch, sleeper=slept.append)
    assert report["reason"] == script.NOT_FOUND_ONLY
    assert report["complete"] is False, "an absent model still fails the set"
    assert len(calls) == 1 and slept == []


# -------------------------------------------------------------- the record

def test_progress_is_written_after_every_attempt(tmp_path):
    """A job killed at its budget still leaves the recovery measurement."""
    manifest = tmp_path / "alphafold_manifest.json"
    progress = tmp_path / "progress.json"
    seen = []

    endpoint = _Endpoint(manifest, total=10, per_attempt=[2, 3, 5])

    def watching(argv):
        # Each attempt sees the record as the previous attempt left it, so the
        # first sees no file at all and the Nth sees N-1 rows.
        seen.append(len(json.loads(progress.read_text())["attempts"])
                    if progress.exists() else None)
        return endpoint(argv)

    script.assemble([], manifest=manifest, attempts=5, sleep_s=1, progress=progress,
                    fetch=watching, sleeper=lambda s: None)
    assert seen == [None, 1, 2], "the record lagged the attempts it is meant to hold"
    final = json.loads(progress.read_text())
    assert final["reason"] == script.COMPLETE
    assert [row["gained"] for row in final["attempts"]] == [2, 3, 5]
    assert final["ok"] == 10


def test_the_report_names_how_much_each_attempt_gained(tmp_path):
    manifest = tmp_path / "alphafold_manifest.json"
    report = script.assemble([], manifest=manifest, attempts=4, sleep_s=1,
                             fetch=_Endpoint(manifest, total=9, per_attempt=[4, 0, 5]),
                             sleeper=lambda s: None)
    assert [row["gained"] for row in report["attempts"]] == [4, 0, 5]
    assert report["sleep_seconds"] == 1.0


def test_models_already_on_disk_are_counted_before_the_first_attempt(tmp_path):
    """A restored cache must not be reported as this run's gain."""
    manifest = tmp_path / "alphafold_manifest.json"
    _manifest(manifest, ok=6, failed=4)
    endpoint = _Endpoint(manifest, total=10, per_attempt=[0])
    endpoint.ok = 6
    report = script.assemble([], manifest=manifest, attempts=1, sleep_s=1,
                             fetch=endpoint, sleeper=lambda s: None)
    assert report["attempts"][0]["gained"] == 0
    assert report["ok"] == 6


# ------------------------------------------------------------------- CLI

def test_the_cli_passes_the_fetcher_only_the_flags_it_was_given(monkeypatch, tmp_path):
    captured = {}

    def fake_assemble(argv, **kwargs):
        captured["argv"] = list(argv)
        captured["kwargs"] = kwargs
        return {"complete": True, "reason": script.COMPLETE, "attempts": []}

    monkeypatch.setattr(script, "assemble", fake_assemble)
    code = script.main(["alphafold", "--out", str(tmp_path), "--ids-json", "ids.json",
                        "--json-key", "processed", "--attempts", "7", "--sleep", "12"])
    assert code == 0
    assert captured["argv"] == ["alphafold", "--out", str(tmp_path),
                                "--ids-json", "ids.json", "--json-key", "processed"]
    assert captured["kwargs"]["attempts"] == 7
    assert captured["kwargs"]["sleep_s"] == 12.0
    assert captured["kwargs"]["manifest"] == tmp_path / "alphafold_manifest.json"


def test_an_incomplete_set_exits_non_zero(monkeypatch, tmp_path):
    """The screen must not run on a partial set, so this has to fail the step."""
    monkeypatch.setattr(script, "assemble", lambda *a, **k: {
        "complete": False, "reason": script.STALLED, "attempts": []})
    assert script.main(["alphafold", "--out", str(tmp_path)]) == 1
