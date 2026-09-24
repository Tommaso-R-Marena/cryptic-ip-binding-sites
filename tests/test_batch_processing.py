"""Tests for batch processing infrastructure."""

import json
import threading

import pytest

from cryptic_ip.database.batch_processing import (
    AlphaFoldBatchDownloader,
    AnalysisCache,
    ParallelProcessor,
    append_results_to_file,
)
from cryptic_ip.errors import OperationTimeoutError, RecoveryStateError, ValidationError


def _double_value(item):
    return {"uniprot_id": item["uniprot_id"], "value": item["value"] * 2}


class FakeResponse:
    def __init__(self, text):
        self.text = text

    def raise_for_status(self):
        return None


class FakeSession:
    def __init__(self):
        self.calls = []

    def request(self, method, url, timeout=30, **kwargs):
        self.calls.append((url, kwargs.get("params", {})))
        proteome = kwargs.get("params", {}).get("query", "").split(":", 1)[1]
        return FakeResponse(f"{proteome}_A\n{proteome}_B\n{proteome}_C\n")


class _Result:
    def __init__(self, ok, not_found=False):
        self.ok, self.not_found, self.error = ok, not_found, "" if ok else "missing"


def test_batch_downloader_resume(tmp_path):
    session = FakeSession()
    downloader = AlphaFoldBatchDownloader(
        output_dir=tmp_path / "af",
        state_path=tmp_path / "state.json",
        requests_per_second=10,
        session=session,
    )

    calls = []

    def fake_fetch(uniprot_ids):
        calls.extend(uniprot_ids)
        # _C has no AlphaFold model: not found, not a failure.
        return {uid: _Result(ok=not uid.endswith("_C"), not_found=uid.endswith("_C")) for uid in uniprot_ids}

    downloader._fetch_models = fake_fetch

    summary_first = downloader.download_proteomes(["UP000002311"], resume=True)
    assert summary_first["downloaded"] == 2
    assert summary_first["not_found"] == 1
    assert summary_first["failed"] == 0

    summary_second = downloader.download_proteomes(["UP000002311"], resume=True)
    assert summary_second["skipped"] == 2
    # Only the absent accession is asked for again.
    assert calls == ["UP000002311_A", "UP000002311_B", "UP000002311_C", "UP000002311_C"]
    # One streamed request per listing, not one per page of 500.
    assert all("stream" in url for url, _ in session.calls)
    assert session.calls[0][1]["format"] == "list"


def test_batch_downloader_corrupt_state(tmp_path):
    state = tmp_path / "state.json"
    state.write_text("{not-json")
    downloader = AlphaFoldBatchDownloader(
        output_dir=tmp_path / "af",
        state_path=state,
        requests_per_second=10,
        session=FakeSession(),
    )
    with pytest.raises(RecoveryStateError):
        downloader.download_proteomes(["UP000002311"], resume=True)


def test_batch_downloader_config_validation(tmp_path):
    with pytest.raises(ValidationError):
        AlphaFoldBatchDownloader(output_dir=tmp_path / "af", requests_per_second=0.05)


def test_analysis_cache_roundtrip_and_invalidate(tmp_path):
    db_path = tmp_path / "cache.sqlite"
    cache = AnalysisCache(db_path, pipeline_version="v1", pipeline_params={"threshold": 0.6})
    cache.set_cached_result("P12345", "test_org", {"score": 0.9})
    cache.set_cached_results_batch([
        ("P99999", "test_org", {"score": 0.7}),
        ("P11111", "test_org", {"score": 0.2}),
    ])

    assert cache.get_cached_result("P12345", "test_org") == {"score": 0.9}
    assert cache.get_cached_result("P99999", "test_org") == {"score": 0.7}

    csv_path = cache.export_results(tmp_path / "cache.csv", "csv")
    json_path = cache.export_results(tmp_path / "cache.json", "json")
    assert csv_path.exists()
    assert json_path.exists()
    cache.vacuum()
    cache.close()

    other = AnalysisCache(db_path, pipeline_version="v2", pipeline_params={"threshold": 0.7})
    deleted = other.invalidate_outdated_cache()
    assert deleted >= 1
    assert other.get_cached_result("P12345", "test_org") is None
    other.close()


def test_analysis_cache_concurrent_access(tmp_path):
    cache = AnalysisCache(tmp_path / "cache.sqlite", pipeline_version="v1")

    def writer(index: int):
        cache.set_cached_result(f"P{index}", "org", {"score": index})

    threads = [threading.Thread(target=writer, args=(i,)) for i in range(10)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()

    assert cache.get_cached_result("P3", "org") == {"score": 3}
    cache.close()


def test_parallel_processor_and_csv_append(tmp_path):
    items = [{"uniprot_id": f"P{i}", "value": i} for i in range(8)]
    processor = ParallelProcessor(
        analyze_function=_double_value,
        workers=2,
        chunk_size=3,
        checkpoint_path=tmp_path / "checkpoint.json",
    )

    results = processor.run(items, resume=False)
    assert len(results) == len(items)

    output_csv = append_results_to_file(results, tmp_path / "results.csv")
    assert output_csv.exists()

    processor_resume = ParallelProcessor(
        analyze_function=_double_value,
        workers=2,
        chunk_size=2,
        checkpoint_path=tmp_path / "checkpoint.json",
    )
    resumed = processor_resume.run(items, resume=True)
    assert resumed == []


def test_parallel_processor_timeout(tmp_path):
    def slow(item):
        import time

        time.sleep(0.2)
        return {"uniprot_id": item["uniprot_id"]}

    processor = ParallelProcessor(
        analyze_function=slow,
        workers=1,
        chunk_size=1,
        checkpoint_path=tmp_path / "checkpoint.json",
        timeout_per_item_seconds=0.01,
    )

    with pytest.raises(OperationTimeoutError):
        processor.run([{"uniprot_id": "P1"}], resume=False)


def test_parallel_processor_corrupt_checkpoint(tmp_path):
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text("{broken")
    processor = ParallelProcessor(
        analyze_function=_double_value,
        workers=1,
        chunk_size=1,
        checkpoint_path=checkpoint,
    )
    with pytest.raises(RecoveryStateError):
        processor.run([{"uniprot_id": "P1", "value": 1}], resume=True)
