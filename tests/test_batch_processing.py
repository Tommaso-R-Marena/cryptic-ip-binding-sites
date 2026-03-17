"""Tests for batch processing infrastructure."""

from pathlib import Path

from cryptic_ip.database.batch_processing import (
    AlphaFoldBatchDownloader,
    AnalysisCache,
    ParallelProcessor,
    append_results_to_file,
)


def _double_value(item):
    return {"uniprot_id": item["uniprot_id"], "value": item["value"] * 2}


class FakeResponse:
    def __init__(self, payload):
        self._payload = payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        return None


class FakeSession:
    def request(self, method, url, timeout=30, **kwargs):
        params = kwargs.get("params", {})
        query = params.get("query", "")
        proteome = query.split(":", 1)[1]
        payload = {
            "results": [
                {"primaryAccession": f"{proteome}_A"},
                {"primaryAccession": f"{proteome}_B"},
            ]
        }
        return FakeResponse(payload)


def test_batch_downloader_resume(tmp_path):
    downloader = AlphaFoldBatchDownloader(
        output_dir=tmp_path / "af",
        state_path=tmp_path / "state.json",
        requests_per_second=100,
        session=FakeSession(),
    )

    calls = []

    def fake_fetch(uniprot_id):
        calls.append(uniprot_id)
        (tmp_path / "af" / f"{uniprot_id}.pdb").write_text("MODEL")
        return tmp_path / "af" / f"{uniprot_id}.pdb"

    downloader.af_client.fetch_structure = fake_fetch

    summary_first = downloader.download_proteomes(["UP000002311"], resume=True)
    assert summary_first["downloaded"] == 2

    summary_second = downloader.download_proteomes(["UP000002311"], resume=True)
    assert summary_second["skipped"] == 2
    assert len(calls) == 2


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
