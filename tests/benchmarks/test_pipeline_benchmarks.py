"""Performance benchmark tests for pipeline regressions."""

import time

from cryptic_ip.database.batch_processing import AnalysisCache


def test_cache_batch_write_benchmark(tmp_path):
    cache = AnalysisCache(tmp_path / "bench.sqlite", pipeline_version="v1")
    rows = [(f"P{i}", "org", {"score": float(i)}) for i in range(500)]

    started = time.perf_counter()
    written = cache.set_cached_results_batch(rows)
    elapsed = time.perf_counter() - started

    assert written == 500
    assert elapsed < 2.0
    cache.vacuum()
    cache.close()
