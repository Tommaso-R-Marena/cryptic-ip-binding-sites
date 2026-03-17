#!/usr/bin/env python
"""Profile cryptic IP analysis pipeline performance and bottlenecks."""

from __future__ import annotations

import argparse
import json
import time
import tracemalloc
from pathlib import Path

import pandas as pd

from cryptic_ip.analysis.analyzer import ProteinAnalyzer
from cryptic_ip.utils.profiling import TIMING_REGISTRY


def profile_single(pdb_path: str, skip_electrostatics: bool) -> dict:
    analyzer = ProteinAnalyzer(pdb_path, skip_electrostatics=skip_electrostatics)
    tracemalloc.start()
    started = time.perf_counter()
    results = analyzer.run_pipeline()
    elapsed = time.perf_counter() - started
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    analyzer.cleanup()
    return {
        "runtime_seconds": elapsed,
        "peak_memory_mb": peak / (1024 * 1024),
        "pocket_count": 0 if results is None else len(results),
    }


def profile_batch(pdb_paths: list[str], skip_electrostatics: bool) -> dict:
    rows = []
    for pdb_path in pdb_paths:
        row = profile_single(pdb_path, skip_electrostatics=skip_electrostatics)
        row["pdb_path"] = pdb_path
        rows.append(row)
    frame = pd.DataFrame(rows)
    return {
        "mean_runtime_seconds": float(frame["runtime_seconds"].mean()),
        "mean_peak_memory_mb": float(frame["peak_memory_mb"].mean()),
        "throughput_proteins_per_hour": float(3600.0 / frame["runtime_seconds"].mean()),
    }


def build_recommendations(summary: dict) -> list[str]:
    recs = []
    if summary["throughput_proteins_per_hour"] < 1000:
        recs.append("Enable --skip-electrostatics for first-pass screening and run APBS only on top hits.")
    recs.append("Tune chunk size and workers to keep CPU utilization near saturation.")
    recs.append("Use AnalysisCache + batch writes to avoid repeated pocket rescoring.")
    return recs


def main() -> int:
    parser = argparse.ArgumentParser(description="Profile cryptic pipeline runtime and memory")
    parser.add_argument("--pdb", action="append", required=True, help="One or more PDB files")
    parser.add_argument("--skip-electrostatics", action="store_true")
    parser.add_argument("--report", default="profiling_report.json")
    args = parser.parse_args()

    TIMING_REGISTRY.clear()
    single = profile_single(args.pdb[0], skip_electrostatics=args.skip_electrostatics)
    batch = profile_batch(args.pdb, skip_electrostatics=args.skip_electrostatics)
    timings = [record.__dict__ for record in TIMING_REGISTRY.summary()]

    payload = {
        "single_protein": single,
        "batch_summary": batch,
        "function_timings": timings,
        "recommendations": build_recommendations(batch),
    }
    Path(args.report).write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
