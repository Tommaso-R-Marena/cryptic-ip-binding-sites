#!/usr/bin/env python
"""Batch infrastructure entrypoint for large AlphaFold proteome screening."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from cryptic_ip.database.batch_processing import (
    AlphaFoldBatchDownloader,
    AnalysisCache,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="AlphaFold batch download + cache export utility")
    sub = parser.add_subparsers(dest="command", required=True)

    dl = sub.add_parser("download", help="Download AlphaFold structures for UniProt proteomes")
    dl.add_argument(
        "--proteomes",
        nargs="+",
        default=["UP000002311", "UP000005640", "UP000002195"],
        help="UniProt proteome IDs",
    )
    dl.add_argument("--output-dir", default="data/alphafold")
    dl.add_argument("--state-path", default="data/alphafold/download_state.json")
    dl.add_argument("--rps", type=float, default=3.0, help="requests per second")

    cache = sub.add_parser("cache-export", help="Export cached results to csv/json/hdf5")
    cache.add_argument("--db", required=True)
    cache.add_argument("--pipeline-version", required=True)
    cache.add_argument("--pipeline-params", default="{}", help="JSON string for pipeline parameters")
    cache.add_argument("--output", required=True)
    cache.add_argument("--format", choices=["csv", "json", "hdf5"], required=True)

    args = parser.parse_args()

    if args.command == "download":
        downloader = AlphaFoldBatchDownloader(
            output_dir=Path(args.output_dir),
            state_path=Path(args.state_path),
            requests_per_second=args.rps,
        )
        summary = downloader.download_proteomes(args.proteomes, resume=True)
        print(summary)

    if args.command == "cache-export":
        params = json.loads(args.pipeline_params)
        cache_db = AnalysisCache(args.db, pipeline_version=args.pipeline_version, pipeline_params=params)
        try:
            cache_db.export_results(args.output, args.format)
        finally:
            cache_db.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
