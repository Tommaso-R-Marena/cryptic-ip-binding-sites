#!/usr/bin/env python3
"""Workflow 3: batch proteome screen -> comparative stats -> publication figures.

Expected outputs
- `batch_predictions.csv`
- comparative tables (`hit_rates.csv`, `pairwise_tests.csv`)
- figure assets in output directory

Typical runtime
- Download/screen: hours for whole proteomes
- Comparative + figures: <5 minutes

Troubleshooting
- UniProt API throttling: lower request rate or retry with resume.
- Incomplete proteome runs: rerun with checkpoint resume enabled.
- Figure failures: verify all required CSV columns exist.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from cryptic_ip.analysis.comparative_analysis import ComparativeIPAnalysis
from cryptic_ip.database.batch_processing import AlphaFoldBatchDownloader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Screen proteomes and generate comparative outputs.")
    parser.add_argument("--proteome", action="append", required=True, help="UniProt proteome ID")
    parser.add_argument("--output-dir", type=Path, default=Path("results/workflow3"))
    parser.add_argument("--skip-download", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        downloader = AlphaFoldBatchDownloader(output_dir=args.output_dir / "structures")
        downloader.download_proteomes(args.proteome, resume=True)

    predictions = pd.read_csv(args.output_dir / "batch_predictions.csv")
    hits = predictions[["organism", "protein_id", "is_hit"]]

    ip6_map = {org: float(val) for org, val in predictions.groupby("organism")["ip6_uM"].first().items()}
    comp = ComparativeIPAnalysis()
    hit_rates = comp.compute_hit_rates(hits, ip6_map=ip6_map)
    pairwise = comp.pairwise_organism_tests(hit_rates)

    hit_rates.to_csv(args.output_dir / "hit_rates.csv", index=False)
    pairwise.to_csv(args.output_dir / "pairwise_tests.csv", index=False)

    ortholog_summary = pd.read_csv(args.output_dir / "ortholog_summary.csv")
    go_summary = pd.read_csv(args.output_dir / "go_summary.csv")
    comp.generate_figures(hit_rates, ortholog_summary, go_summary, output_dir=str(args.output_dir / "figures"))


if __name__ == "__main__":
    main()
