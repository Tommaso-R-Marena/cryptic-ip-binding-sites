#!/usr/bin/env python
"""Phase 5: Comparative multi-organism cryptic IP-binding analysis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cryptic_ip.analysis.comparative_analysis import ComparativeIPAnalysis


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run comparative analysis across 3 proteomes")
    parser.add_argument("--hits", required=True, help="CSV with columns: organism, protein_id, is_hit")
    parser.add_argument(
        "--orthologs",
        required=True,
        help="CSV with columns: orthogroup, organism, protein_id",
    )
    parser.add_argument(
        "--go-associations-json",
        required=True,
        help=(
            "JSON mapping organism -> GO association CSV path. Each association CSV must have "
            "columns protein_id, go_id"
        ),
    )
    parser.add_argument("--go-obo", required=True, help="Path to go-basic.obo")
    parser.add_argument(
        "--ip6-json",
        default='{"S_cerevisiae": 20, "H_sapiens": 25, "D_discoideum": 520}',
        help="JSON mapping organism -> cellular IP6 concentration in uM",
    )
    parser.add_argument("--output-dir", default="results/comparative_analysis", help="Output folder")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    hits_df = pd.read_csv(args.hits)
    ortholog_df = pd.read_csv(args.orthologs)
    go_assoc_files = json.loads(args.go_associations_json)
    ip6_map = json.loads(args.ip6_json)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    analysis = ComparativeIPAnalysis(confidence_level=0.95)
    result = analysis.run_pipeline(
        hits_df=hits_df,
        ip6_map=ip6_map,
        ortholog_df=ortholog_df,
        go_association_files=go_assoc_files,
        go_obo_path=args.go_obo,
        output_dir=str(output_dir),
    )

    print("Comparative analysis complete.")
    print(result.hit_rate_table.to_string(index=False))
    print(result.correlation_table.to_string(index=False))
    print(f"Outputs written to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
