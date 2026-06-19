#!/usr/bin/env python3
"""Pilot MD validation on top cryptic-pocket candidates (short production runs)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.validation.md_validation import MDSimulationConfig, OpenMMMDValidationPipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates-csv",
        type=Path,
        default=Path("results/publication/gallery/gallery_inputs.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/md_validation/pilot"))
    parser.add_argument("--top-n", type=int, default=5)
    parser.add_argument("--production-ns", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = MDSimulationConfig(production_ns=args.production_ns, equilibration_ns=0.5)
    pipeline = OpenMMMDValidationPipeline(output_dir=args.output_dir, config=config)
    report = pipeline.validate_top_candidates(args.candidates_csv, top_n=args.top_n)
    print(f"MD pilot complete: {len(report)} candidates processed")
    print(f"Report: {args.output_dir / 'md_validation_report.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
