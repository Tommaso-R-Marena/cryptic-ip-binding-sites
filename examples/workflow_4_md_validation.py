#!/usr/bin/env python3
"""Workflow 4: select top screening candidates -> run MD validation -> filter stable hits.

Expected outputs
- `md_validation_report.csv`
- `stable_candidates.csv`
- per-candidate trajectory + visualization scripts

Typical runtime
- Mocked CI mode: <1 minute
- Real OpenMM production runs: hours to days

Troubleshooting
- `OpenMM is required`: install openmm and CUDA toolkit if using GPU.
- `MDTraj is required`: install mdtraj for trajectory analysis.
- No stable candidates: adjust thresholds or increase MD sampling.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from cryptic_ip.validation.md_validation import OpenMMMDValidationPipeline


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run MD validation on top cryptic-site candidates.")
    parser.add_argument("--candidates-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("results/workflow4"))
    parser.add_argument("--top-n", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    pipeline = OpenMMMDValidationPipeline(output_dir=args.output_dir)
    report = pipeline.validate_top_candidates(args.candidates_csv, top_n=args.top_n)
    print(report[["candidate_id", "classification"]].to_string(index=False))


if __name__ == "__main__":
    main()
