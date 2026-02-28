#!/usr/bin/env python3
"""Workflow 1: full analysis of one protein (download -> pockets -> score).

Expected outputs
- Downloaded structure (PDB)
- CSV with ranked pockets and composite scores

Typical runtime
- ~2-5 minutes with fpocket/ProDy/APBS installed

Troubleshooting
- `fpocket not found`: install fpocket and ensure it is on PATH.
- `pdb2pqr/APBS` failures: rerun with `--skip-electrostatics`.
- Empty pocket table: verify the structure file is complete and not truncated.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from cryptic_ip.analysis.analyzer import ProteinAnalyzer
from cryptic_ip.database.alphafold_client import AlphaFoldClient


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Analyze one UniProt protein for cryptic IP pockets.")
    parser.add_argument("--uniprot-id", default="P78563", help="UniProt accession (default: ADAR2)")
    parser.add_argument("--output-dir", type=Path, default=Path("results/workflow1"))
    parser.add_argument("--skip-electrostatics", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    client = AlphaFoldClient(cache_dir=args.output_dir / "structures")
    structure_path = client.fetch_structure(args.uniprot_id)

    analyzer = ProteinAnalyzer(str(structure_path), work_dir=str(args.output_dir / "analysis"))
    analyzer.detect_pockets()

    if not args.skip_electrostatics:
        analyzer.calculate_electrostatics()

    scored = analyzer.score_all_pockets()
    scored.to_csv(args.output_dir / f"{args.uniprot_id}_pocket_scores.csv", index=False)
    print(scored.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
