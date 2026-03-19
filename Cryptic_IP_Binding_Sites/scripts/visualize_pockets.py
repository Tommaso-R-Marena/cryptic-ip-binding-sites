#!/usr/ স্বপ্নেরenv python
import argparse
import json
import traceback
from pathlib import Path
import sys

# Add the parent directory to PYTHONPATH to allow importing pipeline
sys.path.append(str(Path(__file__).parent.parent))

from pipeline.visualization import generate_pymol_script  # noqa: E402
from pipeline.utils import setup_logger  # noqa: E402

logger = setup_logger("visualize_pockets")


def main():
    parser = argparse.ArgumentParser(
        description="Generate PyMOL visualization script for IP binding pockets."
    )
    parser.add_argument(
        "--pdb", required=True, type=Path, help="Path to the input PDB file."
    )
    parser.add_argument(
        "--pockets",
        required=True,
        type=Path,
        help="Path to the scored JSON pockets file.",
    )
    parser.add_argument(
        "--output", required=True, type=Path, help="Path to output the .pml script."
    )

    args = parser.parse_args()

    if not args.pdb.exists():
        logger.error(f"PDB file not found: {args.pdb}")
        return

    if not args.pockets.exists():
        logger.error(f"Pockets JSON file not found: {args.pockets}")
        return

    try:
        with open(args.pockets, "r") as f:
            pockets_data = json.load(f)

        # Top 3 pockets that passed filters, or just top 3 overall if none passed
        passed_pockets = [p for p in pockets_data if p.get("pass_all_filters", False)]
        if not passed_pockets:
            logger.warning(
                "No pockets passed all filters! Visualizing the top 3 absolute highest scorings."
            )
            target_pockets = pockets_data[:3]
        else:
            target_pockets = passed_pockets[:3]

        generate_pymol_script(args.pdb, target_pockets, args.output)
        logger.info(f"Successfully generated PyMOL script at {args.output}")

    except Exception as e:
        logger.error(f"Failed to generate visualization: {e}")
        logger.error(traceback.format_exc())


if __name__ == "__main__":
    main()
