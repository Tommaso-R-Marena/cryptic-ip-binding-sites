import sys
import argparse
from pathlib import Path
from pipeline.utils import load_config, setup_logger

logger = setup_logger(__name__)


def qc_proteome(organism: str) -> bool:
    """Runs Quality Control checks on downloaded proteome."""
    config = load_config()

    if organism not in config["organisms"]:
        logger.error(f"Unknown organism: {organism}")
        return False

    expected_count = config["organisms"][organism]["expected_file_count"]
    data_dir = Path(f"data/structures/{organism}")

    if not data_dir.exists():
        logger.error(f"Data directory {data_dir} does not exist.")
        return False

    pdb_files = list(data_dir.glob("*.pdb"))
    actual_count = len(pdb_files)

    if actual_count != expected_count:
        logger.error(
            f"QC FAILED: Expected {expected_count} PDB files, found {actual_count}"
        )
        return False

    zero_byte = 0
    for pdb in pdb_files:
        if pdb.stat().st_size == 0:
            zero_byte += 1

    if zero_byte > 0:
        logger.error(f"QC FAILED: Found {zero_byte} zero-byte PDB files.")
        return False

    logger.info(
        f"QC PASSED: Validated {actual_count} non-empty PDB files for {organism}."
    )
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism", required=True)
    args = parser.parse_args()

    if not qc_proteome(args.organism):
        sys.exit(1)
