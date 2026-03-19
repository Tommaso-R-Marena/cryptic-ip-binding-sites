import os
import sys
import argparse
from pathlib import Path
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from pipeline.utils import setup_logger, load_config
from pipeline.output_formatter import (
    init_master_csv,
    append_to_master_csv,
    write_protein_json,
)
from validation.validation_report import process_protein

logger = setup_logger(__name__)


def worker(pdb_path: Path, out_dir: Path) -> dict:
    """Worker function for multiprocessing."""
    try:
        pockets = process_protein(pdb_path, out_dir)
        return {"status": "success", "pdb_path": pdb_path, "pockets": pockets}
    except Exception as e:
        logger.error(f"Error processing {pdb_path}: {e}")
        return {"status": "error", "pdb_path": pdb_path, "error": str(e)}


def screen_proteome(organism: str, workers: int, batch_size: int = 500):
    config = load_config()
    if organism not in config["organisms"]:
        logger.error(f"Unknown organism: {organism}")
        sys.exit(1)

    data_dir = Path(f"data/structures/{organism}")
    out_dir = Path(f"data/results/{organism}")
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"{organism}_master.csv"

    init_master_csv(csv_path)

    pdb_files = list(data_dir.glob("*.pdb"))
    if not pdb_files:
        logger.error(
            f"No PDB files found in {data_dir}. Did you run fetch_{organism}.yml?"
        )
        sys.exit(1)

    logger.info(
        f"Starting screen for {len(pdb_files)} proteins using {workers} workers."
    )

    # Process in batches to manage memory
    for i in range(0, len(pdb_files), batch_size):
        batch = pdb_files[i : i + batch_size]
        logger.info(
            f"Processing batch {i//batch_size + 1}/{(len(pdb_files) + batch_size - 1)//batch_size}"
        )

        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(worker, pdb, out_dir): pdb for pdb in batch}

            for future in tqdm(as_completed(futures), total=len(batch)):
                result = future.result()
                if result["status"] == "success":
                    pdb = result["pdb_path"]
                    # parse UniProt ID: Assuming format AF-{ID}-F1-model_v4.pdb or similar
                    parts = pdb.name.split("-")
                    uniprot_id = parts[1] if len(parts) > 1 else pdb.stem

                    pockets = result["pockets"]
                    if pockets:
                        write_protein_json(uniprot_id, organism, pockets, out_dir)
                        append_to_master_csv(
                            csv_path, uniprot_id, organism, "Unknown", pockets
                        )

    logger.info(f"Screening complete. Results written to {csv_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--organism",
        required=True,
        help="Organism to screen (yeast, human, dictyostelium)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() - 1,
        help="Number of concurrent workers",
    )
    parser.add_argument(
        "--batch-size", type=int, default=500, help="Batch size for writing to disk"
    )
    args = parser.parse_args()

    screen_proteome(args.organism, args.workers, args.batch_size)
