import json
import csv
from pathlib import Path
from typing import Dict, Any, List
from .utils import setup_logger, load_config

logger = setup_logger(__name__)


def write_protein_json(
    protein_id: str, organism: str, pockets: List[Dict[str, Any]], out_dir: Path
) -> Path:
    """Writes the detailed pocket JSON for a single protein."""
    out_file = out_dir / f"{protein_id}_results.json"
    data = {"protein_id": protein_id, "organism": organism, "pockets": pockets}
    with open(out_file, "w") as f:
        json.dump(data, f, indent=2)
    return out_file


def init_master_csv(csv_path: Path):
    """Initializes the master CSV file with headers if it doesn't exist."""
    headers = [
        "UniProtID",
        "Organism",
        "ProteinName",
        "PocketID",
        "CompositeScore",
        "PocketVolume",
        "MeanSASA",
        "ElectrostaticPotential",
        "BasicResidueCount",
        "MeanPLDDT",
        "PassAllFilters",
        "ManualReviewFlag",
    ]
    if not csv_path.exists():
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(headers)


def append_to_master_csv(
    csv_path: Path,
    protein_id: str,
    organism: str,
    protein_name: str,
    pockets: List[Dict[str, Any]],
):
    """Appends results for a protein to the master CSV."""
    config = load_config()
    flag_threshold = config["pipeline"]["scoring"]["min_composite_score_for_flag"]

    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        for pkt in pockets:
            row = [
                protein_id,
                organism,
                protein_name,
                pkt["id"],
                pkt["composite_score"],
                pkt.get("volume", 0.0),
                round(pkt.get("mean_sasa", 0.0), 2),
                round(pkt.get("electrostatic_potential", 0.0), 2),
                pkt.get("basic_residue_count", 0),
                round(pkt.get("mean_plddt", 0.0), 2),
                pkt.get("pass_all_filters", False),
                pkt["composite_score"] >= flag_threshold,
            ]
            writer.writerow(row)
