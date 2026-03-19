import math
import json
from pathlib import Path
from typing import Dict, Any, List
from .utils import setup_logger, load_config

logger = setup_logger(__name__)


def calculate_basic_residue_count(
    pocket_center: tuple, pdb_path: Path, distance_cutoff: float = 5.0
) -> int:
    """Counts Arg, Lys, His residues within distance_cutoff of pocket_center."""
    count = 0
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                if res_name in ["ARG", "LYS", "HIS"]:
                    try:
                        # Use C-alpha for distance Check
                        atom_name = line[12:16].strip()
                        if atom_name == "CA":
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())

                            dist = math.sqrt(
                                (x - pocket_center[0]) ** 2
                                + (y - pocket_center[1]) ** 2
                                + (z - pocket_center[2]) ** 2
                            )
                            if dist <= distance_cutoff:
                                count += 1
                    except ValueError:
                        pass
    return count


def compute_composite_score(pocket: Dict[str, Any], config: Dict[str, Any]) -> float:
    """Computes the 0-1 composite score for a pocket."""
    w = config["pipeline"]["scoring"]["weights"]

    # Normalize inputs (heuristics for 0-1 scaling based on expected ranges)
    # Depth: typically 10-30A. Let's cap at 30.
    norm_depth = min(1.0, pocket.get("depth", 0) / 30.0)

    # SASA: lower is better for deeply buried. Cap at 50A2 for "exposed".
    # (1 - normalized_SASA) means 0 SASA = 1.0 score
    raw_sasa = pocket.get("mean_sasa", 50.0)
    norm_sasa = min(1.0, raw_sasa / 50.0)
    inv_sasa = 1.0 - norm_sasa

    # Electrostatic: Positive target. Cap at +15 kT/e
    raw_elec = max(0.0, pocket.get("electrostatic_potential", 0.0))
    norm_elec = min(1.0, raw_elec / 15.0)

    # Basic Residues: target is >= 4. Let's cap at 8.
    raw_basic = pocket.get("basic_residue_count", 0)
    norm_basic = min(1.0, raw_basic / 8.0)

    s = (
        w["depth"] * norm_depth
        + w["sasa"] * inv_sasa
        + w["electrostatic"] * norm_elec
        + w["basic_residues"] * norm_basic
    )

    return round(s, 3)


def score_and_filter_pockets(
    pockets: List[Dict[str, Any]],
    sasa_data: Dict[str, float],
    plddt_path: Path,
    cleaned_pdb_path: Path,
) -> List[Dict[str, Any]]:
    """Applies filters and computes final scores for pockets."""
    config = load_config()
    scored_pockets = []

    # Load pLDDT scores
    with open(plddt_path, "r") as f:
        plddt_data = json.load(f)

    for pkt in pockets:
        # 1. Calc Mean SASA
        pocket_sasa_vals = [
            sasa_data.get(r, 50.0) for r in pkt["residues"] if r in sasa_data
        ]
        pkt["mean_sasa"] = (
            sum(pocket_sasa_vals) / len(pocket_sasa_vals) if pocket_sasa_vals else 50.0
        )

        # 2. Basic Residue Count
        pkt["basic_residue_count"] = calculate_basic_residue_count(
            pkt["center"], cleaned_pdb_path
        )

        # 3. Mean pLDDT
        pocket_plddt_vals = [
            plddt_data.get(r[4:], 0.0) for r in pkt["residues"]
        ]  # strip RES_ prefix to match chain_res
        pkt["mean_plddt"] = (
            sum(pocket_plddt_vals) / len(pocket_plddt_vals)
            if pocket_plddt_vals
            else 0.0
        )

        # Calculate Composite Score
        pkt["composite_score"] = compute_composite_score(pkt, config)

        # Evaluate Filters
        pass_sasa = (
            pkt["mean_sasa"] <= config["pipeline"]["freesasa"]["max_mean_sasa_A2"]
        )
        pass_elec = (
            pkt.get("electrostatic_potential", 0.0)
            >= config["pipeline"]["apbs"]["min_electrostatic_kTe"]
        )
        pass_basic = (
            pkt["basic_residue_count"]
            >= config["pipeline"]["scoring"]["min_basic_residues"]
        )
        pass_plddt = (
            pkt["mean_plddt"]
            >= config["pipeline"]["plddt"]["min_plddt_for_pocket_residues"]
        )

        pkt["pass_all_filters"] = pass_sasa and pass_elec and pass_basic and pass_plddt

        # Depth is already calculated accurately in run_fpocket.py parsing

        scored_pockets.append(pkt)

    # Sort by composite score
    scored_pockets.sort(key=lambda x: x["composite_score"], reverse=True)
    return scored_pockets
