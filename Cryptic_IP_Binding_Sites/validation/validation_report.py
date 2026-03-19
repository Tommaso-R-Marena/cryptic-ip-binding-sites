import os
import sys
import argparse
from pathlib import Path

# Add project root to python path
sys.path.append(str(Path(__file__).parent.parent))

from pipeline.utils import setup_logger, load_config, load_controls
from pipeline.prepare_structure import clean_structure_and_extract_plddt
from pipeline.run_fpocket import run_fpocket
from pipeline.run_freesasa import run_freesasa
from pipeline.run_apbs import run_pdb2pqr, run_apbs
from pipeline.score_pockets import score_and_filter_pockets

logger = setup_logger(__name__)

def process_protein(pdb_path: Path, output_dir: Path) -> list:
    """End-to-end processing of a single PDB."""
    logger.info(f"Processing {pdb_path.name}...")
    
    # 1. Clean & Extract pLDDT
    cleaned_pdb, plddt_json = clean_structure_and_extract_plddt(pdb_path, output_dir)
    
    # 2. Run fpocket
    pockets = run_fpocket(cleaned_pdb, output_dir)
    
    if not pockets:
        logger.warning(f"No valid pockets found for {pdb_path.name}")
        return []
        
    # 3. FreeSASA
    sasa_data = run_freesasa(cleaned_pdb, output_dir)
    
    # 4. APBS (Warning: computationally heavy per pocket)
    pqr_path = run_pdb2pqr(cleaned_pdb, output_dir)
    
    for pkt in pockets:
        if pkt.get("center"):
            val = run_apbs(pqr_path, output_dir, pkt["center"])
            pkt["electrostatic_potential"] = val
            
    # Clean up PQR
    if pqr_path.exists(): pqr_path.unlink()
    
    # 5-7. Score and Filter
    scored = score_and_filter_pockets(pockets, sasa_data, plddt_json, cleaned_pdb)
    
    return scored
    

def run_validation():
    """Runs all positive and negative controls."""
    val_dir = Path("data/structures/validation")
    out_dir = Path("data/results/validation")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    pos_ctls = load_controls("positive")
    neg_ctls = load_controls("negative")
    config = load_config()
    
    results = {}
    all_passed = True
    
    logger.info("Running Positive Controls")
    for name, data in pos_ctls.items():
        pdb_id = data.get("alphafold_id") or data["pdb_id"]
        # In validation dir, files are named by pdb_id.pdb or alphafold_id.pdb
        pdb_path = val_dir / f"{pdb_id}.pdb"
        
        if not pdb_path.exists():
            logger.error(f"Missing validation file: {pdb_path}")
            all_passed = False
            continue
            
        pockets = process_protein(pdb_path, out_dir)
        top_score = pockets[0]["composite_score"] if pockets else 0.0
        results[name] = {"type": "positive", "top_score": top_score, "pockets": len(pockets)}
        
        if name == "adar2":
            # Gating Criterion 1: ADAR2 must score well
            # For validation we expect it in top 3 pockets
            rank = -1
            for i, p in enumerate(pockets):
                 if p["pass_all_filters"]:
                     rank = i + 1
                     break
            if rank == -1 or rank > config["pipeline"]["validation"]["adar2_max_pocket_rank"]:
                 logger.error(f"ADAR2 failed gating! Rank: {rank}")
                 all_passed = False

    logger.info("Running Negative Controls")
    for name, data in neg_ctls.items():
        pdb_path = val_dir / f"{data['pdb_id']}.pdb"
        
        if not pdb_path.exists():
            logger.error(f"Missing validation file: {pdb_path}")
            all_passed = False
            continue
            
        pockets = process_protein(pdb_path, out_dir)
        top_score = pockets[0]["composite_score"] if pockets else 0.0
        results[name] = {"type": "negative", "top_score": top_score, "pockets": len(pockets)}
        
        if top_score > config["pipeline"]["validation"]["negative_control_max_score"]:
            logger.error(f"Negative control {name} scored too high: {top_score}")
            all_passed = False
            
    # Generate HTML report
    html = "<html><body><h1>Validation Report</h1><table border='1'><tr><th>Target</th><th>Type</th><th>Top Score</th></tr>"
    for name, res in results.items():
        html += f"<tr><td>{name}</td><td>{res['type']}</td><td>{res['top_score']}</td></tr>"
    html += f"</table><h2>Overall Status: {'PASSED' if all_passed else 'FAILED'}</h2></body></html>"
    
    with open(out_dir / "validation_report.html", "w") as f:
        f.write(html)
        
    if not all_passed:
        logger.error("NO-GO: Validation criteria not met.")
        sys.exit(1)
    else:
        logger.info("GO: Validation passed successfully.")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", default="validation", choices=["validation"], help="Mode to run in")
    args = parser.parse_args()
    
    if args.mode == "validation":
        run_validation()
