import json
from pathlib import Path
from typing import Dict, Any, Tuple
from .utils import setup_logger

logger = setup_logger(__name__)


def clean_structure_and_extract_plddt(
    pdb_path: Path, output_dir: Path
) -> Tuple[Path, Path]:
    """
    Cleans a PDB structure (removes HETATM except IPs) and extracts pLDDT sidecar.
    Returns paths to the cleaned structure and the pLDDT JSON file.
    """
    if not pdb_path.exists():
        raise FileNotFoundError(f"PDB file not found: {pdb_path}")

    cleaned_pdb_path = output_dir / f"{pdb_path.stem}_cleaned.pdb"
    plddt_path = output_dir / f"{pdb_path.stem}_plddt.json"

    plddt_scores: Dict[str, float] = {}

    logger.debug(f"Cleaning structure {pdb_path}")
    try:
        with open(pdb_path, "r") as fin, open(cleaned_pdb_path, "w") as fout:
            for line in fin:
                if line.startswith("ATOM"):
                    fout.write(line)
                    # Extract pLDDT from B-factor column
                    try:
                        res_id = line[22:26].strip()
                        chain = line[21]
                        b_factor = float(line[60:66].strip())
                        key = f"{chain}_{res_id}"
                        if key not in plddt_scores:
                            plddt_scores[key] = b_factor
                    except ValueError:
                        pass  # Ignore parsing errors for non-standard lines
                elif line.startswith("HETATM"):
                    # Keep IP ligands for validation/crystal structures
                    res_name = line[17:20].strip()
                    if res_name in ["IHP", "IP6", "IP4", "IP3"]:
                        fout.write(line)
                elif line.startswith("TER") or line.startswith("END"):
                    fout.write(line)

        with open(plddt_path, "w") as f_plddt:
            json.dump(plddt_scores, f_plddt, indent=2)

    except Exception as e:
        logger.error(f"Error preparing structure {pdb_path}: {e}")
        raise

    return cleaned_pdb_path, plddt_path
