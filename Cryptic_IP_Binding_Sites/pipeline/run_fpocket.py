import subprocess
import shutil
from pathlib import Path
from typing import Dict, Any, List, Tuple
from .utils import setup_logger, load_config

logger = setup_logger(__name__)


def parse_fpocket_output(
    info_file: Path, config: Dict[str, Any]
) -> List[Dict[str, Any]]:
    """Parses fpocket info.txt output and extracts/filters pockets."""
    pockets = []
    current_pocket: Dict[str, Any] = {}

    min_vol = config["pipeline"]["fpocket"]["min_volume_A3"]
    max_vol = config["pipeline"]["fpocket"]["max_volume_A3"]

    try:
        with open(info_file, "r") as f:
            for line in f:
                line = line.strip()
                if line.startswith("Pocket"):
                    if current_pocket:
                        # Process previous pocket
                        if min_vol <= current_pocket.get("volume", 0) <= max_vol:
                            pockets.append(current_pocket)
                    try:
                        current_pocket = {"id": int(line.split()[1])}
                    except (IndexError, ValueError):
                        current_pocket = {"id": -1}
                elif "Score :" in line:
                    current_pocket["score"] = float(line.split(":")[1].strip())
                elif "Druggability Score :" in line:
                    current_pocket["druggability_score"] = float(
                        line.split(":")[1].strip()
                    )
                elif "Volume :" in line:
                    current_pocket["volume"] = float(line.split(":")[1].strip())
                elif "Total SASA :" in line:
                    current_pocket["sasa"] = float(line.split(":")[1].strip())
                # Extract coordinates could be parsed from pocket PDBs later

        # Handle the last pocket
        if current_pocket and (min_vol <= current_pocket.get("volume", 0) <= max_vol):
            pockets.append(current_pocket)

    except Exception as e:
        logger.error(f"Error parsing fpocket output {info_file}: {e}")

    return pockets


def get_pocket_center_and_depth(
    pocket_pdb: Path,
) -> Tuple[Tuple[float, float, float], float]:
    """Calculates the center of mass of a pocket and its maximum depth from its fpocket PDB."""
    x_sum, y_sum, z_sum, count = 0.0, 0.0, 0.0, 0
    atoms = []
    import math

    with open(pocket_pdb, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    x_sum += x
                    y_sum += y
                    z_sum += z
                    count += 1
                    atoms.append((x, y, z))
                except ValueError:
                    pass
    if count == 0:
        return ((0.0, 0.0, 0.0), 0.0)

    cx = x_sum / count
    cy = y_sum / count
    cz = z_sum / count

    max_dist = 0.0
    for ax, ay, az in atoms:
        dist = math.sqrt((ax - cx) ** 2 + (ay - cy) ** 2 + (az - cz) ** 2)
        if dist > max_dist:
            max_dist = dist

    return ((cx, cy, cz), max_dist)


def get_pocket_residues(pocket_pdb: Path) -> List[str]:
    """Extracts lining residues from a pocket PDB."""
    residues = set()
    with open(pocket_pdb, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                res_id = line[22:26].strip()
                chain = line[21]
                res_name = line[17:20].strip()
                residues.add(f"{res_name}_{chain}_{res_id}")
    return list(residues)


def run_fpocket(pdb_path: Path, output_dir: Path) -> List[Dict[str, Any]]:
    """Runs fpocket wrapper and returns extracted pockets."""
    config = load_config()

    # Check if fpocket exists
    if not shutil.which("fpocket"):
        logger.error("fpocket command not found. Ensure it is installed and in PATH.")
        raise FileNotFoundError("fpocket not found")

    try:
        # fpocket generates output in the same directory as the input PDB
        # We'll copy the PDB to a temporary work dir in output_dir
        work_pdb = output_dir / pdb_path.name
        shutil.copy(pdb_path, work_pdb)

        # Run fpocket
        cmd = ["fpocket", "-f", str(work_pdb)]
        subprocess.run(
            cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300
        )

        # Output directory is created by fpocket: <pdbname>_out
        fpocket_out_dir = output_dir / f"{pdb_path.stem}_out"
        info_file = fpocket_out_dir / f"{pdb_path.stem}_info.txt"
        pocket_pdb_dir = fpocket_out_dir / "pockets"

        if not info_file.exists():
            raise FileNotFoundError(f"fpocket info file not generated: {info_file}")

        pockets = parse_fpocket_output(info_file, config)

        # Augment with centers, depth, and residues
        for pkt in pockets:
            pocket_pdb = pocket_pdb_dir / f"pocket{pkt['id']}_atm.pdb"
            if pocket_pdb.exists():
                center, depth = get_pocket_center_and_depth(pocket_pdb)
                pkt["center"] = center
                pkt["depth"] = depth
                pkt["residues"] = get_pocket_residues(pocket_pdb)
            else:
                pkt["center"] = (0.0, 0.0, 0.0)
                pkt["depth"] = 0.0
                pkt["residues"] = []

        # Clean up
        if fpocket_out_dir.exists():
            shutil.rmtree(fpocket_out_dir)
        if work_pdb.exists():
            work_pdb.unlink()

        return pockets

    except subprocess.TimeoutExpired:
        logger.error(f"fpocket timed out for {pdb_path}")
        return []
    except subprocess.CalledProcessError as e:
        logger.error(f"fpocket execution failed for {pdb_path}: {e.stderr.decode()}")
        return []
    except Exception as e:
        logger.error(f"Unexpected error running fpocket on {pdb_path}: {e}")
        return []
