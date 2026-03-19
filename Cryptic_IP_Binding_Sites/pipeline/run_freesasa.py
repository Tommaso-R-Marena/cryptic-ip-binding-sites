import subprocess
import shutil
from pathlib import Path
from typing import Dict, Any, List
from .utils import setup_logger, load_config

logger = setup_logger(__name__)


def parse_freesasa_output(rsa_file: Path) -> Dict[str, float]:
    """Parses FreeSASA RSA output to extract absolute sidechain SASA for Arg, Lys, His."""
    sasa_data: Dict[str, float] = {}

    with open(rsa_file, "r") as f:
        for line in f:
            if line.startswith("RES"):
                # NACCESS/FreeSASA RSA format:
                # RES NAM C SEQ    ABS_ALL REL_ALL ABS_SIDE REL_SIDE ABS_MAIN REL_MAIN ABS_NONP REL_NONP ABS_POL REL_POL
                # 012345678901234567890123456789012345
                res_name = line[4:7].strip()
                if res_name in ["ARG", "LYS", "HIS"]:
                    chain = line[8:9].strip()
                    res_id = line[9:13].strip()
                    try:
                        abs_side = float(line[29:35].strip())
                        key = f"{res_name}_{chain}_{res_id}"
                        sasa_data[key] = abs_side
                    except ValueError:
                        pass
    return sasa_data


def run_freesasa(pdb_path: Path, output_dir: Path) -> Dict[str, float]:
    """Runs FreeSASA to compute per-residue SASA."""
    config = load_config()
    probe_radius = config["pipeline"]["freesasa"]["probe_radius"]

    if not shutil.which("freesasa"):
        logger.error("freesasa command not found.")
        raise FileNotFoundError("freesasa not found")

    out_rsa = output_dir / f"{pdb_path.stem}.rsa"

    cmd = ["freesasa", "--format=rsa", f"--probe-radius={probe_radius}", str(pdb_path)]

    try:
        with open(out_rsa, "w") as f_out:
            subprocess.run(
                cmd, check=True, stdout=f_out, stderr=subprocess.PIPE, timeout=120
            )

        return parse_freesasa_output(out_rsa)

    except subprocess.TimeoutExpired:
        logger.error(f"FreeSASA timed out for {pdb_path}")
        return {}
    except subprocess.CalledProcessError as e:
        logger.error(f"FreeSASA execution failed: {e.stderr.decode()}")
        return {}
    except Exception as e:
        logger.error(f"Unexpected error in FreeSASA: {e}")
        return {}
