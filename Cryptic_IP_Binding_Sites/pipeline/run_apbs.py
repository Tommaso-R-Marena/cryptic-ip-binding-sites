import subprocess
import shutil
import math
from pathlib import Path
from typing import Tuple
from .utils import setup_logger, load_config

logger = setup_logger(__name__)


def parse_apbs_output(dx_file: Path, target_coord: Tuple[float, float, float]) -> float:
    """Finds the electrostatic potential at the nearest grid point to target_coord exactly."""
    if not dx_file.exists():
        return 0.0

    import numpy as np

    grid = []
    origin = [0.0, 0.0, 0.0]
    delta = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    counts = [0, 0, 0]

    try:
        with open(dx_file, "r") as f:
            lines = f.readlines()

        data_started = False
        for line in lines:
            line = line.strip()
            if line.startswith("object 1 class gridpositions"):
                parts = line.split()
                counts = [int(parts[-3]), int(parts[-2]), int(parts[-1])]
            elif line.startswith("origin"):
                parts = line.split()
                origin = [float(parts[-3]), float(parts[-2]), float(parts[-1])]
            elif line.startswith("delta"):
                parts = line.split()
                dim = -1
                if delta[0] == [0.0, 0.0, 0.0]:
                    dim = 0
                elif delta[1] == [0.0, 0.0, 0.0]:
                    dim = 1
                elif delta[2] == [0.0, 0.0, 0.0]:
                    dim = 2
                if dim != -1:
                    delta[dim] = [float(parts[-3]), float(parts[-2]), float(parts[-1])]
            elif line.startswith("object 3 class array"):
                data_started = True
            elif data_started:
                if line.startswith("attribute") or line.startswith("object"):
                    break  # End of data array
                grid.extend([float(x) for x in line.split()])

        if len(grid) == 0:
            return 0.0

        # The typical DataExplorer flat output iterates Z fastest, then Y, then X.
        # So we reshape to (nx, ny, nz)
        grid_array = np.array(grid).reshape((counts[0], counts[1], counts[2]))

        # Identify the exact grid point
        x_idx = int(round((target_coord[0] - origin[0]) / delta[0][0]))
        y_idx = int(round((target_coord[1] - origin[1]) / delta[1][1]))
        z_idx = int(round((target_coord[2] - origin[2]) / delta[2][2]))

        # Enforce exact bounds
        x_idx = max(0, min(x_idx, counts[0] - 1))
        y_idx = max(0, min(y_idx, counts[1] - 1))
        z_idx = max(0, min(z_idx, counts[2] - 1))

        return float(grid_array[x_idx, y_idx, z_idx])

    except Exception as e:
        logger.error(f"Error parsing APBS dx file {dx_file}: {e}")
        return 0.0


def run_pdb2pqr(pdb_path: Path, output_dir: Path) -> Path:
    """Converts PDB to PQR using pdb2pqr."""
    config = load_config()
    ff = config["pipeline"]["apbs"]["forcefield"]
    ph = config["pipeline"]["apbs"]["pH"]

    if not shutil.which("pdb2pqr"):
        logger.error("pdb2pqr not found.")
        raise FileNotFoundError("pdb2pqr")

    pqr_path = output_dir / f"{pdb_path.stem}.pqr"
    cmd = [
        "pdb2pqr",
        f"--ff={ff}",
        f"--titration-state-method=propka",
        f"--with-ph={ph}",
        "--nodebump",
        str(pdb_path),
        str(pqr_path),
    ]

    subprocess.run(
        cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120
    )
    return pqr_path


def run_apbs(
    pqr_path: Path, output_dir: Path, target_coord: Tuple[float, float, float]
) -> float:
    """Runs APBS and queries the resulting potential map at the target coordinate."""
    if not shutil.which("apbs"):
        logger.error("apbs not found.")
        raise FileNotFoundError("apbs")

    # Generate an APBS input (.in) file
    in_path = output_dir / f"{pqr_path.stem}.in"
    dx_path = output_dir / f"{pqr_path.stem}"  # APBS will append .dx

    # Simple APBS config writing
    apbs_input = f"""
read
    mol pqr {pqr_path.name}
end
elec
    mg-auto
    dime 65 65 65
    cglen 40 40 40
    fglen 20 20 20
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    chgm spl2
    srfm smol
    srad 1.4
    swin 0.3
    sdens 10.0
    temp 298.0
    calcenergy total
    calcforce no
    write pot dx {dx_path.name}
end
quit
"""
    with open(in_path, "w") as f:
        f.write(apbs_input)

    try:
        # Run APBS in the output dir so it finds the PQR easily
        subprocess.run(
            ["apbs", in_path.name],
            cwd=output_dir,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=300,
        )

        actual_dx_path = output_dir / f"{pqr_path.stem}.dx"
        val = parse_apbs_output(actual_dx_path, target_coord)

        # Clean up
        if actual_dx_path.exists():
            actual_dx_path.unlink()
        if in_path.exists():
            in_path.unlink()

        return val

    except subprocess.TimeoutExpired:
        logger.error(f"APBS timed out for {pqr_path}")
        return 0.0
    except subprocess.CalledProcessError as e:
        logger.error(f"APBS failed: {e.stderr.decode()}")
        return 0.0
    except Exception as e:
        logger.error(f"Unexpected error in APBS: {e}")
        return 0.0
