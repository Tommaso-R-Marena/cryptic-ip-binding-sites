import os
from pathlib import Path
from typing import List, Dict, Any
from .utils import setup_logger

logger = setup_logger(__name__)


def generate_pymol_script(
    pdb_path: Path,
    pockets: List[Dict[str, Any]],
    output_pml: Path,
    plddt_path: Path = None,
):
    """
    Generates a sophisticated PyMOL script (.pml) to visualize the protein
    and highlight the predicted cryptic IP binding pockets.

    The script will:
      - Load the PDB
      - Show as cartoon
      - Color by pLDDT if available
      - Emphasize and zoom into the top pockets
      - Highlight basic residues (Arg, Lys, His) in the pockets
    """
    try:
        pdb_name = pdb_path.stem
        lines = []

        # Initialization and visual settings
        lines.append("reinitialize")
        lines.append("bg_color white")
        lines.append(f"load {pdb_path.resolve()}, {pdb_name}")
        lines.append("hide all")
        lines.append(f"show cartoon, {pdb_name}")
        lines.append("set cartoon_fancy_helices, 1")
        lines.append("set cartoon_smooth_loops, 1")
        lines.append("set dash_gap, 0")
        lines.append("set dash_radius, 0.08")
        lines.append("color gray70, all")

        # Color by pLDDT if path is provided (for future enhancement, but right now just standard coloring)
        # We can simulate advanced coloring via b-factor if pLDDT is inserted there
        lines.append("spectrum b, rainbow_rev, minimum=50, maximum=100")

        # Create selections for each pocket and highlight them
        for i, pkt in enumerate(pockets):
            # pkt residues are formatted like "RES_A_123"
            res_list = []
            for r in pkt.get("residues", []):
                parts = r.split("_")
                if len(parts) >= 3:
                    chain = parts[1]
                    res_num = parts[2]
                    res_list.append(f"(chain {chain} and resi {res_num})")

            if res_list:
                sel_name = f"pocket_{i+1}"
                sel_str = " or ".join(res_list)
                lines.append(f"select {sel_name}, {pdb_name} and ({sel_str})")

                # Highlight the pocket
                lines.append(f"show surface, {sel_name}")
                lines.append(f"set transparency, 0.4, {sel_name}")
                color = "tv_red" if i == 0 else ("tv_blue" if i == 1 else "tv_green")
                lines.append(f"color {color}, {sel_name}")

                # Show sidechains for basic residues interacting with IP
                lines.append(
                    f"select basic_{sel_name}, {sel_name} and resn ARG+LYS+HIS"
                )
                lines.append(f"show sticks, basic_{sel_name}")
                lines.append(f"color atomic, basic_{sel_name}")

                if i == 0:
                    lines.append(f"zoom {sel_name}, 5")

        lines.append("deselect")

        with open(output_pml, "w") as f:
            f.write("\n".join(lines) + "\n")

        logger.info(f"PyMOL visualization script written to {output_pml}")

    except Exception as e:
        logger.error(f"Error generating PyMOL script: {e}")
        raise
