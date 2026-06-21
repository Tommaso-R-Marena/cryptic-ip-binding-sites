"""
Parser for fpocket output files.
"""

import re
from pathlib import Path
from typing import List, Dict
import pandas as pd


class FpocketParser:
    """
    Parse fpocket output files to extract pocket properties.
    """
    
    def parse_info_file(self, info_file: Path) -> pd.DataFrame:
        """
        Parse fpocket info file.
        
        Args:
            info_file: Path to _info.txt file
            
        Returns:
            DataFrame with pocket properties
        """
        pockets = []
        current_pocket = {}
        
        with open(info_file, 'r') as f:
            lines = f.readlines()
        
        for line in lines:
            line = line.strip()
            
            # New pocket starts
            if line.startswith('Pocket'):
                if current_pocket:
                    pockets.append(current_pocket)
                pocket_match = re.match(r'Pocket\s+(\d+)', line)
                current_pocket = {'pocket_id': int(pocket_match.group(1))}
            
            # Extract properties
            elif ':' in line and current_pocket:
                key, value = line.split(':', 1)
                key = key.strip().lower().replace(' ', '_').replace('-', '_')
                value = value.strip()
                
                # Try to convert to float
                try:
                    value = float(value)
                except ValueError:
                    pass
                
                current_pocket[key] = value
        
        # Add last pocket
        if current_pocket:
            pockets.append(current_pocket)

        frame = pd.DataFrame(pockets)
        if frame.empty:
            return frame

        pockets_dir = info_file.parent / "pockets"
        centers = []
        residue_sets = []
        for pocket_id in frame["pocket_id"]:
            pocket_file = pockets_dir / f"pocket{int(pocket_id)}_atm.pdb"
            centers.append(self._pocket_centroid(pocket_file))
            residue_sets.append(self._pocket_residue_ids_from_atoms(pocket_file))
        frame["center_x"] = [center[0] for center in centers]
        frame["center_y"] = [center[1] for center in centers]
        frame["center_z"] = [center[2] for center in centers]
        frame["fpocket_residue_ids"] = [",".join(str(r) for r in residues) for residues in residue_sets]

        return frame

    def _pocket_residue_ids_from_atoms(self, pocket_file: Path) -> List[int]:
        """Extract fpocket-reported residue numbers from pocket atom records."""
        if not pocket_file.exists():
            return []

        residue_ids: set[int] = set()
        with pocket_file.open("r", encoding="utf-8") as handle:
            for line in handle:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    residue_ids.add(int(line[22:26].strip()))
                except ValueError:
                    continue
        return sorted(residue_ids)

    def _pocket_centroid(self, pocket_file: Path) -> tuple[float, float, float]:
        """Compute alpha-sphere centroid for a pocket coordinates file."""
        if not pocket_file.exists():
            return (0.0, 0.0, 0.0)

        atoms = self.parse_pocket_atoms(pocket_file)
        if not atoms:
            return (0.0, 0.0, 0.0)

        xs = [atom["x"] for atom in atoms]
        ys = [atom["y"] for atom in atoms]
        zs = [atom["z"] for atom in atoms]
        return (float(sum(xs) / len(xs)), float(sum(ys) / len(ys)), float(sum(zs) / len(zs)))
    
    def parse_pocket_atoms(self, pocket_file: Path) -> List[Dict]:
        """
        Parse individual pocket PDB file to get alpha sphere coordinates.
        
        Args:
            pocket_file: Path to pocket PDB file
            
        Returns:
            List of atom dictionaries
        """
        atoms = []
        
        with open(pocket_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atom = {
                        'x': float(line[30:38]),
                        'y': float(line[38:46]),
                        'z': float(line[46:54]),
                        'atom_type': line[12:16].strip()
                    }
                    atoms.append(atom)
        
        return atoms
