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
        
        return pd.DataFrame(pockets)
    
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
