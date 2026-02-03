"""
Protein structure analysis and pocket detection.
"""

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select
from prody import parsePDB, calcSASA

from .fpocket_parser import FpocketParser
from .scorer import PocketScorer


class ProteinAnalyzer:
    """
    Main analyzer class for detecting and characterizing pockets in protein structures.
    
    Integrates fpocket for pocket detection, FreeSASA for accessibility,
    and APBS for electrostatics.
    """
    
    def __init__(self, pdb_path: str, work_dir: Optional[str] = None):
        """
        Initialize analyzer with a protein structure.
        
        Args:
            pdb_path: Path to PDB file
            work_dir: Working directory for outputs (temp if None)
        """
        self.pdb_path = Path(pdb_path)
        self.work_dir = Path(work_dir) if work_dir else Path(tempfile.mkdtemp())
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # Parse structure
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', str(self.pdb_path))
        
        # Initialize components
        self.fpocket_parser = FpocketParser()
        self.scorer = PocketScorer()
        
        # Storage for results
        self.pockets = None
        self.sasa_data = None
        self.electrostatic_data = None
        
    def detect_pockets(self, min_alpha_sphere: int = 3) -> pd.DataFrame:
        """
        Detect pockets using fpocket.
        
        Args:
            min_alpha_sphere: Minimum number of alpha spheres for a pocket
            
        Returns:
            DataFrame with pocket properties
        """
        # Check if fpocket is available
        try:
            subprocess.run(['fpocket', '-h'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(
                "fpocket not found. Please install from https://github.com/Discngine/fpocket"
            )
        
        # Run fpocket
        output_dir = self.work_dir / f"{self.pdb_path.stem}_out"
        cmd = [
            'fpocket',
            '-f', str(self.pdb_path),
            '-m', str(min_alpha_sphere)
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(self.work_dir))
        
        if result.returncode != 0:
            raise RuntimeError(f"fpocket failed: {result.stderr}")
        
        # Parse fpocket output
        info_file = output_dir / f"{self.pdb_path.stem}_info.txt"
        if not info_file.exists():
            raise RuntimeError(f"fpocket output not found: {info_file}")
        
        self.pockets = self.fpocket_parser.parse_info_file(info_file)
        self.pockets['pdb_path'] = str(self.pdb_path)
        
        return self.pockets
    
    def calculate_sasa(self) -> Dict[int, float]:
        """
        Calculate solvent accessible surface area for all residues.
        
        Returns:
            Dictionary mapping residue number to SASA value
        """
        try:
            # Use ProDy for SASA calculation
            structure = parsePDB(str(self.pdb_path))
            sasa = calcSASA(structure)
            
            # Get per-residue SASA
            residue_sasa = {}
            for residue in structure.iterResidues():
                resnum = residue.getResnum()
                atom_indices = residue.getIndices()
                residue_sasa[resnum] = np.sum(sasa[atom_indices])
            
            self.sasa_data = residue_sasa
            return residue_sasa
            
        except Exception as e:
            raise RuntimeError(f"SASA calculation failed: {e}")
    
    def calculate_electrostatics(self) -> Optional[np.ndarray]:
        """
        Calculate electrostatic potential using APBS.
        
        Returns:
            Electrostatic potential grid (placeholder - requires APBS installation)
        """
        # This is a placeholder - full APBS integration requires
        # pdb2pqr conversion and APBS execution
        print("Warning: Electrostatic calculation requires APBS installation")
        print("See: http://www.poissonboltzmann.org/")
        return None
    
    def get_pocket_residues(self, pocket_id: int, distance_cutoff: float = 5.0) -> List[int]:
        """
        Get residue numbers within distance of pocket center.
        
        Args:
            pocket_id: Pocket identifier
            distance_cutoff: Distance threshold in Angstroms
            
        Returns:
            List of residue numbers
        """
        if self.pockets is None:
            raise ValueError("Run detect_pockets() first")
        
        pocket = self.pockets[self.pockets['pocket_id'] == pocket_id].iloc[0]
        center = np.array([pocket['center_x'], pocket['center_y'], pocket['center_z']])
        
        # Get all CA atoms
        residue_numbers = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if 'CA' in residue:
                        ca_coord = residue['CA'].get_coord()
                        distance = np.linalg.norm(ca_coord - center)
                        if distance <= distance_cutoff:
                            residue_numbers.append(residue.id[1])
        
        return residue_numbers
    
    def count_basic_residues(self, pocket_id: int, distance_cutoff: float = 5.0) -> int:
        """
        Count basic residues (Arg, Lys, His) near pocket.
        
        Args:
            pocket_id: Pocket identifier
            distance_cutoff: Distance threshold in Angstroms
            
        Returns:
            Number of basic residues
        """
        pocket_residues = self.get_pocket_residues(pocket_id, distance_cutoff)
        basic_count = 0
        
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.id[1] in pocket_residues:
                        if residue.get_resname() in ['ARG', 'LYS', 'HIS']:
                            basic_count += 1
        
        return basic_count
    
    def analyze_pocket(self, pocket_id: int) -> Dict:
        """
        Complete analysis of a single pocket.
        
        Args:
            pocket_id: Pocket identifier
            
        Returns:
            Dictionary with all pocket metrics
        """
        if self.pockets is None:
            self.detect_pockets()
        
        if self.sasa_data is None:
            self.calculate_sasa()
        
        pocket = self.pockets[self.pockets['pocket_id'] == pocket_id].iloc[0]
        pocket_residues = self.get_pocket_residues(pocket_id)
        
        # Calculate average SASA for pocket residues
        pocket_sasa = np.mean([self.sasa_data.get(r, 0) for r in pocket_residues])
        
        # Count basic residues
        basic_count = self.count_basic_residues(pocket_id)
        
        return {
            'pocket_id': pocket_id,
            'volume': pocket.get('volume', 0),
            'depth': pocket.get('mean_local_hydrophobic_density', 0),
            'sasa': pocket_sasa,
            'basic_residues': basic_count,
            'residue_count': len(pocket_residues),
            'center': (pocket['center_x'], pocket['center_y'], pocket['center_z'])
        }
    
    def score_all_pockets(self) -> pd.DataFrame:
        """
        Score all detected pockets for cryptic IP binding.
        
        Returns:
            DataFrame with scores for all pockets
        """
        if self.pockets is None:
            self.detect_pockets()
        
        results = []
        for pocket_id in self.pockets['pocket_id']:
            try:
                analysis = self.analyze_pocket(pocket_id)
                score = self.scorer.calculate_composite_score(
                    volume=analysis['volume'],
                    depth=analysis['depth'],
                    sasa=analysis['sasa'],
                    basic_count=analysis['basic_residues']
                )
                analysis['composite_score'] = score
                results.append(analysis)
            except Exception as e:
                print(f"Warning: Failed to analyze pocket {pocket_id}: {e}")
                continue
        
        return pd.DataFrame(results).sort_values('composite_score', ascending=False)
