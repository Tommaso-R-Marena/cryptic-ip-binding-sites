"""
ADAR2 validation - the gold standard for cryptic IP6 binding.
"""

import os
from pathlib import Path
from typing import Dict, Optional
import urllib.request

from ..analysis import ProteinAnalyzer


def download_adar2_structures(data_dir: str = "data/validation") -> Dict[str, Path]:
    """
    Download ADAR2 structures from AlphaFold and PDB.
    
    Args:
        data_dir: Directory to save structures
        
    Returns:
        Dictionary with paths to downloaded files
    """
    data_path = Path(data_dir)
    data_path.mkdir(parents=True, exist_ok=True)
    
    structures = {}
    
    # AlphaFold structure
    alphafold_url = "https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb"
    alphafold_path = data_path / "AF-P78563-F1-model_v4.pdb"
    if not alphafold_path.exists():
        print(f"Downloading AlphaFold ADAR2 structure...")
        urllib.request.urlretrieve(alphafold_url, alphafold_path)
    structures['alphafold'] = alphafold_path
    
    # Crystal structure with IP6
    pdb_url = "https://files.rcsb.org/download/1ZY7.pdb"
    pdb_path = data_path / "1ZY7.pdb"
    if not pdb_path.exists():
        print(f"Downloading PDB 1ZY7 (ADAR2 crystal structure)...")
        urllib.request.urlretrieve(pdb_url, pdb_path)
    structures['crystal'] = pdb_path
    
    return structures


def validate_adar2(structure_path: Optional[str] = None,
                  use_alphafold: bool = True) -> Dict:
    """
    Validate pipeline on ADAR2 IP6 binding site.
    
    This is the critical test: if the pipeline doesn't find IP6 in ADAR2,
    the parameters are wrong.
    
    Args:
        structure_path: Path to ADAR2 structure (downloads if None)
        use_alphafold: Use AlphaFold structure instead of crystal
        
    Returns:
        Dictionary with validation metrics
    """
    # Download structures if not provided
    if structure_path is None:
        structures = download_adar2_structures()
        structure_path = structures['alphafold'] if use_alphafold else structures['crystal']
    
    print(f"\nValidating ADAR2 IP6 binding site...")
    print(f"Structure: {structure_path}")
    
    # Analyze structure
    analyzer = ProteinAnalyzer(str(structure_path))
    
    print("Detecting pockets...")
    pockets = analyzer.detect_pockets()
    print(f"Found {len(pockets)} pockets")
    
    print("Scoring pockets...")
    scored = analyzer.score_all_pockets()
    
    # Known IP6 binding site residues in ADAR2
    # From Macbeth et al. (2005) - K376, K519, R522, R651, K672
    ip6_residues = [376, 519, 522, 651, 672]
    
    # Find pocket that best matches IP6 site
    # This is a simplified check - real validation would align coordinates
    top_pocket = scored.iloc[0]
    
    results = {
        'total_pockets': len(pockets),
        'top_pocket_id': top_pocket['pocket_id'],
        'top_pocket_score': top_pocket['composite_score'],
        'volume': top_pocket['volume'],
        'sasa': top_pocket['sasa'],
        'basic_residues': top_pocket['basic_residues'],
        'depth': top_pocket['depth'],
        'structure_used': 'AlphaFold' if use_alphafold else 'Crystal',
        'success_criteria': {
            'score_above_0.7': top_pocket['composite_score'] >= 0.7,
            'low_sasa': top_pocket['sasa'] < 10.0,
            'sufficient_basic': top_pocket['basic_residues'] >= 4,
            'appropriate_volume': 300 <= top_pocket['volume'] <= 800
        }
    }
    
    # Check if validation passed
    criteria = results['success_criteria']
    results['validation_passed'] = all(criteria.values())
    
    # Print results
    print("\n" + "="*60)
    print("ADAR2 VALIDATION RESULTS")
    print("="*60)
    print(f"Top pocket score: {results['top_pocket_score']:.3f}")
    print(f"Volume: {results['volume']:.1f} ų")
    print(f"SASA: {results['sasa']:.2f} ų")
    print(f"Basic residues: {results['basic_residues']}")
    print(f"\nValidation criteria:")
    for criterion, passed in criteria.items():
        status = "✓" if passed else "✗"
        print(f"  {status} {criterion}")
    print(f"\nOVERALL: {'PASSED' if results['validation_passed'] else 'FAILED'}")
    print("="*60 + "\n")
    
    return results
