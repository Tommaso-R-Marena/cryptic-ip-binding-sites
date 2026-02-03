#!/usr/bin/env python3
"""Validate structural data integrity and metadata."""

import sys
import json
import logging
from pathlib import Path
from typing import Dict, List
import argparse

try:
    from Bio import PDB
except ImportError:
    print("Error: Biopython required. Install with: pip install biopython")
    sys.exit(1)

from cryptic_ip.database.alphafold_client import AlphaFoldClient
from cryptic_ip.database.pdb_client import PDBClient
from cryptic_ip.database.uniprot_client import UniProtClient

logger = logging.getLogger(__name__)


def validate_protein_structure(pdb_file: Path) -> Dict:
    """Validate PDB file structure and extract metrics.
    
    Args:
        pdb_file: Path to PDB file
        
    Returns:
        Dictionary with validation results
    """
    parser = PDB.PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('protein', pdb_file)
        model = structure[0]
        
        # Count residues
        n_residues = sum(1 for _ in model.get_residues() if PDB.is_aa(_))
        
        # Count atoms
        n_atoms = sum(1 for _ in model.get_atoms())
        
        # Check for missing residues
        residue_numbers = [r.id[1] for r in model.get_residues() if PDB.is_aa(r)]
        missing_residues = []
        if residue_numbers:
            for i in range(min(residue_numbers), max(residue_numbers)):
                if i not in residue_numbers:
                    missing_residues.append(i)
        
        # Extract b-factors (pLDDT for AlphaFold)
        b_factors = [atom.bfactor for atom in model.get_atoms()]
        avg_bfactor = sum(b_factors) / len(b_factors) if b_factors else 0.0
        
        return {
            'valid': True,
            'n_residues': n_residues,
            'n_atoms': n_atoms,
            'missing_residues': len(missing_residues),
            'avg_bfactor': round(avg_bfactor, 2),
            'file_size': pdb_file.stat().st_size,
            'error': None
        }
        
    except Exception as e:
        return {
            'valid': False,
            'error': str(e)
        }


def validate_adar2_structures() -> Dict:
    """Validate ADAR2 reference structures.
    
    Downloads and validates:
    - AlphaFold prediction (P78563)
    - Crystal structure (1ZY7)
    
    Returns:
        Validation report
    """
    report = {
        'protein': 'ADAR2',
        'uniprot_id': 'P78563',
        'structures': {}
    }
    
    # Download and validate AlphaFold structure
    logger.info("Validating ADAR2 AlphaFold structure...")
    af_client = AlphaFoldClient()
    try:
        af_structure = af_client.fetch_structure('P78563')
        af_validation = validate_protein_structure(af_structure)
        
        # Get metadata
        af_metadata = af_client.get_metadata('P78563')
        
        report['structures']['alphafold'] = {
            'file': str(af_structure),
            'validation': af_validation,
            'metadata': af_metadata
        }
        
        logger.info(f"✓ AlphaFold: {af_validation['n_residues']} residues, "
                   f"pLDDT={af_validation['avg_bfactor']:.1f}")
    except Exception as e:
        logger.error(f"✗ AlphaFold validation failed: {e}")
        report['structures']['alphafold'] = {'error': str(e)}
    
    # Download and validate crystal structure
    logger.info("Validating ADAR2 crystal structure (1ZY7)...")
    pdb_client = PDBClient()
    try:
        pdb_structure = pdb_client.fetch_structure('1ZY7')
        pdb_validation = validate_protein_structure(pdb_structure)
        
        # Get metadata
        pdb_info = pdb_client.get_entry_info('1ZY7')
        
        report['structures']['crystal'] = {
            'pdb_id': '1ZY7',
            'file': str(pdb_structure),
            'validation': pdb_validation,
            'metadata': pdb_info
        }
        
        logger.info(f"✓ Crystal: {pdb_validation['n_residues']} residues, "
                   f"resolution={pdb_info['resolution']} Å")
    except Exception as e:
        logger.error(f"✗ Crystal structure validation failed: {e}")
        report['structures']['crystal'] = {'error': str(e)}
    
    # Get UniProt annotations
    logger.info("Fetching UniProt annotations...")
    uniprot_client = UniProtClient()
    try:
        protein_info = uniprot_client.get_protein_info('P78563')
        report['uniprot'] = protein_info
        logger.info(f"✓ UniProt: {protein_info['gene_name']} "
                   f"({protein_info['sequence_length']} aa)")
    except Exception as e:
        logger.error(f"✗ UniProt fetch failed: {e}")
        report['uniprot'] = {'error': str(e)}
    
    # Validation summary
    report['summary'] = {
        'alphafold_valid': report['structures'].get('alphafold', {}).get('validation', {}).get('valid', False),
        'crystal_valid': report['structures'].get('crystal', {}).get('validation', {}).get('valid', False),
        'uniprot_valid': 'error' not in report.get('uniprot', {}),
        'all_valid': False
    }
    
    report['summary']['all_valid'] = all([
        report['summary']['alphafold_valid'],
        report['summary']['crystal_valid'],
        report['summary']['uniprot_valid']
    ])
    
    return report


def main():
    parser = argparse.ArgumentParser(
        description="Validate structural data integrity"
    )
    parser.add_argument(
        '--proteins',
        nargs='+',
        help='UniProt IDs to validate (default: P78563 for ADAR2)'
    )
    parser.add_argument(
        '--output',
        type=Path,
        default='validation_report.json',
        help='Output file for validation report'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Verbose output'
    )
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Run validation
    proteins = args.proteins or ['P78563']
    
    print("\n" + "="*60)
    print("  STRUCTURAL DATA VALIDATION")
    print("="*60 + "\n")
    
    all_reports = []
    
    for uniprot_id in proteins:
        print(f"\nValidating {uniprot_id}...")
        
        if uniprot_id == 'P78563':  # ADAR2 - full validation
            report = validate_adar2_structures()
        else:
            # Basic validation for other proteins
            report = {'uniprot_id': uniprot_id, 'error': 'Not implemented'}
        
        all_reports.append(report)
        
        # Print summary
        if 'summary' in report:
            summary = report['summary']
            status = "✓ PASS" if summary['all_valid'] else "✗ FAIL"
            print(f"\n{status}: AlphaFold={summary['alphafold_valid']}, "
                  f"Crystal={summary['crystal_valid']}, "
                  f"UniProt={summary['uniprot_valid']}")
    
    # Save report
    with open(args.output, 'w') as f:
        json.dump({
            'validation_date': pd.Timestamp.now().isoformat(),
            'proteins': all_reports
        }, f, indent=2)
    
    print(f"\nValidation report saved to: {args.output}")
    print("="*60 + "\n")
    
    # Exit with error if any validation failed
    if not all(r.get('summary', {}).get('all_valid', False) for r in all_reports):
        sys.exit(1)


if __name__ == '__main__':
    import pandas as pd
    main()
