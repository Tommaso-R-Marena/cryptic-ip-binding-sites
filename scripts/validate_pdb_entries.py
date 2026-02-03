#!/usr/bin/env python3
"""Validate PDB entries exist and are accessible."""

import sys
import argparse
import logging
from typing import List

from cryptic_ip.database.pdb_client import PDBClient

logger = logging.getLogger(__name__)


def validate_pdb_entry(pdb_id: str, client: PDBClient) -> dict:
    """Validate single PDB entry.
    
    Args:
        pdb_id: PDB identifier
        client: PDB client instance
        
    Returns:
        Validation result dictionary
    """
    try:
        # Fetch entry info
        info = client.get_entry_info(pdb_id)
        
        # Try to download structure
        structure_path = client.fetch_structure(pdb_id)
        
        return {
            'pdb_id': pdb_id,
            'valid': True,
            'title': info['title'],
            'resolution': info['resolution'],
            'method': info['method'],
            'deposition_date': info['deposition_date'],
            'file_size': structure_path.stat().st_size,
            'error': None
        }
    except Exception as e:
        return {
            'pdb_id': pdb_id,
            'valid': False,
            'error': str(e)
        }


def main():
    parser = argparse.ArgumentParser(
        description="Validate PDB entries"
    )
    parser.add_argument(
        '--pdb-ids',
        required=True,
        help='Comma-separated PDB IDs'
    )
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    
    pdb_ids = [pid.strip().upper() for pid in args.pdb_ids.split(',')]
    
    print("\n" + "="*60)
    print("  PDB ENTRY VALIDATION")
    print("="*60 + "\n")
    
    client = PDBClient()
    all_valid = True
    
    for pdb_id in pdb_ids:
        print(f"Validating {pdb_id}...")
        result = validate_pdb_entry(pdb_id, client)
        
        if result['valid']:
            print(f"  ✓ Valid")
            print(f"    Title: {result['title'][:60]}...")
            print(f"    Resolution: {result['resolution']} Å")
            print(f"    Method: {result['method']}")
        else:
            print(f"  ✗ Invalid: {result['error']}")
            all_valid = False
        
        print()
    
    print("="*60 + "\n")
    
    if not all_valid:
        sys.exit(1)


if __name__ == '__main__':
    main()
