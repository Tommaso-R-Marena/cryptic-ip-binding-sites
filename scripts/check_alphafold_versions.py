#!/usr/bin/env python3
"""Check AlphaFold database versions and freshness."""

import sys
import argparse
import logging
from datetime import datetime
import requests

logger = logging.getLogger(__name__)

PROTEOME_IDS = {
    'yeast': 'UP000002311',
    'human': 'UP000005640',
    'dictyostelium': 'UP000002195'
}


def check_proteome_availability(proteome_id: str) -> dict:
    """Check if proteome is available in AlphaFold DB.
    
    Args:
        proteome_id: UniProt proteome ID
        
    Returns:
        Dictionary with availability info
    """
    # Check FTP listing
    ftp_url = f"https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"
    
    try:
        response = requests.get(ftp_url, timeout=10)
        response.raise_for_status()
        
        # Check if proteome tar file exists in listing
        content = response.text
        proteome_pattern = f"{proteome_id}_"
        available = proteome_pattern in content
        
        return {
            'proteome_id': proteome_id,
            'available': available,
            'ftp_url': ftp_url,
            'status': 'OK' if available else 'NOT_FOUND'
        }
    except Exception as e:
        return {
            'proteome_id': proteome_id,
            'available': False,
            'error': str(e),
            'status': 'ERROR'
        }


def main():
    parser = argparse.ArgumentParser(
        description="Check AlphaFold database versions"
    )
    parser.add_argument(
        '--proteomes',
        nargs='+',
        choices=list(PROTEOME_IDS.keys()),
        default=list(PROTEOME_IDS.keys()),
        help='Proteomes to check'
    )
    
    args = parser.parse_args()
    
    logging.basicConfig(level=logging.INFO)
    
    print("\n" + "="*60)
    print("  ALPHAFOLD DATABASE VERSION CHECK")
    print("="*60 + "\n")
    
    all_ok = True
    
    for proteome_name in args.proteomes:
        proteome_id = PROTEOME_IDS[proteome_name]
        print(f"Checking {proteome_name} ({proteome_id})...")
        
        result = check_proteome_availability(proteome_id)
        
        if result['status'] == 'OK':
            print(f"  ✓ Available at: {result['ftp_url']}")
        else:
            print(f"  ✗ {result['status']}: {result.get('error', 'Not found')}")
            all_ok = False
        
        print()
    
    print("="*60 + "\n")
    
    if not all_ok:
        sys.exit(1)


if __name__ == '__main__':
    main()
