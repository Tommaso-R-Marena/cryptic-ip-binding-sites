#!/usr/bin/env python
"""
Phase 2: Download and organize AlphaFold proteomes.

Downloads structures for yeast, human, and Dictyostelium.
"""

import sys
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from cryptic_ip.database import ProteomeDownloader


def main():
    parser = argparse.ArgumentParser(
        description='Download AlphaFold proteomes for screening'
    )
    parser.add_argument(
        '--organisms',
        nargs='+',
        choices=['yeast', 'human', 'dictyostelium', 'all'],
        default=['yeast'],
        help='Organisms to download'
    )
    parser.add_argument(
        '--data-dir',
        default='data/structures',
        help='Data directory for structures'
    )
    
    args = parser.parse_args()
    
    # Expand 'all' option
    if 'all' in args.organisms:
        organisms = ['yeast', 'human', 'dictyostelium']
    else:
        organisms = args.organisms
    
    print("="*70)
    print("PHASE 2: PROTEOME DATABASE CONSTRUCTION")
    print("="*70)
    
    downloader = ProteomeDownloader(args.data_dir)
    
    # Show what will be downloaded
    print("\nProteomes to download:")
    total_size = 0
    for org in organisms:
        info = downloader.get_info(org)
        print(f"  - {info['organism']}: {info['proteins']} proteins, ~{info['size_gb']} GB")
        total_size += info['size_gb']
    
    print(f"\nTotal estimated size: ~{total_size} GB")
    print("\nNote: Downloads may take several hours depending on connection speed\n")
    
    # Confirm
    response = input("Proceed with downloads? [y/N]: ")
    if response.lower() != 'y':
        print("Download cancelled")
        return 0
    
    # Download each proteome
    for organism in organisms:
        try:
            print(f"\n{'='*70}")
            proteome_dir = downloader.download_proteome(organism)
            print(f"\n✓ {organism} complete: {proteome_dir}")
        except Exception as e:
            print(f"\n✗ Failed to download {organism}: {e}")
            continue
    
    print("\n" + "="*70)
    print("PHASE 2 COMPLETE")
    print("="*70)
    print("\nProteomes downloaded and ready for Phase 3 screening\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
