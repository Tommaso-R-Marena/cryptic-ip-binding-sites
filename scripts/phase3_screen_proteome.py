#!/usr/bin/env python
"""
Phase 3: Screen proteome for cryptic IP binding sites.

Runs validated pipeline across entire proteome.
"""

import sys
import argparse
from pathlib import Path
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent.parent))

from cryptic_ip.analysis import ProteinAnalyzer
from cryptic_ip.database import ProteomeManager


def main():
    parser = argparse.ArgumentParser(
        description='Screen proteome for cryptic IP binding sites'
    )
    parser.add_argument(
        'proteome_dir',
        type=str,
        help='Path to proteome structures directory'
    )
    parser.add_argument(
        '--output',
        '-o',
        default='screening_results.csv',
        help='Output CSV file'
    )
    parser.add_argument(
        '--threshold',
        '-t',
        type=float,
        default=0.60,
        help='Minimum composite score threshold'
    )
    parser.add_argument(
        '--max-structures',
        '-n',
        type=int,
        help='Maximum number of structures to process (for testing)'
    )
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from previous run'
    )
    
    args = parser.parse_args()
    
    print("="*70)
    print("PHASE 3: PROTEOME-WIDE SCREENING")
    print("="*70)
    print(f"\nProteome: {args.proteome_dir}")
    print(f"Score threshold: {args.threshold}")
    print(f"Output: {args.output}\n")
    
    # Initialize manager
    manager = ProteomeManager(args.proteome_dir)
    catalog = manager.build_catalog()
    
    if args.max_structures:
        catalog = catalog.head(args.max_structures)
        print(f"Processing first {args.max_structures} structures (test mode)\n")
    
    # Check for existing results
    processed = set()
    if args.resume and Path(args.output).exists():
        existing = pd.read_csv(args.output)
        processed = set(existing['uniprot_id'].unique())
        print(f"Resuming: {len(processed)} structures already processed\n")
    
    print(f"Screening {len(catalog)} structures...\n")
    
    all_results = []
    error_count = 0
    
    for idx, row in tqdm(catalog.iterrows(), total=len(catalog), desc="Screening"):
        uniprot_id = row['uniprot_id']
        
        # Skip if already processed
        if uniprot_id in processed:
            continue
        
        try:
            pdb_path = row['filepath']
            
            # Analyze structure
            analyzer = ProteinAnalyzer(pdb_path)
            scored = analyzer.score_all_pockets()
            
            # Add metadata
            scored['uniprot_id'] = uniprot_id
            scored['protein_file'] = row['filename']
            
            # Filter by threshold
            candidates = scored[scored['composite_score'] >= args.threshold]
            
            if len(candidates) > 0:
                all_results.append(candidates)
                
        except Exception as e:
            error_count += 1
            tqdm.write(f"Error processing {uniprot_id}: {e}")
            continue
    
    # Combine and save results
    print("\nProcessing complete!")
    print(f"Errors: {error_count}\n")
    
    if all_results:
        final_results = pd.concat(all_results, ignore_index=True)
        final_results = final_results.sort_values('composite_score', ascending=False)
        
        # Save results
        final_results.to_csv(args.output, index=False)
        
        # Print summary
        print("="*70)
        print("SCREENING RESULTS")
        print("="*70)
        print(f"\nTotal candidate sites found: {len(final_results)}")
        print(f"Unique proteins with candidates: {final_results['uniprot_id'].nunique()}")
        print(f"\nScore distribution:")
        print(f"  High confidence (>0.75): {len(final_results[final_results['composite_score'] >= 0.75])}")
        print(f"  Moderate (0.60-0.75): {len(final_results[(final_results['composite_score'] >= 0.60) & (final_results['composite_score'] < 0.75)])}")
        print(f"\nTop 10 candidates:")
        print(final_results[['uniprot_id', 'pocket_id', 'composite_score', 'volume', 'sasa', 'basic_residues']].head(10).to_string())
        print(f"\nFull results saved to: {args.output}")
        print("="*70 + "\n")
    else:
        print("No candidates found above threshold\n")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
