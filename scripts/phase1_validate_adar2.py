#!/usr/bin/env python
"""
Phase 1: Validate pipeline on ADAR2.

This is the critical test - if the pipeline doesn't find IP6 in ADAR2,
the parameters need adjustment.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from cryptic_ip.validation import validate_adar2, ValidationSuite


def main():
    print("="*70)
    print("PHASE 1: PIPELINE VALIDATION")
    print("="*70)
    print("\nTesting on ADAR2 - the gold standard for cryptic IP6 binding")
    print("If this fails, pipeline parameters need adjustment\n")
    
    # Validate ADAR2
    print("\n1. Testing AlphaFold structure...")
    alphafold_results = validate_adar2(use_alphafold=True)
    
    print("\n2. Testing crystal structure (1ZY7)...")
    crystal_results = validate_adar2(use_alphafold=False)
    
    # Check if both passed
    if alphafold_results['validation_passed'] and crystal_results['validation_passed']:
        print("\n" + "="*70)
        print("✓ PHASE 1 VALIDATION: PASSED")
        print("="*70)
        print("\nPipeline successfully identifies cryptic IP6 site in ADAR2")
        print("Ready to proceed to Phase 2 (database construction)\n")
        return 0
    else:
        print("\n" + "="*70)
        print("✗ PHASE 1 VALIDATION: FAILED")
        print("="*70)
        print("\nPipeline failed to correctly identify ADAR2 IP6 site")
        print("Parameters need adjustment before proceeding\n")
        return 1


if __name__ == '__main__':
    sys.exit(main())
