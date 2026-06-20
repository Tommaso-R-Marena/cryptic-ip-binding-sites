#!/usr/bin/env python
"""
Phase 1: Validate pipeline on ADAR2.

The crystal structure (1ZY7) with bound IHP is the tier-1 gold standard.
AlphaFold apo is run for comparison but does not gate CI/CD.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from cryptic_ip.validation import validate_adar2


def main():
    print("=" * 70)
    print("PHASE 1: PIPELINE VALIDATION")
    print("=" * 70)
    print("\nTesting on ADAR2 - the gold standard for cryptic IP6 binding")
    print("Crystal structure (1ZY7) is required to pass; AlphaFold is informational.\n")

    print("\n1. Testing AlphaFold structure (informational)...")
    alphafold_results = validate_adar2(use_alphafold=True, use_electrostatics=False)

    print("\n2. Testing crystal structure (1ZY7)...")
    crystal_results = validate_adar2(use_alphafold=False, use_electrostatics=False)

    if crystal_results["validation_passed"]:
        print("\n" + "=" * 70)
        print("✓ PHASE 1 VALIDATION: PASSED")
        print("=" * 70)
        print("\nPipeline successfully identifies cryptic IP6 site in ADAR2 crystal")
        if not alphafold_results["validation_passed"]:
            print("Note: AlphaFold apo structure did not pass all criteria (expected without ligand).")
        print("Ready to proceed to Phase 2 (database construction)\n")
        return 0

    print("\n" + "=" * 70)
    print("✗ PHASE 1 VALIDATION: FAILED")
    print("=" * 70)
    print("\nPipeline failed to correctly identify ADAR2 IP6 site in crystal structure")
    print("Parameters need adjustment before proceeding\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
