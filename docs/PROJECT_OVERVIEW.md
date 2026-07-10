# Project Overview

## Scientific motivation

Inositol phosphates (IPs), especially IP6, can bind protein pockets that are **buried from solvent** and stabilize fold or assembly — distinct from canonical surface PH-domain signaling. ADAR2 (PDB 1ZY7) is the reference cryptic IP6 site: ligand SASA ≈ 0 Ų with multiple coordinating basic residues.

This repository implements a **computational screening pipeline** to:

1. Detect geometric pockets (fpocket)  
2. Score burial and basic-residue complementarity  
3. Validate against tiered structural controls  
4. Screen AlphaFold proteomes (yeast, human, *Dictyostelium*)  
5. Compare hit rates with cellular IP6 concentration hypotheses  

## Experimental plan

| Phase | Goal | Status |
|-------|------|--------|
| 1 | ADAR2 + negative control validation | Complete (tier-1 CI gate) |
| 2 | RCSB IP-binding benchmark dataset | Complete (136 structures) |
| 3 | Proteome-wide AlphaFold screens | Yeast pilot only |
| 4 | ML vs threshold benchmark | Complete (exploratory) |
| 5 | Comparative statistics + figures | Framework ready |
| 6 | Manuscript + Zenodo archive | In progress |

## Key hypothesis

Organisms with higher intracellular IP6 may retain more cryptic IP-binding sites in their proteomes — testable after full three-organism screens via `scripts/phase5_comparative_analysis.py`.

## Further reading

- [METHODS.md](METHODS.md) — computational methodology  
- [VALIDATION.md](VALIDATION.md) — control design and pass criteria  
- [PUBLICATION_PACKAGE.md](PUBLICATION_PACKAGE.md) — one-command manuscript artifacts  
- [manuscript/README.md](../manuscript/README.md) — submission checklist  
