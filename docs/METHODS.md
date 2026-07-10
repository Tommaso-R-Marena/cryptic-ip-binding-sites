# Methods

This document summarizes the computational methods used in the cryptic IP-binding site pipeline. Auto-generated run-specific text is written to `results/publication/METHODS_AUTO.md` when you execute the publication package.

## Overview

The pipeline identifies buried inositol phosphate (IP) binding pockets in protein structures where IPs may act as structural cofactors rather than canonical signaling ligands. It combines:

1. **fpocket** — geometric cavity detection  
2. **BioPython SASA** — residue burial (Shrake–Rupley)  
3. **Optional APBS** — pocket electrostatic potential  
4. **Composite scoring** — depth, burial, basic residues, volume, electrostatics  
5. **Tiered validation** — ADAR2 positive vs PH-domain negative controls  

## Pocket detection and scoring

| Step | Tool | Key parameters |
|------|------|----------------|
| Pocket detection | fpocket | min α-spheres = 3 |
| SASA | BioPython ShrakeRupley | per-residue, pocket mean |
| Electrostatics | pdb2pqr + APBS | pH 7.4 (optional) |
| Composite score | `PocketScorer` | SASA weight 40%, depth 25% |

**Proteome candidate filters** (default): composite ≥ 0.75, SASA ≤ 10 Ų, ≥ 4 basic residues, volume 300–800 Ų, pocket pLDDT ≥ 70.

## Validation design

### Tier-1 (Phase 1 gate)

| Control | PDB | Role |
|---------|-----|------|
| ADAR2 | 1ZY7 | Cryptic IP6 positive |
| PLCδ1 PH | 1MAI | Surface IP3 negative |

**Pass criterion:** both controls pass individually **and** tier-1 score separation > 0.50.

### Tier-2 (extended panel)

| Control | PDB | Notes |
|---------|-----|-------|
| Pds5B | 5HDT | Often crystal artifact (surface IP6) |
| HDAC1 | 5ICN | Semi-cryptic interface |
| Btk PH | 1BWN | Surface negative (decoy mode) |

Burial classes: `cryptic`, `semi_cryptic`, `surface`, `crystal_artifact`.

## Machine learning benchmark

- **Labels:** pockets within 8 Å of annotated ligand atoms  
- **Splits:** GroupKFold by PDB ID  
- **Metrics:** ROC AUC, PR AUC, MCC  
- **Deployment:** threshold scoring is primary; ML is exploratory given label scarcity  

## RCSB validation dataset

Built with `scripts/build_ip_validation_dataset.py`:

- Query: PDB structures with IP3–IP6 ligands  
- Burial labels from ligand SASA: Cryptic (≤5 Ų), Semi-cryptic (5–50), Surface (>50)  

## Reproducibility

```bash
python scripts/run_publication_package.py --output-dir results/publication
```

Outputs: `provenance_manifest.jsonld`, pinned `requirements.txt`, git commit hash.

See also: [VALIDATION.md](VALIDATION.md), [REPRODUCIBILITY.md](REPRODUCIBILITY.md), [PUBLICATION_PACKAGE.md](PUBLICATION_PACKAGE.md).
