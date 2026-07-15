# Methods (auto-generated)

## Structural input and pocket detection

Protein structures were obtained from the RCSB Protein Data Bank (PDB) for validation benchmarks and from the AlphaFold Database for proteome screening. Pockets were detected with **fpocket** (default minimum alpha-sphere count = 3). Pocket-lining residues were assigned from fpocket atom records; solvent-accessible surface area (SASA) was computed per residue with BioPython Shrake–Rupley.

## Cryptic-site scoring

Each pocket received a composite score (0–1) weighting pocket depth (25%), burial/SASA (40%), basic residue count (20%), volume fit for IP3–IP6 (10%), and optional electrostatic potential from APBS (5%). Proteome candidates were filtered at composite score ≥ 0.75, pocket SASA ≤ 10 Ų, ≥ 4 basic residues (Arg/Lys/His), volume 300–800 Ų, and pocket mean pLDDT ≥ 70 (AlphaFold models).

## Validation controls

**Tier-1 gate (Phase 1):** ADAR2/IP6 (PDB 1ZY7) positive control vs PLCδ1 PH/IP3 (PDB 1MAI) negative control. Passing requires both controls to pass individual criteria and tier-1 score separation > 0.50. Scores use burial-aware cryptic-likeness (ligand SASA + pocket composite).

**Tier-2 controls:** Pds5B (5HDT), HDAC1 (5ICN), Btk PH (1BWN). Crystal artifacts (ligand SASA > 50 Ų) are flagged and excluded from positive-pass claims.

## Machine learning benchmark

Pocket-level labels were assigned within 8 Å of annotated ligand atoms. Models (random forest, XGBoost) were evaluated with grouped cross-validation by PDB ID to prevent structure leakage. Metrics: ROC AUC, PR AUC, Matthews correlation coefficient (MCC).

## RCSB validation dataset

Structures containing inositol phosphate ligands (IP3–IP6) were retrieved from PDB and classified by ligand SASA: Cryptic (≤ 5 Ų), Semi-cryptic (5–50 Ų), Surface (> 50 Ų).

## Reproducibility

- Random seed: 42
- Score threshold: 0.75
- Pipeline commit: 94776adc6a5f96d975312c11de60c4a1590aea4a
- AlphaFold DB snapshot: 2024-01-15
- PDB fetch date: 2024-01-15

### Software versions

- biopython: 1.86
- numpy: 2.2.6
- pandas: 2.3.3
- python: 3.12.3
- scikit-learn: 1.7.2
- scipy: 1.15.3
- xgboost: 3.2.0

Floating-point note: compare scores with tolerance 1e-6 across BLAS implementations.