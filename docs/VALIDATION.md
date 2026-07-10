# Validation Protocol

This document describes how to interpret validation results and what must pass before publication claims.

## Phase 1 tier-1 gate (required)

Run locally or in CI:

```bash
python scripts/phase1_validate_adar2.py
# or
cryptic-ip validate-suite
```

**Required for proteome screening and manuscript tier-1 claims:**

| Metric | Threshold |
|--------|-----------|
| ADAR2 control | pass |
| PLCδ1 PH control | pass |
| Tier-1 separation (ADAR2 − PLCδ1) | > 0.50 |

CI enforces this via `.github/workflows/tier1-validation.yml`.

## Full control benchmark

`ValidationSuite.run_full_validation()` scores all tier-1 and tier-2 controls.

| Tier | Positives | Negatives |
|------|-----------|-----------|
| 1 | ADAR2 | PLCδ1 PH |
| 2 | Pds5B, HDAC1 | Btk PH |

**Interpretation rules:**

1. **Tier-1 gate passing** is necessary and sufficient for Phase 2/3 screening.  
2. **Tier-2 failures** (e.g. Pds5B crystal artifact) do **not** invalidate tier-1 if flagged as `crystal_artifact`.  
3. **Full separation > 0.30** is a secondary benchmark across all controls.  

## Scoring scale

`cryptic_likeness` (validation) and `composite_score` (screening) use aligned 0–1 scales where **higher = more cryptic/buried**.

Burial-aware validation score blends:

- Ligand SASA (when ligand present in crystal)  
- Pocket composite score (depth, SASA, basic residues, volume)  

## ML validation

Do **not** claim ML outperforms thresholds unless held-out ROC AUC exceeds threshold scoring on the same split **and** positive class n > 10 structures.

Current RCSB set has severe imbalance (~1 cryptic structure); treat ML as supplementary.

## Proteome screening validation

Before reporting hit rates:

1. Re-run yeast pilot with strict filters (≥0.75, SASA ≤10, ≥4 basic).  
2. Compare hit rate to literature estimates (~0.2–1% for cryptic sites).  
3. Manually inspect top 10 hits.  

## Publication checklist

- [ ] Tier-1 gate passes (CI green)  
- [ ] `results/publication/RESULTS_SUMMARY.md` reflects tier-2 caveats  
- [ ] Figure legends cite measured n and AUC  
- [ ] Supplementary tables exported (`supplementary/SUPPLEMENTARY_INDEX.md`)  
- [ ] Full proteome screens complete (yeast → human → dicty) for Figure 2 organism claims  

See [METHODS.md](METHODS.md) and [PUBLICATION_PACKAGE.md](PUBLICATION_PACKAGE.md).
