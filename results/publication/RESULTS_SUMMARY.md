# Cryptic IP Binding Sites — Publication Results Summary

Generated: 2026-06-19T18:23:47.372412+00:00

## Structural control benchmark

- Positive controls passed: False
- Negative controls passed: True
- Mean positive score: 0.519
- Mean negative score: 0.456
- Score separation: 0.063
- Tier-1 separation (ADAR2 vs PLCδ1): 0.584
- Phase 1 gate passed: True

### Positive controls

| protein | control_type | tier | pdb_id | ip_type | score | pocket_score | expected | passed | total_pockets | sasa | basic_residues | description |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ADAR2 | positive | 1 | 1ZY7 | IP6 | 0.8557630936719473 | 0.6394077341798682 | 0.75 | True | 76 | 0.0 | 8 | Gold standard - completely buried IP6 |
| Pds5B | positive | 2 | 5HDT | IP6 | 0.212117565563816 | 0.53029391390954 | 0.65 | False | 313 | 269.70226099239414 | 10 | Cohesin regulator with buried IP6 |
| HDAC1 | positive | 2 | 5ICN | IP4 | 0.4902588915652315 | 0.7610180000000001 | 0.6 | False | 68 | 36.06112313260764 | 11 | Histone deacetylase with IP4 at interface |

### Negative controls

| protein | control_type | tier | pdb_id | ip_type | score | pocket_score | expected | passed | total_pockets | sasa | basic_residues | description |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PLCd1_PH | negative | 1 | 1MAI | IP3 | 0.2716676314980122 | 0.6791690787450305 | 0.3 | True | 11 | 99.11890698167906 | 6 | Classic surface-exposed PH domain |
| Btk_PH | negative | 2 | 1BWN | IP4 | 0.6401115586051882 | 0.7238878205582052 | 0.35 | True | 52 | 23.708267721357046 | 9 | Kinase PH domain with surface-accessible IP4 |

## ML classifier benchmark (held-out test set)

| method | test_roc_auc | test_pr_auc |
| --- | --- | --- |
| ML (random_forest) | 0.5814401814401815 | 0.0267667096893359 |
| Threshold scoring | 0.4470925470925471 | 0.013889916435016 |

## Validation dataset

Source: `data/validation/ip_binding_validation_dataset.csv`

## Figure outputs

- `figures/publication/Figure1_Overview.pdf`
- `figures/publication/Figure2_Comparative_Proteomics.pdf`
- `figures/publication/Figure3_Top_Candidates.pdf`

## Interpretation

Phase 1 gate **passed**: ADAR2 and PLCδ1 PH controls both passed, with tier-1 score separation (0.584) above the 0.50 threshold. Burial-aware validation scores separate known cryptic positives from canonical surface PH-domain negatives.
Full control benchmark separation (0.063) is below the 0.30 target; extended tier-2 controls should be reviewed before publication claims.
On the RCSB-derived validation set, ML scoring (ROC AUC 0.58) outperforms fixed thresholds (ROC AUC 0.45).
