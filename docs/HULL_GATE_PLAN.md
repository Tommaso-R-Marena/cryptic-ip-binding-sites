# The screen's hull-depth gate: plan

Pre-registered on 2026-09-24, before the proteome screen's shards were re-aggregated with
any gate other than the current one. Changes are dated amendments in separate files,
written before the results they could affect are read.

## Why

The calibrated screen (`CALIBRATED_CRITERIA` in `cryptic_ip/analysis/proteome_stats.py`)
requires a pocket to sit at least 10 Å inside the protein's convex hull. The threshold
was set between ADAR2 and four PH-domain controls: one positive, four negatives. Two
later results contradict the choice:

- in the transfer set (914 buried sites in 174 families), adding hull depth to a
  learned model changes ROC-AUC by +0.002, and outside the nucleotide-fold super-group,
  under strict grouping, it lowers it by 0.086;
- on the benchmark's temporal holdout, hull depth lowered the rule score (0.934 → 0.890).

This plan decides whether to keep, relax or remove the gate, by a rule fixed here.

## Data

- Pockets: every `*_pockets_part*.csv.gz` shard of the proteome screen, run 35935291031
  (yeast, human, Dictyostelium; all 39 descriptors, pLDDT, hull depth). The composite
  score is recomputed with the current `PocketScorer`, as `proteome_screen.py aggregate`
  does.
- Proteins: the learned screen's `proteins.csv.gz` (run 35992676292) supplies, for every
  scored protein, the `seen` flag (MMseqs2 homologue among benchmark proteins), the
  MMseqs2 cluster (30 % identity, 50 % coverage) and the UniProt IP annotation
  (`known_ip_annotation`). Only **unseen** proteins are evaluated.

## Arms

Identical except for the hull-depth gate; every other calibrated criterion stays
(score ≥ 0.54, no SASA gate, ≥ 4 basic residues, cavity volume window, mean pLDDT ≥ 70):

- **gate10**: hull depth ≥ 10 Å (current);
- **gate5**: hull depth ≥ 5 Å (relaxed);
- **none**: no hull-depth gate.

**Protein ranking under an arm:** the highest composite score among the protein's
pockets that pass every non-score gate of the arm; a protein with no such pocket ranks
below all others (ties). **Hits** under an arm: proteins with a pocket passing every
gate, including the score gate.

## Measures (pooled over the three proteomes)

- **Primary:** paired ROC-AUC differences of the protein rankings for annotated IP
  binders among unseen proteins: Δnone = AUC(none) − AUC(gate10) and
  Δrelax = AUC(gate5) − AUC(gate10). 95 % intervals from 2,000 bootstrap resamples of
  MMseqs2 clusters (seed 20260929), applied to all arms together.
- **Secondary:** recall of annotated binders among the top K unseen proteins of each
  arm's ranking, with K = the number of unseen proteins called hits under gate10; the
  paired recall differences with the same resampling. Also the number of hits, and
  recall and precision among hits, per arm and per proteome.

## Decision (fixed before any re-aggregation)

A gate stays only if it demonstrably helps. With a non-inferiority margin of 0.01
ROC-AUC:

1. **Remove** if the lower bound of Δnone is at least −0.01 and recall@K(none) is at
   least recall@K(gate10) (point estimates).
2. Otherwise **relax** (to 5 Å) if the lower bound of Δrelax is at least −0.01 and
   recall@K(gate5) is at least recall@K(gate10).
3. Otherwise **keep** the 10 Å gate.

The decision is implemented as a code change to `CALIBRATED_CRITERIA` (and the tests and
documentation that state its value) in the same pull request, citing this plan.

## Execution and outputs

`.github/workflows/hull-gate.yml`. Results printed between `BEGIN_HULL_GATE_JSON` and
`END_HULL_GATE_JSON` and extracted into `results/hull_gate/`.
