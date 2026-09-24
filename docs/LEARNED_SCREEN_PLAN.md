# Learned proteome ranking: plan

Pre-registered on 2026-09-24. At that time no proteome pocket had been scored with a
learned model, and the three-proteome screen (run 35935291031) had not finished.

## Why

The benchmark's learned `ip_site` model generalises to protein families released after
its training data. On the temporal holdout it scores ROC-AUC 0.871 [0.782, 0.932],
against 0.730 for the rule-based score. The proteome screen ranks proteins with the
rule-based score. This plan ranks every screened AlphaFold model with the learned
model instead, and tests the ranking against an external truth before naming any
candidate.

## Model

- **Training data.** The model is locked on the **whole** inositol phosphate benchmark
  table: development and holdout rows, task `ip_site`, the full 39-descriptor arm, from
  run 35949588200. The table is identified by its SHA-256.
- **Model selection.** Family, hyperparameters, threshold and Platt calibration come from
  grouped inner cross-validation over strict homology groups
  (`protocol.select_candidate`, 10 draws per family, seed 0).
- **Performance estimate.** Nothing about the model is tuned after this point. Its
  expected performance on new families is the holdout estimate above.

## Scores

- **Pocket score.** The model's calibrated probability for each pocket.
- **Protein score.** The highest pocket score among the protein's confident pockets,
  meaning pockets with a mean pLDDT of at least 70. The rule-based protein score is
  the highest `composite_score` among the same pockets.
- **Unscoreable proteins.** A protein with no confident pocket gets no score and is
  counted as such.

## Truth, and the leak it must avoid

- **Positives.** Proteins whose UniProt binding-site or function annotation names an
  inositol phosphate (`known_ip_annotation`: tris-, tetrakis-, pentakis- or
  hexakisphosphate, phytate, InsP3-8).
- **Negatives.** All other screened proteins.
- **The leak.** A protein homologous to a benchmark training protein would be ranked
  from memory, not from physics.
- **The exclusion.** Every screened protein with an MMseqs2 hit against the UniProt
  sequences of the benchmark's entries is excluded from the evaluation. A hit means
  ≥ 30 % identity over ≥ 50 % of the shorter sequence, with E ≤ 1e-3, the benchmark's
  own grouping criterion. Full-length UniProt sequences are used, which is more
  conservative than PDB chains. Benchmark entries without an accession are counted.

## Evaluation (per proteome, and pooled)

**L1.** Among unseen proteins, the learned protein score ranks annotated binders
above other proteins. The measure is ROC-AUC with a 95 % interval from resampling
MMseqs2 clusters (30 % identity, 50 % coverage) of the proteome's proteins, so that
paralogues are not independent. The result is **supported** if the pooled interval's
lower bound exceeds 0.5.

**L2.** The learned ranking beats the rule-based ranking on the same proteins. The
measure is the paired ROC-AUC difference, with the same resampling. It is
**supported** if the pooled interval's lower bound exceeds 0.

L1 and L2 are Holm-corrected together. If an organism has fewer than 5 unseen
annotated binders, it is reported but gets no per-organism decision.

## Candidates

Candidates are the 25 top-ranked unseen proteins in each proteome that are **not**
annotated binders. Each is reported with:

- its top pocket's score;
- the residues lining that pocket;
- its hull depth and pLDDT;
- the fraction of annotated binders among unseen proteins at the same rank or above,
  as an empirical estimate of precision at that depth.

That estimate is a lower bound, because unannotated true binders count as negatives.
A candidate is a hypothesis for experiment, not a finding. It is labelled that way
wherever it appears.
