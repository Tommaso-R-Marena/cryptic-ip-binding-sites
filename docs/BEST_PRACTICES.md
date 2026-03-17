# Best Practices

## ML classifier vs threshold scoring

- Use **threshold scoring** when you need transparent, mechanistic ranking and immediate interpretability.
- Use the **ML classifier** when you have curated positives/negatives from your own domain and need better precision-recall tradeoffs.
- Recommended: run threshold scoring first, then apply ML rescoring to top pockets.

## Interpreting composite scores

- Treat the score as a **prioritization metric**, not ground truth.
- Inspect contributing features (SASA, depth, basic residues, electrostatics, volume).
- High score + weak biological context should be deprioritized vs moderate score + strong orthology/function evidence.

## Statistical power for proteome size

- Small proteomes (<7k proteins): wider confidence intervals; use exact tests and report Wilson CIs.
- Mid-size proteomes (7k–15k): balanced exploratory + confirmatory design.
- Large proteomes (>15k): pre-register thresholds to avoid multiple-testing drift.

## MD simulation length recommendations

- **Screening stage:** 2–5 ns per candidate for quick triage.
- **Refinement stage:** 20–50 ns for top hits.
- **Mechanistic claims/publication figures:** 100+ ns with replicate seeds and convergence checks.

## Publication checklist

- Report software versions (pipeline, fpocket, APBS, OpenMM).
- Include positive/negative controls (ADAR2 and PH-domain negatives).
- Provide thresholds and sensitivity analyses.
- Share candidate tables + scripts for figure regeneration.
- Deposit trajectories (or representative subsets) and analysis notebooks.
