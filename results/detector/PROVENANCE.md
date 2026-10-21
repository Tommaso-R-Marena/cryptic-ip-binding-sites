# Detector benchmark: provenance

**Run.** `scripts/detect_benchmark.py run`, seed 20261012, commit `d38fc0f`, on the session container (single core).

**Plan.** `docs/DETECTOR_PLAN.md` as amended by `docs/DETECTOR_PLAN_AMENDMENT_1.md`, which pre-registers this benchmark as a method validation and states what it may not be used for.

Every figure below is read from `benchmark.json` by script. Nothing is transcribed.

## What this is

48 measurements across 7 suites, 0 failing. Synthetic fixtures with valid covalent geometry: idealised helix barrels, chosen sequences, no energy minimisation.

**It says nothing about whether the detector finds real inositol-phosphate sites.** That is Q1-Q5 in the plan, asked of real structures, and none has run.

## Suites

| code | title | measurements | failures | seconds |
|---|---|---|---|---|
| B1 | Determinism | 6 | 0 | 0.9 |
| B2 | Invariance | 2 | 0 | 0.1 |
| B3 | Robustness to coordinate noise | 8 | 0 | 0.9 |
| B4 | Discrimination on known construction | 5 | 0 | 0.7 |
| B5 | Cost and scaling | 12 | 0 | 3.8 |
| B6 | Statistical recovery | 9 | 0 | 1.5 |
| B7 | Head to head against the deployed scorer | 6 | 0 | 1.6 |

## The figures worth quoting

| quantity | value |
|---|---|
| field spectrum rotation discrepancy | 2.02e-03 |
| field spectrum translation discrepancy | 6.1e-16 |
| field discrepancy at 0.1 A coordinate noise | 0.014 |
| field discrepancy at 0.5 A coordinate noise | 0.138 |
| cavity volume change at 0.5 A noise | 0.215 |
| basic minus acidic well depth | 3.43 kT/e |
| cavity detection scaling exponent | 1.02 |
| field panel scaling exponent | 0.032 |
| conformal coverage at alpha = 0.1 | 0.903 |
| minimum detectable AUC at 5 homology groups | 0.87 |
| monotonicity violations in 2000 perturbations | 0 |
| deployed scorer separation margin | 0.052 |
| monotone score separation margin | 0.095 |

## Two results worth a sentence each

**The field descriptor's cost does not grow with the protein.** Its scaling exponent is 0.032 in d log t / d log n, measured over barrels from 360 to 1,260 atoms whose radius grows with their helix count so the structure really does get bigger. That is what the 5-Debye-length screening cutoff was for, and it means the 2,352-residue tail of the candidate set costs what a small protein costs.

**Nothing below AUC 0.87 is detectable at this dataset's size.** The buried-site positives occupy 5 strict homology groups with one holding 87% of them. The previously reported held-out AUC of 0.499 is therefore a statement about power, not about biology, and Q4's decision rule requires this figure to be published before any AUC is.

## What the benchmark found wrong with itself

Two things, both recorded because they are the reason to have it:

1. **The ground-truth labels were wrong.** Two poly-alanine channels were labelled positive because they are buried. They have a well depth of 0.00 kT/e, zero basic residues and a shape correlation of 0.00 - buried, and not inositol-phosphate sites. Both scorers duly reported a negative separation margin. Burial is tested on its own in B4, where it belongs, and `GROUND_TRUTH_REASON` now publishes the reason for every label.
2. **The cost suite timed its own warm-up.** The smallest structure measured slower than the next size up, because the first call pays for scipy's ndimage import. The suite now warms up and takes the best of repeats.

## Limits

* The fixtures are geometrically valid and biologically arbitrary.
* The cost figures are from barrels, not candidates. Study Q's first stage is a timing probe on real structures for that reason.
* No parameter in `docs/DETECTOR_PLAN.md` may be tuned on these numbers, and none has been.
* The benchmark is a regression gate. A failing measurement means behaviour changed; the response is to fix the code or amend the plan, never to loosen the bound.

