# Detector benchmark

**all suites pass** - 48 measurement(s) across 7 suite(s), seed 20261012.

These are engineering measurements on synthetic fixtures. They say nothing
about whether the detector finds real binding sites; that is Q1-Q5 in
docs/DETECTOR_PLAN.md and none of those has run.

## B1 - Determinism

*Inputs support: positions. A result that is not reproducible cannot be a provenance claim, and this project's reports are extracted from logs by script on the assumption that a rerun gives the same numbers.*

| measurement | value | expected | ok |
|---|---|---|---|
| cavity map reproducible (acidic_lined) | 1 | identical | yes |
| cavity map reproducible (basic_lined) | 1 | identical | yes |
| cavity map reproducible (narrow) | 1 | identical | yes |
| cavity map reproducible (trp_plugged) | 1 | identical | yes |
| cavity map reproducible (wide_open) | 1 | identical | yes |
| field spectrum reproducible | 1 | bit-identical | yes |


## B2 - Invariance

*Inputs support: positions. The field descriptor's whole purpose is to need no alignment, which is what study N's rigid transplantation needed and got wrong.*

| measurement | value | expected | ok |
|---|---|---|---|
| rotation discrepancy | 0.002018 | < 0.01 | yes |
| translation discrepancy | 6.063e-16 | < 1e-9 | yes |

- rotation discrepancy: largest change in the power spectrum, as a fraction of its own scale, over three rotations about the pocket centre
- translation discrepancy: exact, because translating the system and its centre together changes nothing

## B3 - Robustness to coordinate noise

*Inputs support: positions. Every candidate is a predicted model, so sensitivity to small coordinate error is sensitivity to the predictor rather than to the protein.*

| measurement | value | expected | ok |
|---|---|---|---|
| field discrepancy at 0.1 A noise | 0.01449 | grows smoothly | - |
| cavity volume change at 0.1 A noise | 0.007828 fraction | grows smoothly | - |
| field discrepancy at 0.25 A noise | 0.0358 | grows smoothly | - |
| cavity volume change at 0.25 A noise | 0.01826 fraction | grows smoothly | - |
| field discrepancy at 0.5 A noise | 0.1384 | grows smoothly | - |
| cavity volume change at 0.5 A noise | 0.2146 fraction | grows smoothly | - |
| field discrepancy at 1.0 A noise | 0.2378 | grows smoothly | - |
| cavity volume change at 1.0 A noise | 0.4175 fraction | grows smoothly | - |


## B4 - Discrimination on known construction

*Inputs support: geometry. The fixtures differ in exactly one property each, so a detector that cannot separate them is broken rather than merely weak. Separating them is necessary and nowhere near sufficient for the biological questions.*

| measurement | value | expected | ok |
|---|---|---|---|
| basic minus acidic well depth | 3.43 kT/e | > 0, strongly | yes |
| acidic well is repulsive | -1.935 kT/e | < 0 | yes |
| acidic pocket refused as a polyanion site | 1 | no_polyanion_site | yes |
| narrow barrel is more enclosed | 0.3522 directions | > 0 | yes |
| veto fires on a repulsive well | 1 | field_depth_kt vetoed | yes |

- basic minus acidic well depth: a pocket lined with arginine and lysine must present a more favourable well to a -9 polyanion than one lined with aspartate and glutamate
- acidic pocket refused as a polyanion site: The field is unfavourable for a polyanion at all (-1.93 kT/e), so no comparison between polyanions is meaningful.
- veto fires on a repulsive well: burial and model confidence must not outvote the physics

## B5 - Cost and scaling

*Inputs support: positions. The candidates run to 2,352 residues and the cavity grid grows with the bounding box, so the scaling exponent matters more than the absolute figure. The barrel radius grows with the helix count, so the protein really does get bigger rather than just denser.*

| measurement | value | expected | ok |
|---|---|---|---|
| cavity detection, 360 atoms | 34.53 ms | - | - |
| field panel, 360 atoms | 131 ms | - | - |
| cavity detection, 540 atoms | 42.76 ms | - | - |
| field panel, 540 atoms | 130 ms | - | - |
| cavity detection, 720 atoms | 64.56 ms | - | - |
| field panel, 720 atoms | 135.3 ms | - | - |
| cavity detection, 900 atoms | 83.28 ms | - | - |
| field panel, 900 atoms | 133.5 ms | - | - |
| cavity detection, 1260 atoms | 116.7 ms | - | - |
| field panel, 1260 atoms | 0 ms | no pocket | - |
| cavity detection scaling exponent | 1.016 d log t / d log n | < 2 | yes |
| field panel scaling exponent | 0.03156 d log t / d log n | < 2 | yes |

- field panel, 1260 atoms: no enclosed cavity at this barrel radius, so nothing to describe
- cavity detection scaling exponent: fitted over the measured sizes; above 2 the method would not reach the long tail of the candidate set
- field panel scaling exponent: fitted over the measured sizes; above 2 the method would not reach the long tail of the candidate set

## B6 - Statistical recovery

*Inputs support: statistics. Synthetic populations with a chosen label frequency and a chosen class structure, so the estimator can be compared with the truth rather than with another estimate.*

| measurement | value | expected | ok |
|---|---|---|---|
| label frequency recovered at c=0.2 | 0.9761 | 0.92 to 1.08 | yes |
| label frequency recovered at c=0.5 | 0.9773 | 0.92 to 1.08 | yes |
| label frequency recovered at c=0.9 | 0.9788 | 0.92 to 1.08 | yes |
| conformal coverage at alpha=0.1 | 0.903 | 0.85 to 0.95 | yes |
| minimum detectable AUC at 5 groups | 0.8724 | > 0.8 | yes |
| minimum detectable AUC at 40 groups | 0.6518 | < 0.7 | yes |
| group AUC collapses copies | 1 | 1.0 on one positive group | yes |
| monotonicity violations in 2000 perturbations | 0 | 0 | yes |
| elastic network rigid-body modes | 6 | exactly 6 | yes |

- label frequency recovered at c=0.2: the mean estimator gives 0.70 of the truth, biased low by the class overlap
- label frequency recovered at c=0.5: the mean estimator gives 0.70 of the truth, biased low by the class overlap
- label frequency recovered at c=0.9: the mean estimator gives 0.70 of the truth, biased low by the class overlap
- conformal coverage at alpha=0.1: the guarantee is the point of the method, so it is measured rather than cited
- group AUC collapses copies: fifty copies of one protein are one observation, not fifty
- elastic network rigid-body modes: three translations and three rotations; any other count is a wrong Hessian

## B7 - Head to head against the deployed scorer

*Inputs support: geometry. Both scorers are handed the same measurements from the same fixtures, so the comparison is between the methods and not between two measurement paths. The cases are the ones where being wrong means broken, not the ones where being right means useful.*

| measurement | value | expected | ok |
|---|---|---|---|
| deployed composite: worst positive minus best negative | 0.05158 | > 0 | yes |
| monotone score: worst positive minus best negative | 0.09458 | > 0 | yes |
| constructed negatives the deployed scorer admits | 0 | (reported) | - |
| constructed negatives the monotone score admits | 0 | 0 | yes |
| constructed positives the monotone score misses | 0 | 0 | yes |
| abstains when burial and confidence are both unmeasured | 1 | abstain | yes |

- deployed composite: worst positive minus best negative: a negative margin means the raw ranking puts a constructed negative above a constructed positive
- monotone score: worst positive minus best negative: a negative margin means the raw ranking puts a constructed negative above a constructed positive
- constructed negatives the deployed scorer admits: none, at this threshold
- constructed negatives the monotone score admits: none
- constructed positives the monotone score misses: none
- abstains when burial and confidence are both unmeasured: interval 0.459 to 0.803, missing relative_sasa, cavity_openness, mean_plddt
