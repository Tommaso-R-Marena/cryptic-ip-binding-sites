# Does the sequence predict pocket formation, type, or function?

**Run.** `scripts/pocket_sequence.py analyse`, seed 20261020, 2000 shuffles per window.

**Pre-registered** in `docs/POCKET_SEQUENCE_PLAN.md`, dated before any sequence was read from any structure for this study.

## What is being tested, and what is not

A window of sequence factorises the way a lining of positions does — **composition** (which letters) **⊕ order** (in what arrangement) — and the null for order is permutation of the window's own letters. That is the same quotient `architecture.py` takes on positions.

**Every primary test below is on `composition_z`, never on a raw measure.** A raw difference is already explained: these sites are basic and pocket linings are hydrophobic. A non-zero z says the *arrangement* carries something the letters alone do not, which is the only form of the question this study calls open.

2001 `local` pocket windows and 8407 `dense` ones, of which 2001 carry an architecture class; 9 annotated-site windows against 45 length-matched controls from the same chains. 3219 windows carried an order null.

**Amendment 1** (`docs/POCKET_SEQUENCE_PLAN_AMENDMENT_1.md`, dated after the first run and before any of what it adds was read) re-specifies the formation arm, which the pre-registered rule could only run on ten windows out of two thousand, and adds a second window definition. Of what it adds, **only `formation_heldout` may become a finding** — the rule restricted to the resampling groups the pre-registered rule selected nothing from, so disjoint from everything already looked at. Every arm named `*_dense`, and `formation_matched`, is **exploratory by declaration**: no Holm family, and `is_finding` cannot become true for it however the row reads.

The three arms: **formation** is fpocket's best-ranked pockets against its worst within the same protein (rank ≤ 5), which asks whether the sequence anticipates a geometric judgement that reads no sequence; **type** is pockets the architecture module calls a helical face against those it calls none, which asks whether the letters anticipate the positions; **function** is an annotated inositol-phosphate contact shell against length-matched windows of its own chain.

## Primary

| arm | statistic | windows | groups | mean z (site) | mean z (control) | difference | 95% interval | Holm | could have detected | finding |
|---|---|---|---|---|---|---|---|---|---|---|
| formation | `helix_charge_periodicity` | 10 | 9 | -0.33 | -0.27 | -0.068 | [-0.919, +0.758] | 1.0000 | 1.278 | no |
| formation | `helix_hydropathy_periodicity` | 10 | 9 | +1.00 | -0.31 | +1.310 | [+0.383, +2.249] | 0.0800 | 1.417 | no |
| formation_heldout | `helix_charge_periodicity` | 400 | 63 | -0.10 | -0.15 | +0.053 | [-0.124, +0.226] | 1.0000 | 0.252 | no |
| formation_heldout | `helix_hydropathy_periodicity` | 400 | 63 | +0.27 | +0.26 | +0.014 | [-0.198, +0.215] | 1.0000 | 0.300 | no |
| function | `helix_charge_periodicity` | 9 | 4 | +0.11 | +0.08 | +0.029 | [-0.521, +0.603] | 1.0000 | 0.912 | **no evidence** — 4 groups |
| function | `helix_hydropathy_periodicity` | 9 | 4 | +0.65 | +0.93 | -0.284 | [-1.320, +0.716] | 1.0000 | 1.692 | **no evidence** — 4 groups |
| type | `helix_charge_periodicity` | 10 | 1 | +1.38 | -0.77 | +2.143 | [+2.143, +2.143] | 1.0000 | inf | **no evidence** — 1 groups |
| type | `helix_hydropathy_periodicity` | 10 | 1 | +0.47 | -0.89 | +1.356 | [+1.356, +1.356] | 1.0000 | inf | **no evidence** — 1 groups |

**0 of 8** primary tests is a finding. A finding needs a Holm-corrected p below 0.05, an interval excluding zero, five resampling groups, an observed difference larger than the minimum detectable one, and matched window lengths — all five, fixed by the plan before any of this ran.

## Secondary and exploratory, reported and never promoted

| arm | statistic | windows | groups | mean z (site) | mean z (control) | difference | 95% interval | could have detected | family |
|---|---|---|---|---|---|---|---|---|---|
| formation | `dominant_hydropathy_period_power` | 10 | 9 | +1.57 | +0.15 | +1.420 | [+0.247, +2.487] | 1.774 | secondary |
| formation | `dominant_period_power` | 10 | 9 | +0.01 | -0.13 | +0.137 | [-0.839, +1.020] | 1.396 | secondary |
| formation | `entropy_rate_k2` | 10 | 9 | -0.12 | +0.28 | -0.401 | [-1.195, +0.380] | 1.192 | secondary |
| formation | `lz_normalised` | 10 | 9 | -0.27 | -0.15 | -0.116 | [-1.045, +0.678] | 1.294 | secondary |
| formation | `strand_charge_periodicity` | 10 | 9 | +0.04 | +0.00 | +0.033 | [-0.703, +0.691] | 1.053 | secondary |
| formation | `strand_hydropathy_periodicity` | 10 | 9 | -0.36 | -0.19 | -0.172 | [-0.773, +0.478] | 0.953 | secondary |
| formation_dense | `dominant_hydropathy_period_power` | 375 | 75 | +0.20 | +0.08 | +0.125 | [-0.048, +0.300] | 0.253 | exploratory |
| formation_dense | `dominant_period_power` | 375 | 75 | +0.01 | -0.16 | +0.167 | [+0.012, +0.316] | 0.218 | exploratory |
| formation_dense | `entropy_rate_k2` | 375 | 75 | -0.09 | -0.24 | +0.151 | [-0.045, +0.333] | 0.275 | exploratory |
| formation_dense | `helix_charge_periodicity` | 375 | 75 | -0.18 | -0.15 | -0.028 | [-0.170, +0.118] | 0.208 | exploratory |
| formation_dense | `helix_hydropathy_periodicity` | 375 | 75 | +0.34 | +0.26 | +0.079 | [-0.121, +0.270] | 0.280 | exploratory |
| formation_dense | `lz_normalised` | 375 | 75 | -0.04 | -0.12 | +0.085 | [-0.081, +0.256] | 0.233 | exploratory |
| formation_dense | `strand_charge_periodicity` | 375 | 75 | +0.02 | -0.09 | +0.106 | [-0.059, +0.277] | 0.250 | exploratory |
| formation_dense | `strand_hydropathy_periodicity` | 375 | 75 | -0.17 | -0.16 | -0.001 | [-0.170, +0.170] | 0.251 | exploratory |
| formation_heldout | `dominant_hydropathy_period_power` | 400 | 63 | +0.12 | -0.01 | +0.135 | [-0.078, +0.340] | 0.287 | secondary |
| formation_heldout | `dominant_period_power` | 400 | 63 | -0.01 | +0.02 | -0.037 | [-0.203, +0.124] | 0.239 | secondary |
| formation_heldout | `entropy_rate_k2` | 400 | 63 | -0.05 | -0.23 | +0.172 | [-0.015, +0.387] | 0.301 | secondary |
| formation_heldout | `lz_normalised` | 400 | 63 | -0.07 | -0.04 | -0.034 | [-0.207, +0.147] | 0.263 | secondary |
| formation_heldout | `strand_charge_periodicity` | 400 | 63 | +0.02 | +0.08 | -0.060 | [-0.245, +0.122] | 0.272 | secondary |
| formation_heldout | `strand_hydropathy_periodicity` | 400 | 63 | -0.05 | -0.14 | +0.088 | [-0.075, +0.254] | 0.228 | secondary |
| formation_matched | `dominant_hydropathy_period_power` | 400 | 70 | +0.22 | +0.11 | +0.103 | [-0.090, +0.311] | 0.283 | exploratory |
| formation_matched | `dominant_period_power` | 400 | 70 | -0.04 | +0.04 | -0.077 | [-0.245, +0.076] | 0.240 | exploratory |
| formation_matched | `entropy_rate_k2` | 400 | 70 | -0.05 | -0.24 | +0.193 | [+0.001, +0.414] | 0.300 | exploratory |
| formation_matched | `helix_charge_periodicity` | 400 | 70 | -0.12 | -0.09 | -0.032 | [-0.200, +0.135] | 0.252 | exploratory |
| formation_matched | `helix_hydropathy_periodicity` | 400 | 70 | +0.36 | +0.26 | +0.098 | [-0.114, +0.313] | 0.310 | exploratory |
| formation_matched | `lz_normalised` | 400 | 70 | -0.05 | -0.18 | +0.136 | [-0.024, +0.293] | 0.237 | exploratory |
| formation_matched | `strand_charge_periodicity` | 400 | 70 | +0.06 | +0.04 | +0.013 | [-0.147, +0.165] | 0.231 | exploratory |
| formation_matched | `strand_hydropathy_periodicity` | 400 | 70 | -0.09 | -0.08 | -0.016 | [-0.172, +0.131] | 0.219 | exploratory |
| formation_matched_dense | `dominant_hydropathy_period_power` | 400 | 75 | +0.12 | +0.00 | +0.123 | [-0.069, +0.317] | 0.267 | exploratory |
| formation_matched_dense | `dominant_period_power` | 400 | 75 | -0.01 | -0.03 | +0.017 | [-0.134, +0.169] | 0.228 | exploratory |
| formation_matched_dense | `entropy_rate_k2` | 400 | 75 | -0.15 | -0.32 | +0.165 | [-0.038, +0.367] | 0.302 | exploratory |
| formation_matched_dense | `helix_charge_periodicity` | 400 | 75 | -0.06 | -0.11 | +0.056 | [-0.123, +0.240] | 0.261 | exploratory |
| formation_matched_dense | `helix_hydropathy_periodicity` | 400 | 75 | +0.34 | +0.26 | +0.081 | [-0.158, +0.299] | 0.321 | exploratory |
| formation_matched_dense | `lz_normalised` | 400 | 75 | +0.01 | -0.08 | +0.093 | [-0.077, +0.256] | 0.241 | exploratory |
| formation_matched_dense | `strand_charge_periodicity` | 400 | 75 | -0.00 | -0.05 | +0.046 | [-0.105, +0.201] | 0.227 | exploratory |
| formation_matched_dense | `strand_hydropathy_periodicity` | 400 | 75 | -0.18 | -0.18 | -0.000 | [-0.154, +0.163] | 0.233 | exploratory |
| function | `dominant_hydropathy_period_power` | 9 | 4 | -0.05 | +0.37 | -0.419 | [-1.483, +0.578] | 1.747 | secondary |
| function | `dominant_period_power` | 9 | 4 | -0.03 | -0.01 | -0.025 | [-0.794, +1.087] | 1.602 | secondary |
| function | `entropy_rate_k2` | 9 | 4 | -0.25 | -0.21 | -0.042 | [-0.768, +0.861] | 1.329 | secondary |
| function | `lz_normalised` | 9 | 4 | +0.62 | -0.07 | +0.686 | [+0.065, +1.403] | 1.151 | secondary |
| function | `strand_charge_periodicity` | 9 | 4 | +0.19 | +0.19 | +0.002 | [-0.778, +0.760] | 1.263 | secondary |
| function | `strand_hydropathy_periodicity` | 9 | 4 | -0.30 | -0.54 | +0.244 | [-0.523, +1.315] | 1.508 | secondary |
| function_dense | `dominant_hydropathy_period_power` | 129 | 16 | +0.40 | -0.08 | +0.475 | [-0.181, +1.260] | 1.055 | exploratory |
| function_dense | `dominant_period_power` | 129 | 16 | -0.13 | +0.00 | -0.133 | [-0.519, +0.249] | 0.558 | exploratory |
| function_dense | `entropy_rate_k2` | 129 | 16 | -0.02 | +0.06 | -0.082 | [-0.452, +0.254] | 0.536 | exploratory |
| function_dense | `helix_charge_periodicity` | 129 | 16 | +0.19 | +0.13 | +0.061 | [-0.363, +0.478] | 0.596 | exploratory |
| function_dense | `helix_hydropathy_periodicity` | 129 | 16 | +0.43 | +0.24 | +0.191 | [-0.393, +0.848] | 0.934 | exploratory |
| function_dense | `lz_normalised` | 129 | 16 | -0.16 | -0.13 | -0.032 | [-0.352, +0.297] | 0.479 | exploratory |
| function_dense | `strand_charge_periodicity` | 129 | 16 | -0.09 | +0.00 | -0.094 | [-0.479, +0.296] | 0.600 | exploratory |
| function_dense | `strand_hydropathy_periodicity` | 129 | 16 | +0.05 | -0.09 | +0.134 | [-0.349, +0.644] | 0.725 | exploratory |
| type | `dominant_hydropathy_period_power` | 10 | 1 | -0.79 | +0.68 | -1.467 | [-1.467, -1.467] | inf | secondary |
| type | `dominant_period_power` | 10 | 1 | +0.73 | +0.60 | +0.133 | [+0.133, +0.133] | inf | secondary |
| type | `entropy_rate_k2` | 10 | 1 | +0.67 | -1.51 | +2.177 | [+2.177, +2.177] | inf | secondary |
| type | `lz_normalised` | 10 | 1 | -1.58 | -1.74 | +0.161 | [+0.161, +0.161] | inf | secondary |
| type | `strand_charge_periodicity` | 10 | 1 | +1.21 | -0.99 | +2.196 | [+2.196, +2.196] | inf | secondary |
| type | `strand_hydropathy_periodicity` | 10 | 1 | +0.09 | -0.75 | +0.840 | [+0.840, +0.840] | inf | secondary |
| type_dense | `dominant_hydropathy_period_power` | 400 | 56 | +0.17 | +0.22 | -0.047 | [-0.321, +0.209] | 0.390 | exploratory |
| type_dense | `dominant_period_power` | 400 | 56 | +0.03 | -0.07 | +0.103 | [-0.140, +0.341] | 0.339 | exploratory |
| type_dense | `entropy_rate_k2` | 400 | 56 | +0.04 | -0.32 | +0.357 | [+0.107, +0.609] | 0.359 | exploratory |
| type_dense | `helix_charge_periodicity` | 400 | 56 | -0.06 | -0.09 | +0.032 | [-0.206, +0.282] | 0.356 | exploratory |
| type_dense | `helix_hydropathy_periodicity` | 400 | 56 | +0.40 | +0.47 | -0.064 | [-0.391, +0.214] | 0.428 | exploratory |
| type_dense | `lz_normalised` | 400 | 56 | -0.02 | -0.17 | +0.156 | [-0.062, +0.350] | 0.294 | exploratory |
| type_dense | `strand_charge_periodicity` | 400 | 56 | -0.02 | -0.07 | +0.047 | [-0.157, +0.260] | 0.297 | exploratory |
| type_dense | `strand_hydropathy_periodicity` | 400 | 56 | -0.22 | -0.23 | +0.012 | [-0.190, +0.233] | 0.310 | exploratory |

## Why a lining produced no window

Amendment 1 section 6. A refused lining is not a pocket without sequence signal; it is a pocket whose positions the window definition cannot carry a period over.

| source | reason | linings |
|---|---|---|
| `pockets/local` | span over the cap | 5443 |
| `pockets/local` | under the floor | 963 |
| `sites/local` | span over the cap | 114 |
| `sites/local` | under the floor | 6 |

## The held-out arm

The pre-registered rule selected windows from 9 resampling groups, and those are exactly the groups whose values were read before amendment 1 was written. `formation_heldout` runs the amendment's rule with all of them removed, so it is a confirmatory test on data no decision in the amendment could have been tuned to. Its Holm family is its own two primary statistics, and the plan's five decision rules apply to it unchanged.

**The held-out claim is only as good as that reproduction.** **The pre-registered arm did not reproduce the disclosed look** (not reproduced): amendment 1 section 1 records 10 windows across 9 groups at a difference of 1.317, and this run gives 10 across 9 at 1.3098553290247468. The pocket table is therefore not the one that was looked at, so the held-out arm cannot be assumed held out from anything and **must be read as exploratory** until the difference is explained.

## What this cannot conclude

* fpocket's rank is a geometric judgement of a predicted structure. The formation arm predicts fpocket, not nature.
* A crystal contact shell is one deposited conformation. The function arm predicts where an inositol phosphate sat in that crystal, not affinity and not specificity.
* Every arm is within-protein, so none of this is a sequence-only predictor of anything, and none of it generalises to sequences with no structure.
* The transfer ceiling in `results/architecture/ARCHITECTURE.md` stands over the periodicity statistics here as it does there.
