# Does a pocket lining's composition carry anything — and can this pipeline tell?

**Run.** `scripts/pocket_composition.py analyse`, seed 20261025, 2000 subset draws per lining.

**Pre-registered** in `docs/COMPOSITION_PLAN.md`, dated before any composition statistic had been computed on any real lining in this project.

## The gate

**The positive control fired.** Annotated inositol-phosphate contact shells against same-size subsets of their own chains, on `basic_fraction`: mean z **+4.28** across 13 homology groups on 95 shells, interval [+3.590, +4.997], p 0.0009995, against a minimum detectable 1.055.

That a contact shell is arginine-rich is textbook and **is not offered here as a finding**. What it establishes is the thing study Q9 could not: that these linings, in these structures, through this machinery, can produce a detection at all. Q9's nulls are therefore statements about **order**, not about the pipeline's sensitivity.

4267 candidate-pocket linings, 95 annotated shells and 5088 pockets of the same deposited structures carried a null. The window is the **lining itself**, gap-closed: every statistic here is a function of the multiset of its letters, so closing the gaps cannot change one, and `cryptic_ip/detect/signature.py` is the fence that keeps periodicity statistics off that same window.

## Primary

| arm | statistic | linings | groups | mean z (site) | mean z (control) | difference | 95% interval | Holm | could have detected | finding |
|---|---|---|---|---|---|---|---|---|---|---|
| formation | `basic_fraction` | 400 | 74 | +0.07 | +0.27 | -0.198 | [-0.366, -0.039] | 0.0810 | 0.246 | no |
| formation | `composition_divergence` | 400 | 74 | +0.44 | +0.32 | +0.117 | [-0.083, +0.321] | 0.4908 | 0.285 | no |
| formation | `nonbasic_divergence` | 400 | 74 | +0.39 | +0.38 | +0.013 | [-0.155, +0.188] | 0.8836 | 0.251 | no |
| function | `basic_fraction` | 95 | 10 | +3.92 | -0.07 | +3.994 | [+3.239, +4.747] | 0.0030 | 1.137 | **yes** |
| function | `composition_divergence` | 95 | 10 | +2.21 | +0.11 | +2.097 | [+1.401, +2.845] | 0.0030 | 1.100 | **yes** |
| function | `nonbasic_divergence` | 95 | 10 | +2.59 | +0.07 | +2.519 | [+1.658, +3.358] | 0.0030 | 1.289 | **yes** |
| type | `basic_fraction` | 400 | 54 | +0.11 | +0.32 | -0.207 | [-0.473, +0.056] | 0.2519 | 0.372 | no |
| type | `composition_divergence` | 400 | 54 | +0.59 | +0.27 | +0.317 | [+0.106, +0.534] | 0.0135 | 0.310 | **yes** |
| type | `nonbasic_divergence` | 400 | 54 | +0.58 | +0.41 | +0.172 | [-0.034, +0.399] | 0.2519 | 0.321 | no |

**4 of 9** primary tests are findings. A finding needs the positive control to have fired, a Holm-corrected p below 0.05, an interval excluding zero, five resampling groups, an observed difference above the minimum detectable one, and matched lining sizes — all six, fixed before this ran.

## The control arm's other statistics

Reported because they say what a shell's composition looks like against its own chain, and never promoted: this arm is an instrument check, not a test.

| statistic | shells | groups | mean z | 95% interval | p |
|---|---|---|---|---|---|
| `acidic_fraction` | 95 | 13 | -0.51 | [-0.891, -0.031] | 0.06347 |
| `aromatic_fraction` | 95 | 13 | +0.08 | [-0.300, +0.458] | 0.6857 |
| `basic_fraction` | 95 | 13 | +4.28 | [+3.590, +4.997] | 0.0009995 |
| `composition_divergence` | 95 | 13 | +2.63 | [+1.866, +3.400] | 0.0009995 |
| `glycine_fraction` | 95 | 13 | -0.08 | [-0.500, +0.397] | 0.7541 |
| `mean_hydropathy` | 95 | 13 | -2.58 | [-3.029, -2.079] | 0.0009995 |
| `net_charge_per_residue` | 95 | 13 | +3.19 | [+2.588, +3.818] | 0.0009995 |
| `nonbasic_divergence` | 95 | 13 | +3.07 | [+2.130, +4.016] | 0.0009995 |
| `proline_fraction` | 95 | 13 | -0.78 | [-0.872, -0.687] | 0.0009995 |

## Secondary, reported and never promoted

| arm | statistic | linings | groups | difference | 95% interval | could have detected |
|---|---|---|---|---|---|---|
| formation | `acidic_fraction` | 400 | 74 | -0.203 | [-0.355, -0.053] | 0.220 |
| formation | `aromatic_fraction` | 400 | 74 | +0.335 | [+0.158, +0.514] | 0.246 |
| formation | `glycine_fraction` | 400 | 74 | -0.172 | [-0.316, -0.021] | 0.209 |
| formation | `mean_hydropathy` | 400 | 74 | +0.650 | [+0.391, +0.896] | 0.368 |
| formation | `net_charge_per_residue` | 400 | 74 | +0.010 | [-0.119, +0.137] | 0.183 |
| formation | `proline_fraction` | 400 | 74 | -0.194 | [-0.341, -0.048] | 0.220 |
| function | `acidic_fraction` | 95 | 10 | -0.427 | [-0.935, +0.184] | 0.876 |
| function | `aromatic_fraction` | 95 | 10 | -0.046 | [-0.893, +0.712] | 1.215 |
| function | `glycine_fraction` | 95 | 10 | +0.531 | [-0.045, +1.184] | 0.926 |
| function | `mean_hydropathy` | 95 | 10 | -2.443 | [-3.205, -1.717] | 1.136 |
| function | `net_charge_per_residue` | 95 | 10 | +2.955 | [+2.366, +3.559] | 0.909 |
| function | `proline_fraction` | 95 | 10 | -0.589 | [-0.857, -0.322] | 0.403 |
| type | `acidic_fraction` | 400 | 54 | -0.180 | [-0.381, +0.043] | 0.304 |
| type | `aromatic_fraction` | 400 | 54 | +0.256 | [+0.041, +0.493] | 0.322 |
| type | `glycine_fraction` | 400 | 54 | -0.163 | [-0.307, -0.020] | 0.204 |
| type | `mean_hydropathy` | 400 | 54 | +0.456 | [+0.149, +0.760] | 0.439 |
| type | `net_charge_per_residue` | 400 | 54 | -0.026 | [-0.207, +0.144] | 0.253 |
| type | `proline_fraction` | 400 | 54 | -0.165 | [-0.342, +0.007] | 0.250 |

## Why a lining produced no window

| source | reason | linings |
|---|---|---|
| annotated shells | under the floor | 34 |
| candidate pockets | under the floor | 4140 |
| deposited pockets | under the floor | 4537 |
| deposited pockets | position not in the model | 29 |

## What this cannot conclude

* **A compositional difference is not a discovery here.** That these sites are basic is known; the control arm exists to test the instrument and no row of it is a finding about inositol phosphate.
* Composition and order do not exhaust a sequence. Conservation, coevolution, linear motifs and learned representations are outside this study and outside study Q9, and a null in both says nothing about them.
* Every arm is within-protein or within-chain. **Nothing here is a sequence-only predictor of anything**, and this study does not attempt one.
* The candidate arms read a predicted model and the function arm reads one deposited conformation, so neither is a statement about solution behaviour.
