# Study Q10: is the gate calibrated?

**Run.** `scripts/pocket_composition.py validate`, seed 20261025, 600 subset draws, 600 bootstrap resamples.

The positive control decides how every other row of study Q10 — and every null of study Q9 — may be read. So it is calibrated on synthetic pools **before any real lining reaches it**, and this page says whether it opens on a signal that is there and stays shut on one that is not.

Pools: 18 copies, chains of 400 residues, linings of 14, and in the planted pool 8 of those residues made arginine or lysine. The noise pool's linings come from the same chains by the same rule, so a non-zero mean z in it would be bias rather than signal.

| pool | mean z | p | 95% interval | gate |
|---|---|---|---|---|
| planted basic patch | **+5.99** | 0.008319 | [+5.700, +6.333] | **opens** |
| matched noise, seed 20261125 | +0.05 | 0.8769 | [-0.403, +0.491] | stays shut |
| matched noise, seed 20261126 | +0.07 | 0.7687 | [-0.314, +0.476] | stays shut |
| matched noise, seed 20261127 | +0.09 | 0.6689 | [-0.285, +0.476] | stays shut |
| matched noise, seed 20261128 | +0.36 | 0.2479 | [-0.160, +0.933] | stays shut |
| matched noise, seed 20261129 | +0.21 | 0.3877 | [-0.130, +0.577] | stays shut |
| matched noise, seed 20261130 | -0.14 | 0.5258 | [-0.523, +0.245] | stays shut |
| matched noise, seed 20261131 | -0.07 | 0.7421 | [-0.408, +0.255] | stays shut |
| matched noise, seed 20261132 | -0.17 | 0.5141 | [-0.607, +0.300] | stays shut |
| matched noise, seed 20261133 | +0.18 | 0.5973 | [-0.515, +0.876] | stays shut |
| matched noise, seed 20261134 | -0.22 | 0.3161 | [-0.589, +0.139] | stays shut |
| matched noise, seed 20261135 | -0.05 | 0.8586 | [-0.482, +0.323] | stays shut |
| matched noise, seed 20261136 | +0.27 | 0.2113 | [-0.075, +0.618] | stays shut |

Across the 12 noise pools the mean of the mean-z values is **+0.0470**, against an expected zero, with a spread of 0.1846; **0** of them opened the gate.

| check | result |
|---|---|
| `a_wrong_direction_rejection_is_refused` | pass |
| `no_noise_pool_opens_the_gate` | pass |
| `noise_is_centred_on_zero` | pass |
| `planted_pool_opens_the_gate` | pass |
| `planted_separates_from_every_noise_pool` | pass |

**What a noise pool with a significant p but a negative mean z would be.** The plan fixes the control's direction in advance — a contact shell is *more* basic than subsets of its own chain — so a rejection in the other direction does not open the gate. Over a dozen pools one such pool is ordinary, and the check `a_wrong_direction_rejection_is_refused` is what stops it counting.

**This says nothing about inositol phosphate.** It is a statement about the instrument, on synthetic sequence, and it is the statement study Q9 never made about its own.
