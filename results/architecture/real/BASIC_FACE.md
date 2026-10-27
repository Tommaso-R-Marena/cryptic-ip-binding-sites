# Do a pocket's basic residues sit on a face?

**Run.** `scripts/basic_face.py run`, seed 20261013, on `results/triage/triaged_candidates.csv`. Every figure is computed by that command in the same pass.

**Pre-registered** in `docs/ARCHITECTURE_PLAN_AMENDMENT_2.md`, dated before any basic-residue subset was measured.

## The design, and why it is better than amendment 1's

Amendment 1 compared real linings with random position sets and said plainly that the comparison is nearly tautological. Here each pocket's basic residues are compared with **random subsets of that same pocket's lining, of the same size**. Protein, fold, lining and subset size are held fixed by construction, so a difference can only be which residues were chosen.

The per-pocket statistic is the quantile of the subset among 50 matched draws, uniform on [0, 1] under the null that basic residues sit at random within the lining. **Below 0.5 is more structured than chance.**

99 subset measurements from 75 pockets.

Excluded, by cause:

- 20 — basic_positions: fewer than 6 residues
- 1 — basic_positions: unmeasurable subset or lining
- 30 — conserved_positions: fewer than 6 residues

## Result

An **unresolved** half is not scored as a p of 1. The architecture module separates "cannot be tested" from "no structure", and an earlier version of this script collapsed them by substituting 1.0 — which put half of all register quantiles at exactly 1.000 and none at 0.000, because an unresolved subset tied with its equally unresolved draws. The quantile is also a **mid-p** quantile, splitting ties, because permutation p-values are discrete and counting ties as "at or below" pushes every quantile upward.

| subset | half | pockets | ties | mean quantile | difference from 0.5 | 95% interval | p | Holm | finding |
|---|---|---|---|---|---|---|---|---|---|
| `basic_positions` | topology | 54 | 0.12 | 0.435 | -0.065 | [-0.151, +0.020] | 0.1564 | 0.6257 | no |
| `basic_positions` | register | 54 | 0.21 | 0.536 | +0.036 | [-0.049, +0.116] | 0.3923 | 1.0000 | no |
| `conserved_positions` | topology | 45 | 0.13 | 0.488 | -0.012 | [-0.108, +0.077] | 0.8131 | 1.0000 | no |
| `conserved_positions` | register | 45 | 0.21 | 0.527 | +0.027 | [-0.059, +0.111] | 0.5407 | 1.0000 | no |

Holm across 4 tests: two subsets times two halves, fixed by the amendment before any of it was read.

## How to read a null here

The transfer measurement in `results/architecture/ARCHITECTURE.md` stands over all of this: at the φ/ψ scatter real structures show, the method detects about **a third** of what it detects on idealised geometry. A null in this table is therefore weak evidence of absence, and the honest statement is that the design could not see an effect of the size it could see on ideal geometry.

## What this cannot conclude

* These are candidates, not annotated sites.
* Not evidence of binding: a basic face is a suggestive feature, not an interaction.
* Nothing about whether these are real pockets — the design conditions on the lining and therefore cannot ask that.
