# Study Q7 on real pocket linings

**Run.** `scripts/real_architecture.py run`, seed 20261012, on `results/architecture/real/architecture_inputs.jsonl`. Every figure is computed by that command and written in the same pass.

**Pre-registered** in `docs/ARCHITECTURE_PLAN.md` as amended by `docs/ARCHITECTURE_PLAN_AMENDMENT_1.md`, which fixed what is computed here before any of it was read.

**These are candidates, not annotated sites.** `annotated` is False for every row, so nothing below is evidence that annotated IP sites have a sequence-space architecture. It is a description of real candidate linings beside a null.

75 of 75 linings measurable; size 7–75 positions, median 15.

## Real linings against decoys matched on size and span

| arm | topology | register | **either** |
|---|---|---|---|
| real (75) | 0.31 | 0.11 | **0.33** |
| decoy (375) | 0.02 | 0.03 | **0.04** |

Each decoy has the same number of positions and the same span as the real lining it was drawn for, so the comparison is not confounded by either.

## Architecture classes

| class | real | decoy |
|---|---|---|
| `distributed` | 48 | 361 |
| `mosaic` | 23 | 6 |
| `segmental` | 4 | 8 |

Verdicts: clear 15, none 50, weak 10. `unresolved` is reported separately from `none`, as the plan requires: they are different claims.

Secondary-structure elements called: helix_face 17, none 52, periodic_other 4, strand_face 2

## The one pre-registered split

Pockets a competitor polyanion already explains, against the rest. Group means over homology clusters, Holm across the features, five-group floor.

| feature | arms | difference | 95% interval | p | Holm | evidence | finding |
|---|---|---|---|---|---|---|---|
| `either_p` | False vs True | -0.015 | [-0.198, +0.159] | 0.8661 | 1.0000 | yes | no |
| `occupancy` | False vs True | +0.004 | [-0.043, +0.051] | 0.8726 | 1.0000 | yes | no |
| `gap_entropy` | False vs True | +0.050 | [-0.024, +0.117] | 0.2139 | 0.8556 | yes | no |
| `segments` | False vs True | -0.683 | [-1.855, +0.475] | 0.2859 | 0.8576 | yes | no |

## What this shows, and the control it does not have

The decoys control size and span by construction, so the gap between the arms is not either of those. What it is, though, is **nearly tautological**: a random set of positions in a span is not a protein, and *any* real pocket lining - IP site or not - comes from residues that fold together and will therefore look structured against that null. This result establishes that the instrument registers real linings and gives the rate at which it does. It does **not** establish anything specific to inositol-phosphate sites.

The control that would is **other pockets from the same proteins**: the same folds, the same screen, the same models, differing only in being the pocket the screen did not pick. No committed table carries them, so that comparison is named here as missing rather than approximated.

Three further limits:

* Not evidence about annotated IP sites: every lining here is a candidate, and the screen that chose them already favoured buried, basic pockets, so the selection and the signal cannot be separated in this data.
* Not a predictor. Separating classes here is not prediction until it is tested on groups held out from fitting, which this does not do.
* Lining size is reported above because a difference that tracks it is not a finding, and the decoys control for it by construction.
