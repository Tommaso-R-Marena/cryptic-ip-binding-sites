# Architecture plan, amendment 1: the first real linings

Dated 2026-10-12, **after `docs/ARCHITECTURE_PLAN.md` and before any real
lining has been passed to `describe_architecture`**. Written because the study
is about to read real data for the first time and the parent plan was written
when none was reachable.

## What changed

The parent plan assumed the lining positions would come from the IP benchmark
table, which lives in a CI artefact this container cannot fetch. Three tables
of real results are **committed in the repository** and were found by looking
rather than by fetching:

| table | rows | what it carries |
|---|---|---|
| `results/specificity/candidates.csv` | 75 | `top_pocket_residues` — real pocket lining positions — with `cluster`, `plddt_mean`, `annotated`, `explained` |
| `results/redocking/census.csv` | 662 | real IHP copies: `pdb_id`, `chain`, `resseq`, `homology_group_strict`, `relative_sasa`, `n_contact_chains`, `interface` |
| `results/candidates/candidates_summary.csv` | 1 | `coordinating_residues` for the top candidate |

So the study can run on real linings now, with one decisive restriction.

## What can and cannot be asked of this data

**Architecture (Q7) needs only positions, and positions are what these tables
have.** Every statistic in `cryptic_ip/detect/architecture.py` is a function of
the lining position set alone — the gap multiset and the order of the gaps —
so Q7 is runnable offline on real data today.

**Sequence (Q6) needs letters, and there are none here.** No committed table
carries a chain sequence, and the container cannot reach UniProt or RCSB. Q6
therefore stays as `docs/SEQUENCE_PLAN.md` has it: the inputs builder emits a
sequence-inputs file only when sequences are supplied, and the study job still
refuses to run without them.

**These are candidates, not annotated sites.** `annotated` is False for all 75
rows: they are proteome-screen candidates, so this is **not** the site-versus-
control comparison the parent plan pre-registers. It is a descriptive first
look, and nothing it produces is evidence that annotated IP sites have a
sequence-space architecture.

## What will be computed, fixed before any result is read

1. **Class distribution.** `architecture_class` and `signature_verdict` over the
   75 real linings, with the counts of `unresolved` reported separately from
   `none`, as the parent plan requires.
2. **Against the measured null.** The same quantities on 500 random linings gave
   3.0% structured at p ≤ 0.05 and classes 449 distributed / 44 segmental /
   7 mosaic / 0 phase-locked. Real linings will be compared with that baseline
   and with a **matched decoy arm** drawn to the same size and span as each
   real lining, which is the comparison that controls for both.
3. **Which half carries it**, since the halves were shown to swap with pocket
   shape. Reported separately, never pooled.
4. **One stratified comparison**, and only one: `explained` (30 True, 45 False)
   — pockets a competitor polyanion already accounts for, against the rest.
   Group-resampled over `cluster`, Holm across the features, with the standing
   five-group floor.

## Decision rules, unchanged from the parent plan

* A difference that tracks lining size is not a finding; size is reported as a
  covariate and arms are checked for imbalance on it.
* A null is reported with the detection rate measured at that lining size, so
  "nothing" stays separate from "could not have seen anything".
* `unresolved` is never pooled with `none`.
* **Nothing here is evidence about annotated IP sites.** The strongest available
  claim is a description of what candidate pocket linings look like in the
  invariant space, beside a null.

## The builder, and what it refuses

`scripts/build_lining_inputs.py` turns these tables into study inputs. Two
refusals are built in rather than documented:

* **Chain-blind positions are refused for multi-chain sources.** The repo's
  pocket residue numbers are chain-blind by their own docstring. For an
  AlphaFold model that is harmless because there is one chain; for a PDB entry
  it is not, and the builder will not emit a window for a multi-chain source
  without an explicit chain column.
* **No sequence, no sequence inputs.** The architecture file is emitted
  whenever positions are present; the sequence file is emitted only when
  sequences are supplied, and never with invented or placeholder letters.
