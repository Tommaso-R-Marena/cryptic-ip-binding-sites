# Transfer plan, Amendment 1: without the dominant structural group

**Date and status.** 2026-09-24, about 10:00 UTC. Written after the transfer run's
Prepare step (run 35972729955) and before any evaluation, comparison, permutation
or external-test (T3) output was read.

**Why a separate file.** `docs/TRANSFER_PLAN.md` is a trigger path of the running
transfer workflow, and editing it would cancel the run. This amendment is part of
that plan and will be linked from it once the run finishes.

## What Prepare showed

These are counts only; no outcome is known.

| task | development positives | positive groups (sequence / strict) | largest group's share (sequence / strict) |
|---|---|---|---|
| `cryptic_ip_site` | 914 | 174 / 51 | 0.336 / **0.607** |
| `burial` | 914 | 174 / 51 | 0.336 / **0.607** |
| `ip_site` | 2,762 | 390 / 96 | 0.264 / 0.591 |

- **What happened.** Under the strict grouping (sequence links plus Foldseek links at
  TM-score ≥ 0.5, joined as connected components), one group holds 61 % of the
  buried class-ligand positives. Linking structural neighbours transitively merges
  a large structural superfamily into a single component.
- **Primary decision.** Under the plan's 40 % rule, **T1 and T2 are not evaluable**.
  That is the primary decision and is not changed.

## The single secondary analysis (T1b, T2b)

1. **Remove.** Every row of the one strict group with the most `cryptic_ip_site`
   development positives is removed, from development and holdout alike. No other
   group is removed, and no other variant will be tried. The group's identity and
   its entries are reported.
2. **Recompute.** The development positives, positive groups and largest-group
   shares are recomputed on what remains, for both groupings.
3. **Evaluability.** If a largest share still exceeds 0.40 under either grouping,
   T1b/T2b are **not evaluable**, and nothing further is done.
4. **Otherwise, evaluate:**
   - the same paired runs as the plan: sequence repeats 0-2, strict repeat 0,
     the holdout on repeat 0, and 10 permutations per task;
   - decided by the same code (`scripts/transfer.py decide` logic: the benchmark's
     rule with the ten-permutation control).
5. **Reporting.** The result is reported as *hull depth for buried phosphate-dense
   sites outside the dominant structural group*, a narrower claim than T1/T2, and
   always next to the primary "not evaluable".
6. **Exclusions.** `ip_site` and T3 are unaffected.
