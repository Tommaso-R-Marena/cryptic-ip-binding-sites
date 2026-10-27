# Architecture plan, amendment 5: a rank where the threshold could not fire

Dated 2026-10-19, **before any fixture arm, decoy or real lining has been
measured with the rule this amendment defines**.

One disclosure, because the record is worth more than the appearance of one:
the implementation was exercised once on a single synthetic twelve-position
lining as a smoke test that the code path runs, before this file was written.
No arm of the probe, no decoy and no real lining had been scored, and the
comparison this amendment pre-registers is untouched by it.

## What amendment 4 measured, and the hole it left

Amendment 4's pruning rule — keep the window sizes whose permutation null mean
is below `SATURATED_NULL = 0.9` — **never fired**. On the barrel slabs the
size-6 null mean is 0.88, just under the threshold that amendment 3 had fixed
against a 0.90 measured on a different lining. So the pruned arm computed the
same statistic as the unpruned one, both detected nothing, and the question
amendment 4 actually wanted to ask was never put:

> does a local statistic restricted to the window lengths that are **not**
> saturated recover power at realistic jitter?

Moving the threshold now would be the tuning both amendments forbade. The
threshold stays at 0.9 and is not touched.

## The rule, which is the same argument with a rank instead of a cutoff

Amendment 4 established that selecting window sizes from the **null** is free:
the resultant's modulus is translation invariant, so a window's null
distribution depends only on the gap multiset, which the register null holds
fixed; a rule stated in terms of it picks the same window set for the
observation and for every permutation, and a permutation test whose statistic
is the same function on every draw is exactly valid.

That argument never required a threshold. It covers **any** function of the
gap multiset, so it covers a rank:

> **`best`: keep the single window size with the lowest pilot null mean.**

A rank cannot fail to fire. It selects on exactly the quantity the saturation
diagnostic measures — the room left between the null and 1, which is where
power can live — and it uses the same pilot sample, with its own seed, so the
selection still cannot be a statistic of the null sample it will be judged
against.

Nothing else changes. `WINDOW_SIZES` stays `(6, 8, 10)`, `SATURATED_NULL` stays
0.9 and is still reported, the null is still the gap-multiset permutation null,
and `MIN_ORDERINGS` still makes a lining with too few orderings *unresolved*
rather than null.

## What will be measured, and the decision rule

The same three geometry arms, the same fixtures, the same seeds, paired, with
`best` as a fourth column beside pooled, `all` and `unsaturated`.

* **Primary:** the paired difference in detection rate at the realistic jitter
  arm, `best` minus pooled, by exact McNemar on the discordant pairs.
* **Secondary:** the ideal and bulge arms; the decoy false-positive rate; which
  size the rank selects, per arm, since a rule that always picks 10 is a
  long-window statistic by another name and should be described as one.

**Decision rule, inherited unchanged.** Adopted as a reported descriptor only
if it detects strictly more than the pooled statistic at the realistic arm
*and* its decoy false-positive rate is at or below 0.05. More detections at a
higher false-positive cost is a measured trade, not an adoption. No more
detections is a negative and is recorded as one.

**The prediction, recorded so it can be wrong.** `best` will beat `all` and
`unsaturated` — both of which are dragged down by the saturated short window —
and will still lose to pooled, because the size the rank selects will be 10,
and a 10-position window of a 14-position lining is the pooled statistic minus
four positions and plus a maximum over five offsets. If that is what happens,
the conclusion is the one amendment 3 pointed at and this amendment closes:
**locality is not what limits transfer**, the band scan's flexibility on short
runs is, and no windowed variant of this family is worth trying next.

If instead `best` beats pooled at the realistic arm, the conclusion is the
opposite and far more interesting — the pooled statistic's global-phase
assumption *is* what costs the transfer — and the decision in `docs/DETECTOR.md`
not to wire architecture descriptors into the scorer would have to be revisited
on a held-out set before anything was changed.

## What this cannot conclude

* A fixture is not a protein, and recovering power on jittered ideal helices
  would show the statistic tolerates drift, not that it works on real
  backbones.
* The rank is defined relative to its pilot seed, which is reported with every
  result, and two seeds can rank adjacent sizes differently when their null
  means nearly coincide.
* Three arms of thirty fixtures each. The McNemar is exact, but a difference of
  one or two linings is not a finding at this size, and the amendment does not
  promote one into a finding.
