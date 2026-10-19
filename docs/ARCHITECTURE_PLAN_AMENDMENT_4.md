# Architecture plan, amendment 4: selecting window sizes from the null, which is free

Dated 2026-10-13, **before any saturation-pruned statistic has been computed on
any lining**. Written after amendment 3's measurement was read, and the reason
it exists is that amendment 3 failed for a reason its own diagnostic names
exactly.

## What amendment 3 measured

The windowed statistic, with window sizes fixed in advance at 6, 8 and 10
lining positions, detected **nothing** — 0 of the fixtures at every geometry
arm, against 0.83 for the pooled statistic on ideal geometry and 0.50 at
realistic jitter. The declared saturation diagnostic says why, and it says it
in one number: the **mean of the permutation null at window size 6 is 0.90**.
A random reordering of the same gaps already reaches a coherence of 0.9 over
six positions, so no observation can look unusual against it. The statistic
takes a maximum over all three sizes, so the saturated short window sets the
null for the whole statistic and the longer windows — null means 0.81 and 0.74
— cannot rescue it.

This is the same failure as the per-element statistic, for the same reason: a
band scan over a short run of positions fits almost any spacing. Amendment 3
anticipated it and fixed the window sizes anyway, because tuning them against
the result is precisely the move that invalidates the test.

## The correction, and why it costs nothing

Tuning window sizes against the *observation* would invalidate the p-value.
Selecting them from the **null** does not, and the reason is exact rather than
approximate:

> Under the register null the gap multiset is held fixed and only its order
> varies. The resultant's modulus is invariant to translation, so the null
> distribution of a window's maximum resultant depends **only on the gap
> multiset** — not on the order, and therefore not on the observation. A rule
> that picks window sizes by their null mean is consequently a fixed function
> of the quantity the null holds constant, so it selects the *same* window set
> for the observation and for every permutation. A permutation test whose
> statistic is the same function on every draw is exactly valid.

So the pruned statistic is not a tuned statistic; it is a different statistic
whose definition happens to be computed from the invariant. The selection
cannot look at the observed ordering because the information it uses does not
contain it.

One technical point, which is the whole of the care needed. The null mean has
to be *estimated* when the orderings are too many to enumerate, and an
estimate computed from the same draws that form the null is a statistic of
those draws. To keep the argument exact, the selection is estimated from a
**pilot** permutation sample with its own seed, drawn before and independently
of the null sample, and the pilot seed is part of the statistic's definition.
With exact enumeration no pilot is needed and the selection is deterministic.

## The statistic, fixed before it is computed

1. For each window size in `WINDOW_SIZES`, estimate the null mean of that
   size's maximum resultant from `PILOT_DRAWS` permutations of the gap
   multiset, using `pilot_seed`.
2. **Keep the sizes whose pilot null mean is below `SATURATED_NULL = 0.9`**,
   the threshold already defined and already reported in amendment 3. Nothing
   else about the statistic changes: the maximum is over the kept sizes and
   every band probe, and the identical maximum is taken on every draw of the
   same null.
3. If **no** size survives, the statistic is **unresolved** — not a p of 1, and
   not a silent fallback to the saturated set. A lining whose every window
   length saturates is a lining this question cannot be asked of.
4. `SATURATED_NULL` is not to be moved. It was fixed at 0.9 in amendment 3,
   before any pruned statistic existed, and the pruning inherits it.

## What will be measured, and the decision rule

Exactly amendment 3's design, with the pruned statistic as a third arm beside
the pooled and unpruned ones, on the same fixtures at the same seeds:

* **Primary:** the paired difference in detection rate at the realistic jitter
  arm, pruned minus pooled, by exact McNemar on the discordant pairs.
* **Secondary:** the same on the ideal and bulge arms; the decoy
  false-positive rate; which sizes survive pruning, per arm; and how often the
  statistic is unresolved because nothing survived.

**Decision rule, unchanged from amendment 3 and fixed now.** Adopted as a
reported descriptor only if it detects strictly more than the pooled statistic
at the realistic arm *and* its decoy false-positive rate is at or below 0.05.
Detecting more at a higher false-positive cost is reported as a measured trade
and not adopted. Detecting no more is a negative and is recorded as one.

**A specific prediction, recorded so it can be wrong.** Pruning removes the
saturated short windows, so the pruned statistic should behave like the
long-window statistic: better than the unpruned windowed one, and most likely
still **worse than pooled**, because a 10-position window of a 14-position
lining is nearly the pooled statistic with fewer positions and the same scan.
If that is what happens, the conclusion to draw is the one amendment 3 already
pointed at — the limit on transfer is the band scan's flexibility on short
runs, not the global-phase assumption — and no further windowed variant is
worth trying.

## What this cannot conclude

* The validity argument is about the **p-value**, not about power. A valid test
  with no power is still useless, and the measurement is the only thing that
  can say which this is.
* The pilot is a sample. Two pilot seeds can select different window sets for
  the same lining when a size sits near the threshold; the statistic is
  therefore defined relative to its pilot seed, which is reported with it.
* Nothing here changes what the detector scores, and this amendment does not
  pre-authorise revisiting that decision.
