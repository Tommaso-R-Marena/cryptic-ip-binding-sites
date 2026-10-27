# Architecture plan, amendment 3: attacking the transfer ceiling with a local statistic

Dated 2026-10-13, **before any windowed statistic has been computed on any
lining, fixture or real**. Written after amendment 2's analysis was read and
before this one is run.

## The limit this addresses, and why it is the right one to attack

The transfer measurement is the largest number in this line and the worst one:
with each residue's φ and ψ perturbed — leaving bond lengths ideal to 1e-15, so
that only realism changes — detection of a planted face falls from **0.73 on
ideal geometry to 0.33 at the 8° scatter real structures show**, and to
**0.20** once bulges interrupt the face. Every fixture number in
`results/architecture/ARCHITECTURE.md` is therefore an upper bound about twice
the realistic value, and the real-lining rate in
`results/architecture/real/REAL.md` sits exactly at that realistic 0.33.

Two thirds of the method's measured sensitivity is an artefact of ideal
fixtures. No other change available here moves a number that large.

## The mechanism, stated so that it can be wrong

The pooled resultant is the Neyman–Pearson statistic for a face of **one**
period and **one** unknown phase, and that is a property of the derivation, not
of the code: marginalising a single nuisance phase over a uniform prior gives
`I₀(κ·n·R̄)`. A real helix with 8° of φ/ψ scatter does not have one phase.
The rise per residue drifts, so phase accumulates error along the chain, and a
bulge resets it outright. The pooled statistic asks every residue of the lining
to agree with every other, and over a 100-residue span that is a stronger
demand than the structure can meet.

**The claim to be tested:** the face is still locally phase-locked. If phase
drifts slowly, a resultant taken over a *window* of consecutive lining
positions still concentrates, while the pooled one does not. A statistic that
requires only piecewise phase lock should therefore recover power at realistic
geometry, and should recover nothing on a decoy, because a decoy has no face to
be locally coherent about.

There is a documented reason to expect this to fail, and it is the reason this
amendment is needed rather than obvious. The same derivation already predicted
that a **per-element** statistic would beat the pooled one on multi-element
linings, and it did not: 0 of 12 against 10 of 12, because a band scan over a
run of two or three residues fits almost any spacing, so both the statistic and
its permutation null saturate. A windowed statistic is the same idea at a
larger window, so the measurement to make is whether there is a window size
large enough not to saturate and small enough to tolerate drift. **That
interval may be empty, and if it is, this amendment records a second negative
of the same kind.**

## The statistic, fixed before it is computed

For a lining of sorted positions `p₁ < … < pₙ`:

1. A **window** is a run of `k` consecutive positions *by index*, for every
   offset, with `k` taken from `WINDOW_SIZES = (6, 8, 10)` and `k ≤ n`. Windows
   are indexed in the lining, not in the sequence, so a window spans whatever
   stretch of chain those residues occupy.
2. The windowed statistic is the **maximum over every window and every band
   probe** of the window's mean resultant length. The maximum over window
   sizes is *inside* the statistic, so there is no selection to correct for
   afterwards.
3. The null is the **same gap-multiset permutation null** already used for
   register — it holds topology, n, span, occupancy and gap entropy exactly
   fixed — with the identical maximum taken on every permutation. A maximum
   over more things has a larger null, and computing the null on the same
   maximum is the whole of what makes the p-value calibrated. Exact
   enumeration below `MAX_EXACT_ORDERINGS` orderings, as now.
4. `MIN_ORDERINGS` still applies: below it the windowed half is **unresolved**,
   which is not a p of 1.

No other parameter is introduced. `WINDOW_SIZES` is fixed here, in advance,
and is not to be tuned against the result; the lower end is the existing
`MIN_POSITIONS`, and the upper end is where a window stops being local on a
15-residue median real lining.

## What will be measured, and the decision rule

The comparison is **paired by fixture**: the same planted faces, the same
seeds, the same three geometry arms already in the transfer probe.

| arm | geometry |
|---|---|
| ideal | exact 3.6-residue helices |
| realistic | φ/ψ jitter 8° |
| worst | bulges + φ/ψ jitter 8° |

Primary estimand: **the paired difference in detection rate at the realistic
arm**, windowed minus pooled, with detection meaning the half's permutation p
at or below `CLEAR_P` after Holm within the windowed statistic's own family.
The windowed statistic is reported as its **own test with its own Holm family**
and is not folded into the existing four, so neither detector's published
numbers change.

Secondary, reported whatever the primary shows:

* the same paired difference on the ideal and worst arms;
* the **false-positive rate on decoys** — linings with no planted face, matched
  in size and span — which is the cost side of the trade;
* the saturation diagnostic that killed the per-element statistic: the **mean
  of the permutation null** for each window size, since a null mean near 1 is
  saturation and makes the test powerless by construction.

**Decision rule, fixed now.** The windowed statistic is adopted as a reported
descriptor only if, at the realistic arm, it detects strictly more than the
pooled statistic on the paired fixtures *and* its decoy false-positive rate is
at or below 0.05. If it detects more but costs more false positives, it is
reported as a measured trade and not adopted. If it detects no more, that is a
negative and it is recorded as one, beside the per-element negative, with the
null means that explain it.

Nothing here changes what the detector scores. The decision recorded in
`docs/DETECTOR.md` — that architecture descriptors are not features of the
score — stands unless the realistic-arm detection rate rises enough to revisit
it, and this amendment does not pre-authorise that revision.

## What this cannot conclude

* A fixture is not a protein. Recovering power on jittered ideal helices shows
  the statistic tolerates drift, not that it works on real backbones.
* The jitter model perturbs φ and ψ independently and identically, which real
  backbones do not: real scatter is correlated along the chain and concentrated
  at specific positions. The arm is a stress test, not a sample of nature.
* A bulge removes periodicity rather than hiding it. No statistic can recover
  a signal that is not there, and the worst arm is included to show where that
  floor is, not in the expectation of beating it.
