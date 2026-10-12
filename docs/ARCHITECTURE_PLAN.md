# Lining-architecture plan: pre-registration for study Q7

Dated 2026-10-12, **before any annotated lining has been read by any function
in `cryptic_ip/detect/architecture.py`**. The module exists, is unit-tested,
and has been run on hand-built position sets and on idealised helix barrels
from `cryptic_ip.detect.geometry_fixtures` whose lining is known from the 3-D
construction. Every number below is a property of the method or of an
idealised fixture, measured without looking at any outcome, and
`results/architecture/ARCHITECTURE.md` is the committed run that produced
them. Changes after today are dated amendments written before the results they
could affect are read.

This plan covers study Q7 only. It does not amend studies A–P, Q1–Q5, or Q6
(`docs/SEQUENCE_PLAN.md`), which asks a different question about the same
windows: Q6 is about the *letters* at a site, Q7 is about the *positions*.

## The question

Is there something deep about sequences and pocket formation — and about
pocket *type* — or is a lining's sequence-space structure just a restatement of
its size?

## The framework, which is a small theorem

A pocket is the statement that residues far apart along the chain are close
together in space. So the object is the **set of lining positions**, and it
decomposes exactly.

**The complete H₀ invariant of a finite set of integer positions, under the gap
(Vietoris–Rips) filtration, is its sorted multiset of gaps.** Each position is
born at scale 0 and dies at the gap that merges it into its neighbour, so the
barcode *is* the gap multiset. Verified against explicit single-linkage on 300
random position sets: zero disagreements.

Two consequences, and they are the whole design:

1. Computing persistent homology on a lining buys nothing beyond the gaps.
2. The gap multiset is therefore **exactly the topological half** of a lining's
   structure, and everything else lives in the **order** of the gaps. That
   leftover is real: over the 70 distinct orderings of the gap multiset
   {3,4,3,4,3,4,3,4}, helical phase coherence runs from 0.572 to 0.902 — same
   topology, same n, same span, same occupancy, same gap entropy, and a 0.330
   spread in coherence.

So a lining is **topology ⊕ register**, and each half gets the null that holds
the other fixed:

| half | statistic | null | what it holds fixed |
|---|---|---|---|
| topology | fraction of gaps on one face (3–4 helical, 2 extended) | positions resampled uniformly within the same span | n, span |
| register | phase coherence over a band | the gap multiset permuted | n, span, the entire gap multiset |

The register null is not a choice among many. It is the quotient by the
topological invariant, which makes it the only null that varies the order and
nothing else. The parallel with `docs/SEQUENCE_PLAN.md` is exact: there,
letters ⊕ order with a composition-preserving shuffle; here, gaps ⊕ order with
a gap-preserving shuffle.

## Which half detects a face — measured, not assumed

**Register does not detect secondary-structure faces, and this is the single
most important thing measured before the study.** A helical face whose gaps are
drawn at random from {3, 4} has *no register*: permuting those gaps produces
another helical face, so the ordering carries nothing and the register test
correctly returns p = 1.00. Its periodicity lives entirely in the gap multiset,
where the topology test finds it at p = 0.0027.

Reading face-ness off a coherence statistic would therefore miss every
irregularly spaced face — which is most of them. The topology test is primary
for the *type*; register is a finer, secondary coordinate meaning "order beyond
what the gaps already force".

## Four artefacts, each measured before the design was fixed

1. **An analytic threshold on a scanned coherence is useless.** The Rayleigh
   95th percentile at n = 12 is 0.500 and period-independent. Its true
   false-positive rate is 0.019 at period 3.5 and 0.102 at period 2 — integer
   positions on a rational period cannot give uniform phases — and for the
   *maximum over a band*, which is what a scan actually computes, 0.500 admits
   **26%** of random linings over the helical band and **92%** over a full
   2–12 scan. Only resampling nulls on the same scanned statistic are used.
2. **A borrowed statistic that does not survive the lattice.** The axial
   (doubled-angle) resultant, standard in directional statistics for data with
   two modes 180° apart — which is what a β-strand's two faces are — is
   **identically 1** at period 2 for *any* integer position set (max deviation
   0.0e+00 over 2000 sets), and at period 4 for any set of uniform parity. A
   strand band containing period 2 would call every lining coherent. It is not
   used, and the reason is recorded rather than the choice.
3. **The scan needs a derived probe density.** The coherence peak half-width is
   `0.58 · ω²/span`, with a spread of 0.03 in that constant across a sixteenfold
   range of span. A fixed thirteen-probe scan of the helical band therefore
   loses up to **0.287** in resultant on a perfect face — 0.871 read back as
   0.584, a missed face. Probe density follows from the law.
4. **The permutation null has a floor.** With k distinct orderings no
   permutation p can fall below about 1/k: exactly 0.0667 at 15 orderings,
   0.0286 at 35, 0.0079 at 252. Significance at 0.05 is unreachable below 20,
   while a sampled z off such a null returned +2.60. Linings below
   `MIN_ORDERINGS = 20` are **unresolved**, which is a different claim from "no
   register", and the verdict is taken from the p-value, never the z.

## Multiplicity, measured

Four tests are asked of every lining — two face widths, two bands. Each alone
is calibrated or conservative (α = 0.020, 0.027, 0.057, 0.033 at nominal 0.05
over 300 random linings). **Taking the best of the four gives 0.123.** Holm
over the tests that ran restores 0.030 at 0.05 and 0.000 at 0.01. Every verdict
reads the Holm-corrected p; the raw p-values are reported alongside and are
never used for a claim.

## The spectrum and the categories

Two axes, both fixed here. The **run count** at a face-sized gap, which the
barcode gives exactly, and whether the gap multiset is a face, which the
topology test decides. Register is reported as a third coordinate and is
deliberately not a classification axis, for the reason measured above.

| class | runs | face-like | reading |
|---|---|---|---|
| **phase-locked** | 1–2 | yes | one secondary-structure element lines the pocket |
| **mosaic** | ≥3 | yes | several elements each contribute a face |
| **segmental** | 1–2 | no | a contiguous stretch or loop, no periodic face |
| **distributed** | ≥3 | no | scattered positions, interface-like |
| **unresolved** | — | — | neither half can be tested on this lining |

The graded scale, which is the yes / sort-of / no answer, is the Holm-corrected
p of the stronger half: **clear** (p ≤ 0.01), **weak** (p ≤ 0.05), **none**, or
**unresolved**. On 300 random linings the classes came out 449 distributed, 44
segmental, 7 mosaic, 0 phase-locked, and the graded verdict was `none` for 485
of 500 — so the categories are not a partition that assigns structure to
everything.

Boundaries and their reasons, all fixed before any data: a face gap of 4
because a helix repeats every 3.6 residues so the next same-face residue is 3
or 4 along; three runs because one element contributes one run and two runs is
still a hairpin; p ≤ 0.05 and ≤ 0.01 as the project's existing conventions;
20 orderings from the attainable-p floor above. **Two runs sitting in
`phase_locked` rather than `mosaic` is the boundary's weakest point**, and it is
recorded as such rather than defended.

## Power, and the design consequence

On the barrel's 42-residue face, detection against the fraction of the face a
lining captures. **Each half of the decomposition is reported separately and the
headline is the better of the two**, because the halves detect different pocket
shapes and one half alone is not the detection rate of the method:

| sampling | half | 100% | 80% | 60% | 45% | 33% (n=14) |
|---|---|---|---|---|---|---|
| spatial slab | topology | 1.00 | 1.00 | 1.00 | 1.00 | 0.10 |
| spatial slab | register | 1.00 | 1.00 | 1.00 | 1.00 | 0.73 |
| spatial slab | **either** | 1.00 | 1.00 | 1.00 | 1.00 | **0.77** |
| random subset | **either** | 1.00 | 1.00 | 1.00 | 0.97 | **0.40** |

against a false-positive rate of **0.05** on random linings of 14 positions, so
none of it is bought with false positives.

**This corrects an earlier version of this plan**, which reported the topology
half alone and so put the 14-position spatial rate at 0.10 rather than 0.77 — an
eightfold understatement that made the method look far weaker than it is. The
error was reading a sub-component where the composite was the right statistic;
it is recorded here rather than quietly fixed, and the probe now reports all
three rates so the same mistake cannot recur.

The corrected reading: detection is essentially complete down to about 19
positions and degrades gracefully rather than collapsing below that. Annotated
IP-site linings here are 10–20 residues, so **most of that range is usable**,
and lining size remains a reported covariate rather than a gate. The
element-extended lining stays the pre-registered primary unit, now as a power
*improvement* rather than a necessity. Pre-registered consequence: the
primary unit for Q7 is the **element-extended lining** — each lining residue
together with the contiguous run it belongs to in the local window — with the
bare lining reported as a secondary analysis whose expected result is null for
reasons of power rather than of biology. An MDE-style statement accompanies
every null: the fraction of the face captured, beside the detection rate
measured above at that fraction.

## The falsifiable prediction

At a fixed lining size of 14 residues on the same structure:

| shape | elements | helix-gap fraction | topology | register | **either** |
|---|---|---|---|---|---|
| ring, across all elements | 6.0 | 0.338 | 0.07 | 0.73 | **0.73** |
| column, along 2 elements | 2.0 | 0.615 | 1.00 | 1.00 | **1.00** |

Same structure, same lining size, a real difference — and a difference in
**degree**, not the presence-versus-absence an earlier topology-only reading
implied. The mechanism is visible in the halves, and it is the more interesting
result: a column is a long run on one element and carries a **topology** signal;
a ring is short runs separated by linkers in a repeating arrangement and carries
a **register** signal. The halves swap over with shape.

> **Prediction Q7-P1 (revised, and weaker than its first form).** At equal
> lining size the signature is **stronger** in elongated element-aligned pockets
> than in ring-shaped ones, and **the half that carries it differs between
> them** — topology for columns, register for rings.

Falsifiable in two directions: if detection is uniform across shapes, or if the
same half carries both, the framework is wrong. The second clause is the
sharper test and it was not available before the halves were reported
separately. Shape is assigned from existing structural descriptors before the
split and recorded per pocket in the inputs, never derived inside the analysis.

## Q7-H: architecture homology, and the controls it must beat

A second question the same invariants can be asked: are two pockets lined in
the *same way*? `cryptic_ip/detect/homology.py` answers it with a metric that
respects the decomposition — Wasserstein distance between H₀ persistence
diagrams for the topological half (the canonical TDA metric, applicable only
because the barcode *is* the gap multiset), and the correlation distance
between coherence spectra for the register half, which is exactly invariant to
sliding a lining along the chain. Both halves are metrics, so any non-negative
combination is one.

**This is architecture homology and it is a weaker claim than sequence or
structural homology.** It must never be used for the project's statistical
grouping: the independence `homology_group_strict` supplies comes from sequence
and structure comparison, and a metric on lining spacing supplies none of it.

Every retrieval runs beside two deliberately stupid controls — a distance that
knows only lining size, and one that knows only span. On two fold families
matched on lining size (18 positions each), residues per element, linker length
and with overlapping spans (`results/architecture/HOMOLOGY.md`):

| metric | P@1 | cluster purity | permutation p |
|---|---|---|---|
| size-only *(control)* | 0.50 | 0.53 | 0.49 |
| span-only *(control)* | 0.53 | 0.50 | 0.40 |
| topology | 0.97 | 0.67 | 0.0005 |
| **register** | **1.00** | **1.00** | **0.0005** |
| architecture (combined) | 0.97 | 0.67 | 0.0005 |

Two pre-registered readings. **The halves have opposite utility on the two
tasks** — topology detects a face within one lining and register cannot, while
register separates folds and topology does so less well — which is why both are
always reported and never collapsed. And **the combined metric is worse than
its better half**, because its stated scale of 4 residues per unit of register
distance leaves topology dominating about ten to one; that constant is a
documented convention and is deliberately **not** retuned against this result.

The graded claim fails its control and is recorded as a negative: distance
tracks how many elements a fold has at Spearman +0.95, but the span-only
control reaches +0.95 too, so **the graded ordering is explained by span and is
not a finding**. Family separation survives every control; the ordering does
not.

## Decision rules

* A signature is claimed for a stratum only with a Holm-corrected p ≤ 0.05
  **and** five or more homology groups on each side **and** the arms not
  confounded on lining size or on fraction-of-face captured.
* **A result that tracks lining size is not a finding.** Size is reported as a
  covariate in every comparison, and a difference that disappears when arms are
  matched on it is reported as size-explained.
* A null is reported with the detection rate measured above at the relevant
  fraction-of-face, so "we found nothing" stays separate from "we could not
  have".
* `unresolved` is never pooled with `none`. They are different claims and the
  counts are reported separately.
* No boundary moves after a result. Any change is a dated amendment.

## What this study cannot conclude

* Nothing causal. "Leads to pocket formation" is not testable from a
  cross-section; the strongest claim is association between a sequence-space
  invariant and a structural property, with homology groups resampled.
* Nothing from the barrel about real proteins. It is one idealised fold with
  one exact periodicity, so its power and shape numbers are an **upper bound**
  on what an irregular real protein gives. The identities, degeneracies and the
  peak-width law are general; the detection rates are not.
* Nothing about prediction. A descriptor that separates classes here is not a
  predictor until tested on groups held out from its fitting, which this study
  does not do.
* Nothing about Kolmogorov complexity, which remains uncomputable, and nothing
  here is named after it.
