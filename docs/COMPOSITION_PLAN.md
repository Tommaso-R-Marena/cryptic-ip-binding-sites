# Study Q10: the composition half of the decomposition, and the positive control this project has never had

**Dated 2026-10-25**, written and committed **before any composition statistic has
been computed on any real lining in this project.** Study Q9's artefacts are in
hand and were read; nothing below is a function of them except where section 7
says so explicitly and says what was read.

## 1. Why this study exists

A window of sequence factorises as **composition** (which letters) **⊕ order**
(in what arrangement). Study Q9 took the quotient by composition and asked
whether order carries anything. The answer was no for pocket formation and for
architectural type, and unresolved for inositol-phosphate binding.

**Q9 could not have answered the other half, and did not try.** Its null
permutes a window's own letters, so composition is held exactly fixed and is
invisible to every test it ran. That was deliberate — these sites are known to
be basic and pocket linings are known to be hydrophobic, so a raw compositional
difference is already explained and would have been a finding about nothing.

This study takes the other quotient. It is **not** interesting as biology: that
an inositol-phosphate contact shell is arginine-rich is textbook, and nothing
here claims otherwise. It is interesting for one reason:

> **Q9's three nulls have no real-data positive control.** The only thing
> demonstrating that its machinery can detect a signal is a synthetic planted
> one. Nothing shows that *these* linings, in *these* structures, through *this*
> pipeline, can produce a detection at all.

So the point of Q10 is to put a known-true signal through as much of the same
pipeline as possible and see whether it survives. Two outcomes, both worth
having:

* **The positive control fires.** The instrument is demonstrated on real data,
  and Q9's nulls become statements about order rather than statements about the
  pipeline. The other arms then ask what composition carries beyond basicity.
* **The positive control does not fire.** Then the pipeline cannot detect a
  signal that is certainly present, Q9's nulls are uninterpretable, and that is
  a far more important result than anything Q9 reported. **This plan commits to
  reporting it that way**, and section 6 makes it a code fence rather than a
  promise.

## 2. The window is the lining, and that is legitimate here and nowhere else

Q9 scored the **span** a lining covers — the contiguous stretch from its first
residue to its last, capped at forty. That window contains every residue
*between* the lining's residues as well, so a twelve-residue shell spread over a
long span is scored as a forty-residue window holding a handful of it. Signal is
diluted by whatever the occupancy happens to be.

Q10 scores the **lining** — the gap-closed concatenation of the lining's own
residues, and nothing else. This is the window `cryptic_ip/detect/signature.py`
already calls `lining`, and that module already fences which statistics may use
it: a positional statistic may not, because closing the gaps destroys what a
period means, and `POSITIONAL_STATISTICS` exists precisely to stop one being
computed there by accident.

**Every statistic in this study is a function of the multiset of letters alone**
— invariant to their order — so the gap closing cannot affect any of them. The
existing fence is the reason this is sound, and section 5 adds the matching
fence from the other side.

Consequences, fixed here before they are measured: no span cap and no occupancy
threshold apply, so every lining with enough residues is usable, where Q9's
window admitted 2001 of 8407 pockets and 9 of 129 shells.

## 3. The null is a same-size subset of the same chain

For each lining, draw `SUBSET_DRAWS = 2000` random residue subsets of the same
size from the same chain, compute the statistic on each, and report the
observed value's z and two-sided permutation p against that distribution.

This holds the chain, the organism, the fold and the lining's size fixed and
randomises only **which** residues. It is the compositional counterpart of Q9's
order null, and it answers "is this lining compositionally unusual among
same-size residue sets of its own chain".

It is a weaker control than another pocket of the same protein, and the
architecture line has already said in print that a random position set is not a
protein. That criticism was about *shape*, where it is decisive. For
composition it is not the same objection — the composition of a random subset is
the chain's background composition, which is the right reference for "unusual
for this chain". The stronger control is used as well, in arm Q10c.

## 4. Statistics

All are order-invariant functions of a window's letters.

**Primary (three):**

* `basic_fraction` — the fraction of R and K. **The positive control's
  statistic.** Direction pre-stated: an annotated inositol-phosphate contact
  shell is **more** basic than same-size subsets of its own chain, so a positive
  z.
* `composition_divergence` — the Jensen–Shannon divergence, in bits, between the
  lining's amino-acid distribution and its chain's background distribution. The
  omnibus: it fires on any compositional difference, in any direction, without
  naming one in advance.
* `nonbasic_divergence` — the same divergence computed over the eighteen
  residues that are not R or K, renormalised. This is the "what else, besides
  being basic" statistic, and it is the only primary one whose outcome is not
  already known.

**Secondary, reported in full and never promoted:** `acidic_fraction`,
`aromatic_fraction`, `mean_hydropathy`, `glycine_fraction`, `proline_fraction`,
`net_charge_per_residue`.

## 5. Two fences, one on each side of the decomposition

* **`composition_z` must refuse an order-invariant statistic.** A statistic that
  does not read order has zero variance under an order-permuting null, so its z
  is 0/0 — and the current code returns a silent `0.0` for it, which reads as
  "no ordering signal" when the truth is "the question is undefined". Study Q9
  never hit this because none of its statistics was order-invariant. Q10 adds
  nine that are, so the guard is added with them.
* **No statistic may be registered in both families.** A name appearing in both
  `sequence.STATISTICS` and the new composition registry is a contradiction, and
  a test fails on it, in the spirit of the `POSITIONAL_STATISTICS` fence.

## 6. Arms, and the gate

* **Q10-pc — the positive control.** Annotated inositol-phosphate contact
  shells, `basic_fraction`, subset-z against zero, grouped by
  `homology_group_strict`. **It is a control and not a test**: it is in no Holm
  family, and it gates everything else.
* **Q10a — formation.** A protein's better-ranked pockets against its
  worse-ranked ones, matched on lining size to within two residues, on subset-z.
  Same rank-split rule as `docs/POCKET_SEQUENCE_PLAN_AMENDMENT_1.md` section 3,
  because the same rank/size collinearity applies.
* **Q10b — type.** Pockets the architecture module calls `helix_face` against
  those it calls `none`, matched on lining size. Reported; its own Holm family.
* **Q10c — function.** Annotated shells against **the other pockets of the same
  deposited structure**, matched on lining size — the within-protein design of
  `docs/WITHIN_PROTEIN_PLAN.md`, which is the hard version and the one that
  would mean something.

**The gate, which is the whole point of the study.** If Q10-pc does not clear
every one of the plan's floors in the pre-stated direction, then **no other arm
in this study or in study Q9 is read as a null result**, and every row of this
study is reported as *instrument not demonstrated* rather than as evidence of
absence. This is enforced in code: `is_finding` is forced false on every arm and
the report says so in its own first paragraph. A plan that only promised this
would not be worth writing.

## 7. What was read from study Q9 before this was written

Q9's artefacts were read in full, so this plan is not blind to them. What
informed the design, named so it cannot be mistaken for a prediction:

* that its primary tests returned +0.014, +0.032, −0.064, +0.061 and +0.191
  against minimum detectable differences of 0.300, 0.356, 0.428, 0.596 and
  0.934;
* that its window admitted 2001 of 8407 pockets and 9 of 129 contact shells,
  114 of those 129 refused for span alone.

The second is why section 2 changes the window. The first is why section 1
argues the pipeline needs a positive control: **a design that returns five
nulls in a row has either found a real absence or has no sensitivity, and
nothing in Q9 distinguishes those.** No composition statistic of any real lining
has been computed, by me or by any committed code in this repository, before
this file.

## 8. Decision rules, inherited unchanged

* A finding needs a Holm-corrected p below 0.05 across the primary family, a
  bootstrap interval excluding zero, at least **five** resampling groups, an
  observed difference exceeding the minimum detectable difference, and matched
  window sizes.
* Intervals resample whole homology groups (`homology_group_strict`) or MMseqs2
  30% clusters; a stratum with fewer than five groups is reported as no
  evidence.
* Holm runs across the primary family: three statistics on Q10a and Q10c, six
  tests. Q10b has its own family. Q10-pc is in neither.
* The direction is pre-stated only for the positive control. The other arms are
  two-sided.

## 9. What this cannot conclude, whatever it shows

* **A compositional difference is not a discovery.** That these sites are basic
  is known, and the positive control is in this study to test the instrument,
  not the biology. No row of it will be presented as a finding about inositol
  phosphate.
* A lining's composition is read from a predicted model for the candidate arms
  and from one deposited conformation for the function arm, so neither is a
  statement about solution behaviour.
* Every arm is within-protein or within-chain. **Nothing here is a sequence-only
  predictor of anything**, and this study does not attempt one.
* Composition and order do not exhaust a sequence. Conservation, coevolution,
  linear motifs and learned representations are outside both this study and
  Q9, and a null in both says nothing about them.
