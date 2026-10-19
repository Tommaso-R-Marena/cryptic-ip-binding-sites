# Lining architecture: what the instruments can and cannot see

**Run.** `scripts/architecture_probe.py run`, seed 20261012, on the session container. Every figure below is computed by that command and written into this file in the same pass, so nothing is transcribed.

Synthetic and analytic only. No annotated sequence is read; `docs/ARCHITECTURE_PLAN.md` pre-registers the study that would use this.

## 1. Two identities that collapse

The H0 persistence barcode of a lining equals its sorted gap multiset: 0 disagreements against explicit single-linkage in 300 random position sets. So persistent homology on a lining buys nothing beyond the gaps, and the gaps are therefore *exactly* the topological half of a lining's structure.

## 2. One borrowed statistic that does not survive the lattice

The axial (doubled-angle) resultant, standard in directional statistics for data with two modes 180 degrees apart, is **identically 1** at period 2 for any integer position set: max deviation 0.0e+00 over 2000 arbitrary sets, and 0.0e+00 at period 4 for uniform parity. A strand band containing period 2 would call every lining axially coherent, so it is not used.

## 3. An analytic threshold on a scanned coherence is useless

At n = 12 the Rayleigh 95th percentile is 0.500, independent of period. Its true false-positive rate on random linings:

| statistic | probes | mean R | true alpha at the analytic cutoff |
|---|---|---|---|
| single period 3.6 | 1 | 0.234 | **0.021** |
| helix band | 13 | 0.433 | **0.233** |
| strand band | 6 | 0.392 | **0.192** |
| full scan 2-12 | 101 | 0.537 | **0.714** |

A full 2-12 scan at the nominal cutoff admits most random linings. Only a resampling null on the same scan is honest.

## 4. The peak-width law, and what a fixed scan loses

Half-width = 0.58 x w^2 / span, with a spread of only 0.038 across a sixteen-fold range of span:

| span | measured half-width | ratio to w^2/span |
|---|---|---|
| 20 | 0.3606 | 0.556 |
| 40 | 0.1818 | 0.561 |
| 80 | 0.0940 | 0.581 |
| 160 | 0.0481 | 0.594 |
| 320 | 0.0240 | 0.593 |

A fixed thirteen-probe scan therefore loses up to **0.287** in resultant on a perfect face, which is a missed face rather than rounding. Probe density is derived from the law.

## 5. The permutation null has a floor

| gap multiset | orderings | smallest attainable p |
|---|---|---|
| [3, 1, 3, 3, 1, 3] | 15 | 0.0667 |
| [1, 1, 1, 8, 1, 1, 1] | 7 | 0.2857 |
| [4, 3, 4, 3, 4, 3, 4] | 35 | 0.0286 |
| [3, 4, 3, 4, 3, 4, 3, 4] | 70 | 0.0286 |
| [4, 3, 4, 3, 4, 3, 4, 3, 4, 3] | 252 | 0.0079 |

Significance at 0.05 is unreachable below about 20 orderings, so a lining under 20 is **unresolved** rather than "no register". A z-score off such a null happily returns +2 and means nothing.

## 6. Calibration

On 300 random linings, each test alone is calibrated or conservative:

| test | fraction runnable | alpha at 0.05 |
|---|---|---|
| helix_topology | 1.00 | 0.020 |
| strand_topology | 1.00 | 0.027 |
| helix_register | 1.00 | 0.057 |
| strand_register | 1.00 | 0.033 |

But the family of four is not: taking the best of them gives **0.123** at a nominal 0.05. Holm over the tests that ran restores it to **0.030** at 0.05 and 0.000 at 0.01, so every verdict reads the corrected p.

## 7. Recovery of known architectures

| lining | n | runs | occ | dom w | topology p | register p | signature | element | class |
|---|---|---|---|---|---|---|---|---|---|
| barrel, no linkers | 42 | 1 | 0.80 | 3.73 | 0.0027 | 0.0027 | clear | helix_face | phase_locked |
| barrel + 10-residue linkers | 42 | 6 | 0.74 | 3.58 | 0.0027 | 0.0027 | clear | helix_face | mosaic |
| three helices + linkers | 21 | 3 | 0.72 | 3.57 | 0.0040 | 0.0027 | clear | helix_face | mosaic |
| helix face, exact 3,4 spacing | 12 | 1 | 0.67 | 3.51 | 0.0027 | 0.0065 | clear | helix_face | phase_locked |
| helix face, jittered 3/4 | 13 | 1 | 0.67 | 3.61 | 0.0027 | 1.0000 | clear | helix_face | phase_locked |
| strand face, every 2nd | 14 | 1 | 0.80 | 2.00 | 0.0013 | -- | clear | strand_face | phase_locked |
| contiguous stretch of 14 | 14 | 1 | 1.00 | 9.79 | -- | -- | unresolved | none | unresolved |
| hairpin, two runs of 7 | 14 | 2 | 0.73 | 10.26 | 1.0000 | -- | none | none | segmental |
| random 14 in 200 | 14 | 9 | 0.50 | 11.05 | 0.5889 | 1.0000 | none | none | distributed |
| random 20 in 400 | 20 | 17 | 0.51 | 8.80 | 0.6216 | 0.4344 | none | none | distributed |

The row that matters most is the jittered helical face. Its gaps are drawn at random from {3, 4}, so permuting them gives another helical face and **register finds nothing** - correctly. Its periodicity lives entirely in the gap multiset, which the topology test sees. Reading face-ness off register would miss every irregularly spaced face, which is most of them.

## 8. Power: how much of a lining must be seen

The barrel has 42 face residues. Subsampling them, with each half of the decomposition reported separately and the headline being the better of the two:

| sampling | fraction of face | n | topology | register | **either** |
|---|---|---|---|---|---|
| random | 1.00 | 42 | 1.00 | 1.00 | **1.00** |
| random | 0.80 | 34 | 0.97 | 1.00 | **1.00** |
| random | 0.60 | 25 | 0.87 | 1.00 | **1.00** |
| random | 0.45 | 19 | 0.37 | 0.97 | **0.97** |
| random | 0.33 | 14 | 0.03 | 0.37 | **0.40** |
| spatial | 1.00 | 42 | 1.00 | 1.00 | **1.00** |
| spatial | 0.80 | 34 | 1.00 | 1.00 | **1.00** |
| spatial | 0.60 | 25 | 1.00 | 1.00 | **1.00** |
| spatial | 0.45 | 19 | 1.00 | 1.00 | **1.00** |
| spatial | 0.33 | 14 | 0.10 | 0.73 | **0.77** |

Against a false-positive rate of **0.05** on random linings of 14 positions, so none of this is bought with false positives.

**The two halves detect different things, and reading either one alone understates the method.** An earlier version of this probe reported the topology half only - on the reasoning that topology is what detects a face, which is true of a single element and wrong as a detection rate - and so understated a 14-position slab eightfold.

And at one fixed lining size, pocket **shape** still matters, though less than the topology half alone suggested:

| shape | lining size | elements spanned | helix-gap fraction | topology | register | **either** |
|---|---|---|---|---|---|---|
| ring across all elements | 14 | 6.0 | 0.338 | 0.07 | 0.73 | **0.73** |
| column along 2 elements | 14 | 2.0 | 0.615 | 1.00 | 1.00 | **1.00** |

Same structure, same lining size, and a real difference - 1.00 against 0.73 - but a difference in degree rather than the presence-versus-absence the topology-only measurement implied. The mechanism is visible in the two halves: a column lining is a long run on one element, which is a **topology** signal (1.00 against 0.07); a ring lining is short runs separated by linkers in a repeating arrangement, which is a **register** signal (0.73 against 0.07 for its own topology).

So the falsifiable prediction stands in a weakened and more precise form: the signature should be **stronger** in elongated element-aligned pockets than in ring-shaped ones at equal lining size, and the half that carries it should differ between them. It is no longer a claim that ring pockets carry no signature.


## 9. Does any of this transfer? The sharpest test, and it is a hard negative

Every number above comes from barrels built out of exact 3.6-residue helices, a periodicity no real protein has. `build_backbone` now perturbs each residue's phi and psi and can insert bulges that shift the register of everything after them. Both leave bond lengths and angles **exactly** ideal - only torsions move - so the structures stay valid and only the realism changes.

| condition | face residues | topology | register | **either** |
|---|---|---|---|---|
| ideal | 42 | 0.17 | 0.70 | **0.73** |
| jitter 4 deg | 45 | 0.07 | 0.50 | **0.53** |
| jitter 8 deg (realistic) | 46 | 0.07 | 0.27 | **0.33** |
| jitter 12 deg | 45 | 0.10 | 0.20 | **0.23** |
| jitter 20 deg | 46 | 0.27 | 0.00 | **0.27** |
| one bulge per helix | 51 | 0.03 | 0.30 | **0.30** |
| bulges + jitter 8 deg | 47 | 0.17 | 0.03 | **0.20** |

**Detection falls from 0.73 on the idealised barrel to 0.33 at a realistic 8 degrees of phi/psi scatter, and to 0.20 with bulges + jitter 8 deg.** The documentation has said throughout that the barrel numbers are an upper bound on what an irregular protein would give. This measures the gap rather than asserting it, and the gap is large: roughly a halving at realistic jitter and a three- to fourfold fall once bulges interrupt the face.

That is the single most important limit on this method and it is not a fixable one - a bulge genuinely destroys the periodicity the statistic looks for, because the periodicity is no longer there to find. It bounds what study Q7 can expect from real structures, and it is why the real-lining analysis in `results/architecture/real/REAL.md` must be read beside it.


## 10. Can a *local* statistic recover what jitter takes away?

Pre-registered in `docs/ARCHITECTURE_PLAN_AMENDMENT_3.md`, and as each answer came back in amendments 4 and 5 - every one dated before the statistic it defines had been measured on any fixture, decoy or real lining.

The reasoning is the derivation that justifies the pooled resultant in the first place: it is Neyman-Pearson optimal for a face of one period and **one unknown phase**, and a jittered helix does not have one phase - the rise per residue drifts and a bulge resets it outright. The pooled statistic therefore asks every residue of the lining to agree with every other, which over a hundred-residue span is a stronger demand than the structure can meet. A **local** statistic asks only for piecewise phase lock: the maximum resultant over windows of [6, 8, 10] consecutive lining positions, with the identical maximum taken on every draw of the same gap-permutation null, so nothing is left to correct afterwards.

`pruned` and `best` differ from it only in **which window sizes they keep**, and they choose from the null rather than from the observation - `pruned` by the saturation threshold, `best` by the lowest null mean. That selection is free, and it is the one subtle point here: the resultant's modulus is translation invariant, so a window's null distribution depends only on the **gap multiset**, which the register null holds fixed. A rule stated in terms of it picks the same window set for the observation and for every permutation, and a permutation test whose statistic is the same function on every draw is exactly valid. The estimate comes from a 200-draw pilot with its own seed, so it cannot be a statistic of the null sample it will be judged against.

Paired by fixture - every statistic sees the same lining at the same seed - so each comparison is an exact McNemar on the discordant pairs. All three variants are **register** statistics, so the column they are judged against is the pooled *register* half alone and not the `either` rate quoted elsewhere in this report: comparing a register statistic with a best-of-both rate would flatter or damn it for the wrong reason.

| arm | pooled register | local | pruned | best | McNemar (local) | McNemar (pruned) | McNemar (best) |
|---|---|---|---|---|---|---|---|
| ideal | 0.80 | 0.00 | 0.00 | 0.10 | 1.2e-07 | 1.2e-07 | 9.5e-07 |
| jitter 8 deg (realistic) | 0.30 | 0.00 | 0.00 | 0.07 | 0.0039 | 0.0039 | 0.0391 |
| bulges + jitter 8 deg | 0.10 | 0.00 | 0.00 | 0.10 | 0.2500 | 0.2500 | 1.0000 |

On 30 decoys matched in size and span, with no planted face: pooled 0.00, local 0.07, pruned 0.07, best 0.03, at 0.05.

Read that beside the detection columns rather than on its own: a statistic that detects nothing on a planted face and still calls the occasional decoy is not trading sensitivity for specificity, it is weaker on both counts. At this decoy count the rate is two or three linings either way and carries no interval worth quoting - it is here because the decision rule named it in advance, and it is not what decides the question.

**Window-size null means, which is the diagnostic that explains all of it.**

| arm | size 6 | size 8 | size 10 | size the rank chose |
|---|---|---|---|---|
| ideal | 0.884 | 0.796 | 0.717 | 10: 1.00 |
| jitter 8 deg (realistic) | 0.883 | 0.791 | 0.709 | 10: 1.00 |
| bulges + jitter 8 deg | 0.849 | 0.749 | 0.667 | 10: 1.00 |

**local**: not adopted, **pruned**: not adopted, **best**: not adopted, by the rule fixed in advance: a variant is adopted only if it detects strictly more than the pooled statistic at the realistic arm and its decoy false-positive rate is at or below 0.05; otherwise it is recorded as a measured negative.

The null mean **rises as the window shrinks** - 0.88 at size 6 against 0.72 at size 10 - because a band scan over a handful of positions fits almost any spacing, which is the same failure as the per-element statistic in section 8. Since the statistic maximises over the sizes it keeps, the shortest window it keeps sets its null, and that is why the three variants separate at all.

**The 0.9 threshold never fired**, which is itself a result: every window size was kept on every lining, so `pruned` computed the same statistic as `local` and the two columns agree by construction. The threshold was fixed in amendment 3 against a size-6 null mean measured at 0.90 on a different lining; on the barrel slabs it sits just *below*. That is a lesson about fixing a cutoff against one observed number - pre-registration kept the test honest and also left it with nothing to do - and moving the threshold afterwards would be exactly the tuning the amendment forbade. Amendment 5 answered it with a **rank** instead, which cannot fail to fire and rests on the same validity argument: any function of the gap multiset alone is free.

The rank selected window size 10 on every lining of every arm, so `best` is a long-window statistic described by its selection rule rather than by its number.

The predictions, recorded before the run: amendment 4 predicted that dropping the saturated windows would beat the unpruned statistic and still lose to pooled; amendment 5 predicted that the rank would select the longest window on every lining, beat both of the other variants and still lose to pooled, because a 10-position window of a 14-position lining is the pooled statistic minus four positions plus a maximum over five offsets.

What this establishes is worth more than the negative. **Locality is not what limits transfer.** At the realistic arm the pooled register half detects 0.30 against local 0.00, pruned 0.00, best 0.07. The ordering within the windowed variants is the predicted one and confirms the mechanism - the longer the window a variant keeps, the better it does, because the short windows carry a saturated null - but the best of them is the pooled statistic minus a few positions plus a maximum over offsets, and what that maximum buys is less than what the lost positions cost. So the pooled statistic is not losing power on jittered helices because it demands global phase lock. The limit is the band scan's flexibility on short runs of positions, a property of the scan rather than of the geometry, and that bounds every windowed variant of this family rather than one implementation of it.

## What this cannot be used for

* It is not evidence about real binding sites. Every lining here is synthetic or derived from an idealised barrel, and no annotated sequence has been read.
* The barrel is one fold with one periodicity. The peak-width law and the degeneracies are general; the power and shape numbers are properties of this fixture and are an upper bound on what a real, irregular protein would give.
* The class boundaries - three runs, |z| and p cutoffs, a face gap of 4 - are fixed from geometry and from the statistics above, not from any outcome, and they have never been adjusted against a result.
