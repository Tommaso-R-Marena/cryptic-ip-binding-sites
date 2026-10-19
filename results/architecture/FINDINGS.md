# The lining-architecture line: what it establishes

This page is **generated** by `scripts/architecture_findings.py` from the committed artefacts listed below. It recomputes nothing and no number on it is retyped; a figure absent from an artefact is absent from the page.

| section | artefact | command |
|---|---|---|
| probes | `results/architecture/architecture.json` | `scripts/architecture_probe.py run` |
| homology | `results/architecture/homology.json` | `scripts/homology_probe.py run` |
| real | `results/architecture/real/real.json` | `scripts/real_architecture.py run` |
| paired | `results/architecture/real/basic_face.json` | `scripts/basic_face.py run` |

The pre-registered plans are `docs/ARCHITECTURE_PLAN.md` and its five dated amendments; the module-level account is the `architecture.py` and `homology.py` sections of `docs/DETECTOR.md`; figures are in `results/architecture/figures/`.

## 1. What is settled, and not by a p-value

**The H₀ barcode of a lining *is* its sorted gap multiset.** For a finite set of integer positions under the gap filtration, the complete H₀ persistence invariant carries exactly the information in the sorted gaps — 0 disagreements in 300 comparisons against explicit single-linkage. The consequence is the useful part: persistent homology buys nothing here beyond the gaps, so the honest decomposition of a lining is

> **topology** (which gaps occur) **⊕ register** (in which order they occur),

and each half gets the null that holds the other fixed. The topological null resamples positions uniformly in the same span, holding size and span; the register null permutes the gap multiset, holding topology *entirely*. The second is the quotient by the invariant the first cannot control.

**The resultant is the Neyman–Pearson-optimal statistic** for a von Mises face of unknown phase: marginalising the phase gives a likelihood ratio in `I₀(κ·n·R̄)`, monotone in the resultant, so no function of the positions can beat it at that alternative. Verified against direct numerical marginalisation to nine decimal places.

**The axial (doubled-angle) resultant is unusable here**, which is a proof and not a tuning note: it is *identically* 1 at period 2 for any integer set — max deviation 0.0e+00 — so it cannot distinguish anything. It was dropped.

## 2. Four ways the instrument lied, each measured and fixed

| artefact | what was measured | fix |
|---|---|---|
| Analytic Rayleigh threshold | 0.500 admits a large fraction of *random* linings at n = 12, span 60 | permutation nulls, never the analytic threshold |
| Band-scan resolution | a fixed 13-probe grid loses up to 0.287 in resultant at wide span | probe density derived from the measured peak half-width law, 0.58·ω²/span, flat to 0.04 across a 16× span range |
| Permutation floor | a lining with few distinct gap orderings cannot reach a small p at all | `MIN_ORDERINGS = 20`; below it the register half is **unresolved**, which is not a p of 1 |
| Family multiplicity | best-of-four over two face widths × two bands gives a false-positive rate of 0.123 at a nominal 0.05 | Holm across the family restores 0.030 |

Calibration is on 300 random linings: after Holm, 291 none, 9 weak, and none *clear*.

## 3. What it detects, and the ceiling on all of it

On planted faces in helix barrels with ideal covalent geometry, at lining size 14:

| planted shape | topology | register | either |
|---|---|---|---|
| ring across all elements | 7% | 73% | 73% |
| column along 2 elements | 100% | 100% | 100% |
| *decoy (no planted face)* | 2% | 2% | 5% |

The decoy arm is the row that makes the others readable, and reporting all three columns is a correction: an earlier version of this table read `topology_p` alone and so understated the shapes register detects — by more than a factor of seven on the ring.

**How much of a face is enough.** Diluting a planted face into a larger lining: **random** dilution holds to 45% of the face (19 of 42 positions, 97% detected) and falls to 40% at 33%; **spatial** dilution holds to 45% of the face (19 of 42 positions, 100% detected) and falls to 77% at 33%. Topology is the half that fails first: at the edge of the random curve it is already down to 37% while register is at 97% — the same split the shape table shows.

**The ceiling.** Those are barrels of exact 3.6-residue helices. Perturbing each residue's φ and ψ leaves bond lengths ideal to 1e-15, so only realism changes, and detection falls from **73%** to **33%** at the 8° scatter real structures show, and to **20%** under *bulges + jitter 8 deg*. This is the single largest limit on the method and it is not a tuning problem: a bulge removes the periodicity rather than hiding it. **Every fixture number above should be read as an upper bound about twice the realistic value.**

## 3a. The ceiling was attacked, and it held

Pre-registered in `docs/ARCHITECTURE_PLAN_AMENDMENT_3.md`, and after its answer came back in `docs/ARCHITECTURE_PLAN_AMENDMENT_4.md`, each dated before the statistic it describes had been computed on anything.

The pooled resultant is NP-optimal for a face of one period and **one** unknown phase. A jittered helix does not have one phase: the rise per residue drifts and a bulge resets it. So the obvious repair is a statistic that asks only for *piecewise* phase lock — the maximum resultant over windows of [6, 8, 10] consecutive lining positions, with the same maximum taken on every draw of the same gap-permutation null. Two further variants keep only some of those window sizes — the unsaturated ones, and the single one with the lowest null mean — and **selecting them is free**: the resultant's modulus is translation invariant, so a window's null depends only on the gap multiset, which the register null holds fixed, and a rule stated in terms of the null therefore picks the same windows for the observation and for every permutation, leaving the permutation p-value exact.

| arm | pooled register | windowed | pruned | best window | McNemar (best) |
|---|---|---|---|---|---|
| ideal | 0.80 | 0.00 | 0.00 | 0.10 | 9.5e-07 |
| jitter 8 deg (realistic) | 0.30 | 0.00 | 0.00 | 0.07 | 3.9e-02 |
| bulges + jitter 8 deg | 0.10 | 0.00 | 0.00 | 0.10 | 1.0e+00 |

**No variant was adopted**, by the rule fixed before any of them was computed. The diagnostic says why in one number: the permutation null's own *mean* at window size 6 is 0.883 out of a possible 1, so a random reordering of the same gaps is already nearly as coherent as a lining can be over that many positions, and there is almost nothing left above it for an observation to be unusual in. Since the statistic maximises over the sizes it keeps, the shortest window it keeps sets the null for the whole test. It is the same failure as the per-element statistic, for the same reason: a band scan over a handful of positions fits almost any spacing. False positives are not where this was decided — on decoys the pooled statistic calls 0% and the best-window one 3%.

**What the failure establishes is worth more than the attempt.** Locality is not what limits transfer. The pooled statistic is not losing power on jittered helices because it demands global phase lock, because a statistic that demands only local phase lock does no better, and nor does the long-window variant the rank selects — which it selects on every lining of every arm, making it the pooled statistic minus four positions plus a maximum over five offsets, a trade the measurement says is a loss. The limit is the band scan's flexibility on short runs of positions — a property of the scan, not of the geometry — and that bounds every windowed variant of this family rather than one implementation of it.

## 4. Two pockets lined the same way

A metric that respects the decomposition instead of ignoring it: Wasserstein-1 between H₀ persistence diagrams for topology — the canonical TDA metric, usable here *only* because the barcode is the gap multiset — and correlation distance between coherence spectra for register, which is exactly invariant to sliding a lining along the chain. Identity, symmetry and the triangle inequality are tested for both halves.

Retrieving fold family (helix barrel, face gaps 3–4, against β sheet, gap 2) with lining size 18, residues per element 18 and linker 10 held equal and spans overlapping:

| metric | P@1 | cluster purity | permutation p | |
|---|---|---|---|---|
| `size_only` | 0.50 | 0.53 | 0.4918 | control |
| `span_only` | 0.53 | 0.50 | 0.4033 | control |
| `topology` | 0.97 | 0.67 | 0.0005 |  |
| `register` | 1.00 | 1.00 | 0.0005 |  |
| `architecture` | 0.97 | 0.67 | 0.0005 |  |

**`register` separates the two families** where the size-only and span-only controls sit at chance. Three things follow, two of them limits, and both limits are left standing rather than tuned away:

* **The combined metric is worse than its better half** — purity 0.67 against 1.00 — because its stated 4-residues-per-unit scale leaves topology dominating about ten to one. The constant has deliberately not been retuned against the result that exposed it.
* **The halves have complementary blind spots, not a ranking.** Within one lining topology finds tight, uniform faces that leave register with nothing to test — the panel's strand face, every 2nd residue, has exactly one gap ordering, so register is unresolved while topology gives p = 0.001 — and register finds the rings and diluted faces whose gap multiset is unremarkable. Across linings register separates folds and topology does so less well. That is the reason to keep both halves and never collapse them into one number.
* **The graded relationship is not a finding.** Distance tracks element count at Spearman +0.95 in the barrel family — but the span-only control reaches +0.95 on the same pairs, so span explains it.

And the standing caveat: this is *architecture* homology. It does not show two proteins are evolutionarily related, it does not replace Foldseek or MMseqs2, and it must never be used for the project's statistical grouping.

## 5. On real linings

`scripts/build_lining_inputs.py` turns the project's committed tables into study inputs. Q7 runs on 75 real candidate-pocket linings (median size 15, range 7–75), against 375 decoys matched to each lining's size *and* span.

| | topology | register | either |
|---|---|---|---|
| real linings | 31% | 11% | 33% |
| matched decoys | 2% | 3% | 4% |

`mosaic` is the class that moves: 23 of 75 real linings against 6 of 375 decoys. 17 are called helical faces and 2 strand faces.

The real-lining rate (33%) is the same figure the jittered fixtures give. That is a consistency observation, not a tested prediction.

The one pre-registered split — pockets another polyanion explains against the rest — is null on every feature after Holm (0 findings of 4 features).

**The control this comparison does not have, stated where it is read.** A random position set is not a protein, so *any* real pocket lining would beat that null. The table shows the instrument registers real linings; it shows nothing specific to inositol-phosphate sites. The control that would is other pockets from the same proteins, and no committed table carries residue lists for them — `results/ml_training/validation_pocket_features.csv` has 12,190 pockets across 136 structures but descriptors only.

## 6. The design whose control is not tautological

Section 5's weakness is fixable in one direction without new data: compare each pocket's **basic residues** with random subsets **of that same pocket's lining, of the same size**. Protein, fold, lining and subset size are then held fixed by construction, and a difference can only be which residues were chosen. Pre-registered in `docs/ARCHITECTURE_PLAN_AMENDMENT_2.md` before any basic-residue subset was measured.

The per-pocket statistic is the subset's mid-p quantile among 50 matched draws, uniform on [0, 1] under the null. **Below 0.5 is more structured than chance.** 99 subset measurements from 75 pockets.

| subset | half | pockets | mean quantile | difference | 95% interval | Holm | finding |
|---|---|---|---|---|---|---|---|
| `basic_positions` | topology | 54 | 0.435 | -0.065 | [-0.151, +0.020] | 0.6257 | no |
| `basic_positions` | register | 54 | 0.536 | +0.036 | [-0.049, +0.116] | 1.0000 | no |
| `conserved_positions` | topology | 45 | 0.488 | -0.012 | [-0.108, +0.077] | 1.0000 | no |
| `conserved_positions` | register | 45 | 0.527 | +0.027 | [-0.059, +0.111] | 1.0000 | no |

**0 of 4** pre-registered tests clears both the Holm threshold and its own power floor.

Two corrections are baked into that table, both of which changed the answer. An earlier version scored an **unresolved** register half as a p of 1 — collapsing the exact distinction `architecture.py` insists on, *cannot be tested* against *no structure* — which put half of all register quantiles at exactly 1.000 and none at 0.000. And permutation p-values are discrete, so counting ties as "at or below" pushes every quantile upward; the quantile here splits ties.

Section 3's ceiling stands over all of this: at realistic φ/ψ scatter the method detects about a third of what it detects on ideal geometry, so a null here is weak evidence of absence. The honest statement is that this design could not see an effect of the size it could see on ideal geometry.

## 7. What this line is, and is not, for

These descriptors are **deliberately not** features of the detector's score, and `docs/DETECTOR.md` records the decision with its reason: a feature earns a place in the scorer by carrying signal on the inputs the scorer sees, which are real structures with real backbone scatter, and two thirds of the measured sensitivity here is an artefact of ideal fixtures. Section 6 adds that there is no real-data effect to carry in either. Section 3a is the attempt to lift that ceiling and the measurement that says this family of statistics cannot.

So the line is an instrument with its proofs and its calibration attached, reported beside the score and never folded into it. Nothing on this page is evidence that a pocket binds an inositol phosphate.
