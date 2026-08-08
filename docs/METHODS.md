# Methods

Computational methods for the cryptic inositol phosphate (IP) binding site
pipeline. Run-specific text is written to `results/publication/METHODS_AUTO.md`
when the publication package executes.

## Overview

The pipeline identifies buried IP binding pockets — sites where an inositol
phosphate acts as a structural cofactor rather than a diffusible signalling
ligand. It has four stages:

1. **Ground truth** — every inositol phosphate structure in the PDB, with
   per-ligand-copy burial measurements (`scripts/build_ip_validation_dataset.py`).
2. **Description** — 40 physically interpretable pocket descriptors
   (`cryptic_ip.analysis.features`).
3. **Labelling** — pockets matched to observed ligands by atom-level overlap
   (`cryptic_ip.analysis.labeling`).
4. **Modelling** — nested, grouped, calibrated cross-validation against an
   interpretable rule-based baseline (`cryptic_ip.analysis.ml_classifier`).

---

## 1. Burial measurement

### Per ligand copy, size-normalised

Burial is measured **for each ligand copy separately** and normalised by the
ligand's own surface:

```
relative_sasa = SASA(copy inside the complex) / SASA(same copy in isolation)
```

This ratio is dimensionless and independent of both ligand size and copy number.
Absolute, copy-summed SASA — used in earlier versions — is neither: a structure
with six InsP6 copies reported roughly six times the SASA of one with a single
copy, so a crystallisation artefact determined the burial class.

Four complementary measurements are reported, because no single number separates
a buried cofactor site from a deep surface groove:

| Measurement | What it answers |
|---|---|
| `relative_sasa` | How much of the ligand's surface can solvent still reach? |
| `relative_phosphate_sasa` | The same, for the phosphate groups a polyanion site must sequester |
| `burial_depth` | How far below the molecular surface does the site sit? |
| `enclosure` | What fraction of directions out of the site are blocked by protein? |

### Classification thresholds

Applied to the larger of the whole-ligand and phosphate ratios — a conservative
choice, so a ligand buried to the ring but with phosphates in solvent is not
counted as sequestered.

| Class | Relative SASA |
|---|---|
| `cryptic` | ≤ 0.12 |
| `semi_cryptic` | ≤ 0.25 |
| `surface` | > 0.25 |
| `crystal_artifact` | fewer than 8 protein heavy-atom contacts within 4.5 Å |

**These boundaries are calibrated on measurements, not on the literature
description.** `scripts/calibrate_controls.py` measures the deposited controls:

| Control | PDB | Rel. SASA | Rel. P-SASA | Depth | Enclosure | Basic | Site volume |
|---|---|---|---|---|---|---|---|
| ADAR2 | 1ZY7 | **0.093** | 0.089 | 5.76 Å | 0.941 | 8 | 1525 Å³ |
| Btk PH | 1BWN | 0.253 | 0.266 | 5.98 Å | 0.668 | 3 | 491 Å³ |
| PLCδ1 PH | 1MAI | 0.373 | 0.348 | 4.68 Å | 0.598 | 5 | 1532 Å³ |
| HDAC1 | 5ICN | 0.436 | 0.444 | — | 0.629 | — | — |
| Pds5B | 5HDT | 0.466 | 0.460 | 5.63 Å | 0.738 | 8 | 895 Å³ |

HDAC1's depth, basic-residue count and site volume are shown as pending: the
values previously in this table described the `6A0` site and are not carried
over, since that component is not an inositol phosphate.

The boundary was initially set to 0.05, from the description of the ADAR2 InsP6
as encapsulated with only an 8.4 × 4.6 Å window (Macbeth et al., *Science*
309:1534, 2005). The measurement puts it at 0.093: a narrow window still exposes
a measurable fraction of a 36-atom ligand.

The panel leaves a gap — the sequestered ligand at 0.093, exposed ones at
0.253–0.466 — and the boundary of **0.12** sits inside it. Pds5B, annotated in
the literature as a frequent crystallisation artefact, measures as the most
exposed of all five, consistent with that annotation.

**The buried side of the panel is a single structure.** HDAC1 was previously
recorded as buried at 0.093, which would have given two. That measurement was
made on chemical component `6A0`, which carries no phosphate and is therefore not
an inositol phosphate at all; it was admitted only because the old identifier
whitelist listed it, and it was selected because burial is measured on the *most
buried* matching copy. Identifying ligands from coordinates (below) excludes it,
and the most buried genuinely phosphorylated copy in 5ICN is a solvent-exposed
InsP6 at 0.436. HDAC1 is an exposed control.

That correction is worth stating plainly rather than absorbing quietly: the
boundary separating the positive class now rests on **one** sequestered
structure. It should not be treated as settled. `scripts/burial_survey.py`
re-measures every deposited inositol phosphate complex in the bundled set for
exactly this reason, and reports whether the distribution is bimodal at all —
if it is not, the cryptic/surface split is a modelling choice rather than a
discovered boundary, and every downstream metric should be read in that light.

### The whole deposited set: burial is continuous, not two classes

Five controls cannot say whether a cryptic/surface dichotomy exists at all.
`scripts/burial_survey.py` re-measures every inositol phosphate complex in the
bundled dataset with the per-copy definition and reports the distribution.
**135 of 136 entries were measured**; one contained no phosphorylated inositol.
That the coordinate-based ligand test resolved 99 % of deposited entries is
itself the strongest evidence that it is not too strict.

| q01 | q05 | q10 | q25 | q50 | q75 | q90 | q95 | q99 |
|---|---|---|---|---|---|---|---|---|
| 0.070 | 0.091 | 0.126 | 0.216 | 0.338 | 0.553 | 0.773 | 0.816 | 0.896 |

Two independent boundary estimates were computed, deliberately, because either
one alone would be easy to over-read:

| estimate | value |
|---|---|
| density minimum (trough between modes) | **0.138** |
| Otsu (maximum between-class variance) | **0.463** |
| configured `CRYPTIC_RELATIVE_SASA_MAX` | 0.120 |

**They disagree, and that disagreement is the result.** A genuine two-class
structure would place both estimates in the same gap. Instead the density
estimator finds a shallow trough near the low tail while Otsu — which splits to
maximise between-class variance — cuts near the middle of a broad spread, which
is what it does when there is one wide mode rather than two. The quantiles agree:
burial runs smoothly from 0.07 to 0.90 with no chasm anywhere.

The conclusion is therefore mixed, and both halves matter:

* **The configured boundary is well placed.** 0.12 sits within one histogram bin
  (0.025) of the density minimum at 0.138, and selects 14 of 135 entries as
  cryptic. It was calibrated on ADAR2 alone and independently lands on the
  population's sparsest region, which is meaningful corroboration.
* **But it is a threshold on a continuum, not a natural class boundary.**
  Sequestration is a matter of degree across deposited inositol phosphate
  complexes. The positive class is *defined* by this cutoff rather than
  discovered in the data, and every downstream metric — AUROC, precision@k,
  enrichment — is conditioned on that choice. Reporting sensitivity of the model
  to the boundary is more informative than reporting performance at one value.

Two limits on this survey. The bundled set is essentially one chemistry (134
InsP6, 1 InsP5), so it calibrates InsP6 sequestration and says nothing about
whether InsP3/InsP4 sites distribute the same way. And 6 entries fall under the
crystal-artefact rule (fewer than 8 protein contacts), which is a floor on how
clean any deposited-structure ground truth can be.

### Ligands are identified from coordinates, not from a list of codes

Which molecules count as inositol phosphates decides both what burial is measured
on and which pockets become positive training labels, so it cannot rest on a
hand-written list of PDB chemical component identifiers. Such a list has two
defects: a site whose ligand code is absent is invisible rather than negative,
and a code can name the wrong chemistry — the previous list contained `INS`
(*myo*-inositol), which carries no phosphate at all and is recorded elsewhere in
this codebase as an unphosphorylated negative reference.

`cryptic_ip.analysis.inositol_detection` decides from the atoms instead:

| Criterion | Test |
| --- | --- |
| Inositol core | six carbons in a ring at C–C bonding distance (≤ 1.75 Å) |
| Hexahydroxylation | an oxygen within 1.65 Å of ≥ 5 of the 6 ring carbons |
| Phosphorylation | phosphorus within 1.90 Å of one of those oxygens |
| Series | `InsP{n}` from the phosphorus count, as for formulae |

This mirrors the rule the database module applies to reported formulae, so both
halves of the pipeline now define an inositol phosphate the same way. It needs no
network access and no vocabulary, so a regioisomer, a pyrophosphate or a deoxy
analogue is recognised on its structure rather than on whether anyone typed its
code. Exact name matching remains available for callers measuring one named
component.

### Burial depth does not discriminate on real structures

The panel shows something the synthetic benchmark cannot: **depth is
uninformative on deposited structures.** The values are completely interleaved,
and the *largest* depth in the panel belongs to a surface negative:

| | buried | exposed |
|---|---|---|
| depth (Å) | ADAR2 5.76 | Btk **5.98**, PLCδ1 4.68, Pds5B 5.63 |
| enclosure | ADAR2 0.941 | Btk 0.668, PLCδ1 0.598, HDAC1 0.629, Pds5B 0.738 |

No threshold on depth isolates the buried control: its 5.76 Å sits strictly
inside the 4.68-5.98 Å spread of the exposed ones, so any cutoff admitting it
admits a surface site too. On idealised synthetic spheres
the same measure separates them perfectly (22 Å against 5 Å).

The reason is the definition. Depth is the distance to the *nearest* solvent-
exposed atom — a minimum over a large set. A real protein surface is irregular
enough that some exposed atom lies within a few Å of almost any interior point,
so the minimum saturates. An idealised sphere has no such irregularity, which is
exactly why the synthetic benchmark could not have revealed this, and it is a
concrete limit on what synthetic validation can establish.

Enclosure separates the panel cleanly (buried 0.941 against maximum exposed
0.738) because it integrates over directions rather than taking a
minimum. Enclosure is therefore the reliable burial discriminator on real data,
and depth should be read as a supporting descriptor rather than a criterion.

The rule-based scorer still assigns depth 22 % of its weight, which on this
evidence buys nothing. Rebalancing toward enclosure is the indicated change, but
it is deferred: five controls are too few to fit weights on, and the scorer is a
baseline rather than the deployed model. The finding is pinned in
`tests/test_control_calibration.py` so it cannot be quietly forgotten.

### Implementation details that matter

- **SASA** is a vectorised Shrake–Rupley integration (`cryptic_ip.analysis.geometry`)
  with a 1.4 Å probe and 512 sample points per atom by default. A local
  implementation is required because the occluding atom set and the measured atom
  set must be chosen independently — that is what makes the in-complex vs
  in-isolation ratio possible.
- **The surface is defined on the holo structure.** Removing the ligand first
  would let the probe enter the vacated cavity, so a fully enclosed site would
  report near-zero depth.
- **An atom counts as surface at ≥ 5 Å² SASA.** A near-zero threshold makes depth
  fragile: atoms bordering a narrow crevice pick up a fraction of an Å², and a
  ligand at the centre of a 22 Å sphere then reports ~5 Å instead of ~22 Å.
- **Sibling ligand copies are excluded from the measurement context**, so crystal
  packing between copies cannot make a surface ligand look buried.
- **Phosphate groups are identified by element and P–O bond distance**, not atom
  names. Name-prefix matching both over-matches (bridging `PA`/`PB` atoms) and
  under-matches (unconventional names).
- **Waters and hydrogens are excluded**, so structures modelled with and without
  them are comparable.

---

## 2. Pocket descriptors

`cryptic_ip.analysis.features` computes 40 descriptors in seven blocks. Full
per-descriptor rationale: `feature_documentation()`.

| Block | Descriptors |
|---|---|
| Geometry | fpocket volume, convex-hull volume, alpha-sphere count and density, radius of gyration, asphericity, extent |
| Burial | burial depth, enclosure, buried-residue fraction, mean relative SASA |
| Accessibility | lining-residue SASA (mean, median, min, max, total) |
| Composition | basic / acidic / aromatic / hydroxyl / polar counts and fractions, Kyte–Doolittle hydropathy |
| Charge geometry | coordinating nitrogen count, their distance distribution and dispersion, net formal charge, charge density, charge balance |
| Electrostatics | screened Coulomb potential (always available), APBS potential (optional) |
| Confidence | pLDDT mean, minimum, fraction ≥ 70 |

Two corrections relative to earlier versions:

- **`pocket_depth` now means depth.** It was previously populated from fpocket's
  *mean local hydrophobic density*, a composition statistic unrelated to depth,
  while the genuine geometric depth was computed and discarded. fpocket's value
  is still reported, under its accurate name.
- **Electrostatics are always populated.** The APBS descriptor was `NaN` for
  every row whenever APBS was unavailable, which was the common case. A
  Debye–Hückel screened Coulomb sum over formal side-chain charges now provides
  an always-computable surrogate in the same kT/e units. It is a continuum
  approximation — no explicit solvent, no titration shifts, no dielectric
  boundary — and is not a substitute for a Poisson–Boltzmann calculation.

All descriptors use chain-aware residue keys `(model, chain, resseq, icode)`.
Indexing by residue number alone collides across chains, which inflated basic
residue counts by the number of chains in a homomer.

---

## 3. Pocket labelling

A pocket is labelled by the **fraction of ligand heavy atoms** within 4.0 Å of
its alpha spheres:

| Overlap | Label | Used in training |
|---|---|---|
| ≥ 30 % | positive | yes |
| ≤ 5 % | negative | yes |
| between | ambiguous | **no — excluded** |

Excluding the ambiguous band rather than calling it negative avoids injecting
label noise exactly where the decision boundary lies.

Centroid-to-centroid distance, used previously, is unsuitable: InsP6 spans about
11 Å, so its centroid can lie more than 8 Å from the centre of the pocket that
holds it, while an unrelated neighbouring pocket can fall inside 8 Å.

**Detector recall is reported.** A structure containing a ligand that yields no
matching pocket is a pocket-detection failure, and it bounds every downstream
result: a site fpocket never proposes cannot be scored. The count appears in
`results/ml_training/labeling_summary.json`.

**Negatives come from two sources**: non-site pockets in IP-binding proteins, and
all pockets of matched high-resolution structures containing no inositol
phosphate (`--n-decoys`). Without the second source every negative comes from an
IP-binding protein, which makes the benchmark easier than a proteome screen.

---

## 4. Model training and evaluation

### Protocol

| Aspect | Choice | Reason |
|---|---|---|
| Outer loop | 5-fold `StratifiedGroupKFold` | Unbiased performance estimate |
| Inner loop | 3-fold grouped randomised search | Hyperparameters never see the outer fold |
| Grouping | UniProt accession, falling back to PDB ID | The PDB holds many entries per protein; splitting by entry leaks |
| Selection metric | Average precision | Positives are rare; AUROC is insensitive to precision in that regime |
| Imbalance | Class weights | Oversampling duplicates rows from one protein across a split boundary |
| Calibration | Group-disjoint held-out slice per fold | Calibrating on fitted data maps overconfidence onto overconfidence |
| Missing values | Median impute **with indicator** | A missing APBS value is itself informative |
| Threshold | Chosen on out-of-fold predictions | The threshold is a fitted parameter |

Candidate models: L2 logistic regression (the interpretable floor), random
forest, extremely randomised trees, histogram gradient boosting, and XGBoost when
installed. Every candidate is trained on identical splits with identical seeds,
so comparisons are paired.

### Reported metrics

- **Discrimination** — AUROC, average precision, MCC, F1, balanced accuracy.
- **Calibration** — Brier score, expected calibration error.
- **Screening utility** — enrichment factor and precision at the top 1 %, 5 %
  and 10 % of the ranked list. A proteome screen inspects the head of a ranked
  list, so these predict saved effort better than accuracy at a threshold.
- **Uncertainty** — bootstrap confidence intervals resampled over **groups**.
  Resampling pockets treats several pockets from one protein as independent
  observations and understates uncertainty.
- **Model comparison** — DeLong's test on the shared out-of-fold predictions.
  Overlapping bootstrap intervals do not establish that two models are
  indistinguishable when both score the same samples.

### Baseline comparison

The learned model is compared with the rule-based `PocketScorer` on the same
rows, with a DeLong test. If the transparent weighted score matches the learned
model, the learned model adds complexity without adding information, and the
report says so.

---

## 5. Rule-based score

`PocketScorer` encodes the published criteria directly, as smooth monotone
functions of each measurement:

| Component | Weight | Midpoint |
|---|---|---|
| Mean lining SASA | 0.25 | 20 Å² |
| Burial depth | 0.22 | 12 Å |
| Basic residue count | 0.20 | 3.5 |
| Enclosure | 0.13 | 0.75 |
| Electrostatic potential | 0.10 | 3 kT/e |
| Cavity volume | 0.10 | 300–800 Å³ plateau |

The components were previously step functions — 4 basic residues scored 0.8 and
3 scored 0.4 — which made the score unstable under measurement noise for no
reason grounded in the biology. `CalibratedRuleScorer` maps the composite onto a
probability by logistic regression when labels are available; the raw composite
is a weighted average of bounded components and must not be read as one.

---

## 6. Validation controls

### Tier-1 gate

| Control | PDB | Role |
|---|---|---|
| ADAR2 | 1ZY7 | Cryptic InsP6 positive |
| PLCδ1 PH | 1MAI | Surface InsP3 negative |

Pass criterion: both controls pass individually **and** tier-1 separation > 0.50.

**Site selection.** The validators select the pocket to grade by *ligand-atom
overlap*, the same criterion used for training labels. They previously selected
the pocket whose centre was nearest the ligand centroid, which picks the wrong
pocket on ADAR2: the pocket with 100 % ligand overlap scores 0.644, while the
nearest-centre pocket scores 0.432. The gate was grading a pocket that does not
contain the inositol phosphate.

**Burial criteria use the relative measure.** The control pass criteria compare
relative SASA against the calibrated boundary rather than raw Å², since an
absolute cutoff means different things for InsP3 and InsP6 and scales with the
number of ligand copies in the crystal.

### Tier-2 panel

| Control | PDB | Notes |
|---|---|---|
| Pds5B | 5HDT | Often a surface/crystal artefact |
| HDAC1 | 5ICN | Measures as surface (0.436); its InsP6 is solvent-exposed |
| Btk PH | 1BWN | Surface negative (decoy mode) |

### Offline self-check

`cryptic-ip self-check` builds synthetic structures whose burial follows from
their construction — a ligand at the centre of a closed shell *must* measure as
cryptic — and verifies the pipeline recovers it. It needs no network access. It
validates the machinery, not the biology.

---

## 7. Statistics

- Proportions: Wilson score intervals for screening summaries, exact
  Clopper–Pearson intervals for confirmatory claims.
- Multiple testing: Benjamini–Hochberg FDR for screening, Holm–Bonferroni
  family-wise control for confirmatory comparisons.
- Group comparisons: assumption checks **select** the test — Welch's t-test when
  Shapiro–Wilk normality holds, Mann–Whitney U otherwise. Earlier versions
  reported the assumption check and then applied a t-test regardless.
- Effect sizes: Hedges' g (small-sample corrected) for the parametric case,
  Cliff's delta for the non-parametric case.
- Enrichment: two-sided permutation tests with FDR correction.

---

## 8. Reproducibility

```bash
python scripts/build_ip_validation_dataset.py --n-decoys 500   # requires RCSB access
python scripts/extract_pocket_features.py --jobs 8
python scripts/train_ml_classifier.py
```

Every stage writes provenance: the ligand registry and the queries that produced
it, per-file SHA-256 checksums and retrieval timestamps, API request statistics,
software versions, seeds, and label counts. Re-running is free where responses
are cached, and interrupted runs resume.

See also: [DATA_ACQUISITION.md](DATA_ACQUISITION.md), [VALIDATION.md](VALIDATION.md),
[REPRODUCIBILITY.md](REPRODUCIBILITY.md), [PUBLICATION_PACKAGE.md](PUBLICATION_PACKAGE.md).
