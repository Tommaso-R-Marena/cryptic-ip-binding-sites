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
| ADAR2 | 1ZY7 | **0.093** | 0.089 | 5.76 Å | 0.941 | 8 | 1488 Å³ |
| Btk PH | 1BWN | 0.253 | 0.266 | 5.98 Å | 0.668 | 3 | 500 Å³ |
| PLCδ1 PH | 1MAI | 0.373 | 0.348 | 4.68 Å | 0.598 | 5 | 1534 Å³ |
| HDAC1 | 5ICN | 0.436 | 0.444 | 4.39 Å | 0.629 | 4 | 939 Å³ |
| Pds5B | 5HDT | 0.466 | 0.460 | 5.63 Å | 0.738 | 8 | 909 Å³ |

HDAC1's row is measured at its InsP6 site. Every value in this table comes
from a single calibration run under a single scorer, so composites and site
volumes are comparable across rows.

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

That correction is worth stating plainly rather than absorbing quietly: within
this five-structure panel the boundary rests on **one** sequestered example.
That is why the panel is not the only evidence for it — the population survey
below measures all 136 deposited complexes, and independently places the
boundary where the controls put it.

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

### How much rides on the cutoff

Because burial is continuous, the positive class is whatever the boundary says it
is. The survey therefore reports class size across candidate cutoffs rather than
at one:

| boundary | 0.05 | 0.08 | 0.10 | **0.12** | 0.15 | 0.20 | 0.25 | 0.30 |
|---|---|---|---|---|---|---|---|---|
| positives | 0 | 2 | 11 | **14** | 16 | 27 | 40 | 58 |
| fraction | 0.0 % | 1.5 % | 8.1 % | **10.4 %** | 11.9 % | 20.0 % | 29.6 % | 43.0 % |

This is the strongest evidence for the configured value, and it is independent of
the density estimate. **The class size is flat where the boundary sits.** Moving
the cutoff across 0.10-0.15 changes the positive count by 5 entries; the same
0.05 shift at 0.20-0.25 changes it by 13, and at 0.25-0.30 by 18. The boundary is
therefore in the least sensitive part of the curve — nearly three times less
sensitive than the region immediately above it — which is what makes a threshold
on a continuum defensible even though the continuum has no natural break.

It also shows the cost of the original choice. The boundary was first set to 0.05
from the literature description of the ADAR2 site; on the deposited set that
cutoff yields **zero** positives, so the pipeline would have had no positive class
at all.

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

### Burial depth discriminates weakly; enclosure strongly

This section previously said depth was **uninformative** on deposited
structures. That was an overstatement from five structures, and the population
survey corrects it.

On the five-structure panel the values do interleave, and the *largest* depth
belongs to a surface control:

| | buried | exposed |
|---|---|---|
| depth (Å) | ADAR2 5.76 | Btk **5.98**, PLCδ1 4.68, HDAC1 4.39, Pds5B 5.63 |
| enclosure | ADAR2 0.941 | Btk 0.668, PLCδ1 0.598, HDAC1 0.629, Pds5B 0.738 |

No threshold on depth isolates ADAR2 there. But a panel that small cannot say
whether depth carries information, only that it does not separate these five.
PR #40's manuscript draft had reached the opposite conclusion from an even
smaller comparison — ADAR2 at 6.55 Å against about 2 Å for two PH domains, "a
~3× difference". Both claims rested on three to five structures.

`scripts/burial_survey.py` settles it on all 129 eligible deposited entries
(crystal artefacts excluded), scoring each descriptor against relative SASA:

| descriptor | Spearman ρ vs relative SASA | AUROC, buried vs exposed | 95 % CI |
|---|---|---|---|
| burial depth | **+0.02** | **0.72** | 0.58–0.84 |
| enclosure | −0.81 | 0.99 | 0.97–1.00 |

Buried is relative SASA ≤ 0.12 (14 entries), exposed > 0.25 (89); the
ambiguous middle band is left out rather than forced into either class.

Both earlier claims were wrong, in opposite directions:

* **Depth is not uninformative.** It separates the buried extreme from the
  exposed one with AUROC 0.72, and the interval excludes chance.
* **But it does not discriminate the way PR #40 described.** 0.72 is weak, and
  ρ ≈ 0 means depth has essentially no monotone relationship with burial across
  the continuum — it tells the extremes apart somewhat without tracking the
  degree of burial in between. A "~3× difference" between one buried and two
  surface structures does not generalise to a clean separation.

The likely reason is the definition. Depth is the distance to the *nearest*
solvent-exposed atom — a minimum over a large set — and a real protein surface
is irregular enough that some exposed atom lies within a few Å of most interior
points, so the minimum compresses. An idealised sphere has no such
irregularity, which is why the synthetic benchmark, where depth separates
perfectly (22 Å against 5 Å), could not have revealed this.

Two caveats on reading the table:

* **Enclosure's 0.99 is partly by construction.** Enclosure and relative SASA
  both measure how completely protein surrounds the same ligand copy, so their
  agreement is two views of one property rather than independent validation.
  Depth is the more independent measure, which makes its weaker showing more
  informative, not less.
* **This tests depth measured at the ligand centroid.** PR #40's figure used
  the fpocket *pocket centre*, a proxy for the same location. The ligand
  centroid is the more direct measurement, but the survey does not test the
  pocket-centre operationalisation itself.

The rule-based scorer assigns depth 22 % of its weight and enclosure 13 %.
Against these numbers that allocation is backwards, and rebalancing toward
enclosure is the indicated change. It is left to the repository owner: the
scorer is a baseline rather than the deployed model. The panel observation and
the population result are both pinned in tests, so neither can be quietly
dropped.

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

### Pre-registered benchmark: results

The protocol above is replaced by `docs/ANALYSIS_PLAN.md`, run end to end in CI
by `.github/workflows/benchmark.yml`. The run reported here is 35949588200 at
commit 72f1ae4. The report is `results/benchmark/benchmark_report.json`, and the
page rendered from it is `results/benchmark/report.html`.

**Data.**

- 367 entries (409 qualified; the rest were over the atom limit or failed extraction,
  and each is counted in the run's accounting).
- 78,920 pockets from fpocket on ligand-free structures.
- Two homology groupings:
  - sequence: MMseqs2, ≥ 30 % identity over ≥ 50 % of the shorter chain; 85 groups;
  - strict: sequence links plus Foldseek links with TM-score ≥ 0.5; 40 groups.
- Temporal holdout: the latest 20 % of strict groups by first release date. That is
  21 entries (7,312 pockets), never seen by any fitting, selection or threshold.

| task | development positives (sequence / strict groups) | largest group's share | holdout positives |
|---|---|---|---|
| `ip_site` | 765 (70 / 29) | 7 % / 40 % | 55 (8 groups) |
| `cryptic_ip_site` | 60 (6 / 5) | 87 % / 87 % | 0 |
| `burial` | 60 (6 / 5) | 87 % / 87 % | 0 |

**Finding the inositol phosphate site (`ip_site`).** Nested grouped
cross-validation jointly selects the model family and its hyperparameters in the
inner loop, and fits the threshold and calibration there too. It uses 39 descriptors,
with no B-factors and no ligand. Every interval resamples whole homology groups.

| evaluation | learned model ROC-AUC | PR-AUC | rule-based ROC-AUC | learned − rule ROC-AUC |
|---|---|---|---|---|
| CV, sequence groups (3 repeats) | 0.951 [0.932, 0.966] | 0.358 [0.260, 0.489] | 0.879 | +0.072 [0.046, 0.093] |
| CV, strict groups | 0.929 [0.907, 0.954] | 0.187 [0.166, 0.590] | 0.879 | +0.051 [0.022, 0.069] |
| **temporal holdout (locked model)** | **0.871 [0.782, 0.932]** | 0.136 [0.105, 0.250] | 0.730 | **+0.141 [0.033, 0.200]** |
| shuffled labels (control) | 0.491 [0.477, 0.509] | 0.010 | 0.493 | – |

The learned model generalises to protein families released after everything it
was trained on, and there it beats the rule-based score by more than it does in
cross-validation. The gap from 0.95 (cross-validation) to 0.87 (holdout) is the
cost of meeting new families. The shuffled-label control is at chance.

**The hull-depth hypotheses (H1, H2) are not evaluable.** For both buried-site
tasks, one homology group holds 87 % of the 60 positives. By the rule fixed before
the run, no grouped cross-validation estimate is then evidence, and the holdout
contains no buried site. The point estimates are not interpreted:

- H1 (paired ROC-AUC with minus without hull depth): −0.087 [−0.256, 0.398] with
  sequence groups, 0.308 [−0.001, 0.369] with strict groups;
- H2: 0.063 [−0.212, 0.146] and −0.332 [−0.410, 0.048].

The PDB holds about six independent families with a buried inositol phosphate.
That is a property of the deposited record, not of this pipeline, and no
re-analysis of PDB data alone can answer H1 or H2. Four `cryptic_ip_site` runs
also failed to fit, on folds holding a single positive. The fitting guard that
prevents this was added after the run (commit 19d193d), and it cannot change these
decisions.

For `ip_site`, hull depth adds a little in cross-validation: 0.951 against 0.931
without it (sequence), and 0.929 against 0.924 (strict). On the holdout it is
0.871 against 0.861. The plan did not pre-register this comparison, so it is
reported as a description, not a test.

**Controls.**

- `ip_site` and `cryptic_ip_site` pass their shuffled-label controls.
- `burial` fails: 0.613 [0.506, 0.682]. Diagnostic D1 in the plan separates a chance
  excursion from a leak with 30 further shuffles (run 35965940997,
  `results/null/permutation_null.md`). The verdict is **chance**: across 30 shuffles the
  pooled ROC-AUC is 0.497 ± 0.051 (2.5–97.5 percentiles 0.406–0.600), and only one
  shuffle exceeds 0.60. That one is the full run's own shuffle (repeat 0), which
  reproduces 0.613 exactly. Ten `cryptic_ip_site` shuffles give 0.498 ± 0.049.
  The control's criterion is too strict: its interval resamples groups within **one**
  shuffle and so ignores the between-shuffle standard deviation of about 0.05, which is
  as large as the interval's own half-width. A single failure is therefore not evidence
  of a leak. As the plan fixed in advance, H2 stays not evaluable.

**Which single descriptors carry the signal** (`docs/EXPLORATION_PLAN.md`, with
every analysis in `results/exploration/ledger.jsonl`). Each descriptor was signed
by a direction predicted from physics before any value was computed. Each was scored
on development pockets with strict-group intervals, and Benjamini–Hochberg was applied
across all 82 tests.

- **Electrostatics comes first.** The Debye-screened Coulomb potential at the pocket
  centre ranks the true site first in 43 % of structures, against 4 % by chance
  (ROC-AUC 0.855). Counts of basic residues and basic nitrogens follow (0.85).
- **No zero-parameter descriptor beats the rule-based score** (0.879).
- **The confirmatory set was chosen by rule and all three members replicate** on the
  temporal holdout (Holm across the set):

  | descriptor | holdout ROC-AUC [95 % CI] | Holm p |
  |---|---|---|
  | `n_basic_nitrogens` | 0.782 [0.660, 0.847] | < 0.001 |
  | `electropositive_enclosure` | 0.766 [0.625, 0.833] | 0.012 |
  | `n_strong_basic_residues` | 0.692 [0.583, 0.830] | 0.012 |

- **This is the textbook determinant of phosphate binding, now shown out of family.**
  It is not a new descriptor.
- **Four priors failed**: `positive_charge_density`, `alpha_sphere_density`,
  `burial_depth` and `basic_nitrogen_dispersion` all rank IP sites *below* other
  pockets. Each is a ratio or a nearest-distance, and each falls with pocket size,
  while fpocket's pockets on IP sites are large. So on this data the four are
  mostly size in disguise. Testing that needs a size-adjusted analysis, which is
  not yet pre-registered.

### Transfer test: buried phosphate-dense ligand sites

**Plan and runs.**

- **Plan.** `docs/TRANSFER_PLAN.md`, with `docs/TRANSFER_PLAN_AMENDMENT_1.md`. The
  amendment was written after Prepare and before any outcome was read.
- **Main run.** Transfer run 35972729955.
- **Secondary run.** Transfer-secondary run 35986927155 (T1b/T2b).

**Why.** The inositol phosphate benchmark could not evaluate its hull-depth
hypotheses: its buried positives come from about six families. So the same
physical question was asked of a broader class of ligands, defined by a formula
rule:

- at least 2 phosphorus atoms;
- at least 0.07 phosphorus per heavy atom;
- no inositol phosphate;
- not lipid-linked.

Thirty-five components passed. Among them are the nucleotide di- and
triphosphates and their analogues, PRPP, fructose-1,6-bisphosphate, isoprenoid
pyrophosphates and pyrophosphate.

**Data.**

- **Sample.** X-ray entries at 2.5 Å or better that hold a class ligand: 8,434, of
  which 38 also held an inositol phosphate and were removed. From the rest, a seeded
  sample of 1,200 was drawn; 1,195 were measured, giving 101,338 pockets.
- **Pipeline.** Descriptors, labels, homology groups and the temporal holdout come
  from the benchmark's code, unchanged.
- **Joint groups.** Homology groups were computed jointly with the inositol
  benchmark's entries.
- **Scale.** There are 914 buried (cryptic) class-ligand pockets in 174 sequence
  families (51 strict groups). The inositol set has 60 in about 6.

**The structural super-group.**

- Under strict grouping (Foldseek TM ≥ 0.5 links, joined transitively), one group,
  G:10JT, holds 61 % of the buried positives. It spans 688 of the 1,195 entries,
  across the nucleotide-binding folds: Ras-family GTPases (KRAS, HRAS),
  heterotrimeric G proteins, protein kinases (CDK2), HSP90, ATP synthase, myosin,
  carbamoyl-phosphate synthetase and others.
- Joining structural neighbours as connected components chains these folds together.
  That makes the strict grouping very conservative for broad ligand classes.

**Results.** All intervals are 95 % intervals from resampling whole homology groups.

| task | learned ROC-AUC, sequence CV | strict CV | temporal holdout (8 families) | rule-based ROC-AUC |
|---|---|---|---|---|
| buried class-ligand site (`cryptic_ip_site`) | 0.915 [0.901, 0.943] | 0.941 [0.927, 0.970] | **0.983 [0.955, 1.000]** | 0.676 (CV), 0.867 (holdout) |
| buried vs surface (`burial`) | 0.830 [0.799, 0.882] | 0.788 [0.750, 0.871] | 0.900 [0.754, 0.995] | 0.508 (CV), 0.806 (holdout) |
| any class-ligand site (`ip_site`) | 0.905 [0.878, 0.925] | 0.892 [0.845, 0.916] | 0.917 [0.863, 0.962] | 0.681 (CV), 0.761 (holdout) |

| hypothesis | hull depth added (paired ROC-AUC), sequence | strict | holdout | decision |
|---|---|---|---|---|
| T1 (all entries) | +0.002 [−0.008, 0.023] | +0.009 [−0.019, 0.015] | +0.004 [−0.003, 0.010] | **not evaluable** (largest strict group 61 %) |
| T2 (all entries) | +0.001 [−0.006, 0.011] | +0.005 [−0.017, 0.010] | +0.002 [−0.036, 0.050] | **not evaluable** |
| T1b (without G:10JT) | −0.010 [−0.021, 0.002] | **−0.086 [−0.131, −0.047]** | −0.000 | **inconclusive** |
| T2b (without G:10JT) | −0.004 [−0.015, 0.008] | **−0.021 [−0.039, −0.006]** | +0.007 [−0.026, 0.036] | **inconclusive** |

T1b/T2b use the remaining 507 entries: 359 buried positives in 96 sequence and 50 strict
groups, with largest shares of 9 % and 20 %.

**Controls.** The ten-permutation controls are at chance for every tested task:

- `cryptic_ip_site`: 0.504 ± 0.013;
- `burial`: 0.498 ± 0.018;
- T1b: 0.507 ± 0.015;
- T2b: 0.499 ± 0.028.

The untested `ip_site` task's single-permutation control gave 0.511 [0.5004, 0.522].
With about 99,000 pockets its interval is narrow. This is the single-permutation
miscalibration that diagnostic D1 identified, and it is reported, not re-decided.

**T3 (external, descriptive).**

- **Training set.** The model was trained on this dataset only, after removing every
  entry that shares a strict group with an inositol phosphate benchmark entry. That
  left 279 buried positives in 46 groups.
- **Buried IP sites.** On the inositol benchmark's buried IP sites it scores ROC-AUC
  **0.974 [0.839, 0.993]**, against 0.923 for the rule-based score.
- **Any IP site.** On all IP sites it scores 0.842 against 0.868.
- **Status.** The buried IP positives lie in only 5 groups, so this is descriptive
  and supports no claim of transfer.

**What this settles.**

1. **Recognition is learnable.** Buried phosphate-dense sites can be recognised by a
   learned model in families it has never seen, far above the hand-built score. That
   score was tuned on IP6 and is near chance at telling buried from surface sites
   outside inositol phosphates.
2. **Hull depth adds nothing.** Across 914 buried sites in 174 families, adding hull
   depth to a learned model changes ROC-AUC by +0.002. Outside the super-group, under
   strict grouping, it lowers it by 0.086.
   - The pre-registered decisions are "not evaluable" and "inconclusive", not
     "refuted", because the sequence-grouping intervals are not inside ±0.01.
   - The evidence points one way: hull depth does not help, and the screen's
     hull-depth gate should be reconsidered.

### Proteome screen and learned ranking

**Scope.** The screen covers the yeast, human and Dictyostelium AlphaFold proteomes: 37,384
proteins scored, from run 35935291031. The rule-based digest is in
`results/proteome_screen/DIGEST_run35935291031.md`.

**Known binders under the rule-based score.** All six rank near the top of the human
proteome. The rank percentiles are:

| protein | rank percentile |
|---|---|
| ADAR1 | 99.7 |
| ADAR2 | 99.0 |
| HDAC3 | 97.1 |
| HDAC1 | 96.9 |
| PDS5B | 93.7 |
| ADAT1 | 92.8 |

ADAR2 anchored the score threshold, so its rank is not independent evidence.

**The rule-based score's top hits are explained.** They fall into two classes:

- **Enzymes and carriers of other phosphate-dense anions:**
  - BPGM and its phosphoglycerate-mutase relatives (GPM2, GPM3, gpmA);
  - isopentenyl-diphosphate isomerase (IDI1, ipi);
  - aconitase (ACO2) and sulfite oxidase (SUOX);
  - mitochondrial carriers (SLC25A16, CTP1, mcfR).
- **Kelch β-propellers**, whose central channel is buried and basic: HCFC1, HCFC2, KEL1,
  KLHDC10 and ATRNL1.

**Learned ranking** (`docs/LEARNED_SCREEN_PLAN.md`, run 35992676292; report in
`results/learned_screen/LEARNED_SCREEN.md`).

- **Model.** The `ip_site` model was locked on the whole benchmark table (extra trees,
  selected by grouped inner CV).
- **Protein score.** Each protein's score is its best confident pocket (pLDDT ≥ 70).
- **Exclusion.** 2,226 proteins with an MMseqs2 homologue among the benchmark's UniProt
  sequences were excluded (≥ 30 % identity, ≥ 50 % coverage of the shorter sequence,
  E ≤ 1e-3).
- **Truth.** UniProt inositol phosphate annotations.
- **Intervals.** Resampled over proteome sequence clusters.

| organism | unseen proteins | annotated binders | learned ROC-AUC | rule ROC-AUC | learned − rule |
|---|---|---|---|---|---|
| human | 18,668 | 27 | 0.805 [0.709, 0.878] | 0.734 [0.655, 0.807] | +0.071 [−0.029, 0.162] |
| Dictyostelium | 11,117 | 4 | 0.900 [0.781, 0.963] | 0.943 [0.877, 0.997] | −0.043 [−0.131, 0.034] |
| yeast | 5,373 | 4 | 0.877 [0.584, 0.999] | 0.863 [0.630, 0.975] | +0.014 [−0.296, 0.303] |
| **pooled** | **35,158** | **35** | **0.832 [0.760, 0.890]** | 0.779 [0.716, 0.840] | +0.053 [−0.026, 0.127] |

- **Decisions (Holm).** L1 is **supported**: the learned model ranks annotated IP binders
  above other proteins across whole proteomes, on proteins homologous to nothing it was
  trained on. L2 is **not supported**: it does not significantly beat the rule-based score.
- **Precision of the candidates.**
  - The base rate of annotated binders among unseen human proteins is 0.14 %.
  - Among the top 27 human proteins it is 7.4 % (2 of 27): roughly 50-fold enrichment.
    Yeast is about 100-fold. Dictyostelium's top 25 contains no annotated binder.
  - These precisions are lower bounds, because unannotated true binders count as
    negatives. A candidate is roughly a 1-in-10 to 1-in-15 hypothesis, not a finding.
- **What the candidates are.** Most are explained by other phosphate- or sulfate-rich
  ligands:
  - PAPS-dependent sulfotransferases: CHST1, CHST4, HS3ST1, HS3ST5;
  - UDP-sugar glycosyltransferases: EXT1, GYS1, GSY2;
  - nucleotide and anion carriers: AAC3, YHM2, SLC25A27, mcfA/U/O;
  - ATP-binding motors and pumps: myosins, P-type ATPases, TMEM94, ATP8B4;
  - sugar-phosphate enzymes: PFK2, G6PD/ZWF1;
  - prenyl pyrophosphate and polyphosphate enzymes: COQ1, ppkA.

  WD40 and other β-propellers recur: COPA/COP1, WDR46, CSTF1, SEMA4A. So do known
  phosphoinositide-headgroup binders: ASAP1 (PH domain), MTMR11 (myotubularin), RLBP1
  (CRAL-TRIO).
- **The one coherent unexplained lead.** Two α-arrestins rank high: human ARRDC2 (rank 5)
  and yeast ART5 (rank 27). β-arrestins, which share the arrestin fold, are established
  IP6 binders, and α-arrestins are too distant in sequence to be excluded as homologues.
  Whether their top pocket corresponds to the β-arrestin IP6 site is the first thing to
  check. It does not (next section).

### Redocking benchmark

**Plan and run.** `docs/REDOCKING_PLAN.md`, redocking run 36020459246.

- The report is attempt 2. Dock alphafold 21 was re-run once: its first attempt docked
  all 9 copies but lost the artifact upload to a connection reset.
- Results are in `results/redocking/`: `redocking.json`, `copies.csv` (one row per
  selected copy), `census.csv`, `REPORT.md` and `report.html`. All were extracted from
  the Report job's log with `scripts/extract_log_block.py`.
- `configuration_audit.csv` holds the audit of the configuration exclusions (run
  36029338442).

**Census.**

- 662 IP copies were found in 367 benchmark entries; 462 copies are eligible.
- Exclusions: 167 whose configuration differs from the CCD template (the audit shows the
  difference is real, with a median absolute chiral volume of 2.36 Å³ at the differing
  centres, and 153 of the 167 are IHP), 32 crystal artefacts, and 1 copy that fails the
  CCD template on valence.
- 272 copies were selected, one per entry and burial class; 269 form the primary set, and
  the 3 incomplete copies are a flagged stratum.

**Accounting.**

- No copy went unreached: every shard finished inside its time budget.
- 19 of the 269 primary copies (7.1 %, in 5 strict groups) failed receptor preparation
  and have no docking outcome. 14 are PDB2PQR exceptions, 11 of them hydrogen-placement
  errors. The other 5 are large cryo-EM assemblies (7T3P, 7T3Q, 7T3R, 9YKY, 9YLI, each
  130,000–140,000 atoms) whose PQR and PDB outputs disagree in atom count, which is a
  limit of the receptor reader at that size.
- 17 of the 19 are surface copies, so the estimates below are conditional on a receptor
  that can be prepared.
- Repairing the reader would have re-run the whole workflow and the plan has no rule for
  it, so these failures are reported rather than fixed.
- 250 copies in 31 strict groups have an outcome.

**Decisions (Holm across R1–R4).**

| | question | estimate (group estimand) | Holm p | decision |
|---|---|---|---|---|
| R1 | protocol reliability | 0.110 [0.046, 0.193] | < 0.001 | **unreliable** |
| R2 | success(cryptic) − success(surface) | 0.299 [−0.056, 0.942] | – | **not evaluable** (4 cryptic groups) |
| R3 | AlphaFold cross-docking success | 0.045 [0.000, 0.122] | < 0.001 | **not trustworthy** |
| R4 | Vina score separates true site from decoy (ROC-AUC) | 0.744 [0.658, 0.873] | < 0.001 | **discriminates** |

**Outcomes on the primary set** (250 copies, 31 groups; per copy and per group).

| outcome | per copy | per group |
|---|---|---|
| top-pose success at 2 Å | 0.096 [0.056, 0.141] | 0.110 [0.046, 0.193] |
| best of 20 poses at 2 Å | 0.295 [0.236, 0.419] | 0.324 [0.221, 0.441] |
| success at 1 Å | 0.012 [0.002, 0.025] | 0.016 [0.001, 0.039] |
| success at 3 Å | 0.199 [0.135, 0.317] | 0.202 [0.115, 0.302] |
| phosphorus-only success at 2 Å | 0.100 [0.059, 0.143] | 0.111 [0.046, 0.194] |
| Spearman (score vs RMSD) within runs | 0.072 [0.013, 0.134] | 0.058 [−0.025, 0.146] |

**Where the failures are.** Of the 242 copies whose seed-1 top pose misses 2 Å, 230 are
scoring failures: the minimised crystal pose scores worse than the top docked pose, so
the search reached a near-native pose region and the function preferred something else.
Only 12 are sampling failures. The best-of-20 rate, three times the top-pose rate, says
the same thing.

**Strata (top-pose success, group estimand).**

- Burial: surface 0.044 [0.010, 0.084]; semi-cryptic 0.175 [0.085, 0.285]; cryptic 0.343
  [0.000, 0.750] in 4 groups, which is not evidence.
- Interface 0.041 [0.000, 0.103] against single-chain 0.124 [0.050, 0.223].
- Metal within 3 Å 0.202 [0.067, 0.402] against no metal 0.079 [0.029, 0.142].
- Method: X-ray 0.117 [0.045, 0.210]; cryo-EM 0.035 [0.000, 0.091] in 6 groups.
- Species: InsP6 0.114 [0.047, 0.199]; InsP3 0.049; InsP4 0.052; InsP5 0.000 in 4 groups,
  which is not evidence.
- Resolution shows no trend: 0.077, 0.140, 0.072 and 0.077 from ≤ 2.0 Å to > 3.0 Å.

**Controls.**

- *Seed noise.* The three seeds agree on success for 0.818 [0.699, 0.915] of copies, and
  the top score's seed-to-seed standard deviation is 0.222 [0.185, 0.265] kcal/mol.
- *Scoring functions and protonation, paired against the primary seed 1.* AD4 is better
  by 0.058 [0.009, 0.135]; Vinardo by 0.014 [0.000, 0.034]; the fully deprotonated
  ligand by 0.013 [−0.003, 0.035]; keeping metals by 0.051 [−0.060, 0.149] on the 42
  copies with a metal. Every arm is near the floor: AD4 reaches 0.100 [0.025, 0.195].
- *Site finding.* fpocket's top 3 pockets contain a true IP site for 0.330 [0.197, 0.472]
  of copies, and 39 % of structures have any positive pocket among them.
- *AlphaFold cross-docking.* 233 copies had a usable model. Success is 0.045
  [0.000, 0.122], and 0.000 [0.000, 0.001] on surface copies. Paired against the crystal
  receptor on the same copies the difference is −0.002 [−0.029, 0.016]: the model is not
  worse than the crystal, because both are near zero.

**What this settles.**

- The protocol recovers crystal IP poses in about a tenth of cases, so it is not a
  reliable tool for placing an IP ligand, and R1's "unreliable" label is what the data
  support.
- The failure is in ranking, not in search. That is the prediction the electrostatic
  re-ranking study (`docs/RERANK_PLAN.md`) was written to test, and it was pre-registered
  before this report was read.
- The score still separates a true site from a decoy pocket (R4). Site *detection*,
  which is what the screen does, is a different and easier task than pose prediction.
- R2 cannot be decided: the cryptic class holds 4 strict groups, below the plan's
  minimum of 5. The point estimates run the other way from the usual expectation, with
  buried sites easier than surface ones, but the interval is wide and no claim is made.

### The α-arrestin lead

**Plan and run.** `docs/ARRESTIN_PLAN.md`, arrestin run 36023420103; results in
`results/arrestin/` (`arrestin.json`, `ARRESTIN.md` with one dossier per protein,
`family.csv`, `mapping.json`, `report.html`).

**B1: the family.** α-arrestins are arrestin-fold proteins (Pfam PF00339 or PF02752) that
are not visual or β-arrestins. They were scored against every other unseen protein of
the three pooled proteomes: 35,158 proteins, including 20 unseen α-arrestins in 14
MMseqs2 clusters.

- The learned score gives ROC-AUC 0.757 [0.683, 0.856] (2,000 cluster resamples). The
  5th percentile is 0.694, so the family test is **supported**. The asymptotic
  Mann–Whitney p is 3.4 × 10⁻⁵, which assumes paralogues are independent.
- Per organism (descriptive): yeast 0.820 (6 clusters) and *Dictyostelium* 0.828
  (5 clusters) are both supported. Human is *not evaluable*: its 6 unseen α-arrestins
  fall in 3 clusters.
- The human signal is ARRDC2 alone. It ranks 11th of 35,158 in the pooled ranking;
  TXNIP and ARRDC1, 3, 4 and 5 all rank below 14,000th. ART5 ranks 135th.
- Positive control: the four classic arrestins are all *seen*, homologous to benchmark
  proteins, and rank at the 90.6–99.1th percentile of all scored proteins.

**B2–B3: the site.**

- *References.* The benchmark holds 7 classic-arrestin entries (1ZSH, 5TV1, 7F1W, 7F1X,
  7JTB, 7JXA, 7MOR) with 24 distinct IP sites. They all fall in a single strict
  homology group.
- *ARRDC2.* It aligns best to 5TV1 chain A (TM-score 0.641, normalised by the
  reference). Three of the four residues of site IHP_A_401 map (226→211, 227→212,
  332→290).
- *ART5.* It aligns best to 7F1W chain D (TM-score 0.678). Only one of the five residues
  of IHP_D_501 maps (171→306).
- *Overlap.* In both proteins the mapped site shares no residue with the top learned
  pocket (Jaccard 0.0, threshold 0.20), so criterion 1 fails.
- *Other α-arrestins.* No other α-arrestin's top pocket overlaps a mapped site. Their
  TM-scores against the references are 0.60–0.69; the classic arrestins' own models
  score 0.94–0.96.

**B4: conservation.**

- *ARRDC2.* Its UniRef50 cluster gives 119 homologues, but none of its mapped site
  residues is K, R or H. It is therefore **not conserved** (at least 3 basic positions
  are needed).
- *ART5.* Its cluster has 8 homologues, so conservation is *not evaluable* (at least 10
  are needed).

**B5: docking.**

- *Tasks.* There were 49 docking tasks with no failures.
- *Protocol-validity gate.* IP6 redocked into the 24 crystal arrestin sites had a mean
  top-pose success of **0.00** at 2 Å, against the 0.5 required. Every copy was
  configuration-eligible and surface-bound. The top poses lay 5–15 Å from the crystal
  ligand. The protocol is **not valid** on arrestin sites, so criteria 3 (convergence)
  and 4 (scores) are *not evaluable* and count as failing.
- *Descriptive scores (Vina, kcal/mol).* These are recorded because they were computed,
  not as evidence.

  | | ARRDC2 | ART5 |
  |---|---|---|
  | IP6 at the lead site | −4.28 | −6.58 |
  | IP6 at the top pocket | −5.73 | −6.24 |
  | IP6 at the five random negative pockets | −4.31 to −5.80 | −4.25 to −6.31 |
  | ATP at the lead site | −6.01 | −8.08 |
  | Largest-cluster fraction of the 15 top poses | 0.13 | 0.20 |

  - For ARRDC2, the lead-site score is weaker than every negative pocket and than the
    weakest positive control (−4.35; AlphaFold ARRB1/ARRB2/SAG controls −4.35 to −5.81).
  - For ART5, the lead-site score beats its negatives, but its poses do not converge
    (0.5 is required).
  - At both lead sites, ATP scores better than IP6.

**B6: verdict.**

| protein | 1 overlap | 2 conservation | 3 convergence | 4 scores | verdict |
|---|---|---|---|---|---|
| ARRDC2 (Q8TBH0) | fails (Jaccard 0.0) | fails (no basic site residue) | not evaluable | not evaluable | **not supported** |
| ART5 (P53244) | fails (Jaccard 0.0) | not evaluable (8 homologues) | not evaluable | not evaluable | **not supported** |

**Study E (short MD) was not run.** The task made it conditional on A and B. B's leads
failed on the structural criteria (overlap, conservation), which do not depend on
docking. Simulating a docked pose from a protocol that reproduces none of 24 crystal
arrestin–IP poses would have had no defensible starting structure.

### IP versus other polyanion sites (specificity)

**Plan and run.** `docs/SPECIFICITY_PLAN.md`, specificity run 36020458866; results in
`results/specificity/` (extracted from the run log).

**Question.** The learned `ip_site` model recognises buried polyanion sites in general.
Can the 39 descriptors tell an inositol phosphate (IP) site from a site for another
phosphate-dense ligand?

**Data.**

- **Label 1:** pockets on an IP copy in the IP benchmark table (run 35949588200).
- **Label 0:** pockets on a class ligand in the transfer table (run 35972729955).
- **Excluded:** pockets on neither ligand. The question is conditional on a polyanion
  site.
- **Groups:** the joint MMseqs2 and sequence-plus-Foldseek groups computed over both
  datasets.
- **Holdout:** the latest 20 % of joint strict groups.
- **Size:** 3,638 pockets in 1,522 entries. Development holds 734 IP pockets (64
  sequence / 21 strict groups) and 2,763 other-polyanion pockets (391 / 97).
- **Protocol:** the benchmark's runner, unchanged.

**Variants.**

- S1: all rows.
- S1b: without the joint strict group holding the most other-polyanion pockets, G:10JT
  (the nucleotide-binding folds).
- S1r: IP rows restricted to X-ray entries at 2.5 Å or better, the transfer set's rule.
  This checks that the model is not separating the two datasets instead of the
  chemistry.

| variant | sequence CV ROC-AUC | strict CV ROC-AUC | holdout ROC-AUC | largest IP group (seq / strict) | 10-permutation mean | decision |
|---|---|---|---|---|---|---|
| S1 | 0.862 [0.802, 0.905] | 0.818 [0.699, 0.886] | 0.843 [0.767, 0.931] | 25 % / **49 %** | 0.502 ± 0.013 | **not evaluable** (40 % rule, strict) |
| S1b | 0.830 [0.731, 0.902] | 0.757 [0.580, 0.878] | 0.742 [0.624, 0.844] (8 IP groups) | 14 % / 17 % | 0.493 ± 0.022 | **learnable** |
| S1r | 0.812 [0.704, 0.899] | 0.859 [0.680, 0.925] | 0.773 [0.693, 0.948] (3 IP groups: underpowered) | 15 % / 35 % | 0.506 ± 0.017 | **learnable** |

- **S1.** Joining the two datasets merges IP families with nucleotide-binding folds into
  one joint strict group holding 49 % of the IP pockets. By the plan's 40 % rule, S1 is
  therefore not evaluable, although every interval sits well above 0.5.
- **S1b and S1r.** Both clear 0.5 under both groupings, and S1b also on a powered
  holdout.
- **The gate.** The plan's gate opened on S1b: S1 was not evaluable only because of
  the 40 % rule, S1b is learnable, and S1r's point estimate exceeds 0.5.

**Re-ranking (conditional, run).**

- **Model.** The specificity model was locked on S1b (extra trees). Each pocket's
  combined score is P(site) × P(IP | site). Each protein's score is its best confident
  pocket.
- **Unseen proteins.** Proteins with an MMseqs2 homologue among either dataset's UniProt
  sequences were excluded: 4,300 proteins. That left 33,084 unseen proteins, 34 of
  them annotated IP binders.

| ranking | ROC-AUC [95 %] (cluster bootstrap) | Holm p | decision |
|---|---|---|---|
| L1 learned `ip_site` (same proteins) | 0.844 [0.771, 0.902] | – | – |
| L3a combined | 0.828 [0.756, 0.886] | < 0.001 | **supported** |
| L3b combined − L1 | −0.015 [−0.044, 0.015] | 0.321 | **not supported** |

**What this settles.**

- **Specificity is learnable at the pocket level.** Among polyanion sites, the
  descriptors separate IP sites from nucleotide and sugar-phosphate sites in families
  held out by homology, and also on comparable X-ray data.
- **It does not improve the proteome ranking.** L3b's interval is centred just below 0,
  so this is not "no effect of any size". A candidate list re-weighted by P(IP | site)
  should not be read as more trustworthy than the learned ranking.
- **Candidates.** 75 candidates (the top 25 per proteome); 30 of them are "explained" by
  the plan's keyword rule. The rule is crude. It misses "mitochondrial substrate
  carrier" names (mcfA, mcfU) and matches some broad keywords. Its counts are an
  approximate triage, not a classification.
- **ARRDC2 again.** It moves to rank 1 among unseen human proteins under the combined
  score, with P(IP | site) = 0.90.

### The screen's hull-depth gate

**Plan and run.** `docs/HULL_GATE_PLAN.md`, hull-gate run 36020459316; results in
`results/hull_gate/`.

**Arms.** The calibrated screen was re-aggregated from its 800 shards (2,845,762
pockets) with the hull-depth gate at 10 Å (current), at 5 Å, and removed. Every other
calibrated criterion was unchanged. Proteins were ranked by their best composite score
among pockets passing the arm's non-score gates.

**Evaluation.** 35,158 unseen proteins, 35 annotated IP binders, and 2,000 resamples of
20,073 MMseqs2 clusters.

| arm | ROC-AUC [95 %] | recall@546 | hits (annotated) |
|---|---|---|---|
| 10 Å gate | 0.701 [0.613, 0.789] | 0.114 | 546 (4) |
| 5 Å gate | 0.738 [0.645, 0.817] | 0.114 | 554 (4) |
| no gate | 0.747 [0.663, 0.821] | 0.114 | 554 (4) |

**Decision.**

- Removing the gate changes ROC-AUC by +0.046 [−0.015, 0.101], and relaxing it to 5 Å
  by +0.037 [−0.018, 0.091].
- Neither lower bound reaches the non-inferiority margin of −0.01.
- By the plan's rule the 10 Å gate is **kept**. `CALIBRATED_CRITERIA` is unchanged and
  cites the decision.
- The point estimates favour removing the gate, and the recall of the top 546 proteins
  is identical in all arms. The data do not show that the gate helps. They also do not
  show, at the pre-registered margin, that it can be dropped without cost.

### Results on the deposited set (superseded)

> **Superseded: these numbers overstate performance.** An audit found that the
> protocol below leaks. Cross-validation groups fell back to the PDB entry
> (the grouping read a column the dataset did not have), so the same protein -
> ADAR2 has 12 entries - and its homologues sat on both sides of folds; the
> `plddt_*` descriptors are crystallographic B-factors here and carry the removed
> ligand's ordering back in; the decision threshold and the best of five models
> were chosen on the out-of-fold predictions they were reported with; and the
> DeLong tests treated pockets clustered within proteins as independent. The
> pre-registered benchmark (`docs/ANALYSIS_PLAN.md`, `.github/workflows/benchmark.yml`)
> corrects each of these, on every inositol phosphate complex in the PDB rather
> than 136. Its results are in the section above.

Trained in CI (`.github/workflows/train-real-data.yml`) on the 136 measurable
RCSB entries; nested 5×3 grouped CV, intervals bootstrapped over proteins.
Descriptors are computed on **apo** structures - the ligand and every
non-polymer atom removed - because a descriptor computed with the ligand
present lets the ligand occlude its own pocket, a signal no AlphaFold target
can show. `ip_site_holo` repeats the first task with the ligand left in place,
on the same code, to measure that leak.

| task | positives (structures) | best model | ROC-AUC [95 % CI] | PR-AUC [95 % CI] | rule-based ROC / PR | DeLong ML − rule |
|---|---|---|---|---|---|---|
| ip_site (apo) | 307 (118) | extra trees | 0.975 [0.960–0.986] | 0.786 [0.724–0.840] | 0.934 / 0.497 | +0.040, p = 3.9×10⁻⁸ |
| ip_site (holo) | 309 (118) | extra trees | 0.982 [0.969–0.991] | 0.832 [0.777–0.876] | 0.950 / 0.634 | +0.032, p = 1.2×10⁻⁶ |
| cryptic_ip_site (apo) | 32 (13) | extra trees | 0.937 [0.869–0.992] | 0.864 [0.747–0.957] | 0.953 / 0.528 | −0.016, p = 0.47 |

- **The leak was real but small for the learned model** (ROC-AUC 0.982 → 0.975,
  PR-AUC 0.832 → 0.786) and larger for the rule-based score (PR-AUC 0.634 →
  0.497), whose SASA and depth terms read the ligand's occlusion directly.
- **On `ip_site` the learned model beats the rule-based score** decisively.
- **On `cryptic_ip_site` it does not**: 32 positive pockets from 13 structures
  give an interval 0.12 wide, and the difference from the rule-based score is
  not significant. Across the five candidate models on the same folds the
  AUROC runs from 0.937 to 0.991, itself a sign of how little data there is.
- **What these tasks measure.** `ip_site` asks which pocket in an inositol
  phosphate-binding protein holds the ligand; it does not ask whether the site
  is buried. `cryptic_ip_site` labels only buried sites positive, but surface
  inositol phosphate pockets become *ambiguous* (481 of them) rather than
  negatives, so the model is never asked to tell a buried site from a surface
  one - the distinction the screen depends on. Both tasks are measured inside
  proteins already known to bind an inositol phosphate, so neither number is a
  proteome-screen precision.

---

## 5. Rule-based score

`PocketScorer` encodes the published criteria directly, as smooth monotone
functions of each measurement:

| Component | Weight | Midpoint |
|---|---|---|
| Mean lining SASA | 0.25 | 20 Å² |
| Depth (distance to the nearest solvent-exposed atom) | 0.22 | 12 Å |
| Basic residue count | 0.20 | 3.5 |
| Enclosure | 0.13 | 0.75 |
| Electrostatic potential | 0.10 | 3 kT/e |
| Cavity volume | 0.10 | 300–1600 Å³ plateau (cavity, not ligand, volume) |

**Depth to the convex hull was tested as the depth input, and rejected.**
The depth component reads the distance to the nearest solvent-exposed atom. On
an apo structure - every AlphaFold model - that collapses at enclosed sites: the
walls of an empty cavity are themselves exposed, so ADAR2's InsP6 site reads
4.5 Å deep on its model, no deeper than the PH-domain surface sites (section
6a). Depth to the convex hull ignores internal cavities (ADAR2 17.3 Å; PH
domains 5.7-8.2 Å), and scoring depth from it, with a ramp set from the plan's
own thresholds (midpoint 11.5 Å between "deeper than 15 Å" and "shallower than
8 Å"), widens ADAR2's margin over the PH domains from 0.07 to 0.26 on the
calibration panel. On the deposited benchmark, which that choice never saw, it
made the score worse. The training workflow scores the same out-of-fold
pockets with both depth measures (run on commit 7012ff3):

| task | rule-based ROC-AUC, burial / hull depth | PR-AUC, burial / hull | DeLong hull − burial |
|---|---|---|---|
| ip_site (apo) | 0.934 / 0.890 | 0.496 / 0.389 | −0.044, p = 1×10⁻¹³ |
| ip_site (holo) | 0.951 / 0.910 | 0.634 / 0.507 | −0.041, p = 4×10⁻¹⁶ |
| cryptic_ip_site (apo) | 0.953 / 0.954 | 0.537 / 0.385 | +0.001, p = 0.89 |

A fixed ramp on hull depth penalises the many genuine sites that sit shallow in
a large protein, and a protein's size, not only its site's burial, sets how
deep a point can be. Burial depth therefore remains the default;
`ScoringParameters(depth_measure="hull")` selects the alternative. Hull depth is
kept where it earns its place: as a learned descriptor (with it, the best
cryptic-site model reaches ROC-AUC 0.991, 0.972-1.000, against 0.937 before -
32 positives in 13 structures, so this is suggestive rather than settled), and
as a gate in the screen's calibrated hit definition, which applies it as a
threshold rather than as a ramp inside the score.

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

**The gate scores the deposited holo structures, with the ligand present.**
That is not how a proteome target is scored, and it matters: the bound ligand
occludes the residues lining its own pocket. Scored as the screen scores - the
ligand's chain alone, every non-polymer atom removed - ADAR2's site scores 0.572
and the highest PH-domain site 0.491, a separation of **0.081**, not 0.581; on
the AlphaFold models it is 0.069 (section 6a). The gate's 0.581 measures the
pipeline on an input it never sees in a screen.

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

### 6a. The project plan's Phase 1 criteria, measured

`scripts/phase1_criteria.py` (CI: `.github/workflows/phase1-criteria.yml`)
measures every criterion of the plan's sections 5-6 on the full panel - the
plan's positives ADAR2, Pds5B, HDAC1 and HDAC3, and its four PH-domain
negatives PLCδ1, Btk, DAPP1 and Grp1 - on both the crystal structure and the
control's AlphaFold model. Each is evaluated as a proteome target: reduced to
the ligand's chain, every non-polymer atom removed, pockets detected and scored
as the screen scores them. The holo ligand only says which pocket is the site.
On the AlphaFold model the ligand is placed by superposing the binding region,
residues paired by sequence alignment rather than number.

| site (AlphaFold model) | rank | score | hull depth | depth | lining SASA | basic ≤5 Å (crystal) | APBS kT/e | RMSD to crystal |
|---|---|---|---|---|---|---|---|---|
| ADAR2 | 1/82 | 0.575 | 17.3 Å | 4.5 Å | 30.7 | 8 | +62 | 0.48 Å |
| Pds5B | 45/183 | 0.327 | 38.3 | 2.6 | 79.0 | 8 | +10 | 0.24 |
| HDAC1 | 2/52 | 0.510 | 17.8 | 1.0 | 48.7 | 4 | +40 | 0.44 |
| HDAC3 | 4/63 | 0.385 | 7.1 | 1.1 | 31.2 | 5 | +19 | 0.22 |
| PLCδ1 PH | 3/81 | 0.478 | 8.2 | 2.7 | 82.3 | 6 | +11 | 0.91 |
| Btk PH (1BWN) | 1/90 | 0.506 | 7.8 | 3.8 | 37.1 | 7 | +17 | 1.76 |
| DAPP1 PH | 2/37 | 0.446 | 5.7 | 3.2 | 61.4 | 6 | +18 | 0.32 |
| Grp1 PH | 2/53 | 0.484 | 7.9 | 3.0 | 55.5 | 7 | +21 | 0.23 |

1BTK, the identifier the plan gives for Btk, contains no inositol phosphate;
1BWN is the Btk PH domain with Ins(1,3,4,5)P4.

What holds:

- **ADAR2's site is the top-ranked pocket** in both the crystal structure and
  the AlphaFold model, and the model reproduces the binding region to 0.48 Å.
  Every control's model is within the plan's 2 Å.
- Basic residues (≥ 6 within 5 Å) and a positive potential are met at ADAR2.

What does not:

- **Rank within a protein does not separate buried from surface sites.** Every
  PH-domain site ranks 1-3 of its protein's pockets: in a small domain the
  inositol phosphate site is simply the main pocket. The plan's "bottom half"
  criterion fails for all four negatives.
- **Basic residues and potential do not separate them either.** They are
  properties of phosphate-binding sites, buried or not: the PH-domain sites
  carry 6-7 basic residues and +11 to +21 kT/e.
- **The composite score barely does.** 0.575 for ADAR2 against 0.446-0.506 for
  the PH domains, so no positive-negative overlap for ADAR2, but a 0.07 margin,
  and Pds5B and HDAC3 fall below every negative. None reaches the plan's 0.7.
- **Lining SASA fails on apo sites.** The plan's < 5 Å² target is a property of
  the complex (ADAR2's basic coordinating residues measure 10.6 Å² with the
  ligand present, 42.5 without). With the ligand removed they are exposed.

The separating measurement is **depth to the convex hull**: 13.0 Å (crystal)
and 17.3 Å (model) for ADAR2 against 5.7-8.4 Å for every PH-domain site. The
pipeline's burial depth - distance to the nearest solvent-exposed atom -
collapses on apo structures, because the walls of the empty cavity are
themselves exposed: ADAR2's fully enclosed site reads 3-4.5 Å. Hull depth is
size-dependent (a small domain cannot contain a deep point; Pds5B's elongated
HEAT-repeat model inflates it), which the screen's analysis has to account for.

**Consequence for the screen.** The plan's strict filter - score ≥ 0.75 and
lining SASA ≤ 10 Å² - rejects ADAR2's own site on two gates. The proteome
screen therefore reports the plan's definition and a calibrated one (score ≥
0.54, hull depth ≥ 10 Å, no SASA gate) side by side. The score threshold is the
midpoint between ADAR2 (0.575) and the highest PH-domain site (0.506). The
calibrated thresholds rest on one positive and four negatives, and a test fails
if a scorer change stops them separating the controls. The screen stores every
pocket's descriptors and rescores them at aggregation, so a change to the
scorer or the criteria is applied to a finished screen without re-screening.

**Three APBS defects were found doing this.** The DX reader never parsed real
APBS output (it matched a header line APBS does not write), would have read the
map with x and z transposed, and the fixed 60 Å fine grid did not contain large
proteins, so off-centre pockets returned the boundary value. The grid is now
sized to the molecule and focused on the site, and sampling outside it raises.

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
- Group bootstrap, measured calibration (`tests/test_statistical_calibration.py`,
  run with `-m slow`). These are simulations of a null ROC-AUC with families correlated
  within a group, 300 replicates each. The one-sided false-positive rate of "95 %
  interval above 0.5" is:
  - 2.3 % with 30 groups and 3.3 % with 80 groups, against a nominal 2.5 %;
  - 6–8 % with 5–9 groups.

  The docking plans' 5-group minimum therefore marks intervals from fewer than 5
  groups as not evidence, but it does not make an interval from 5–9 groups exact. A
  decision that rests on 5–9 groups is read with that caveat.
- Cross-fitted re-ranking (study F), in the same simulations:
  - No "improves" decision in 40 null data sets.
  - With a partial electrostatic signal (half of the native poses marked), "improves"
    in 10, 19 and 20 of 20 data sets at small, medium and large effects.
  - The out-of-fold gain never exceeds the in-sample optimum, so choosing the weight
    does not inflate it.

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
