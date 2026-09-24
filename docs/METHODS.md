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
  excursion from a leak with 30 further shuffles (run 35965940997; result pending).

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
