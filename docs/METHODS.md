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
| `cryptic` | ≤ 0.05 |
| `semi_cryptic` | ≤ 0.25 |
| `surface` | > 0.25 |
| `crystal_artifact` | fewer than 8 protein heavy-atom contacts within 4.5 Å |

The boundaries are anchored on the paradigm cases: the ADAR2 InsP6 is described
as completely encapsulated with only a narrow window to the exterior
(Macbeth et al., *Science* 309:1534, 2005), while surface signalling sites such
as the PLCδ1 PH domain leave roughly half the ligand exposed.

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

### Tier-2 panel

| Control | PDB | Notes |
|---|---|---|
| Pds5B | 5HDT | Often a surface/crystal artefact |
| HDAC1 | 5ICN | Semi-cryptic interface site |
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
