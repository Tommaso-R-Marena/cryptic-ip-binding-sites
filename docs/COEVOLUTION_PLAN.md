# Does buried IP-site frequency track IP concentration? A comparative test: plan

Pre-registered on 2026-09-26, before any per-organism hit rate at an absolute threshold
was computed.

## Why this test, and why it has not been run

The project's founding document names this as the headline comparative question:
*Dictyostelium discoideum* holds ~520 µM IP6 against ~15–50 µM in human and ~15–25 µM in
yeast, and it makes the pyrophosphates in quantity too (IP7 ~60 µM, IP8 ~180 µM, rising
25-fold in starvation). If protein architecture co-evolved with inositol-phosphate
availability, the organism swimming in IPs should bury more of them.

All three proteomes have been screened for a year of project time and **this comparison
has never been made**, because every downstream study took candidates as a **rank cutoff
per organism** (the top 25–27 of each), which fixes the count by construction and can
say nothing about rates. A rate needs one absolute threshold applied to all three
proteomes. That is what this study computes.

**What was read first.** The screen's outputs are extensively read: study C's ROC-AUCs,
the 75 candidates and their organisms, and that 34 annotated binders sit among 33,084
unseen proteins. **No per-organism hit rate at an absolute threshold has been computed or
seen.** The quantity this study tests is new.

## The confound, stated before the result

A raw hit-rate difference between these three proteomes is close to uninterpretable, and
the naive comparison is **not** this study's test. At least four things differ between
the organisms besides IP biology:

1. **AlphaFold confidence.** Model quality is not equal across proteomes, and a pocket
   found in a low-confidence region is noise.
2. **Amino-acid composition.** *Dictyostelium* is famously Asn/Gln-rich with long
   low-complexity runs; basic-residue content and pocket statistics follow composition.
3. **Protein length**, which sets how many pockets a protein has at all.
4. **Training-set bias.** The learned model was fitted on PDB-derived data dominated by
   human and other model organisms, so it may simply transfer better to human.

So the **primary test is the matched comparison (C2), not the raw one (C1)**, and the
learned and rule-only scores are both reported (C4) because only the learned one carries
confound 4.

## Data

- **Proteins.** Every scored protein in the three screened proteomes, from study C's
  `proteins_combined.csv.gz` — the same table studies H and I used, carrying
  `combined`, `learned_score`, `plddt_mean`, `hull_depth`, `organism_key`, `annotated`,
  `seen` and the MMseqs2 30 % `cluster`.
- **Grouping for every interval.** MMseqs2 clusters, resampled whole, within organism.
- **Composition.** The fraction of K, R and H in each protein's sequence, and its length,
  from the same UniProt sequences the other studies use.

## The threshold, fixed here

A hit rate needs a threshold that does not depend on the organisms being compared, and
must not be chosen after seeing the rates.

- **Anchor.** The threshold is the `combined` score at a **1 % false-positive rate among
  non-annotated, unseen proteins pooled across all three proteomes**. Pooling the anchor
  leaves the per-organism rates free to differ; it does not force them equal.
- **Sensitivity.** The whole curve is reported at anchors 0.1 %, 0.5 %, 1 %, 2 % and 5 %.
  The decision rests on the 1 % anchor alone; the curve is shown so the answer cannot be
  a property of one cut.
- **pLDDT floor.** A protein is eligible only if its top pocket's mean pLDDT is ≥ 70,
  the founding document's own filter. Proteins below it are excluded from every arm and
  the excluded count is reported per organism.

## Decisions

- **C1, the raw comparison (descriptive, explicitly not the test).** Hit rate per
  organism, and the difference *Dictyostelium* − mean(human, yeast), with cluster
  bootstrap intervals. Reported first precisely so the reader can see how much of it
  the matching removes.
- **C2, primary.** The same difference computed within **matched strata**: proteins are
  binned by mean pLDDT (width 5), length (log₂ width 0.5) and basic fraction (width
  0.02), and only bins containing proteins from *Dictyostelium* and from at least one
  comparator contribute; the difference is the bin-size-weighted mean of the
  within-bin differences, resampling clusters (2,000 resamples, seed 20261005).
  - **higher in Dictyostelium** if the 95 % lower bound > 0.
  - **lower in Dictyostelium** if the upper bound < 0.
  - **no material difference** if the interval lies inside ±0.002, i.e. 0.2 percentage
    points on a hit rate the founding document expects to be 0.2–0.9 %.
  - **inconclusive** otherwise. **Not evaluable** below 5 contributing clusters.
- **C3, where any difference lives (secondary).** The same difference within each pLDDT
  decile. A difference that exists only in the lowest deciles is a model-quality
  artefact and is reported as one.
- **C4, robustness to the learned model (secondary).** C2 recomputed with the rule-only
  score in place of `combined`. The learned model carries a training-distribution bias
  across organisms; the rule does not. **If C2 and C4 disagree in direction, the primary
  is reported as not robust** and no co-evolution claim is made.
- **C5, the pyrophosphate prediction (secondary).** *Dictyostelium*'s excess is largest
  in IP7 and IP8, not IP6. If co-evolution is real and driven by the pyrophosphates, the
  difference should be larger among **deeply buried** pockets (hull depth in the top
  quartile pooled across organisms), which are the ones big enough and enclosed enough to
  hold a bulkier ligand. Reported as C2 restricted to that quartile. This is a directional
  prediction stated before the data, not a post-hoc slice.
- **Holm** across C2 and C4, the two confirmatory tests. C1, C3 and C5 are not corrected.

## What this can and cannot support

- It can say whether the buried-IP-pocket rate differs between these proteomes once model
  confidence, length and composition are held fixed. That is the founding document's
  headline question, asked in the one form in which the answer means something.
- A null is a real answer and is reported as one. With three organisms there are **two
  independent contrasts and no replication**, so even a clean positive is a correlation
  across three points, not a demonstration of co-evolution. This limit is stated in the
  result, not only here.
- It cannot establish that any individual pocket binds an IP. Studies C and H already
  bound what can be claimed there, and nothing in this study changes it.
- It cannot separate IP concentration from everything else that differs between an
  amoebozoan, a fungus and an animal. Phylogeny is not controlled and cannot be with
  n = 3.

## Execution and outputs

- **Workflow.** `.github/workflows/coevolution.yml`: one analysis job, then a report.
- **Code.** `scripts/coevolution.py`.
- **Results.** Printed between `BEGIN_COEVOLUTION_JSON` and `END_COEVOLUTION_JSON`,
  extracted into `results/coevolution/`.
