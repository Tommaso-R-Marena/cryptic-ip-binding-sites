# Pre-registered analysis plan: the expanded inositol phosphate benchmark

This plan is committed **before** any result on the expanded benchmark exists.
Every threshold, grouping rule, test and decision rule below is fixed here; the
analysis code implements it, and the report states each deviation, if any, with
its reason. Results obtained on the earlier 136-entry set motivated these
questions and are not part of the evidence for them.

## 1. Why a new analysis

The 136-entry benchmark answered some questions and left two open:

- whether depth to the convex hull, as a **learned descriptor**, improves the
  identification of buried (cryptic) inositol phosphate sites - the earlier
  estimate (ROC-AUC 0.937 → 0.991) rested on 32 positive pockets in 13
  structures and compared two runs that differed in more than that descriptor;
- whether the pipeline can tell a **buried** inositol phosphate site from a
  **surface** one, which is the distinction a proteome screen depends on and
  which no task measured.

An audit of the evaluation found leaks that inflate every earlier number. They
are corrected here, and each correction is enforced by an assertion or a test:

| defect | effect | correction |
|---|---|---|
| cross-validation groups fell back to the PDB entry because the grouping read a column the dataset did not have | the same protein (ADAR2: 12 entries) on both sides of a fold | homology groups (section 3), and no fallback |
| multi-protein entries formed their own group key | a protein alone and in complex in different groups | groups are connected components over shared homology |
| homologous proteins (ADAR1/ADAR2, HDAC1/2/3, PH domains) in different groups | homologue memorised across folds | sequence and structure homology grouping |
| `plddt_*` descriptors are crystallographic B-factors on PDB entries | residues ordered by the removed ligand have low B-factors: the ligand leaks back; the same column means pLDDT on AlphaFold models | excluded from every model |
| decision threshold chosen on the pooled out-of-fold predictions it is then scored on | optimistic MCC / F1 | chosen inside each training fold |
| best of five model families chosen by, and reported with, the same outer-fold score | winner's-curse selection bias | family chosen inside the inner loop, jointly with hyperparameters |
| DeLong tests treat pockets as independent | pockets cluster within proteins: p-values too small | paired bootstrap over homology groups |

## 2. Data

**Inclusion.** Every PDB entry containing at least one phosphorylated inositol
(non-lipid) component, discovered from the chemical component dictionary as
`scripts/build_ip_validation_dataset.py` does, from X-ray diffraction or
electron microscopy at a resolution of **3.5 Å or better**. NMR entries are
excluded (no resolution; ensemble coordinates).

**Exclusions, each counted and reported:** entries whose coordinates cannot be
fetched or parsed; entries with more than **60,000 protein heavy atoms**
(pocket detection on whole large assemblies is not comparable with the rest);
entries in which fpocket finds no pocket.

**Labels** (unchanged from `cryptic_ip.analysis.labeling`): a pocket is
positive when at least 30 % of an inositol phosphate copy's heavy atoms lie
within 4.0 Å of its alpha spheres, negative at 5 % or less, and ambiguous in
between; ambiguous pockets are excluded. Each copy's burial class comes from
its own relative SASA (cryptic ≤ 0.12, semi-cryptic ≤ 0.25, surface > 0.25).
A copy with fewer than 8 protein heavy-atom contacts is a crystal artefact: a
pocket on it is excluded, being neither a binding site nor evidence of absence.
A copy whose burial could not be measured labels a pocket positive for
`ip_site` only. Pockets are detected and described once; every task's labels are
derived from that one record (the overlap and the burial class of the copy each
pocket touches), so the tasks cannot drift apart. Labels are computed on the deposited structure; **descriptors
on the ligand-free copy** (every non-polymer atom removed).

**Descriptors:** the 43 of `FEATURE_NAMES`, minus the three `plddt_*`
(B-factors here) and `electrostatic_potential` (APBS is not run on the
benchmark), leaving **39**. A test fails if an excluded descriptor reaches a
model.

## 3. Grouping: no protein, homologue or shared fold on both sides of a split

Every protein chain of every entry takes part, not only the chains touching
the ligand: a negative pocket memorised in one fold is as much a leak as a
positive one.

- **Primary grouping (sequence).** All-against-all MMseqs2 search of every
  protein chain. Two entries are linked when any pair of their chains aligns at
  **≥ 30 % identity over ≥ 50 % of the shorter chain, E ≤ 10⁻³**. Groups are
  the connected components of that graph.
- **Strict grouping (structure), a sensitivity analysis.** The primary links,
  plus links from an all-against-all Foldseek search: any chain pair with
  **E ≤ 10⁻³ and alignment TM-score ≥ 0.5 over ≥ 50 % of the shorter chain**.
  This joins homologues too remote for sequence search - the PH domains, for
  example - into one group.

Every split, every bootstrap resample and the holdout respect the groups. An
assertion checks, for every split, that no group appears on both sides. The
size distribution of the groups is reported; if one component holds more than
40 % of a task's positive pockets under a grouping, that grouping cannot support
a five-fold evaluation and is reported as such rather than forced.

## 4. Tasks

| task | positives | negatives | excluded |
|---|---|---|---|
| **ip_site** | pockets on any inositol phosphate copy | pockets touching none | ambiguous overlap |
| **cryptic_ip_site** | pockets on a cryptic copy | pockets touching none | pockets on semi-cryptic or surface copies; ambiguous overlap |
| **burial** (new) | pockets on a cryptic copy | pockets on a surface copy | everything else |

`burial` is the screen's question - buried or surface, among genuine inositol
phosphate sites - and the one the earlier tasks never asked.

## 5. Evaluation protocol

- **Holdout (temporal, group-disjoint).** Groups of the **strict** grouping
  are ordered by the release date of their earliest entry; the latest **20 %**
  of groups form the holdout, so no holdout protein has a sequence or structural
  homologue in the development set. It is chosen from dates and groups alone,
  before any model is fitted, and is used exactly once, by the locked model.
- **Development: nested, grouped, stratified cross-validation**, 5 outer folds
  × 3 inner folds, **repeated 3 times** with different fold seeds under the
  primary grouping and once under the strict grouping. The inner loop chooses
  the model family (L2 logistic regression, extra trees, histogram gradient
  boosting - a linear floor, a bagged and a boosted ensemble) **jointly** with
  its hyperparameters (10 random draws per family), by average precision. The
  chosen configuration's inner out-of-fold predictions then fix, for that outer
  fold, the MCC-optimal decision threshold and a Platt calibrator; neither ever
  sees the outer test fold. Discrimination is measured on the uncalibrated
  scores, calibration (Brier, ECE) on the calibrated ones.
- **Rule-based baseline:** `PocketScorer` with default parameters; no fitting;
  its threshold chosen on the training folds in the same way.
- **Locked model:** the protocol above refitted on the whole development set,
  then scored once on the holdout.
- **Permutation control:** each task is run once more with its labels permuted
  at random across pockets, through the identical code path; its ROC-AUC 95 %
  interval must contain 0.5. If it does not, labels reach the model by some
  route other than the descriptors, and no result of that task is reported as
  evidence.

## 6. Hypotheses and tests

Primary (Holm-corrected across the two):

- **H1.** Adding `hull_depth` to the descriptor set improves identification of
  cryptic sites (`cryptic_ip_site`).
- **H2.** Adding `hull_depth` improves buried-versus-surface discrimination
  (`burial`).

Each is tested as a **paired ablation**: the full protocol run twice on
identical folds, seeds and hyperparameter draws, with and without
`hull_depth` (39 and 38 descriptors). The effect is the difference in out-of-fold ROC-AUC (secondary:
average precision), averaged over the three repeats; its 95 % interval and
two-sided p-value come from **2,000 bootstrap resamples of homology groups**,
applied to both arms together.

Secondary, reported with intervals and not corrected: learned model versus
rule-based score on each task; `burial` discrimination itself (ROC-AUC against
0.5); every result under the strict grouping; every result on the holdout.

## 7. What counts as definitive

For H1 and H2 separately:

- **Supported** only if the development interval lies entirely above 0 under
  the primary grouping **and** under the strict grouping, the Holm-adjusted
  p-value is below 0.05, **and**, where the holdout holds at least 10 positive
  groups, the holdout difference is positive with its interval above 0.
- **Refuted** if the development interval lies entirely below 0, or entirely
  within ±0.01 ROC-AUC (no effect of any practical size), under both
  groupings.
- **Inconclusive** otherwise - reported as such, with the interval.

A holdout with fewer than 10 positive groups for a task is reported and marked
underpowered; it can then neither support nor refute.

## 8. Reporting

All numbers above are produced by CI from the committed code; the report links
the run. Any deviation from this plan is listed with its reason in the report,
next to the result it affects.
