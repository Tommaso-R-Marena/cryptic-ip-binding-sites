# The sampling ceiling: plan

Pre-registered on 2026-09-25, before any pose at an exhaustiveness above 32 was
generated for this study.

**What was known when this was written.** Studies A and F are complete and recorded
(`results/redocking/`, `results/rerank/`):

- Top-pose success at 2 Å is 0.110 [0.046, 0.193] (group estimand, 250 copies, 31 strict
  homology groups).
- 230 of 242 seed-1 failures are *scoring* failures by study A's crystal-pose test, yet
  64.7 % of the 750 seed runs in study F never produced a pose within 2 Å at all.
- The best-of-list rate in study F — any of up to 40 poses within 2 Å — is
  0.399 [0.286, 0.521]. That is the ceiling on any re-scoring of those lists.
- Electrostatic re-ranking moved top-pose success to 0.147 [0.090, 0.213], a difference
  of +0.054 [−0.001, 0.114]: specific to electrostatics (permutation p = 0.0099) but not
  bounded away from zero.

Every one of those numbers came from runs at **exhaustiveness 32**, the value fixed in
`docs/REDOCKING_PLAN.md`. IP6 has 12 rotatable bonds and a ligand of 36 heavy atoms; 32
is AutoDock Vina's default for far smaller ligands. So the ceiling may be a property of
the search budget rather than of the receptor, the scoring function or the pose problem.
Nothing in this plan changes studies A or F, whose decisions stand as recorded.

## Why it matters

The two readings have opposite consequences.

- If the ceiling rises with budget, then "docking cannot place inositol phosphates" is
  wrong as stated: it is "docking as conventionally configured cannot", and the fix is
  compute.
- If the ceiling does not move, the ceiling is real, the pose problem is not a search
  problem, and study F's partial gain is close to all that re-scoring can deliver.

## Data

- **Copies.** The 250 copies with an outcome in study F (`results/rerank/rerank.json`),
  in 31 strict homology groups. The 19 copies that failed receptor preparation in study
  A fail here for the same reason and stay excluded; this study inherits that limitation
  and does not repair it.
- **Receptor, ligand, box, starting pose, RMSD.** Exactly as in `docs/REDOCKING_PLAN.md`,
  through the same code. Ligand: the CCD template in its primary protonation state.
- **Pose lists.** Up to 40 poses within 10 kcal/mol, as in `docs/RERANK_PLAN.md`.
- **Per pose** the Vina score, the symmetric heavy-atom RMSD to the crystal copy with no
  superposition, and E_el from `cryptic_ip.rescoring.electrostatics`.

## Arms

| arm | exhaustiveness | seeds | copies |
|---|---|---|---|
| E32 | 32 | 1, 2, 3 | all 250 — **reused from study F's pose lists, not re-docked** |
| E128 | 128 | 1, 2, 3 | all 250 |
| E512 | 512 | 1 | a stratified subset of 80 |

The E32 arm is taken from the `rerank-arms-*` artifacts of rerank run 36055694750. Vina's
seeds and starting poses are identical by construction, so E32 needs no new compute; this
is what makes the study affordable.

**The E512 subset** is fixed here, before any E512 pose exists: within each burial class,
copies are ordered by `copy_key` and every k-th is taken so that the subset holds 80
copies with the class proportions of the 250 (about 60 surface, 16 semi-cryptic, 4
cryptic). The selection uses no outcome of any arm. Cryptic copies fall in 4 strict
groups, so the subset's cryptic stratum is not evidence, as in every other study here.

## Metrics

- **Sampling ceiling**: per copy, the fraction of its seed runs in which *any* pose of the
  list is within 2.0 Å. This is the quantity G1 tests.
- **Top-pose success**: per copy, the fraction of its seed runs whose top-ranked pose is
  within 2.0 Å.
- **Re-ranked top-pose success**: the same, ranking by Vina + w·E_el with **w = 0.1
  frozen**. That is the weight four of study F's five folds chose; it is a constant here,
  fitted on no data in this study, so no weight is selected in this study at all.
- **Estimands and intervals.** Per copy and per group (each strict homology group
  weighted equally), 2,000 group-bootstrap resamples of `homology_group_strict`, seed
  20260930, percentile 95 % intervals. Arms are paired on the same copies and resampled
  together.
- **Small strata.** Fewer than 5 strict groups: reported with the words "fewer than 5
  independent groups: this interval is not evidence".

## Decisions

Three primary questions, Holm-corrected together over the group estimand.

- **G1, the ceiling (primary).** Paired difference in the sampling ceiling, E128 − E32.
  - **budget-limited** if the 95 % lower bound > 0 and Holm p < 0.05: more search finds
    near-native poses that exhaustiveness 32 misses.
  - **search-saturated** if the 95 % upper bound < 0.05: quadrupling the budget buys less
    than 5 points of ceiling, so the ceiling is a property of the problem, not the budget.
  - **inconclusive** otherwise.
- **G2, end-to-end success.** Paired difference in top-pose success, E128 − E32. Labels
  **improves** (lower bound > 0, Holm p < 0.05), **no gain** (upper bound < 0.05), else
  **inconclusive**. A ceiling that rises without top-pose success rising is a scoring
  problem that a bigger search has made worse, and the report must say so.
- **G3, search plus electrostatics.** Paired difference between re-ranked top-pose success
  at E128 (w = 0.1 frozen) and Vina top-pose success at E32 — the best protocol these two
  studies can assemble against the one study A registered. Same labels as G2.

**Not evaluable** in place of any label if fewer than 5 strict groups contribute.

## Descriptive, not tested

- **G4, rescue rate.** Among the seed runs whose E32 list held no pose within 2 Å, the
  fraction whose E128 list does; and for the subset, whose E512 list does. Reported with a
  group-bootstrap interval, and marked as conditional on an E32 outcome, so it is a
  description of where the gain lands rather than evidence for it.
- **The E512 ladder.** On the 80-copy subset, ceiling and top-pose success at 32, 128 and
  512 with one seed, to show whether the curve is still climbing at 512 or has flattened.
- Strata: burial class, interface, metal, X-ray or cryo-EM.
- Wall-clock seconds per copy per arm, so the cost of any gain is on the record.

## What the result can and cannot support

- **budget-limited.** The pose-prediction failure reported in study A is, in part, a
  configuration choice. It would not show that the scoring function is adequate: G2
  decides that separately, and a ceiling that rises while top-pose success does not is
  evidence *against* the scoring function.
- **search-saturated.** The ceiling is real. Study F's gain is then close to the most that
  re-scoring can deliver, and the remaining failure is in the pose problem itself: the
  receptor is rigid, the crystal sidechain rotamers were fitted to the ligand present, and
  no amount of search in a rigid box will fix that.
- Neither outcome speaks to affinity. These are pose-recovery rates for a ligand of −9
  formal charge scored by functions without explicit electrostatics, which is what studies
  A and F already said.

## Execution and outputs

- **Workflow.** `.github/workflows/sampling.yml`: 30 docking shards, then a report job.
- **Code.** `scripts/sampling.py` (dock), `scripts/sampling_report.py` (G1–G4).
- **Results.** Printed between `BEGIN_SAMPLING_JSON` and `END_SAMPLING_JSON`, extracted
  into `results/sampling/`.
- Changes to this plan are dated amendments in separate files, written before the results
  they could affect are read.
