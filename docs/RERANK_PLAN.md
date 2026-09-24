# Electrostatic re-ranking of docked inositol phosphate poses: plan

Pre-registered on 2026-09-24, before any pose for this study was generated or scored.

**What was known when this was written:**

- The redocking census (run 36020459246) had been read: 272 copies selected, 269 in the
  primary set.
- The α-arrestin study (`results/arrestin/`) had found that Vina top-pose success on 24
  crystal arrestin–IP sites was 0.00.
- No outcome of the redocking benchmark itself (`docs/REDOCKING_PLAN.md`: R1–R4, the
  controls, any per-copy RMSD) had been read. That run was still docking.

This study changes nothing in the redocking plan. Its decisions stand as registered.
Changes to this plan are dated amendments in separate files, written before the results
they could affect are read.

## Why

AutoDock Vina's scoring function has no explicit electrostatics. Its hydrogen-bond
term is short-ranged and charge-blind. For a ligand carrying −5 to −9 charges, whose
crystal sites are lined with Lys, Arg, His and metals, that is the most obvious missing
physics.

The question is whether adding one screened-Coulomb term, with a single weight fitted
out-of-sample, picks the crystal-like pose from Vina's own pose list more often than
Vina's ranking does. The study separates two failure modes:

- **scoring failure:** a pose within 2 Å of the crystal was generated but not ranked first;
- **sampling failure:** no such pose was generated.

## Data

- **Copies.** The primary set of the redocking census: `redocking-census` artifact of
  run 36020459246, rows with `selected` and `primary_set` true.
- **Receptor and ligand.** Receptor preparation, the CCD ligand, its primary
  protonation state, the box and the starting pose all follow `docs/REDOCKING_PLAN.md`,
  through the same code (`scripts/redocking.py`'s `Context`).
- **Docking.** AutoDock Vina 1.2.7, exhaustiveness 32, seeds 1, 2 and 3.
  - The pose list is widened to at most **40 poses within 10 kcal/mol** of the best, so
    the re-ranker has room to act. Vina's own top pose does not depend on how many poses
    are written out.
  - Each seed's top-pose RMSD is compared with the redocking primary arm's, as a
    reproducibility check (descriptive).
- **Per pose, recorded:**
  - the Vina total score;
  - the symmetric heavy-atom RMSD to the crystal copy, with no superposition
    (`cryptic_ip.docking.rmsd`);
  - the electrostatic energy defined below.

## The electrostatic term (fixed here)

E_el = 332.0637 Σ_i Σ_j q_i q_j exp(−κ r_ij) / (ε(r_ij) r_ij), in kcal/mol:

- i runs over ligand atoms and j over receptor atoms, with r_ij ≤ 12 Å.
- Charges: the partial charges written in the PDBQT files. For the receptor these are
  PDB2PQR/AMBER charges at pH 7.4, with non-polar hydrogens merged. For the ligand they
  are Meeko's Gasteiger charges for the primary protonation state.
- ε(r) = 4r, a distance-dependent dielectric.
- κ = 0.127 Å⁻¹, the Debye screening of 150 mM monovalent salt at 298 K.
- r_ij is floored at 1.5 Å.

Nothing in this term is fitted. The only fitted quantity is the weight w below.

## Re-ranking and cross-fitting

- **Score.** S_w(pose) = Vina(pose) + w · E_el(pose), for w in the grid
  {0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1}. w = 0 is Vina.
- **Folds.** 5 folds of strict homology groups (`homology_group_strict`), assigned by a
  fixed random permutation of the groups (seed 20260929). All copies of a group fall in
  the same fold.
- **Fitting.** For each held-out fold, w is the grid value that maximises the
  group-equal mean top-pose success on the other four folds. Ties go to the smaller w.
  That w re-ranks every pose list in the held-out fold.
- **Outcome per copy.** The out-of-fold top-pose success: the fraction of the three
  seeds whose top pose under S_w is within 2.0 Å.

## Decisions

- **F1, primary (the only hypothesis test; Holm is trivial).** The paired difference in
  top-pose success, re-ranked minus Vina (w = 0), under the group estimand.
  - **Intervals.** 2,000 resamples of strict homology groups (seed 20260929). The whole
    cross-fitting procedure is repeated inside every resample: each drawn group keeps its
    fold, and training and evaluation are weighted by the draw counts. The interval
    therefore includes the uncertainty of choosing w.
  - **Improves** if the 95 % lower bound is > 0.
  - **Worsens** if the 95 % upper bound is < 0.
  - **No detectable difference** otherwise.
  - **Not evaluable** with fewer than 5 groups.
  - The per-copy estimand is reported beside the group estimand.
- **F2, failure decomposition (descriptive).** For Vina and for the re-ranked score, the
  share of seed runs whose top pose fails while a pose within 2 Å exists in the list
  (scoring failure), and the share with no such pose (sampling failure). The sampling
  ceiling (best-of-list success) is reported with its interval.
- **F3, permutation control.** In 20 permutations (seeds 1–20), E_el is shuffled among
  the poses of each seed run and the whole F1 point estimate is recomputed. A true
  signal should vanish: the permutation distribution is reported, with the fraction of
  permutations whose difference reaches the observed one.
- **F4, strata (descriptive).** F1's two arms by burial class (surface, semi-cryptic,
  cryptic) and for the classic-arrestin entries (1ZSH, 5TV1, 7F1W, 7F1X, 7JTB, 7JXA,
  7MOR). A stratum with fewer than 5 groups is marked as not evidence.
- **Also reported:**
  - the distribution of the fitted w across folds and resamples;
  - the Spearman correlation of E_el with RMSD within runs;
  - E_el and the Vina score of the minimised crystal pose.

## What the result can and cannot support

- **Improves.** A single cheap term fixes part of Vina's polyanion blind spot on this
  benchmark. It does not establish that the absolute scores are meaningful. It does not
  apply to poses the search never samples (F2).
- **No detectable difference / worsens.** Electrostatics of this form do not help. The
  failures are then either sampling failures or need a desolvation term that this plan
  does not include.

## Execution and outputs

- **Workflow.** `.github/workflows/rerank.yml`: 30 docking shards, then a report job.
- **Code.** `cryptic_ip/rescoring/` holds the electrostatic term and the cross-fitting.
  `scripts/rerank.py` holds the dock and report commands.
- **Results.** Printed between `BEGIN_RERANK_JSON` and `END_RERANK_JSON`, extracted into
  `results/rerank/`.
