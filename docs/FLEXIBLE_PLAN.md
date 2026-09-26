# The rigid receptor: plan

Pre-registered on 2026-09-26, before any pose from a flexible-receptor run exists.

**What was known when this was written.** Studies A, F and G are complete and recorded
(`results/redocking/`, `results/rerank/`, `results/sampling/`):

- Top-pose success at 2 Å is 0.110 [0.046, 0.193] (group estimand, 250 copies, 31 strict
  homology groups).
- 230 of 242 seed-1 failures are *scoring* failures by study A's crystal-pose test: the
  crystal pose scores worse than the pose Vina ranked first.
- Electrostatic re-ranking moves top-pose success by +0.054 [−0.001, 0.114] in study F
  and replicates at +0.041 [0.004, 0.083] in study G. Specific by permutation
  (p = 0.0099) but small.
- **The search is saturated.** Study G's E128 arm moved success by +0.019 with a 95 %
  upper bound of 0.047, below the 0.05 margin, and the 16× E512 arm was flat. More
  sampling is a closed question.

One structural explanation survives all three studies and has never been tested here:
**the receptor is rigid**. Every run redocks into a crystal structure whose side chains
were fitted *with the ligand present*. That should make redocking easier, not harder, so
a rigid receptor cannot by itself explain a 0.110 success rate — unless the scoring
function is exploiting the pre-formed pocket to place a *wrong* pose that fits the fixed
rotamers better than the true one does. Study A's best-of-list ceiling rises sharply with
burial (surface 0.279, semi-cryptic 0.563, cryptic 0.810), which is the signature of a
pocket whose shape, not the ligand's chemistry, is doing the discriminating.

## The question

Does allowing the pocket's side chains to move change top-pose success — and in which
direction? Both directions are informative and both are reportable:

- **Better** would mean the rigid receptor was the binding constraint, and the correct
  configuration for IP ligands is a flexible one.
- **Worse** would mean the pre-formed pocket was carrying the result, and the apparent
  0.110 is optimistic for any prospective use, where no ligand-fitted rotamers exist.
- **No change** would close the last structural explanation and leave the scoring
  function as the sole remaining suspect, which is where studies F and G already point.

The third outcome is the one I expect. It is pre-registered as a full result, not a
failure.

## Data

- **Copies.** The same 250 copies with an outcome in study F, in 31 strict homology
  groups. The 19 copies that failed receptor preparation in study A stay excluded.
- **Receptor, ligand, box, starting pose, RMSD, protonation.** Exactly as in
  `docs/REDOCKING_PLAN.md`, through the same code, at exhaustiveness 32.
- **Flexible residues.** Every receptor residue with a side-chain heavy atom within 4.0 Å
  of any ligand heavy atom in the crystal copy, excluding GLY, ALA and PRO (no rotatable
  side chain) and excluding CYS in a disulphide. Capped at the 8 residues nearest the
  ligand centroid, ordered by minimum heavy-atom distance, because Vina's flexible
  sampling degrades beyond roughly that many. The cap and the 4.0 Å radius are fixed
  here.
- **Pose lists.** Up to 40 poses within 10 kcal/mol, as in `docs/RERANK_PLAN.md`, so the
  arms are directly comparable to studies F and G.

## Arms

| arm | receptor | seeds | copies |
|---|---|---|---|
| rigid | as studies A/F/G | 1, 2, 3 | all 250 — **reused, not re-docked** |
| flex | 4 Å side chains movable, ≤ 8 | 1, 2, 3 | all 250 |

The rigid arm is taken from the `rerank-arms-*` artifacts of rerank run 36055694750, as
study G did. This is what makes the study affordable.

## Decisions

- **J1, primary.** The paired difference in top-pose success at 2 Å, flex − rigid, over
  copies, resampled by `homology_group_strict` (2,000 resamples, seed 20261003).
  - **better** if the 95 % lower bound > 0.
  - **worse** if the 95 % upper bound < 0.
  - **no material change** if the interval lies inside ±0.05, the same margin study G
    used.
  - **inconclusive** otherwise. **Not evaluable** below 5 groups.
- **J2, the ceiling (secondary).** The best-of-list rate in each arm, and their paired
  difference. A flexible receptor that raises the *ceiling* while leaving the top pose
  alone is a scoring result, not a sampling one, and is reported as such.
- **J3, burial (secondary, pre-specified stratification).** J1 within each of study A's
  three burial classes. A stratum with fewer than 5 groups yields no evidence and is
  reported as not evaluable rather than as a null.
- **Holm** across J1 and J2. J3 is exploratory and is not corrected.

## What this can and cannot support

- It can say whether the rigid receptor is the remaining cause of study A's failure rate.
- It cannot say anything about *cryptic* pocket opening, which is a backbone motion; this
  moves side chains only. A null here does not rule out that a cryptic site needs a
  different receptor conformation entirely, only that rotamer freedom is not the answer.
- It cannot change studies A, F or G, whose decisions stand as recorded.

## Execution and outputs

- **Workflow.** `.github/workflows/flexible.yml`: a sharded docking matrix, then a report.
- **Code.** `scripts/flexible.py`, mirroring `scripts/rerank.py`'s docking path with a
  test pinning the two to identical output on a rigid receptor.
- **Results.** Printed between `BEGIN_FLEXIBLE_JSON` and `END_FLEXIBLE_JSON`, extracted
  into `results/flexible/`.
