# Machine-checked geometry for the template transplant (study N)

## What this is

A Lean 4 / Mathlib development of the geometric lemmas that `scripts/template_fit.py`
relies on, produced by **Harmonic Aristotle** from a specification written in this
project. It audits claims the code already assumed; it proves no new biology and
changes no result in studies A–N.

Toolchain: `leanprover/lean4:v4.28.0`, Mathlib pinned at `v4.28.0` (`lakefile.toml`).
Plain-language summary in `RESULTS.md`; the modules are under `RequestProject/`.

## Why it exists: it caught a false claim in our own code

`scripts/template_fit.py::correspondences` carried this comment:

> *"Pruning on distances is exact for a rigid map: a superposition cannot make two points
> match if the distance between them differs by more than the tolerance, so nothing
> admissible is lost and the enumeration stays small."*

**That is false as the code uses it**, and study N ran with it. The search prunes
correspondences whose pairwise distances disagree by more than `DISTANCE_TOLERANCE = 1.5 Å`,
but accepts on `RMSD_CEILING = 4.0 Å` — a root-*mean*-square over 4 anchors, which bounds
the average residual, not the worst one. A concrete, non-degenerate case: template anchors
`(0,0,0), (6,0,0), (0,7,0), (0,0,8)` against the same set with the fourth anchor moved 5 Å
along *z* gives Kabsch **RMSD 2.087 Å** — well inside the ceiling — while the worst pairwise
distance mismatch is **5.0 Å**, so pruning discards it.

Consequence: study N's search was **not exhaustive**. Reported fit scores are therefore
upper bounds on the true minimal anchor RMSD rather than the minimum the plan claims. It
does not rescue study N (still degenerate N3, still surface placements) and it applies to
both arms, so it is unlikely to flip N2's direction — but the claim was wrong and is
corrected here rather than left standing.

## What was proved

| result | statement | module |
|---|---|---|
| `distance_pruning_sound` | per-anchor residuals ≤ ε ⟹ every pairwise distance agrees within **2ε** | `Pruning.lean` |
| `distance_pruning_tight` | 2ε is attained, so it cannot be improved | `Pruning.lean` |
| `pruning_complete` | under the **Chebyshev** criterion (max residual ≤ ε), pruning at τ = 2ε discards nothing admissible | `Pruning.lean` |
| `residual_le_rmsd_mul_sqrt_k` | for **any** rigid motion, max residual ≤ RMSD·√k | `RMSD.lean` |
| `minimizer_residual_le_of_minRMSD_le` | for the **minimising** motion, max residual ≤ ρ·**√(k−1)** | `RMSD.lean` |
| `minimizer_residual_bound_tight` | √(k−1) is attained, so it is sharp | `RMSD.lean` |
| `smallest_safe_tolerance` | pruning at τ keeps every correspondence with R ≤ ρ **iff** τ ≥ ρ·**√(2k)** | `RMSD.lean` |
| `kabsch_optimal`, `kabsch_exists_optimal` | the Kabsch/Procrustes solution minimises Σ‖A tᵢ + b − qᵢ‖² over rotations, with SVD existence proved rather than assumed | `Kabsch.lean`, `SVD.lean` |
| `reflection_necessary` | without `det A = +1` the optimum can be a reflection: a chiral tetrahedron and its mirror image are matched at cost 0 only by `diag(1,1,−1)` | `Reflection.lean` |

## A correction to the specification we sent

The specification asserted that `ρ·√k` is the best bound on the minimiser's residuals.
**That was wrong.** The sharp bound is `ρ·√(k−1)`, because the optimal translation forces
the residual vectors to sum to zero, and that constraint costs one degree of freedom.
Likewise the smallest safe pruning tolerance is `ρ·√(2k)`, not the `2ρ·√k` obtained by
substituting a per-anchor bound into the 2ε lemma: two residuals cannot both sit at the
maximum at once.

Both corrections were checked independently of Lean, by Monte Carlo over 24,000 random
correspondences for k = 2…8 using this repository's own `kabsch_batch`: the residuals sum
to zero at the minimiser to 1.5 × 10⁻¹⁴, and the ratios `max‖r‖ / (R√(k−1))` and
`max discrepancy / (R√(2k))` both stay ≤ 1 and attain exactly 1.000000.

## What this means for the method

For our settings (k = 4 matched anchors, ρ = 4.0 Å ceiling):

- the minimiser's worst anchor can sit **ρ√3 ≈ 6.93 Å** out;
- sound pruning under the RMSD criterion needs **τ ≥ ρ√8 ≈ 11.31 Å**.

An 11.31 Å tolerance exceeds the span of the anchor constellations themselves, so it prunes
essentially nothing. **The RMSD acceptance criterion cannot be made sound by any useful
tolerance.** The fix is to accept on the **Chebyshev** criterion instead — every matched
anchor within ε — and prune at τ = 2ε, which `pruning_complete` certifies as discarding
nothing admissible. That is also the better scientific choice: a transplant with one anchor
5 Å out is not a good fit merely because the mean absorbs it.

Because that changes the method rather than fixing a bug in it, it belongs to the corrected
study N **pre-registration**, not to a silent edit of a study already run.

## Limits of this verification

- **The build was not reproduced here.** This container has no Lean toolchain and cannot
  fetch Mathlib, so "builds with no `sorry`" is not independently confirmed. What was
  checked: the sources contain no `sorry`, no `axiom` declaration, no `native_decide`, no
  `unsafe`, and no `implemented_by`; the theorem statements were read and say what the
  summary claims, with `smallest_safe_tolerance` stated as an `iff`; and the mathematics
  was re-derived numerically as above. Anyone with Lean 4.28 can run `lake build` here.
- These theorems concern the **search**, not the biology. They say the enumeration misses
  nothing under a stated criterion. They say nothing about whether a pocket binds IP6.

## Correction, 2026-10-01: the template library size

This note originally repeated the plan's figure of 131 copies in 77 entries and 14 strict
homology groups. That was wrong. The filter the code applies keeps copies whose
`symmetry_contact` flag is empty, which is the 91 cryo-EM copies (no unit cell, so the
symmetry-contact test is inapplicable rather than failed), so the library that ran was
**222 copies across 137 PDB entries in 17 strict homology groups**. See
`docs/TEMPLATE_PLAN_AMENDMENT_2.md`. No result changes: the run always used the larger
library, only its description was wrong.
