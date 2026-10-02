# Summary of changes for run cc3ee04c-fbe2-4bba-9048-d0f8c277ceb2
Everything is proved in Lean and the project builds with no `sorry` and no custom axioms. Items 1 and 3 and the Kabsch part hold as stated. In item 2, the claim that `ρ·√k` is the best bound for the minimising `T` is wrong: the sharp bound is `ρ·√(k−1)`. A plain-language write-up is in `RESULTS.md`.

Setup: points are in `EuclideanSpace ℝ (Fin 3)`, and a rigid motion is `x ↦ Ax + b` with `AᵀA = 1` and `det A = 1` (`RigidMotion` in `RequestProject/Basic.lean`).

**1. Soundness of distance pruning (`Pruning.lean`): true, and `2ε` is tight.**
- `distance_pruning_sound` proves the bound.
- `distance_pruning_tight` shows `2ε` is reached for every `ε ≥ 0`: take `t₀ = t₁ = 0`, `q₀ = (−ε,0,0)`, `q₁ = (ε,0,0)` and `T = id`. Both residuals are `ε` and the discrepancy is exactly `2ε`.

**2. The RMSD criterion (`RMSD.lean`). Here `R` is the infimum of `RMSD(T)` over all rigid motions.**
- **`R ≤ ρ` does not bound the minimiser's residuals by `ρ`.** See the `k = 4` example below.
- **`√k` bound, any `T`:** `residual_le_rmsd_mul_sqrt_k` proves `max‖T(tᵢ) − qᵢ‖ ≤ RMSD(T)·√k`. This is sharp only for `T` that are not minimisers (`residual_bound_sqrt_k_tight`).
- **Correction for the minimiser:** if `T` is the minimiser and `R ≤ ρ`, then `max‖·‖ ≤ ρ·√(k−1)` (`minimizer_residual_le_of_minRMSD_le`). The reason is that the best translation makes the residual vectors sum to zero. `minimizer_residual_bound_tight` shows `√(k−1)` is attained for every `k ≥ 1`.
- **Smallest safe tolerance:** for `k ≥ 2` and `ρ ≥ 0`, pruning with tolerance `τ` keeps every correspondence with `R ≤ ρ` exactly when `τ ≥ ρ·√(2k)` (`smallest_safe_tolerance`). So `τ_min = ρ√(2k)`, which is smaller than the `2ρ√k` you get from item 1. For `k = 4`, `ρ = 4` this is `8√2 ≈ 11.31`.
- **Counterexample for `k = 4`, `ρ = 4`, `τ = 1.5` (`rmsd_counterexample_k4`):** take all `tᵢ = 0` and `q = ((−5,0,0), (5,0,0), 0, 0)`.
  - `R ≤ 4`: the identity gives `√12.5 ≈ 3.54` and is a minimiser.
  - The minimiser's residual at anchor 0 is `5 > ρ`.
  - The discrepancy between anchors 0 and 1 is `10`, so pruning discards this correspondence whether `τ = 1.5` is read as an absolute tolerance or as `1.5·ρ = 6`.

**3. Completeness under the Chebyshev criterion: true.** `pruning_complete` shows that any correspondence admitting a `T` with `max‖T(tᵢ) − qᵢ‖ ≤ ε` passes the `2ε` test.

**Kabsch (`Kabsch.lean`, `OrthLemmas.lean`, `SVD.lean`).**
- **Optimality:** `kabsch_optimal` takes an SVD `H = U·diag(s)·Vᵀ` of the cross-covariance (`U`, `V` orthogonal, `s₀ ≥ s₁ ≥ s₂ ≥ 0`). It proves that `A* = V·diag(1,1,det V·det U)·Uᵀ` has determinant `+1`, and that with `b* = q̄ − A*t̄` the pair minimises `Σ‖Atᵢ + b − qᵢ‖²` over all orthogonal `A` with `det A = +1` and all `b`.
- **Proof outline:** centring removes `b`, and the cost becomes a constant minus `2·tr(AH)`. The key step is `Σ Mᵢᵢsᵢ ≤ s₀ + s₁ + det(M)·s₂` for orthogonal `M`. The `det M = −1` case uses the fact that a 3×3 rotation has trace at least `−1`.
- **SVD existence:** I also proved that every real 3×3 matrix has such an SVD (`svd_exists`). So `kabsch_exists_optimal` states the result without assuming an SVD is given.

**`det = +1` is necessary (`Reflection.lean`, `reflection_necessary`).** Take the chiral tetrahedron `t = (0, e₀, e₁, e₂)` and its mirror image `q = (0, e₀, e₁, −e₂)`.
- The reflection `diag(1,1,−1)` matches them with cost 0.
- It is the only orthogonal matrix with cost 0, and its determinant is `−1`.
- Every rigid motion has strictly positive cost.

So minimising over all orthogonal matrices returns a reflection here, not a rigid motion.

All modules are imported by `RequestProject/Main.lean`, and the main results are listed in the Properties table.