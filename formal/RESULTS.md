# Rigid-motion pruning, RMSD, and Kabsch: results

Everything below is proved in Lean 4 / Mathlib. There is no `sorry` and no custom axiom.
Points are in `E3 = EuclideanSpace ℝ (Fin 3)`. A rigid motion is `x ↦ A x + b` with `Aᵀ A = 1` and
`det A = 1` (`RigidMotion` in `RequestProject/Basic.lean`).

## 1. Soundness of distance pruning (`RequestProject/Pruning.lean`)
* `distance_pruning_sound`: if `‖T(tᵢ) − qᵢ‖ ≤ ε` for all `i`, then
  `| ‖tᵢ − tⱼ‖ − ‖qᵢ − qⱼ‖ | ≤ 2ε` for all `i, j`. **True.**
* `distance_pruning_tight`: the bound `2ε` is **tight**. For every `ε ≥ 0` take `t₀ = t₁ = 0`,
  `q₀ = (−ε,0,0)`, `q₁ = (ε,0,0)` and `T = id`: both residuals equal `ε` and the discrepancy is exactly `2ε`.

## 2. An RMSD criterion does not give that pruning (`RequestProject/RMSD.lean`)
`RMSD(T) = sqrt((1/k) Σ ‖T(tᵢ) − qᵢ‖²)` and `R = inf_T RMSD(T)` (`minRMSD`).
* `rmsd_counterexample_k4`: `R ≤ ρ` does **not** imply `maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ ρ` for the minimising `T`
  (example given below).
* `residual_le_rmsd_mul_sqrt_k`: for **any** `T`, `maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ RMSD(T)·√k`.
  `residual_bound_sqrt_k_tight`: this bound is attained by a `T` that is *not* the minimiser.
* **Correction to the claim "the best general bound is ρ√k".** For the *minimising* `T`
  the best bound is `ρ·√(k−1)`, not `ρ·√k`:
  * `minimizer_residual_le_of_minRMSD_le`: if `T` minimises the RMSD and `R ≤ ρ`, then
    `maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ ρ√(k−1)`. Reason: the best translation forces `Σᵢ (T(tᵢ) − qᵢ) = 0`.
  * `minimizer_residual_bound_tight`: the factor `√(k−1)` is attained for every `k ≥ 1`
    (`tᵢ = 0`, `q₀ = (−(k−1),0,0)`, `qⱼ = (1,0,0)` for `j ≠ 0`).
* **Smallest safe pruning tolerance** (`smallest_safe_tolerance`): for `k ≥ 2` and `ρ ≥ 0`,
  pruning with tolerance `τ` keeps every correspondence with `R ≤ ρ` **iff** `τ ≥ ρ·√(2k)`.
  So `τ_min = ρ√(2k)`, smaller than the `2ρ√k` you get by plugging `ε = ρ√k` into item 1.
  (`pruning_sound_of_minRMSD_le` proves soundness and `pruning_tolerance_attained` proves sharpness.)
  For `k = 4`, `ρ = 4` this gives `τ_min = 8√2 ≈ 11.31`.
* **Concrete counterexample, `k = 4`, `ρ = 4`, `τ = 1.5`** (`rmsd_counterexample_k4`):
  `tᵢ = 0` for all four anchors, `q = ((−5,0,0), (5,0,0), (0,0,0), (0,0,0))`.
  * `R ≤ 4`. The identity gives `RMSD = √(50/4) = √12.5 ≈ 3.54`, and the identity is an RMSD minimiser.
  * The minimiser's residual at anchor 0 is `5 > 4 = ρ`.
  * `| ‖t₀ − t₁‖ − ‖q₀ − q₁‖ | = |0 − 10| = 10`. That exceeds `1.5`, and it also exceeds `1.5·ρ = 6`.
    So pruning with `τ = 1.5` (read either as an absolute tolerance or as a multiple of `ρ`) discards
    a correspondence with `R ≤ ρ`.

## 3. Completeness under the Chebyshev criterion (`RequestProject/Pruning.lean`)
* `pruning_complete`: if `∃ T, ∀ i, ‖T(tᵢ) − qᵢ‖ ≤ ε` (`Admissible`), then the correspondence passes
  the test `| ‖tᵢ − tⱼ‖ − ‖qᵢ − qⱼ‖ | ≤ 2ε` for all `i, j` (`PassesPruning (2ε)`). **True.**

## Kabsch / orthogonal Procrustes (`RequestProject/Kabsch.lean`, `SVD.lean`, `Reflection.lean`)
Let `H = Σ (tᵢ − t̄)(qᵢ − q̄)ᵀ = U·diag(s₀,s₁,s₂)·Vᵀ` with `U, V` orthogonal and `s₀ ≥ s₁ ≥ s₂ ≥ 0`.
Set `A* = V·diag(1,1,det V·det U)·Uᵀ` and `b* = q̄ − A* t̄`.
* `kabsch_optimal`: `A*` is orthogonal with `det A* = +1`, and `Σ‖A* tᵢ + b* − qᵢ‖² ≤ Σ‖A tᵢ + b − qᵢ‖²`
  for every orthogonal `A` with `det A = +1` and every `b`. `kabsch_optimal_rigid` states the same
  thing using the `RigidMotion` structure.
  * Proof outline: centring removes `b`. The cost then equals `C − 2·tr(A H)`. With `M = Vᵀ A U`,
    `tr(A H) = Σ Mᵢᵢ sᵢ ≤ s₀ + s₁ + det(M)·s₂`. The `det M = −1` case uses the fact that a 3×3 rotation
    has trace `≥ −1`, via the identity `(1+tr N)(3−tr N) = Σ_{i<j}(Nᵢⱼ − Nⱼᵢ)²` for `N ∈ SO(3)`.
* `svd_exists` / `kabsch_exists_optimal`: every real 3×3 matrix has such an SVD. So the Kabsch
  motion exists for every correspondence and is optimal.
* `reflection_necessary`: the restriction to `det = +1` is necessary. Take the chiral tetrahedron
  `t = (0, e₀, e₁, e₂)` and its mirror image `q = (0, e₀, e₁, −e₂)`:
  * the reflection `P = diag(1,1,−1)` has cost `0`;
  * every orthogonal `A` with cost `0` equals `P`, which has `det = −1`;
  * every rigid motion has strictly positive cost.

  So minimising over all orthogonal matrices returns a reflection, which is not a rigid motion.
