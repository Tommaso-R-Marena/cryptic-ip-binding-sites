module

public import RequestProject.OrthLemmas

/-!
# The Kabsch / orthogonal-Procrustes solution

Given anchors `t q : Fin k → E3`, centre them at their centroids `t̄, q̄`, form the
cross-covariance `H = Σᵢ (tᵢ − t̄)(qᵢ − q̄)ᵀ`, take a singular value decomposition
`H = U · diag(s₀, s₁, s₂) · Vᵀ` (`U, V` orthogonal, `s₀ ≥ s₁ ≥ s₂ ≥ 0`), and put
`d = det V · det U`, `A* = V · diag(1, 1, d) · Uᵀ`, `b* = q̄ − A* t̄`.
Then `A*` is a rotation and `(A*, b*)` minimises `Σᵢ ‖A tᵢ + b − qᵢ‖²` over all rotations `A`
and translations `b`.
-/

@[expose] public section

open Matrix

noncomputable section

variable {k : ℕ}

/-- The centroid `(1/k) Σᵢ pᵢ`. -/
def centroid (p : Fin k → E3) : E3 := (k : ℝ)⁻¹ • ∑ i, p i

/-- The least-squares cost `Σᵢ ‖A tᵢ + b − qᵢ‖²`. -/
def lsCost (t q : Fin k → E3) (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3) : ℝ :=
  ∑ i, ‖mact A (t i) + b - q i‖ ^ 2

/-- The cross-covariance matrix `H = Σᵢ (tᵢ − t̄)(qᵢ − q̄)ᵀ`. -/
def crossCov (t q : Fin k → E3) : Matrix (Fin 3) (Fin 3) ℝ :=
  ∑ i, vecMulVec (WithLp.ofLp (t i - centroid t)) (WithLp.ofLp (q i - centroid q))

/-- The Kabsch rotation built from an SVD `H = U · diag(s) · Vᵀ`:
`A* = V · diag(1, 1, det V · det U) · Uᵀ`. -/
def kabschRot (U V : Matrix (Fin 3) (Fin 3) ℝ) : Matrix (Fin 3) (Fin 3) ℝ :=
  V * diagonal ![1, 1, V.det * U.det] * Uᵀ

/-- The Kabsch translation `b* = q̄ − A* t̄`. -/
def kabschTrans (t q : Fin k → E3) (U V : Matrix (Fin 3) (Fin 3) ℝ) : E3 :=
  centroid q - mact (kabschRot U V) (centroid t)

lemma sum_sub_centroid (p : Fin k → E3) : ∑ i, (p i - centroid p) = 0 := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_univ, Fintype.card_fin, centroid,
    ← Nat.cast_smul_eq_nsmul ℝ, smul_smul, mul_inv_cancel₀ (by exact_mod_cast hk.ne'), one_smul,
    sub_self]

lemma sum_norm_add_sq (u : Fin k → E3) (w : E3) :
    ∑ i, ‖u i + w‖ ^ 2 = ∑ i, ‖u i‖ ^ 2 + 2 * inner ℝ (∑ i, u i) w + k * ‖w‖ ^ 2 := by
  simp_rw [norm_add_sq_real]
  rw [Finset.sum_add_distrib, Finset.sum_add_distrib, sum_inner, Finset.mul_sum]
  simp

/-- Centering: `cost(A, b) = Σ ‖A(tᵢ − t̄) − (qᵢ − q̄)‖² + k ‖A t̄ + b − q̄‖²`. -/
lemma lsCost_center (t q : Fin k → E3) (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3) :
    lsCost t q A b = ∑ i, ‖mact A (t i - centroid t) - (q i - centroid q)‖ ^ 2 +
      k * ‖mact A (centroid t) + b - centroid q‖ ^ 2 := by
  have h : ∀ i, mact A (t i) + b - q i = (mact A (t i - centroid t) - (q i - centroid q)) +
      (mact A (centroid t) + b - centroid q) := by
    intro i; rw [mact_sub]; abel
  have hs : ∑ i, (mact A (t i - centroid t) - (q i - centroid q)) = 0 := by
    rw [Finset.sum_sub_distrib]
    simp only [mact]
    rw [← map_sum, sum_sub_centroid, sum_sub_centroid]
    simp
  simp only [lsCost, h]
  rw [sum_norm_add_sq, hs, inner_zero_left]
  ring

/-- `Σᵢ ⟪A xᵢ, yᵢ⟫ = tr(A · Σᵢ xᵢ yᵢᵀ)`. -/
lemma sum_inner_mact (A : Matrix (Fin 3) (Fin 3) ℝ) (x y : Fin k → E3) :
    ∑ i, inner ℝ (mact A (x i)) (y i) =
      (A * ∑ i, vecMulVec (WithLp.ofLp (x i)) (WithLp.ofLp (y i))).trace := by
  rw [Matrix.mul_sum, trace_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  simp [mact, inner, trace, mul_apply, vecMulVec_apply, mulVec, dotProduct, Fin.sum_univ_three]
  ring

/-- For orthogonal `A`, the centred cost equals `C − 2 tr(A H)` with `C` independent of `A`. -/
lemma centered_cost_eq (t q : Fin k → E3) (A : Matrix (Fin 3) (Fin 3) ℝ) (hA : IsOrth A) :
    ∑ i, ‖mact A (t i - centroid t) - (q i - centroid q)‖ ^ 2 =
      ∑ i, ‖t i - centroid t‖ ^ 2 + ∑ i, ‖q i - centroid q‖ ^ 2 -
        2 * (A * crossCov t q).trace := by
  have hpt : ∀ i, ‖mact A (t i - centroid t) - (q i - centroid q)‖ ^ 2 =
      ‖t i - centroid t‖ ^ 2 - 2 * inner ℝ (mact A (t i - centroid t)) (q i - centroid q) +
        ‖q i - centroid q‖ ^ 2 := by
    intro i; rw [norm_sub_sq_real, norm_mact_of_isOrth hA]
  simp only [hpt]
  rw [crossCov, ← sum_inner_mact, Finset.sum_add_distrib, Finset.sum_sub_distrib, Finset.mul_sum]
  ring

/-- `tr(M · diag s) = Σ Mᵢᵢ sᵢ`. -/
lemma trace_mul_diagonal (M : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ) :
    (M * diagonal s).trace = M 0 0 * s 0 + M 1 1 * s 1 + M 2 2 * s 2 := by
  simp [trace, mul_diagonal, Fin.sum_univ_three]

/-- The Kabsch rotation is a rotation (orthogonal with determinant `+1`). -/
theorem kabschRot_isRot (U V : Matrix (Fin 3) (Fin 3) ℝ) (hU : IsOrth U) (hV : IsOrth V) :
    IsRot (kabschRot U V) := by
  have hd : (V.det * U.det) ^ 2 = 1 := by rw [mul_pow, hV.det_sq, hU.det_sq, one_mul]
  have hD : IsOrth (diagonal ![1, 1, V.det * U.det]) := by
    unfold IsOrth
    rw [diagonal_transpose, diagonal_mul_diagonal, ← diagonal_one]
    congr 1
    ext i
    fin_cases i <;> simp [← sq, hd]
  refine ⟨(hV.mul hD).mul hU.transpose, ?_⟩
  rw [kabschRot, det_mul, det_mul, det_transpose, det_diagonal, Fin.prod_univ_three]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two, Matrix.head_cons,
    Matrix.tail_cons, one_mul]
  linear_combination hd

/-- For a rotation `A` and an SVD `H = U diag(s) Vᵀ`,
`tr(A H) ≤ s₀ + s₁ + det V det U · s₂ = tr(A* H)`. -/
lemma trace_le_kabsch (H U V : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ)
    (hU : IsOrth U) (hV : IsOrth V) (h01 : s 1 ≤ s 0) (h12 : s 2 ≤ s 1) (h2 : 0 ≤ s 2)
    (hH : H = U * diagonal s * Vᵀ) (A : Matrix (Fin 3) (Fin 3) ℝ) (hA : IsRot A) :
    (A * H).trace ≤ (kabschRot U V * H).trace := by
  have key : ∀ B : Matrix (Fin 3) (Fin 3) ℝ,
      (B * H).trace = ((Vᵀ * B * U) * diagonal s).trace := by
    intro B
    rw [hH, show B * (U * diagonal s * Vᵀ) = (B * U * diagonal s) * Vᵀ by
      simp only [Matrix.mul_assoc], trace_mul_comm]
    simp only [Matrix.mul_assoc]
  have hK : Vᵀ * kabschRot U V * U = diagonal ![1, 1, V.det * U.det] := by
    rw [kabschRot, show Vᵀ * (V * diagonal ![1, 1, V.det * U.det] * Uᵀ) * U =
      (Vᵀ * V) * diagonal ![1, 1, V.det * U.det] * (Uᵀ * U) by simp only [Matrix.mul_assoc],
      hV, hU, Matrix.one_mul, Matrix.mul_one]
  have hM : IsOrth (Vᵀ * A * U) := (hV.transpose.mul hA.1).mul hU
  have hdetM : (Vᵀ * A * U).det = V.det * U.det := by
    rw [det_mul, det_mul, det_transpose, hA.2, mul_one]
  rw [key, key, hK, trace_mul_diagonal, trace_mul_diagonal]
  have := hM.diag_weighted_le s h01 h12 h2
  rw [hdetM] at this
  simpa using this

/-- **Kabsch's theorem.** Let `H = U · diag(s) · Vᵀ` be a singular value decomposition of the
cross-covariance matrix (`U, V` orthogonal, `s₀ ≥ s₁ ≥ s₂ ≥ 0`). Then the Kabsch motion
`(A*, b*)` is a rigid motion (`A*` orthogonal with `det A* = +1`), and it minimises
`Σᵢ ‖A tᵢ + b − qᵢ‖²` over all orthogonal `A` with `det A = +1` and all `b`. -/
theorem kabsch_optimal (t q : Fin k → E3) (U V : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ)
    (hU : IsOrth U) (hV : IsOrth V) (h01 : s 1 ≤ s 0) (h12 : s 2 ≤ s 1) (h2 : 0 ≤ s 2)
    (hH : crossCov t q = U * diagonal s * Vᵀ) :
    IsRot (kabschRot U V) ∧
      ∀ (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3), IsRot A →
        lsCost t q (kabschRot U V) (kabschTrans t q U V) ≤ lsCost t q A b := by
  refine ⟨kabschRot_isRot U V hU hV, fun A b hA => ?_⟩
  have hK := kabschRot_isRot U V hU hV
  rw [lsCost_center, lsCost_center, centered_cost_eq t q _ hK.1, centered_cost_eq t q _ hA.1]
  have htr := trace_le_kabsch _ U V s hU hV h01 h12 h2 hH A hA
  have h0 : mact (kabschRot U V) (centroid t) + kabschTrans t q U V - centroid q = 0 := by
    rw [kabschTrans]; abel
  rw [h0, norm_zero]
  nlinarith [sq_nonneg ‖mact A (centroid t) + b - centroid q‖,
    (Nat.cast_nonneg k : (0:ℝ) ≤ k)]

/-- The same statement phrased with the `RigidMotion` structure: the Kabsch motion is a rigid
motion whose least-squares cost is no larger than that of any rigid motion. -/
theorem kabsch_optimal_rigid (t q : Fin k → E3) (U V : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ)
    (hU : IsOrth U) (hV : IsOrth V) (h01 : s 1 ≤ s 0) (h12 : s 2 ≤ s 1) (h2 : 0 ≤ s 2)
    (hH : crossCov t q = U * diagonal s * Vᵀ) :
    ∃ Tstar : RigidMotion, Tstar.A = kabschRot U V ∧ Tstar.b = kabschTrans t q U V ∧
      ∀ T : RigidMotion, ∑ i, ‖Tstar (t i) - q i‖ ^ 2 ≤ ∑ i, ‖T (t i) - q i‖ ^ 2 :=
  ⟨⟨kabschRot U V, kabschTrans t q U V, kabschRot_isRot U V hU hV⟩, rfl, rfl,
    fun T => (kabsch_optimal t q U V s hU hV h01 h12 h2 hH).2 T.A T.b T.isRot⟩

end
