module

public import RequestProject.Basic

/-!
# Elementary facts about `3 × 3` orthogonal matrices

The key inequality used in the proof of Kabsch's theorem: for an orthogonal `M` and
`s₀ ≥ s₁ ≥ s₂ ≥ 0`, `Σ Mᵢᵢ sᵢ ≤ s₀ + s₁ + det(M) s₂`.
-/

@[expose] public section

open Matrix

noncomputable section

lemma IsOrth.mul_transpose {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) : M * Mᵀ = 1 :=
  mul_eq_one_comm.mp hM

lemma IsOrth.det_sq {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) : M.det ^ 2 = 1 := by
  have := congrArg Matrix.det hM
  rw [det_mul, det_transpose, det_one] at this
  rw [sq]; exact this

lemma IsOrth.det_eq {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) :
    M.det = 1 ∨ M.det = -1 := by
  have h := hM.det_sq
  have : (M.det - 1) * (M.det + 1) = 0 := by linear_combination h
  rcases mul_eq_zero.mp this with h1 | h1
  · left; linarith
  · right; linarith

lemma IsOrth.mul {M N : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) (hN : IsOrth N) :
    IsOrth (M * N) := by
  unfold IsOrth at *
  rw [transpose_mul, Matrix.mul_assoc, ← Matrix.mul_assoc Mᵀ, hM, Matrix.one_mul, hN]

lemma IsOrth.transpose {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) : IsOrth Mᵀ := by
  unfold IsOrth; rw [transpose_transpose]; exact hM.mul_transpose

lemma IsOrth.neg {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) : IsOrth (-M) := by
  unfold IsOrth at *; simp [hM]

/-- Diagonal entries of an orthogonal matrix are at most `1`. -/
lemma IsOrth.diag_le_one {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) (i : Fin 3) :
    M i i ≤ 1 := by
  have h := congrFun (congrFun hM i) i
  simp only [mul_apply, transpose_apply, one_apply_eq, Fin.sum_univ_three] at h
  fin_cases i <;> simp at h ⊢ <;> nlinarith [sq_nonneg (M 0 0), sq_nonneg (M 1 0),
    sq_nonneg (M 2 0), sq_nonneg (M 0 1), sq_nonneg (M 1 1), sq_nonneg (M 2 1),
    sq_nonneg (M 0 2), sq_nonneg (M 1 2), sq_nonneg (M 2 2)]

/-- A rotation `N` of `ℝ³` has trace at least `-1`. -/
lemma IsRot.trace_ge {N : Matrix (Fin 3) (Fin 3) ℝ} (hN : IsRot N) :
    -1 ≤ N 0 0 + N 1 1 + N 2 2 := by
  obtain ⟨hO, hdet⟩ := hN
  have hadj : adjugate N = Nᵀ := by
    calc adjugate N = adjugate N * (N * Nᵀ) := by rw [hO.mul_transpose, Matrix.mul_one]
      _ = (adjugate N * N) * Nᵀ := by rw [Matrix.mul_assoc]
      _ = Nᵀ := by rw [adjugate_mul, hdet, one_smul, Matrix.one_mul]
  have c00 := congrFun (congrFun hadj 0) 0
  have c11 := congrFun (congrFun hadj 1) 1
  have c22 := congrFun (congrFun hadj 2) 2
  simp [adjugate_fin_three] at c00 c11 c22
  have hF0 := congrFun (congrFun hO 0) 0
  have hF1 := congrFun (congrFun hO 1) 1
  have hF2 := congrFun (congrFun hO 2) 2
  simp [mul_apply, Fin.sum_univ_three] at hF0 hF1 hF2
  have hd0 := hO.diag_le_one 0
  have hd1 := hO.diag_le_one 1
  have hd2 := hO.diag_le_one 2
  have key : (1 + (N 0 0 + N 1 1 + N 2 2)) * (3 - (N 0 0 + N 1 1 + N 2 2)) =
      (N 0 1 - N 1 0) ^ 2 + (N 0 2 - N 2 0) ^ 2 + (N 1 2 - N 2 1) ^ 2 := by
    linear_combination (-1 : ℝ) * (hF0 + hF1 + hF2) - 2 * (c00 + c11 + c22)
  by_contra hlt
  push_neg at hlt
  have : (1 + (N 0 0 + N 1 1 + N 2 2)) * (3 - (N 0 0 + N 1 1 + N 2 2)) < 0 :=
    mul_neg_of_neg_of_pos (by linarith) (by linarith)
  nlinarith [sq_nonneg (N 0 1 - N 1 0), sq_nonneg (N 0 2 - N 2 0), sq_nonneg (N 1 2 - N 2 1)]

/-- An orthogonal `3 × 3` matrix with determinant `-1` has trace at most `1`. -/
lemma IsOrth.trace_le_of_det_neg {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M)
    (hdet : M.det = -1) : M 0 0 + M 1 1 + M 2 2 ≤ 1 := by
  have hN : IsRot (-M) := ⟨hM.neg, by rw [det_neg]; simp [hdet]; norm_num⟩
  have := hN.trace_ge
  simp at this
  linarith

/-- **Key inequality.** For orthogonal `M` and `s₀ ≥ s₁ ≥ s₂ ≥ 0`,
`M₀₀ s₀ + M₁₁ s₁ + M₂₂ s₂ ≤ s₀ + s₁ + det(M) s₂`. -/
lemma IsOrth.diag_weighted_le {M : Matrix (Fin 3) (Fin 3) ℝ} (hM : IsOrth M) (s : Fin 3 → ℝ)
    (h01 : s 1 ≤ s 0) (h12 : s 2 ≤ s 1) (h2 : 0 ≤ s 2) :
    M 0 0 * s 0 + M 1 1 * s 1 + M 2 2 * s 2 ≤ s 0 + s 1 + M.det * s 2 := by
  have hd0 := hM.diag_le_one 0
  have hd1 := hM.diag_le_one 1
  have hd2 := hM.diag_le_one 2
  rcases hM.det_eq with hdet | hdet
  · rw [hdet]
    nlinarith [mul_nonneg (sub_nonneg.mpr hd0) (h2.trans (h12.trans h01)),
      mul_nonneg (sub_nonneg.mpr hd1) (h2.trans h12), mul_nonneg (sub_nonneg.mpr hd2) h2]
  · rw [hdet]
    have ht := hM.trace_le_of_det_neg hdet
    nlinarith [mul_nonneg (sub_nonneg.mpr h01) (sub_nonneg.mpr hd0),
      mul_nonneg (sub_nonneg.mpr h12) (by linarith : (0:ℝ) ≤ 2 - M 0 0 - M 1 1),
      mul_nonneg h2 (by linarith : (0:ℝ) ≤ 1 - M 0 0 - M 1 1 - M 2 2)]

end
