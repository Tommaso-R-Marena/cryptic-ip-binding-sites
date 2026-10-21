module

public import RequestProject.Kabsch

/-!
# Existence of a singular value decomposition for real `3 × 3` matrices

Every real `3 × 3` matrix `H` can be written `H = U · diag(s₀, s₁, s₂) · Vᵀ` with `U, V`
orthogonal and `s₀ ≥ s₁ ≥ s₂ ≥ 0`. Consequently the Kabsch construction is always available,
and it yields an optimal rigid motion for every correspondence.
-/

@[expose] public section

open Matrix

noncomputable section

lemma inner_E3 (x y : E3) : inner ℝ x y = ∑ a, x a * y a := by
  simp [PiLp.inner_apply, mul_comm]

lemma inner_mact_transpose_mul (H : Matrix (Fin 3) (Fin 3) ℝ) (x y : E3) :
    inner ℝ (mact (Hᵀ * H) x) y = inner ℝ (mact H x) (mact H y) := by
  simp only [inner_E3, mact, Matrix.toLpLin_apply]
  simp [mulVec, dotProduct, mul_apply, Fin.sum_univ_three]
  ring

lemma isOrth_of_onb (e : OrthonormalBasis (Fin 3) ℝ E3) : IsOrth (Matrix.of fun a i => e i a) := by
  unfold IsOrth
  ext i j
  have h := e.inner_eq_ite i j
  rw [inner_E3] at h
  simp only [mul_apply, transpose_apply, of_apply, one_apply]
  rw [← h]

/-- **SVD for real `3 × 3` matrices.** -/
theorem svd_exists (H : Matrix (Fin 3) (Fin 3) ℝ) :
    ∃ (U V : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ), IsOrth U ∧ IsOrth V ∧
      s 1 ≤ s 0 ∧ s 2 ≤ s 1 ∧ 0 ≤ s 2 ∧ H = U * diagonal s * Vᵀ := by
  set T : E3 →ₗ[ℝ] E3 := Matrix.toEuclideanLin (Hᵀ * H) with hT
  have hTsym : T.IsSymmetric := by
    intro x y
    have h1 := inner_mact_transpose_mul H x y
    have h2 := inner_mact_transpose_mul H y x
    simp only [mact] at h1 h2
    rw [hT, h1, real_inner_comm ((toEuclideanLin (Hᵀ * H)) y) x, h2]
    exact real_inner_comm _ _
  have hn : Module.finrank ℝ E3 = 3 := finrank_euclideanSpace_fin
  set v := hTsym.eigenvectorBasis hn with hv
  set lam := hTsym.eigenvalues hn with hlam
  have hTv : ∀ i, mact (Hᵀ * H) (v i) = lam i • v i := by
    intro i
    have := hTsym.apply_eigenvectorBasis hn i
    rw [RCLike.ofReal_real_eq_id, id] at this
    exact this
  have hHv : ∀ i j, inner ℝ (mact H (v i)) (mact H (v j)) = if i = j then lam i else 0 := by
    intro i j
    rw [← inner_mact_transpose_mul, hTv, real_inner_smul_left, v.inner_eq_ite]
    split_ifs <;> simp
  have hlam0 : ∀ i, 0 ≤ lam i := by
    intro i
    have := hHv i i
    rw [if_pos rfl, real_inner_self_eq_norm_sq] at this
    rw [← this]; positivity
  have hanti : Antitone lam := hTsym.eigenvalues_antitone hn
  set s : Fin 3 → ℝ := fun i => Real.sqrt (lam i) with hs
  set w : Fin 3 → E3 := fun i => (s i)⁻¹ • mact H (v i) with hw
  set S : Set (Fin 3) := {i | lam i ≠ 0} with hS
  have hortho : Orthonormal ℝ (S.restrict w) := by
    rw [orthonormal_iff_ite]
    rintro ⟨i, hi⟩ ⟨j, hj⟩
    simp only [Set.restrict_apply, hw, real_inner_smul_left, real_inner_smul_right, hHv,
      Subtype.mk.injEq]
    split_ifs with hij
    · subst hij
      have hsi : s i ^ 2 = lam i := Real.sq_sqrt (hlam0 i)
      have hne : s i ≠ 0 := by
        intro h0; apply hi; rw [← hsi, h0]; ring
      field_simp
      exact hsi.symm
    · simp
  obtain ⟨b, hb⟩ := hortho.exists_orthonormalBasis_extension_of_card_eq
    (by rw [hn, Fintype.card_fin])
  have hHvb : ∀ i, mact H (v i) = s i • b i := by
    intro i
    by_cases hi : lam i = 0
    · have h0 : mact H (v i) = 0 := by
        have := hHv i i
        rw [if_pos rfl, hi, real_inner_self_eq_norm_sq] at this
        exact norm_eq_zero.mp (pow_eq_zero_iff two_ne_zero |>.mp this)
      rw [h0]; simp [hs, hi]
    · rw [hb i hi, hw]
      have hsi : s i ^ 2 = lam i := Real.sq_sqrt (hlam0 i)
      have hne : s i ≠ 0 := by
        intro h0; apply hi; rw [← hsi, h0]; ring
      simp only [smul_smul, mul_inv_cancel₀ hne, one_smul]
  refine ⟨Matrix.of fun a i => b i a, Matrix.of fun a i => v i a, s, isOrth_of_onb b,
    isOrth_of_onb v, ?_, ?_, Real.sqrt_nonneg _, ?_⟩
  · exact Real.sqrt_le_sqrt (hanti (by decide))
  · exact Real.sqrt_le_sqrt (hanti (by decide))
  · have hVV : (Matrix.of fun a i => v i a) * (Matrix.of fun a i => v i a)ᵀ = 1 :=
      (isOrth_of_onb v).mul_transpose
    have hHV : H * (Matrix.of fun a i => v i a) = (Matrix.of fun a i => b i a) * diagonal s := by
      ext a i
      have := congrArg (fun x : E3 => x a) (hHvb i)
      simp only [mact, Matrix.toLpLin_apply, PiLp.toLp_apply, PiLp.smul_apply,
        smul_eq_mul] at this
      have lhs : ∑ x, H a x * (v i) x = s i * b i a := by
        rw [← this]; simp [mulVec, dotProduct]
      simp only [mul_apply, of_apply]
      rw [lhs]
      simp [diagonal_apply, mul_comm]
    calc H = H * ((Matrix.of fun a i => v i a) * (Matrix.of fun a i => v i a)ᵀ) := by
            rw [hVV, Matrix.mul_one]
      _ = (H * (Matrix.of fun a i => v i a)) * (Matrix.of fun a i => v i a)ᵀ := by
            rw [Matrix.mul_assoc]
      _ = _ := by rw [hHV]

/-- **Kabsch, unconditional form.** For every correspondence there is an SVD of the
cross-covariance matrix, and the resulting Kabsch motion is a rigid motion minimising
`Σᵢ ‖A tᵢ + b − qᵢ‖²` over all rotations `A` and translations `b`. -/
theorem kabsch_exists_optimal {k : ℕ} (t q : Fin k → E3) :
    ∃ (U V : Matrix (Fin 3) (Fin 3) ℝ) (s : Fin 3 → ℝ),
      IsOrth U ∧ IsOrth V ∧ s 1 ≤ s 0 ∧ s 2 ≤ s 1 ∧ 0 ≤ s 2 ∧
      crossCov t q = U * diagonal s * Vᵀ ∧
      IsRot (kabschRot U V) ∧
      ∀ (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3), IsRot A →
        lsCost t q (kabschRot U V) (kabschTrans t q U V) ≤ lsCost t q A b := by
  obtain ⟨U, V, s, hU, hV, h01, h12, h2, hH⟩ := svd_exists (crossCov t q)
  exact ⟨U, V, s, hU, hV, h01, h12, h2, hH, kabsch_optimal t q U V s hU hV h01 h12 h2 hH⟩

end
