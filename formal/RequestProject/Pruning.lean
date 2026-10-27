module

public import RequestProject.Basic

/-!
# Items 1 and 3: soundness / completeness of pairwise-distance pruning

A correspondence is `t q : Fin k → E3` (template anchors `t i` matched to query anchors `q i`).
-/

@[expose] public section

open Matrix

noncomputable section

/-- Chebyshev admissibility: some rigid motion moves every template anchor to within `ε` of its
query anchor, i.e. `∃ T, maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ ε`. -/
def Admissible {k : ℕ} (ε : ℝ) (t q : Fin k → E3) : Prop :=
  ∃ T : RigidMotion, ∀ i, ‖T (t i) - q i‖ ≤ ε

/-- The pairwise-distance pruning test with tolerance `τ`: keep the correspondence iff
`| ‖tᵢ − tⱼ‖ − ‖qᵢ − qⱼ‖ | ≤ τ` for all `i, j`. -/
def PassesPruning {k : ℕ} (τ : ℝ) (t q : Fin k → E3) : Prop :=
  ∀ i j, |‖t i - t j‖ - ‖q i - q j‖| ≤ τ

/-- **Item 1 (soundness of distance pruning).** If a rigid motion `T` satisfies
`‖T(tᵢ) − qᵢ‖ ≤ ε` for all `i`, then all pairwise distances agree up to `2ε`. -/
theorem distance_pruning_sound {k : ℕ} (t q : Fin k → E3) (ε : ℝ) (T : RigidMotion)
    (hT : ∀ i, ‖T (t i) - q i‖ ≤ ε) (i j : Fin k) :
    |‖t i - t j‖ - ‖q i - q j‖| ≤ 2 * ε := by
  have h1 := dist_dist_dist_le (T (t i)) (T (t j)) (q i) (q j)
  rw [T.dist_apply, Real.dist_eq, dist_eq_norm, dist_eq_norm, dist_eq_norm,
    dist_eq_norm] at h1
  have := hT i
  have := hT j
  linarith

/-- **Item 1, tightness of `2ε`.** For every `ε ≥ 0` there is a two-anchor correspondence and a
rigid motion (the identity) with all residuals `≤ ε` for which the distance discrepancy equals
`2ε` exactly: `t₀ = t₁ = 0`, `q₀ = (−ε,0,0)`, `q₁ = (ε,0,0)`. -/
theorem distance_pruning_tight (ε : ℝ) (hε : 0 ≤ ε) :
    ∃ (t q : Fin 2 → E3) (T : RigidMotion), (∀ i, ‖T (t i) - q i‖ ≤ ε) ∧
      |‖t 0 - t 1‖ - ‖q 0 - q 1‖| = 2 * ε := by
  refine ⟨fun _ => 0, ![-EuclideanSpace.single 0 ε, EuclideanSpace.single 0 ε],
    RigidMotion.id, ?_, ?_⟩
  · intro i
    fin_cases i <;> simp [RigidMotion.coe_apply, RigidMotion.id, mact_one, abs_of_nonneg hε]
  · have : -EuclideanSpace.single (0 : Fin 3) ε - EuclideanSpace.single 0 ε =
        EuclideanSpace.single 0 (-(2 * ε)) := by
      ext m
      by_cases hm : m = 0
      · simp [hm]; ring
      · simp [hm]
    simp only [Matrix.cons_val_zero, Matrix.cons_val_one, sub_self, norm_zero, this,
      EuclideanSpace.norm_single, Real.norm_eq_abs, abs_neg, zero_sub]
    rw [abs_of_nonneg (by linarith : (0:ℝ) ≤ 2 * ε), abs_of_nonneg (by linarith)]

/-- **Item 3 (completeness under the Chebyshev criterion).** Pruning with tolerance `2ε`
never discards an admissible correspondence. -/
theorem pruning_complete {k : ℕ} (ε : ℝ) (t q : Fin k → E3) (h : Admissible ε t q) :
    PassesPruning (2 * ε) t q := by
  obtain ⟨T, hT⟩ := h
  exact distance_pruning_sound t q ε T hT

end
