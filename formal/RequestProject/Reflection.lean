module

public import RequestProject.Kabsch

/-!
# Why the constraint `det A = +1` is necessary

If one minimises `Σᵢ ‖A tᵢ + b − qᵢ‖²` over *all* orthogonal `A` (reflections allowed), the
optimum can be a reflection, which is not a rigid motion. Witness: a chiral tetrahedron and
its mirror image,
`t = (0, e₀, e₁, e₂)`, `q = (0, e₀, e₁, −e₂)`.
The reflection `P = diag(1, 1, −1)` matches them exactly (cost `0`), it is the *only* orthogonal
matrix achieving cost `0`, and every rigid motion has strictly positive cost.
-/

@[expose] public section

open Matrix

noncomputable section

/-- Template: the chiral tetrahedron `0, e₀, e₁, e₂`. -/
def chiralT : Fin 4 → E3 :=
  ![0, EuclideanSpace.single 0 1, EuclideanSpace.single 1 1, EuclideanSpace.single 2 1]

/-- Query: its mirror image `0, e₀, e₁, −e₂`. -/
def chiralQ : Fin 4 → E3 :=
  ![0, EuclideanSpace.single 0 1, EuclideanSpace.single 1 1, EuclideanSpace.single 2 (-1)]

/-- The reflection `P = diag(1, 1, −1)`. -/
def reflP : Matrix (Fin 3) (Fin 3) ℝ := diagonal ![1, 1, -1]

lemma reflP_isOrth : IsOrth reflP := by
  unfold IsOrth reflP
  rw [diagonal_transpose, diagonal_mul_diagonal, ← diagonal_one]
  congr 1
  ext i
  fin_cases i <;> simp

lemma reflP_det : reflP.det = -1 := by
  simp [reflP, det_diagonal, Fin.prod_univ_three]

lemma mact_single_apply (A : Matrix (Fin 3) (Fin 3) ℝ) (c a : Fin 3) :
    (mact A (EuclideanSpace.single c 1)) a = A a c := by
  simp [mact, Matrix.mulVec_single]

/-- The reflection matches the mirror image exactly. -/
lemma lsCost_reflP : lsCost chiralT chiralQ reflP 0 = 0 := by
  have h : ∀ i, mact reflP (chiralT i) - chiralQ i = 0 := by
    intro i
    ext a
    fin_cases i <;> fin_cases a <;>
      simp [chiralT, chiralQ, reflP, mact, Matrix.mulVec_single]
  simp [lsCost, h]

/-- Any orthogonal-Procrustes solution with zero cost must be the reflection `P`. -/
lemma eq_reflP_of_lsCost_zero (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3)
    (h : lsCost chiralT chiralQ A b = 0) : A = reflP := by
  have h0 : ∀ i, mact A (chiralT i) + b - chiralQ i = 0 := by
    intro i
    have := (Finset.sum_eq_zero_iff_of_nonneg (fun j _ => sq_nonneg
      ‖mact A (chiralT j) + b - chiralQ j‖)).mp h i (Finset.mem_univ i)
    exact norm_eq_zero.mp (pow_eq_zero_iff two_ne_zero |>.mp this)
  have hb : b = 0 := by
    have := h0 0
    simpa [chiralT, chiralQ, mact] using this
  subst hb
  ext a c
  have hc := congrArg (fun v : E3 => v a) (h0 c.succ)
  fin_cases a <;> fin_cases c <;>
    simp [chiralT, chiralQ, reflP, mact_single_apply, sub_eq_zero] at hc ⊢ <;> linarith

/-- **Necessity of `det A = +1`.** For the chiral tetrahedron and its mirror image:
* the reflection `P` (orthogonal, `det P = −1`) achieves cost `0`, the global minimum over all
  orthogonal `A`;
* every orthogonal `A` achieving that minimum equals `P`, so it is not a rigid motion;
* every rigid motion (`det A = +1`) has strictly positive cost.
So minimising over the full orthogonal group returns a reflection here; the restriction to
`det = +1` (the sign correction in Kabsch's construction) is needed to obtain a rigid motion. -/
theorem reflection_necessary :
    IsOrth reflP ∧ reflP.det = -1 ∧ lsCost chiralT chiralQ reflP 0 = 0 ∧
    (∀ (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3), lsCost chiralT chiralQ A b = 0 → A = reflP) ∧
    (∀ (A : Matrix (Fin 3) (Fin 3) ℝ) (b : E3), IsRot A → 0 < lsCost chiralT chiralQ A b) := by
  refine ⟨reflP_isOrth, reflP_det, lsCost_reflP, eq_reflP_of_lsCost_zero, ?_⟩
  intro A b hA
  rcases (Finset.sum_nonneg (fun j _ => sq_nonneg
      ‖mact A (chiralT j) + b - chiralQ j‖) : 0 ≤ lsCost chiralT chiralQ A b).lt_or_eq with
    hlt | heq
  · exact hlt
  · have := eq_reflP_of_lsCost_zero A b heq.symm
    have h2 := hA.2
    rw [this, reflP_det] at h2
    norm_num at h2

end
