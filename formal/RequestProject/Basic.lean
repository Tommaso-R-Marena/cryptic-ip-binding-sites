module

public import Mathlib

/-!
# Rigid motions of `ℝ³`

Points live in the Euclidean space `E3 = EuclideanSpace ℝ (Fin 3)`.
A matrix `A` is *orthogonal* if `Aᵀ * A = 1`, and a *rotation* if moreover `det A = 1`.
A rigid motion is a map `x ↦ A x + b` with `A` a rotation.
-/

@[expose] public section

open Matrix

noncomputable section

/-- Euclidean 3-space. -/
abbrev E3 := EuclideanSpace ℝ (Fin 3)

/-- A real `3 × 3` matrix is orthogonal if `Aᵀ A = 1`. -/
def IsOrth (A : Matrix (Fin 3) (Fin 3) ℝ) : Prop := Aᵀ * A = 1

/-- A real `3 × 3` matrix is a rotation if it is orthogonal with determinant `+1`. -/
def IsRot (A : Matrix (Fin 3) (Fin 3) ℝ) : Prop := IsOrth A ∧ A.det = 1

/-- The action of a matrix on a point of `E3`. -/
def mact (A : Matrix (Fin 3) (Fin 3) ℝ) (x : E3) : E3 := Matrix.toEuclideanLin A x

/-- A rigid motion `x ↦ A x + b` of `ℝ³`, with `A` orthogonal and `det A = +1`. -/
structure RigidMotion where
  /-- the linear (rotation) part -/
  A : Matrix (Fin 3) (Fin 3) ℝ
  /-- the translation part -/
  b : E3
  /-- `A` is orthogonal with determinant `+1` -/
  isRot : IsRot A

namespace RigidMotion

/-- Applying a rigid motion to a point. -/
def apply (T : RigidMotion) (x : E3) : E3 := mact T.A x + T.b

instance : CoeFun RigidMotion (fun _ => E3 → E3) := ⟨RigidMotion.apply⟩

lemma coe_apply (T : RigidMotion) (x : E3) : T x = mact T.A x + T.b := rfl

/-- The identity rigid motion. -/
def id : RigidMotion where
  A := 1
  b := 0
  isRot := ⟨by simp [IsOrth], by simp⟩

instance : Inhabited RigidMotion := ⟨id⟩

end RigidMotion

lemma mact_add (A : Matrix (Fin 3) (Fin 3) ℝ) (x y : E3) :
    mact A (x + y) = mact A x + mact A y := map_add _ _ _

lemma mact_sub (A : Matrix (Fin 3) (Fin 3) ℝ) (x y : E3) :
    mact A (x - y) = mact A x - mact A y := map_sub _ _ _

lemma mact_one (x : E3) : mact 1 x = x := by
  simp [mact]

lemma norm_sq_eq_dot (v : E3) : ‖v‖ ^ 2 = (WithLp.ofLp v) ⬝ᵥ (WithLp.ofLp v) := by
  rw [EuclideanSpace.norm_eq, Real.sq_sqrt (by positivity)]
  simp [dotProduct, sq]

/-- Orthogonal matrices preserve norms. -/
lemma norm_mact_of_isOrth {A : Matrix (Fin 3) (Fin 3) ℝ} (hA : IsOrth A) (x : E3) :
    ‖mact A x‖ = ‖x‖ := by
  have h : ‖mact A x‖ ^ 2 = ‖x‖ ^ 2 := by
    rw [norm_sq_eq_dot, norm_sq_eq_dot]
    show (A *ᵥ WithLp.ofLp x) ⬝ᵥ (A *ᵥ WithLp.ofLp x) = _
    rw [dotProduct_mulVec, ← vecMul_transpose, vecMul_vecMul, hA, vecMul_one]
  exact (pow_left_inj₀ (norm_nonneg _) (norm_nonneg _) two_ne_zero).mp h

/-- Rigid motions preserve distances. -/
lemma RigidMotion.dist_apply (T : RigidMotion) (x y : E3) : dist (T x) (T y) = dist x y := by
  rw [dist_eq_norm, dist_eq_norm, coe_apply, coe_apply, add_sub_add_right_eq_sub, ← mact_sub,
    norm_mact_of_isOrth T.isRot.1]

end
