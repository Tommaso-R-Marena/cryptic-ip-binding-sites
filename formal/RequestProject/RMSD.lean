module

public import RequestProject.Pruning

/-!
# Item 2: an RMSD criterion does not justify `2ρ`-style pruning
-/

@[expose] public section

open Matrix

noncomputable section

variable {k : ℕ}

/-- Sum of squared residuals `Σᵢ ‖T(tᵢ) − qᵢ‖²`. -/
def sumSq (t q : Fin k → E3) (T : RigidMotion) : ℝ := ∑ i, ‖T (t i) - q i‖ ^ 2

/-- `RMSD(T) = sqrt((1/k) Σᵢ ‖T(tᵢ) − qᵢ‖²)`. -/
def rmsd (t q : Fin k → E3) (T : RigidMotion) : ℝ := Real.sqrt (sumSq t q T / k)

/-- `R = min_T RMSD(T)` (formally the infimum over all rigid motions; it is attained, but none
of the results below need that). -/
def minRMSD (t q : Fin k → E3) : ℝ := ⨅ T : RigidMotion, rmsd t q T

lemma sumSq_nonneg (t q : Fin k → E3) (T : RigidMotion) : 0 ≤ sumSq t q T :=
  Finset.sum_nonneg fun _ _ => sq_nonneg _

lemma rmsd_nonneg (t q : Fin k → E3) (T : RigidMotion) : 0 ≤ rmsd t q T := Real.sqrt_nonneg _

lemma minRMSD_le_rmsd (t q : Fin k → E3) (T : RigidMotion) : minRMSD t q ≤ rmsd t q T :=
  ciInf_le ⟨0, by rintro _ ⟨T', rfl⟩; exact rmsd_nonneg _ _ _⟩ T

lemma rmsd_sq (t q : Fin k → E3) (T : RigidMotion) : rmsd t q T ^ 2 = sumSq t q T / k := by
  rw [rmsd, Real.sq_sqrt (div_nonneg (sumSq_nonneg _ _ _) (Nat.cast_nonneg _))]

lemma sumSq_eq_k_mul_rmsd_sq (t q : Fin k → E3) (T : RigidMotion) (hk : 0 < k) :
    sumSq t q T = k * rmsd t q T ^ 2 := by
  rw [rmsd_sq]; field_simp

/-- **Item 2(a): the general bound.** For *any* rigid motion `T`,
`maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ RMSD(T) · √k`. -/
theorem residual_le_rmsd_mul_sqrt_k (t q : Fin k → E3) (T : RigidMotion) (i : Fin k) :
    ‖T (t i) - q i‖ ≤ rmsd t q T * Real.sqrt k := by
  have hk : (0 : ℝ) < k := by exact_mod_cast Fin.pos i
  rw [rmsd, ← Real.sqrt_mul (div_nonneg (sumSq_nonneg _ _ _) hk.le), div_mul_cancel₀ _ hk.ne']
  apply Real.le_sqrt_of_sq_le
  exact Finset.single_le_sum (f := fun j => ‖T (t j) - q j‖ ^ 2)
    (fun _ _ => sq_nonneg _) (Finset.mem_univ i)

/-- **Item 2(a), sharpness for arbitrary `T`.** The factor `√k` cannot be improved for an
arbitrary rigid motion with `RMSD(T) = ρ`: put the whole error on one anchor. -/
theorem residual_bound_sqrt_k_tight (m : ℕ) (ρ : ℝ) (hρ : 0 ≤ ρ) :
    ∃ (t q : Fin (m + 1) → E3) (T : RigidMotion), rmsd t q T = ρ ∧
      ‖T (t 0) - q 0‖ = ρ * Real.sqrt (m + 1 : ℕ) := by
  refine ⟨fun _ => 0, fun j => if j = 0 then EuclideanSpace.single 0 (ρ * Real.sqrt (m + 1 : ℕ))
    else 0, RigidMotion.id, ?_, ?_⟩
  · have hs : sumSq (fun _ => (0 : E3))
        (fun j : Fin (m + 1) => if j = 0 then EuclideanSpace.single 0
          (ρ * Real.sqrt (m + 1 : ℕ)) else 0) RigidMotion.id = ρ ^ 2 * (m + 1 : ℕ) := by
      rw [sumSq, Fin.sum_univ_succ]
      simp [RigidMotion.coe_apply, RigidMotion.id, mact_one, mul_pow]
      left; exact Real.sq_sqrt (by positivity)
    rw [rmsd, hs, mul_div_assoc, div_self (by positivity), mul_one, Real.sqrt_sq hρ]
  · simp [RigidMotion.coe_apply, RigidMotion.id, mact_one, abs_of_nonneg hρ]

lemma sum_norm_sub_sq (r : Fin k → E3) (v : E3) :
    ∑ i, ‖r i - v‖ ^ 2 = ∑ i, ‖r i‖ ^ 2 - 2 * inner ℝ (∑ i, r i) v + k * ‖v‖ ^ 2 := by
  simp_rw [norm_sub_sq_real]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib, sum_inner, Finset.mul_sum]
  simp

lemma rmsd_le_iff (t q : Fin k → E3) (T T' : RigidMotion) (hk : 0 < k) :
    rmsd t q T ≤ rmsd t q T' ↔ sumSq t q T ≤ sumSq t q T' := by
  have hkR : (0 : ℝ) < k := by exact_mod_cast hk
  rw [rmsd, rmsd, Real.sqrt_le_sqrt_iff (div_nonneg (sumSq_nonneg _ _ _) hkR.le),
    div_le_div_iff_of_pos_right hkR]

/-- At a minimiser of the RMSD, the residual vectors sum to zero (optimality in the
translation). -/
lemma sum_residual_eq_zero_of_minimizer (t q : Fin k → E3) (T : RigidMotion)
    (hT : ∀ T' : RigidMotion, rmsd t q T ≤ rmsd t q T') :
    ∑ i, (T (t i) - q i) = 0 := by
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp
  set c := ∑ i, (T (t i) - q i) with hc_def
  let T' : RigidMotion := ⟨T.A, T.b - (1 / (k : ℝ)) • c, T.isRot⟩
  have hres : ∀ i, T' (t i) - q i = (T (t i) - q i) - (1 / (k : ℝ)) • c := by
    intro i; simp only [T', RigidMotion.coe_apply]; abel
  have hkR : (0 : ℝ) < k := by exact_mod_cast hk
  have h1 : sumSq t q T' = sumSq t q T - ‖c‖ ^ 2 / k := by
    simp only [sumSq, hres]
    rw [sum_norm_sub_sq, ← hc_def, real_inner_smul_right, real_inner_self_eq_norm_sq, norm_smul,
      Real.norm_eq_abs, abs_of_pos (by positivity)]
    field_simp
    ring
  have h2 := (rmsd_le_iff t q T T' hk).mp (hT T')
  have h3 : ‖c‖ ^ 2 / k ≤ 0 := by linarith
  have h4 : ‖c‖ ^ 2 ≤ 0 := by
    have := mul_le_mul_of_nonneg_right h3 hkR.le
    rwa [div_mul_cancel₀ _ hkR.ne', zero_mul] at this
  have : ‖c‖ = 0 := by nlinarith [norm_nonneg c]
  exact norm_eq_zero.mp this

/-- **Item 2(b): sharper bound at the minimiser.** If `T` minimises the RMSD, then
`maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ RMSD(T) · √(k − 1)`. In particular the bound `ρ√k` claimed to be best is
*not* sharp for the minimising motion. -/
theorem minimizer_residual_le (t q : Fin k → E3) (T : RigidMotion)
    (hT : ∀ T' : RigidMotion, rmsd t q T ≤ rmsd t q T') (i : Fin k) :
    ‖T (t i) - q i‖ ≤ rmsd t q T * Real.sqrt ((k : ℝ) - 1) := by
  have hk : 0 < k := Fin.pos i
  have hkR : (0 : ℝ) < k := by exact_mod_cast hk
  set r : Fin k → E3 := fun j => T (t j) - q j with hr
  have hsum : ∑ j, r j = 0 := sum_residual_eq_zero_of_minimizer t q T hT
  have hS : sumSq t q T = ∑ j, ‖r j‖ ^ 2 := rfl
  have hsplit : r i + ∑ j ∈ Finset.univ.erase i, r j = 0 := by
    rw [Finset.add_sum_erase _ _ (Finset.mem_univ i), hsum]
  have hri : r i = -∑ j ∈ Finset.univ.erase i, r j := eq_neg_of_add_eq_zero_left hsplit
  have hcard : ((Finset.univ.erase i).card : ℝ) = k - 1 := by
    rw [Finset.card_erase_of_mem (Finset.mem_univ i), Finset.card_univ, Fintype.card_fin,
      Nat.cast_sub hk]
    simp
  have hsq_split : ‖r i‖ ^ 2 + ∑ j ∈ Finset.univ.erase i, ‖r j‖ ^ 2 = sumSq t q T := by
    rw [hS, Finset.add_sum_erase _ (fun j => ‖r j‖ ^ 2) (Finset.mem_univ i)]
  have h1 : ‖r i‖ ≤ ∑ j ∈ Finset.univ.erase i, ‖r j‖ := by
    rw [hri, norm_neg]; exact norm_sum_le _ _
  have h2 : ‖r i‖ ^ 2 ≤ ((k : ℝ) - 1) * ∑ j ∈ Finset.univ.erase i, ‖r j‖ ^ 2 := by
    calc ‖r i‖ ^ 2 ≤ (∑ j ∈ Finset.univ.erase i, ‖r j‖) ^ 2 :=
          pow_le_pow_left₀ (norm_nonneg _) h1 2
      _ ≤ _ := by rw [← hcard]; exact sq_sum_le_card_mul_sum_sq
  have hSk := sumSq_eq_k_mul_rmsd_sq t q T hk
  have h3 : (k : ℝ) * ‖r i‖ ^ 2 ≤ (k : ℝ) * (((k : ℝ) - 1) * rmsd t q T ^ 2) := by
    nlinarith
  have h4 : ‖r i‖ ^ 2 ≤ (rmsd t q T * Real.sqrt ((k : ℝ) - 1)) ^ 2 := by
    have hk1 : (0 : ℝ) ≤ (k : ℝ) - 1 := by
      have : (1 : ℝ) ≤ k := by exact_mod_cast hk
      linarith
    rw [mul_pow, Real.sq_sqrt hk1]
    nlinarith [le_of_mul_le_mul_left h3 hkR]
  exact (pow_le_pow_iff_left₀ (norm_nonneg _)
    (mul_nonneg (rmsd_nonneg _ _ _) (Real.sqrt_nonneg _)) two_ne_zero).mp h4

/-- A minimiser realises the minimal RMSD `R`. -/
lemma rmsd_eq_minRMSD_of_minimizer (t q : Fin k → E3) (T : RigidMotion)
    (hT : ∀ T' : RigidMotion, rmsd t q T ≤ rmsd t q T') : rmsd t q T = minRMSD t q :=
  le_antisymm (le_ciInf hT) (minRMSD_le_rmsd t q T)

/-- **Item 2(b), in terms of `R`.** If `R ≤ ρ` and `T` is the minimising motion then
`maxᵢ ‖T(tᵢ) − qᵢ‖ ≤ ρ · √(k − 1)`. -/
theorem minimizer_residual_le_of_minRMSD_le (t q : Fin k → E3) (T : RigidMotion)
    (hT : ∀ T' : RigidMotion, rmsd t q T ≤ rmsd t q T') (ρ : ℝ) (hR : minRMSD t q ≤ ρ)
    (i : Fin k) : ‖T (t i) - q i‖ ≤ ρ * Real.sqrt ((k : ℝ) - 1) := by
  have h := minimizer_residual_le t q T hT i
  rw [rmsd_eq_minRMSD_of_minimizer t q T hT] at h
  exact h.trans (mul_le_mul_of_nonneg_right hR (Real.sqrt_nonneg _))

/-- If all template anchors coincide at the origin and the query anchors have centroid `0`,
then the identity motion minimises the RMSD. -/
lemma id_minimizer_of_t_zero (t q : Fin k → E3) (ht : ∀ i, t i = 0) (hq : ∑ i, q i = 0) :
    ∀ T' : RigidMotion, rmsd t q RigidMotion.id ≤ rmsd t q T' := by
  intro T'
  rcases Nat.eq_zero_or_pos k with rfl | hk
  · simp [rmsd]
  rw [rmsd_le_iff t q _ _ hk]
  have h0 : ∀ i, T' (t i) - q i = -(q i - T'.b) := by
    intro i; simp [RigidMotion.coe_apply, ht i, mact]
  have h1 : ∀ i, RigidMotion.id (t i) - q i = -(q i) := by
    intro i; simp [RigidMotion.coe_apply, RigidMotion.id, ht i, mact_one]
  simp only [sumSq, h0, h1, norm_neg]
  rw [sum_norm_sub_sq, hq, inner_zero_left]
  nlinarith [sq_nonneg ‖T'.b‖, (Nat.cast_nonneg k : (0:ℝ) ≤ k)]

/-- **Item 2(b), sharpness of `√(k − 1)`.** For every `k = m + 1 ≥ 1` there is a correspondence
whose RMSD-minimising motion (the identity) has a residual equal to `R · √(k − 1)`:
all `tᵢ = 0`, `q₀ = (−m, 0, 0)` and `qⱼ = (1, 0, 0)` for `j ≠ 0`. -/
theorem minimizer_residual_bound_tight (m : ℕ) :
    ∃ t q : Fin (m + 1) → E3,
      (∀ T' : RigidMotion, rmsd t q RigidMotion.id ≤ rmsd t q T') ∧
      ‖RigidMotion.id (t 0) - q 0‖ = rmsd t q RigidMotion.id * Real.sqrt (((m + 1 : ℕ) : ℝ) - 1) := by
  set q : Fin (m + 1) → E3 := fun j =>
    if j = 0 then EuclideanSpace.single 0 (-(m : ℝ)) else EuclideanSpace.single 0 1 with hq
  have hsum : ∑ j, q j = 0 := by
    rw [Fin.sum_univ_succ]
    simp only [hq, if_pos rfl, Fin.succ_ne_zero, if_false, Finset.sum_const, Finset.card_univ,
      Fintype.card_fin]
    ext l
    by_cases hl : l = 0 <;> simp [hl]
  have hid : ∀ j, RigidMotion.id ((fun _ => (0 : E3)) j) - q j = -(q j) := by
    intro j; simp [RigidMotion.coe_apply, RigidMotion.id, mact_one]
  have hS : sumSq (fun _ => (0 : E3)) q RigidMotion.id = (m : ℝ) ^ 2 + m := by
    simp only [sumSq, hid, norm_neg]
    rw [Fin.sum_univ_succ]
    simp [hq, Fin.succ_ne_zero]
  have hr : rmsd (fun _ => (0 : E3)) q RigidMotion.id = Real.sqrt m := by
    rw [rmsd, hS]
    congr 1
    push_cast
    field_simp
  refine ⟨fun _ => 0, q, id_minimizer_of_t_zero _ q (fun _ => rfl) hsum, ?_⟩
  rw [hr, hid, norm_neg]
  push_cast
  rw [add_sub_cancel_right, Real.mul_self_sqrt (Nat.cast_nonneg _)]
  simp [hq]

/-- For any rigid motion `T`, each pairwise discrepancy is at most `RMSD(T) · √(2k)`. -/
lemma discrepancy_le_rmsd_mul (t q : Fin k → E3) (T : RigidMotion) (i j : Fin k) :
    |‖t i - t j‖ - ‖q i - q j‖| ≤ rmsd t q T * Real.sqrt (2 * k) := by
  have hk : 0 < k := Fin.pos i
  have h1 := dist_dist_dist_le (T (t i)) (T (t j)) (q i) (q j)
  rw [T.dist_apply, Real.dist_eq, dist_eq_norm, dist_eq_norm, dist_eq_norm,
    dist_eq_norm] at h1
  have hrhs : rmsd t q T * Real.sqrt (2 * k) = Real.sqrt (2 * sumSq t q T) := by
    rw [sumSq_eq_k_mul_rmsd_sq t q T hk, ← Real.sqrt_sq (rmsd_nonneg t q T),
      ← Real.sqrt_mul (sq_nonneg _), Real.sqrt_sq (rmsd_nonneg t q T)]
    congr 1; ring
  rw [hrhs]
  by_cases hij : i = j
  · subst hij; simp
  have hpair : ‖T (t i) - q i‖ ^ 2 + ‖T (t j) - q j‖ ^ 2 ≤ sumSq t q T := by
    have := Finset.sum_le_sum_of_subset_of_nonneg (f := fun l => ‖T (t l) - q l‖ ^ 2)
      (Finset.subset_univ {i, j}) (fun _ _ _ => sq_nonneg _)
    rwa [Finset.sum_pair hij] at this
  apply h1.trans
  apply Real.le_sqrt_of_sq_le
  nlinarith [sq_nonneg (‖T (t i) - q i‖ - ‖T (t j) - q j‖)]

/-- **Item 2(c): soundness of the RMSD-based tolerance.** If `R ≤ ρ` then every pairwise
distance discrepancy is at most `ρ · √(2k)`. -/
theorem pruning_sound_of_minRMSD_le (t q : Fin k → E3) (ρ : ℝ) (hR : minRMSD t q ≤ ρ) :
    PassesPruning (ρ * Real.sqrt (2 * k)) t q := by
  intro i j
  have hk0 : (0 : ℝ) < k := by exact_mod_cast Fin.pos i
  have hk : (0 : ℝ) < 2 * k := by linarith
  have hs := Real.sqrt_pos.mpr hk
  have h : |‖t i - t j‖ - ‖q i - q j‖| / Real.sqrt (2 * k) ≤ minRMSD t q := by
    apply le_ciInf
    intro T
    rw [div_le_iff₀ hs]
    exact discrepancy_le_rmsd_mul t q T i j
  rw [div_le_iff₀ hs] at h
  exact h.trans (mul_le_mul_of_nonneg_right hR hs.le)

/-- **Item 2(c), sharpness.** For `k ≥ 2` and `ρ ≥ 0` the tolerance `ρ√(2k)` is attained:
there is a correspondence with `R ≤ ρ` whose discrepancy between anchors `0` and `1` is exactly
`ρ√(2k)`. -/
theorem pruning_tolerance_attained (m : ℕ) (ρ : ℝ) (hρ : 0 ≤ ρ) :
    ∃ t q : Fin (m + 2) → E3, minRMSD t q ≤ ρ ∧
      |‖t 0 - t 1‖ - ‖q 0 - q 1‖| = ρ * Real.sqrt (2 * (m + 2 : ℕ)) := by
  set a : ℝ := ρ * Real.sqrt (2 * (m + 2 : ℕ)) / 2 with ha
  have ha0 : 0 ≤ a := by positivity
  set q : Fin (m + 2) → E3 := fun j =>
    if j = 0 then EuclideanSpace.single 0 (-a) else
      if j = 1 then EuclideanSpace.single 0 a else 0 with hq
  have hS : sumSq (fun _ => (0 : E3)) q RigidMotion.id = 2 * a ^ 2 := by
    have hid : ∀ j, RigidMotion.id ((fun _ => (0 : E3)) j) - q j = -(q j) := by
      intro j; simp [RigidMotion.coe_apply, RigidMotion.id, mact_one]
    simp only [sumSq, hid, norm_neg]
    rw [Fin.sum_univ_succ, Fin.sum_univ_succ]
    have h1 : ∀ j : Fin m, (j.succ.succ : Fin (m + 2)) ≠ 1 := by
      intro j h
      have := congrArg Fin.val h
      simp [Fin.val_succ] at this
    simp [hq, Fin.succ_ne_zero, h1]
    ring
  have hr : rmsd (fun _ => (0 : E3)) q RigidMotion.id = ρ := by
    rw [rmsd, hS, ha, div_pow, mul_pow, Real.sq_sqrt (by positivity)]
    convert Real.sqrt_sq hρ using 2
    push_cast
    field_simp
  refine ⟨fun _ => 0, q, hr ▸ minRMSD_le_rmsd _ _ _, ?_⟩
  have hq01 : q 0 - q 1 = EuclideanSpace.single 0 (-(2 * a)) := by
    have : (1 : Fin (m + 2)) ≠ 0 := by simp
    simp only [hq, if_pos rfl, if_neg this]
    ext l
    by_cases hl : l = 0
    · simp [hl]; ring
    · simp [hl]
  rw [hq01, EuclideanSpace.norm_single, Real.norm_eq_abs, abs_neg,
    abs_of_nonneg (by linarith : (0:ℝ) ≤ 2 * a)]
  simp only [sub_self, norm_zero, zero_sub, abs_neg]
  rw [abs_of_nonneg (by linarith : (0:ℝ) ≤ 2 * a), ha]
  ring

/-- **Item 2(c), conclusion.** For `k ≥ 2` and `ρ ≥ 0`, the pruning test with tolerance `τ`
keeps every correspondence with `R ≤ ρ` if and only if `τ ≥ ρ√(2k)`. So the smallest safe
tolerance is exactly `τ = ρ√(2k)`. -/
theorem smallest_safe_tolerance (m : ℕ) (ρ : ℝ) (hρ : 0 ≤ ρ) (τ : ℝ) :
    (∀ t q : Fin (m + 2) → E3, minRMSD t q ≤ ρ → PassesPruning τ t q) ↔
      ρ * Real.sqrt (2 * (m + 2 : ℕ)) ≤ τ := by
  constructor
  · intro h
    obtain ⟨t, q, hR, hd⟩ := pruning_tolerance_attained m ρ hρ
    rw [← hd]
    exact h t q hR 0 1
  · intro hτ t q hR i j
    exact (pruning_sound_of_minRMSD_le t q ρ hR i j).trans hτ

/-- The concrete `k = 4` template: all four template anchors at the origin. -/
def exT : Fin 4 → E3 := fun _ => 0

/-- The concrete `k = 4` query: `(−5,0,0), (5,0,0), 0, 0`. -/
def exQ : Fin 4 → E3 :=
  ![EuclideanSpace.single 0 (-5), EuclideanSpace.single 0 5, 0, 0]

/-- **Item 2(d): concrete counterexample, `k = 4`, `ρ = 4`.**
For `tᵢ = 0` and `q = ((−5,0,0), (5,0,0), 0, 0)`:
* `R ≤ 4` (indeed `R = √12.5 ≈ 3.54`);
* the identity is an RMSD-minimising motion, yet its residual at anchor `0` is `5 > 4 = ρ`;
* the discrepancy `| ‖t₀ − t₁‖ − ‖q₀ − q₁‖ | = 10`, so pruning with `τ = 1.5`
  (and even with `τ = 1.5 · ρ = 6`) discards this correspondence. -/
theorem rmsd_counterexample_k4 :
    minRMSD exT exQ ≤ 4 ∧
    (∀ T' : RigidMotion, rmsd exT exQ RigidMotion.id ≤ rmsd exT exQ T') ∧
    ‖RigidMotion.id (exT 0) - exQ 0‖ = 5 ∧
    |‖exT 0 - exT 1‖ - ‖exQ 0 - exQ 1‖| = 10 ∧
    ¬ PassesPruning 1.5 exT exQ ∧ ¬ PassesPruning (1.5 * 4) exT exQ := by
  have hsum : ∑ i, exQ i = 0 := by
    rw [Fin.sum_univ_four]
    ext l
    by_cases hl : l = 0 <;> simp [exQ, hl]
  have hid : ∀ j, RigidMotion.id (exT j) - exQ j = -(exQ j) := by
    intro j; simp [RigidMotion.coe_apply, RigidMotion.id, mact_one, exT]
  have hS : sumSq exT exQ RigidMotion.id = 50 := by
    simp only [sumSq, hid, norm_neg]
    rw [Fin.sum_univ_four]
    simp [exQ]
    norm_num
  have hr : rmsd exT exQ RigidMotion.id ≤ 4 := by
    rw [rmsd, hS, Real.sqrt_le_left (by norm_num)]
    norm_num
  have h01 : exQ 0 - exQ 1 = EuclideanSpace.single 0 (-10) := by
    ext l
    by_cases hl : l = 0
    · simp [exQ, hl]; norm_num
    · simp [exQ, hl]
  have hd : |‖exT 0 - exT 1‖ - ‖exQ 0 - exQ 1‖| = 10 := by
    rw [h01]
    simp [exT]
  refine ⟨(minRMSD_le_rmsd _ _ _).trans hr, id_minimizer_of_t_zero exT exQ (fun _ => rfl) hsum,
    ?_, hd, ?_, ?_⟩
  · rw [hid, norm_neg]; simp [exQ]
  · intro h; have := h 0 1; rw [hd] at this; norm_num at this
  · intro h; have := h 0 1; rw [hd] at this; norm_num at this

end
