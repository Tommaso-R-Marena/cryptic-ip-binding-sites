# The next-generation pocket detector

Code in `cryptic_ip/detect/`, driven by `scripts/detect.py`, pre-registered in
`docs/DETECTOR_PLAN.md`, and cluster-ready as study Q in
`hpc/studies/detector.yml`.

**Nothing here has been run on real candidate structures.** Everything below is
either a property of the method, measured on synthetic input, or a design
decision with its reason. The five questions that would make any of it a
finding are Q1-Q5 in the plan, and none has run.

---

## 1. What problem each piece solves

The deployed screen's own results constrain the design. Each row is a finding
from studies A-P and what it rules out.

| finding | consequence | what replaces it |
|---|---|---|
| burial is a continuum; two boundary estimators disagree (0.138 vs 0.463) and the positive class is *defined* by the 0.12 cut | any classifier against that boundary is circular | a monotone score over continuous terms, with no class boundary |
| the five-gate cascade blocks ADAT1, PDS5B, HDAC1 and HDAC3 — all annotated binders | a conjunction of five gates has the union of five false-negative rates | one score, one physical veto, and abstention for missing input |
| held-out ML ROC-AUC 0.499 against a 0.556 baseline | not a capacity problem: unlabelled pockets were used as negatives, and positives sit in 5 homology groups | positive-unlabelled learning, group-held-out, with the minimum detectable effect computed first |
| ΔAUC −0.015 for IP-specific re-weighting; 30 of 75 candidates explained by another polyanion | an absolute "is this a polyanion pocket" has the answer yes everywhere | a comparative margin against a competitor panel, with `polyanion_generic` as a real outcome |
| every candidate is a single apo conformer | a cryptic site is closed in the apo state, so burial measures the wrong variable | persistence and crypticity across an elastic-network conformer ensemble |

## 2. The six modules

### `cavity.py` — detection in process

LIGSITE's protein-solvent-protein count on a 0.8 Å grid: mark the free points,
count the directions in which protein is met on both sides, keep points
enclosed in at least 5 of 7. Nothing novel, deliberately — the novelty is what
runs on top, and a new cavity detector underneath it would be a second thing to
validate.

Two consequences: `fpocket` stops being a hard requirement (the analyze path
currently raises without it), and hundreds of conformers per protein become
affordable. **It is not fpocket and its volumes are not interchangeable with
fpocket's** — alpha spheres on a Voronoi tessellation against enclosed grid
points. Q5 measures the agreement; no threshold calibrated on one may be
applied to the other.

Verified against geometries whose answer is known by construction: a closed
shell gives exactly one cavity at the origin with 7/7 enclosure, cutting a hole
in it drops the enclosure to 6.21/7, and a flat slab gives none.

### `ensemble.py` — crypticity, the signal that was missing

Displace along the softest elastic-network modes (both directions, ±1 and ±2 Å),
detect cavities in every conformer, match them across conformers by position,
and report per site:

* **persistence** — the fraction of conformers in which it exists;
* **crypticity** — `1 − apo_volume / best_volume`, so a site absent in the apo
  and open elsewhere scores 1.

A site with high persistence and low crypticity is an ordinary open pocket; one
with middling persistence and high crypticity is what this project set out to
find. Verified: an explicit ensemble whose apo has no cavity and which opens to
484 Å³ gives persistence 0.60, crypticity 1.00, classification `cryptic`.

**The sharpest limit:** atoms move rigidly with their residue, so this cannot
open a side-chain cryptic site. That is the reason Q1 is question one and is
allowed to come back negative.

### `openability.py` — cost to open, not whether open

An anisotropic network model Hessian on Cα atoms, its lowest non-trivial modes,
projected onto the pocket lining's outward directions. The **opening overlap**
is how much of a mode is the opening motion; weighted by 1/λ, because a mode's
thermal amplitude goes as the inverse of its eigenvalue.

Correctness rests on the null space: a free elastic network has exactly six
zero modes. The tests check that, plus that a rigid translation and a rigid
rotation both cost zero energy, and that a disconnected network is refused
rather than analysed.

**The cost is in `gamma * A^2`, not kT.** The spring constant has never been
calibrated against a measured opening free energy. The field is named
`opening_cost_gamma_a2` for that reason — it was `opening_cost_kt` until the
mismatch was noticed, and this pipeline has already read a hydrophobic density
as a depth in angstrom once.

### `field.py` — rotation-invariant electrostatics

The pocket's potential, from formal charges under Debye-Hückel screening at
physiological ionic strength, expanded on six concentric shells in real
spherical harmonics. The per-shell, per-degree **power** is kept, and that is
exactly invariant to rotation about the pocket centre — a rotation mixes the
harmonics within a degree orthogonally and leaves the sum of their squares
alone.

This is what removes the alignment that study N needed. Rigid transplantation of
131 crystal IP6 sites gave AUC 0.599 and placements that were not buried
(median buried fraction 0.000); a rotation-invariant descriptor has no
orientation to get wrong.

Measured properties, which is how the constants were chosen:

| quantity | value | at |
|---|---|---|
| rotation discrepancy | 2.4e-2 | 192 quadrature points |
| | 2.7e-3 | **768 (default)** |
| | 6.0e-4 | 3072 |
| screening-cutoff error | 3.0e-4 | 5 Debye lengths (default) |
| cost, 78,920 pockets, 1 core | 3.0 h | 768 points |

The l=4 harmonics were hand-derived at first and were wrong. Every lower degree
was exactly orthonormal, so the obvious checks passed; the only symptom was a
rotation invariance stuck at 130% error however fine the grid. **An invariance
that does not improve with more sampling is a basis error, not a quadrature
error** — that distinction found it, and a Gram-matrix test now pins it.

### `specificity.py` — comparative, not absolute

IP6 against a panel of ATP, PAPS, UDP-glucose, sulfate, phosphate, IP3 and a
heparin tetrasaccharide, reporting one of four outcomes with a pre-registered
margin of 0.15:

* `ip6_selective` — IP6 fits better than every competitor by the margin
* `polyanion_generic` — a competitor is within the margin; it is a polyanion
  pocket and this method cannot say which polyanion
* `competitor_preferred` — a competitor fits better
* `no_polyanion_site` — the well is unfavourable, so no comparison is meaningful

The middle outcome is the one the deployed screen cannot express and is probably
the right answer for a large share of the 75 candidates. Reporting it is not a
weaker claim than a candidate list; it is the claim the evidence supports. The
summary names *which* competitor explains how many pockets, because "30 of 75
are explained by another ligand" is weaker than knowing whether they are ATP
sites or sulfate sites.

A margin sweep is published next to the result, so a reader can see how much
the answer depends on the threshold — a check nobody can run on a tuned one.

### `monotone.py` — abstention and one veto

Every feature carries a direction declared in advance, a monotone ramp, a
non-negative weight, and a reason. Monotonicity is checked over 2,000 random
perturbations rather than trusted.

Three things it does that the gate cascade cannot:

1. **A missing input is an interval, not a failure.** A crystal structure has no
   pLDDT; a pipeline without APBS has no electrostatic term. The score is
   computed over the range the missing feature could take, and **abstains** when
   that interval straddles the threshold. A protein is never rejected for a
   measurement nobody made.
2. **Alternative measurements substitute.** `relative_sasa` (measured on a
   ligand copy) and `cavity_openness` (measured by the cavity detector) are the
   same property by different means, under separate names because they are
   separate measurements. Having one drops the other from the score rather than
   widening it.
3. **One physical veto.** `field_depth_kt >= 0`: a well that repels the ligand
   cannot hold it, however buried the pocket. This was added after an acidic
   synthetic decoy — perfectly enclosed, zero basic residues, −3.6 kT/e —
   scored 0.554 against the 0.55 threshold and was called a site. Burial and
   model confidence were outvoting the physics. An unmeasured well *abstains*
   rather than clearing the veto.

It also answers the question the cascade could not: **which measurement is
holding this pocket back**, which is what an experimentalist needs.

### `pu.py` — learning without negatives, and honest power

Labelled positives plus an unlabelled pool, with the label frequency `c`
estimated per homology group and held out whole. Then
`P(positive | x) = g(x) / c`, clipped and visibly saturated rather than
exceeding one.

Elkan and Noto's mean estimator of `c` is **biased low by 30% at every true
value**, measured with the estimator isolated from any classifier by handing it
the exact `g(x)`:

| true c | mean | q90 | q95 | q99 | max |
|---|---|---|---|---|---|
| 0.2 | 0.70 | 0.96 | 0.98 | 0.99 | 1.00 |
| 0.5 | 0.69 | 0.96 | 0.98 | 0.99 | 1.00 |
| 0.9 | 0.70 | 0.96 | 0.98 | 0.99 | 1.00 |

(each entry is the estimate divided by the truth). The bias is the class
overlap, and it matters in the dangerous direction — the posterior divides by
`c`, so a `c` that is too low inflates P(positive). The default is the 95th
percentile within each group; max and q99 rest on one or two rows, and with 134
entries in 5 groups there are not enough labelled positives to spend on an
extreme quantile.

Also here: split-conformal prediction sets (which can contain *both* labels,
the method saying it does not know), a leakage audit that runs on every fit,
and the **minimum detectable AUC**:

| positive groups | min detectable AUC |
|---|---|
| 2 | 1.00 |
| 3 | 0.98 |
| **5** | **0.87** |
| 10 | 0.77 |
| 40 | 0.65 |

This project's buried-site positives occupy 5 strict homology groups, one
holding 87% of them. **Nothing below AUC 0.87 is detectable at that size**, so
the previously reported 0.499 is a statement about power, not about biology.
`python scripts/detect.py power` prints this, and Q4's decision rule requires
it to be published before any AUC.

## 3. Running it

```bash
python scripts/detect.py spec                      # the specification, auditable
python scripts/detect.py power                     # what this dataset could detect
python scripts/detect.py describe model.pdb        # one structure, full detector
python scripts/detect.py describe model.pdb --no-ensemble   # skip most of the cost
python scripts/detect.py shard --structures DIR --shard 0 --shards 272 --out part.json
```

On the cluster, study Q in `hpc/studies/detector.yml`. It is **not** plan-gated,
because `docs/DETECTOR_PLAN.md` exists — unlike study M and study P, which
`rockfish.py submit` still refuses.

Measured cost: the full detector with a 13-conformer ensemble took **4.0 s** on
a 1,202-atom synthetic structure at 384 quadrature points. The candidates are
much larger (median 436 residues, maximum 2,352) and the cavity grid grows with
the bounding box, so the real per-structure cost is expected to be one to two
orders of magnitude higher on the long tail. That is an expectation, not a
measurement, which is why study Q's first stage is a timing probe on real
structures rather than a production sweep — the same discipline that let a
90-minute CPU probe kill a co-folding study before it wasted far more.

## 4. What it cannot do

* No output is calibrated against a measured affinity. None is a prediction.
* The opening cost is ordinal at best and is not in kT.
* The electrostatic field is analytic, in uniform bulk-water permittivity. A
  buried pocket's lower local permittivity is not modelled, so well depths are
  *underestimated* for the most buried sites — the direction that loses
  candidates rather than inventing them.
* The ligand templates are spherically averaged. An inositol phosphate is not.
* Elastic-network modes cannot open a side-chain cryptic site.
* One chemistry dominates the ground truth: 134 InsP6 against 1 InsP5. IP3 is
  in the competitor panel, not the target, for that reason.
* No wet-lab validation. The output is a set of hypotheses.
