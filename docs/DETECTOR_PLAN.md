# Detector plan: pre-registration for studies Q1-Q5

Dated 2026-10-10, **before any candidate, control or binder has been scored by
any method in `cryptic_ip/detect/`**. The code exists and is unit-tested; it
has been run on synthetic structures from `cryptic_ip.testing.synthetic` and on
nothing else. No figure in this document is a result, and the two numbers that
appear in it — the minimum detectable AUC and the measured descriptor
tolerances — are properties of the method and the dataset's size, computed
without looking at any outcome.

This plan covers the next-generation detector in `cryptic_ip/detect/` and
`scripts/detect.py`. It does not amend studies A-P; those stand as reported.

## Why a new detector at all

Four findings from this project constrain the design, and each one rules
something out.

1. **Burial is a continuum.** The 0.12 relative-SASA cutoff was calibrated on
   five controls and tested against 135 deposited entries; two boundary
   estimators disagree (density minimum 0.138, Otsu 0.463) and the quantiles
   run smoothly from 0.070 to 0.896. The positive class is *defined* by the
   cutoff rather than discovered, so every downstream AUROC is conditioned on
   it. A detector that classifies against that boundary inherits the
   circularity.
2. **The cascade of hard gates misses known binders.** ADAR2 and ADAR1 pass;
   ADAT1, PDS5B and HDAC1 are blocked at the score gate and HDAC3 at the
   volume gate. All four are annotated binders.
3. **Machine learning did not beat the rules**, held-out ROC-AUC 0.499 against
   0.556. The stated cause was class imbalance. The real cause is that
   unlabelled pockets were used as negatives, and that the positives occupy
   five strict homology groups with one holding 87% of them.
4. **Specificity is unresolved.** Re-weighting by P(IP | site) gave ΔAUC
   −0.015 [−0.044, 0.015], and 30 of 75 candidates are already flagged as
   explained by another phosphate- or sulfate-rich ligand.

And one more, which is the reason the detector's centre of gravity moves:

5. **A cryptic site is closed in the apo structure.** Every candidate is a
   single apo AlphaFold conformer, and the screen measures how buried a pocket
   is in it. The method is looking for the thing it is named after in the one
   state where it is hardest to see.

## What is new, and what it is replacing

| component | replaces | the claim being made |
|---|---|---|
| `detect/cavity.py` | the hard dependency on fpocket | a LIGSITE-style grid detector in process, so hundreds of conformers are affordable |
| `detect/ensemble.py` | single-conformer burial | persistence and **crypticity** across elastic-network conformers |
| `detect/openability.py` | "is the pocket open" | "what would it cost to open it", from low-frequency mode overlap |
| `detect/field.py` | the ≥4 basic residues gate | rotation-invariant electrostatic power spectrum and shell complementarity |
| `detect/specificity.py` | an absolute candidate call | a comparative margin against a competitor panel, with `polyanion_generic` as an outcome |
| `detect/monotone.py` | the five-gate cascade | a monotone additive score that abstains on missing input, with one physical veto |
| `detect/pu.py` | supervised classification | positive-unlabelled learning with group-held-out calibration and conformal sets |

Three measurements of the methods themselves, taken before any study data and
recorded so the parameters are not free:

* the field power spectrum is rotation-invariant to a discrepancy of
  **2.7 × 10⁻³** at the default 768 quadrature points, falling to 6.0 × 10⁻⁴ at
  3072, which is how 768 was chosen;
* the screening cutoff at 5 Debye lengths costs **3.0 × 10⁻⁴**, an order of
  magnitude below the quadrature error, which is why it is the default;
* the Elkan–Noto mean estimator of the label frequency is biased low by
  **30%** at every true value, independent of the classifier, so the default is
  the 95th percentile within each homology group.

## The questions, with their decision rules fixed now

Each question names its estimand, its interval method, and what would make it
answer "no". Intervals are 2000-resample bootstraps over whole
`homology_group_strict` groups, never over copies. Holm correction across the
five primary tests. **A stratum with fewer than 5 groups gives no evidence**
and is reported as not evaluable whichever way it points.

### Q1 — Does crypticity find sites that burial misses?

**Primary test.** Over the annotated binder set, the fraction whose binding
site is recovered by the ensemble's cryptic or semi-cryptic classification,
against the fraction recovered by the deployed burial gate on the apo
conformer. Paired by protein.

**Decision rule.** Q1 answers yes only if the paired difference excludes zero
and is at least **+0.10** in favour of crypticity. A smaller positive
difference is reported as no material change, because the ensemble costs 13×
the detection time and a gain below 10 points does not justify replacing a
method that is already understood.

**What would make it no.** Elastic-network modes are global and side-chain
free. If the known cryptic sites open by a side-chain rotation rather than a
domain motion — which is likely for a large share of them — crypticity will
not recover them and Q1 is negative. That is a real possible outcome and it is
the reason this is question one.

### Q2 — Does the field shape separate IP6 sites from other polyanion sites?

**Primary test.** Group-level ROC-AUC of the IP6-minus-best-competitor margin,
separating the annotated IP6 binders from a matched set of annotated
ATP, PAPS, UDP-sugar and sulfate binders.

**Decision rule.** Q2 answers yes only if the group-level AUC interval excludes
**0.5 from below and its lower bound exceeds 0.65**. Below that the margin is
not usable as a filter and the honest report is that the descriptor cannot do
the one thing it was built for.

**Matching.** Competitor-binder controls are matched on pocket volume, basic
residue count **and** net formal charge within the pocket, because study O's
control arm failed by matching on pLDDT and hull depth while the score
structurally required basic residues. The matching is run and its balance
published **before** any margin is computed.

### Q3 — Is the monotone score with abstention better than the gate cascade?

**Primary test.** On the annotated binder set, the recall of the monotone score
at the pre-registered threshold 0.55, against the recall of the five-gate
cascade, with abstentions counted as misses for the monotone score. Counting
abstentions as misses is deliberately unfavourable to the new method.

**Decision rule.** Q3 answers yes if recall improves and the number of
false positives on the matched control set does not rise by more than 20%.

**How the veto is counted.** A binder ruled out by the electrostatic veto
counts as a miss, like any other. If the veto costs recall on the annotated
binders, that is evidence against the veto and it is reported as such rather
than being relaxed: the point of fixing it in advance is that a condition
which loses true positives has to be allowed to say so.

### Q4 — Does positive-unlabelled learning beat the rule-based score?

**Primary test.** Group-level ROC-AUC of the PU posterior against the
rule-based score, both scored out-of-fold with whole homology groups held out.

**Decision rule.** The minimum detectable AUC at this dataset's size is
computed and published **first**. With buried-site positives in 5 strict
homology groups that figure is **0.87**, so any observed AUC below 0.87 is
reported as *no power to detect*, not as *no signal*. Q4 answers yes only if
the observed group-level AUC exceeds the minimum detectable value. Anything
else is a statement about the PDB's sampling of IP6 complexes.

**This is the question most likely to be unanswerable**, and saying so in
advance is the point. The previously reported 0.499 was read as a modelling
failure; at 5 groups it is indistinguishable from a method that works.

### Q5 — Does the in-process cavity detector agree with fpocket?

**Primary test.** Over the 272 redocking copies, the fraction whose
ligand-centre is within 4 Å of the centre of the largest cavity found by each
detector, and the Spearman correlation of their volumes.

**Decision rule.** Q5 is a calibration study, not a hypothesis test. Its
output is an agreement figure published alongside every result the new
detector produces. Volumes from the two detectors are **not** interchangeable
whatever the correlation turns out to be — fpocket fits alpha spheres to a
Voronoi tessellation and this counts enclosed grid points — and no downstream
threshold calibrated on one may be applied to the other.

## Parameters fixed now, and not to be moved on results

| parameter | value | basis |
|---|---|---|
| decision threshold | 0.55 | the midpoint of the monotone score, chosen before any pocket was scored |
| electrostatic veto | `field_depth_kt >= 0` | a well that repels the ligand cannot hold it, however buried the pocket; added after an acidic synthetic decoy scored 0.554 against the 0.55 threshold and was called a site, with burial and model confidence outvoting the physics. It is the only veto, and an unmeasured well abstains rather than clearing it |
| specificity margin | 0.15 | deliberately large, because the descriptor is unvalidated |
| field shells | 2.0, 3.5, 5.0, 6.5, 8.0, 10.0 Å | brackets an IP6's phosphate radius of 4-5 Å |
| field max degree | 4 | above it the power measures atomic granularity |
| quadrature points | 768 | measured: 2.7e-3 rotation discrepancy |
| screening cutoff | 5 Debye lengths | measured: 3.0e-4, below the quadrature error |
| Debye length | 7.9 Å | 150 mM monovalent salt at 310 K; physiological, not fitted |
| grid spacing | 0.8 Å | finer than LIGSITE's 1.0 because phosphates are 2.5 Å apart |
| enclosure threshold | 5 of 7 | LIGSITE's own recommendation for buried sites |
| ANM cutoff | 15 Å | the long-range convention for Cα networks |
| ensemble | modes 0-2, amplitudes ±1, ±2 Å | 12 conformers plus the apo |
| crypticity boundaries | cryptic ≥ 0.5, semi ≥ 0.25 | stated, not fitted |
| conformal miscoverage | 0.1 | sets contain the truth 90% of the time |
| label frequency | 95th percentile per group | measured: the mean is 30% biased low |
| bootstrap | 2000 resamples of strict groups | the project's standing convention |
| group floor | 5 | the project's standing convention |
| seed | 20261010 | |

## What this detector cannot do, stated before anyone asks

* **It does not predict affinity.** No output is calibrated against a measured
  binding constant, and none should be reported as one.
* **The opening cost is not in kT.** It is in units of the elastic-network
  spring constant times angstrom squared. The spring constant has never been
  calibrated against a measured opening free energy, so the cost is ordinal at
  best. The field carries its units in its name for that reason.
* **The electrostatic field is analytic, not Poisson-Boltzmann.** Formal
  charges under Debye-Hückel screening in uniform bulk-water permittivity. A
  buried pocket has a lower local permittivity, which this does not model; the
  consequence is that well depths are underestimated for the most buried sites,
  which is the direction that loses candidates rather than inventing them.
* **The ligand templates are spherically averaged.** An inositol phosphate is
  not spherically symmetric. The averaging is consistent with a descriptor
  that is itself rotation-averaged within each degree, and it is still an
  approximation of the ligand rather than the ligand.
* **Elastic-network modes cannot open a side-chain cryptic site.** Atoms move
  with their residue, rigidly. This is the sharpest limit on Q1.
* **One chemistry dominates the ground truth**: 134 InsP6 entries against 1
  InsP5. Nothing here speaks to InsP3 or InsP4 sites, and IP3 is in the
  competitor panel rather than the target for that reason.
* **No wet-lab validation.** Everything is computational and the output is a
  set of hypotheses to test.

## Amendment policy

Any change to this plan is a dated amendment in `docs/DETECTOR_PLAN_AMENDMENT_N.md`,
written **before** the affected results are read, stating what changed and why.
That is how studies J and O were handled and it is not optional here.
