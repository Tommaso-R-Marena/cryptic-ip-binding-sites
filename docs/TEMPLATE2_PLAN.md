# Can a real IP6 site be transplanted into a *buried* candidate pocket? Study O: plan

Pre-registered on 2026-10-01, **before any fit under this design has been computed**. The
only numbers quoted below come from `results/redocking/census.csv` (reference crystal
structures, already in the repository) and from study N's published outcome. No candidate
has been scored under this design.

## Why study N needs replacing rather than amending

Study N (`docs/TEMPLATE_PLAN.md`, run 36743030815) passed both its primary tests and was
still uninformative. Three defects, all recorded in `results/template/PROVENANCE.md`:

1. **N3's threshold was degenerate.** The 95th percentile of the control arm landed on the
   censoring floor of −4.0 Å, because at least 71 of 75 controls admitted no placement at
   all. "Above the 95th percentile of controls" reduced to "was not censored".
2. **The placements were not buried.** Median buried fraction **0.000**, zero for 22 of 35.
   A transplant constrained only by "do not clash" finds room on the surface, which is not
   what this project is about.
3. **The search was not exhaustive.** The code claimed distance pruning lost nothing
   admissible. It is false under an RMSD acceptance criterion: RMSD bounds the *mean*
   residual, so pruning at 1.5 Å discards correspondences well inside a 4.0 Å ceiling. The
   machine-checked bounds are in `formal/` — sound pruning there would need **τ ≥ ρ√(2k)**,
   which for k = 4, ρ = 4 Å is **11.31 Å**, wider than the anchor constellations themselves.

A fourth defect was in the control arm, which the fit score makes unscoreable: the score
needs at least 4 basic residues in the pocket, and controls were matched on pLDDT and hull
depth but never on basic-residue count, while the candidates come from a screen that
rewards basic pockets.

Each fix changes the method, so this is a new pre-registration, not an amendment.

## Three changes, each closing one defect

### 1. Chebyshev acceptance, with provably lossless pruning

Acceptance is on the **worst** matched anchor, not the mean: a placement is admissible when
**every** matched anchor lies within **ε = 2.5 Å** of its partner. Pruning then runs at
**τ = 2ε = 5.0 Å**, and `pruning_complete` in `formal/RequestProject/Pruning.lean` proves
that *no* admissible correspondence is discarded. The search becomes exhaustive by
construction rather than by assertion.

This is also the better scientific criterion. A transplant with one anchor 5 Å out is not a
good fit merely because three others absorb it in the average.

The score is the **lowest achievable worst-anchor deviation** over all eligible templates
and correspondences, among placements that are clash-free. Lower is better; as in study N it
is reported negated so that larger means better. Queries with no admissible placement are
censored at ε.

### 2. Burial required of the *placed* ligand, on a threshold taken from crystal IP6 sites

A placement counts only if the transplanted ligand is actually buried, measured the same way
the project measures every crystal copy — `cryptic_ip.validation.burial_metrics.compute_ligand_burial`,
the function that produced the census's own `relative_sasa` column.

**The threshold is fixed here, from the reference structures, before any candidate is
scored.** Among the eligible crystal IHP copies, relative SASA separates the project's
burial classes as:

| class | copies | max relative SASA |
|---|---|---|
| cryptic | 21 | 0.1044 |
| semi-cryptic | 43 | 0.2363 |
| surface | 158 | 0.8664 |

A placed ligand is **buried** when its relative SASA is **≤ 0.2363** — the largest value
observed among the 64 crystal IP6 copies this project classifies as cryptic or semi-cryptic.
The stricter cryptic-only cut (≤ 0.1044) is reported alongside as a sensitivity, and cannot
change a decision.

### 3. Controls matched on basic-residue count

For each candidate, a **pool of 5** pLDDT- and depth-matched low-ranking proteins is drawn
by the existing rule (`scripts/triage.py::matched_controls`, now with `n_per_candidate`),
and the control kept is the pool member whose **pocket basic-residue count** is closest to
the candidate's, ties broken by accession. Counts are measured from the same AlphaFold models
the fit uses, so the arms are matched on the quantity the score structurally requires.

The count used for matching, and the residual imbalance after matching, are both reported.

## The template library, stated correctly this time

`results/redocking/census.csv` restricted to `comp_id == IHP`, `status == eligible`,
`complete == True`, and **not** flagged as a symmetry contact. That is **222 copies across
137 PDB entries in 17 strict homology groups**.

Study N's plan described this filter as `symmetry_contact == False` and reported 131 copies
in 77 entries and 14 groups. **That was wrong**, and it understated the library that
actually ran. The flag is empty for the **91 cryo-EM copies**, which have no unit cell, so
the symmetry-contact test is *inapplicable* to them rather than failed. Including them is
correct — there is no crystal packing to contaminate a cryo-EM site — and it matters here,
because cryo-EM is exactly where this project expects IP6 to act as a structural cofactor in
large assemblies. The record is corrected in `docs/TEMPLATE_PLAN_AMENDMENT_2.md`.

## Arms, and the leakage rules

Four arms, scored identically: **annotated binders** and their basic-count-matched controls
(which make the guard possible), and the **75 candidates** with theirs.

A template is ineligible for a query when it shares a UniProt accession with it, when its
accession lies in the query's MMseqs2 30 % cluster, or when the query is a PDB-derived binder
in the template's strict homology group. Counts removed by each rule are reported.

## Decisions

Intervals are 95 % bootstrap intervals over **2,000 resamples of whole MMseqs2 30 %
clusters**, seed **20261007**. Per-copy and per-group estimands are both reported. **A
stratum with fewer than 5 clusters yields no evidence.** Holm runs across O1 and O2.

**O1 — calibration guard, and it runs first.** AUC of the fit score separating annotated
binders from their basic-count-matched controls.

* **Passes** if the 95 % lower bound exceeds 0.5.
* Otherwise the method cannot recognise a site it is shown, so it cannot recognise one it is
  not. **The candidate arm is not run at all** — the workflow gates it on this verdict, as
  study N's did — and the study reports the failure and stops.

**O2 — primary, only if O1 passes.** Paired difference in fit score, each candidate minus its
own matched control, cluster bootstrap.

* **better in candidates** if the lower bound > 0; **worse** if the upper bound < 0;
* **no difference** if the whole interval lies within ±**0.25 Å**; **inconclusive** otherwise.

**O3 — the shortlist, reported only if O1 passes.** A candidate is **template-compatible**
when it admits a placement that is simultaneously (a) clash-free, (b) within ε = 2.5 Å on
every matched anchor, and (c) buried at relative SASA ≤ 0.2363. These are **absolute**
criteria read off reference structures, not a percentile of the control arm — which is what
made N3 degenerate. The count of candidates and of controls meeting them are both reported,
so the reader can see the contrast rather than infer it.

O3 is descriptive. It carries no p-value, and the vocabulary is fixed here as it was in study
N: such proteins are **template-compatible**, never *predicted binders*.

## What this study still cannot conclude

* Nothing about affinity, kinetics, or whether the protein meets IP6 in a cell.
* Nothing about specificity between inositol phosphates or against other polyanions. Study C
  showed the descriptors cannot make that call and geometry cannot either.
* Nothing that survives the fact that every candidate structure is a **single apo predicted
  conformer**. A pocket that cannot host IP6 rigidly might host it after rearrangement, and
  this design cannot see that.
* A candidate that passes O3 is a protein whose pocket **can** host a real IP6 site in a
  buried position. That is a hypothesis worth an experiment, and it is not a finding of
  binding.
