# Can a real IP6 site be transplanted onto a candidate pocket? Study N: plan

Pre-registered on 2026-09-30, **before any template was extracted from a PDB entry and
before any candidate pocket was fitted**. No number in this document was computed; the
census counts quoted below come from `results/redocking/census.csv`, which studies F and G
already produced and which this study does not extend.

## Why this study exists

Studies A–L answered a great deal about the *ranking* and almost nothing about the
*ligand*. The two facts that set this study up are both negative:

* **Docking cannot score these sites.** Study F: top-pose success 0.110 over 242 crystal
  copies, with 230 of 242 failures attributable to scoring rather than search. Study G
  confirmed it: a 16× search budget moved success by +0.019. Re-ranking by a screened
  Coulomb term recovered +0.054 (study G's permutation p = 0.0099) and no more.
* **The descriptors cannot separate IP from other polyanions** well enough to call a
  proteome hit specific (study C), and the properties that *do* separate candidates from
  matched controls — conserved-basic enrichment (study H, +0.646 [0.485, 0.808]) and
  sequence order at matched pLDDT (study I, −0.093 [−0.149, −0.040]) — are properties of
  the ranking, not evidence that IP6 in particular would bind.

So the open question is narrow and structural: **is a candidate pocket geometrically
capable of hosting an IP6 molecule in the way real IP6 sites do?** That question does not
need a scoring function, which is the one component known to be broken here. It needs only
a rigid superposition and a clash test, both of which are exact.

The method is deliberately old-fashioned: take the basic-residue constellation that
coordinates the phosphates in a real IP6 complex, find the best rigid map onto the
candidate's own basic residues, carry the ligand along, and ask whether it lands without
clashing. This is a *necessary-condition* test. Passing it does not make a protein an IP6
binder. Failing it is strong evidence against.

## What co-folding would have been, and why this is not a substitute

A CPU feasibility probe for Boltz-2 (`scripts/cofold_probe.py`) is running in parallel and
is reserved as **study M**. If co-folding turns out affordable, study M is pre-registered
separately and this study stands on its own terms rather than being folded into it. This
study is not a cheaper approximation of co-folding: it answers a strictly geometric
question and makes no energetic claim.

## Data

**Templates.** `results/redocking/census.csv`, restricted to `comp_id == IHP` (myo-inositol
hexakisphosphate) with `status == eligible`, `complete == True` and
`symmetry_contact == False`. That is 131 copies across 77 PDB entries and 14 strict
homology groups. Coordinates come from the RCSB entry itself; no ligand is rebuilt or
minimised. Templates with fewer than 3 anchor residues (below) are dropped and the
dropped count is reported.

**Queries.** Three arms, all scored by the identical procedure:

1. **Annotated binders** — the proteins flagged `annotated` in study C's
   `proteins_combined.csv.gz`, which is what makes the calibration guard possible.
2. **Candidates** — the 75 proteins in `results/triage/triaged_candidates.csv`.
3. **Matched controls** — one per candidate, drawn by the *existing*
   `scripts/triage.py::matched_controls` rule (same organism, `combined` at or below the
   organism median, pLDDT and hull depth within the established tolerances, each control
   used once). The annotated-binder arm gets its own controls drawn by the same rule.

Query structures are the AlphaFold models the screen itself used; the pocket is the
protein's `top_pocket_residues` from the same table. No pocket is re-detected, so this
study cannot quietly redefine what was screened.

**The binder/control clash, settled in advance.** As in study I, a protein drawn as a
matched control that is also an annotated binder counts as a **control only**, and the
clash count is logged.

## Anchors, fits and the score

**Anchor points.** For a template, the anchor residues are Lys, Arg and His whose
side-chain nitrogens (NZ; NE/NH1/NH2; ND1/NE2) come within **4.0 Å** of any ligand
phosphate oxygen. Each such residue contributes **one** point: the mean position of its
own contacting nitrogens. A template keeps at most the **6** anchors nearest the ligand
centroid, to bound the search.

For a query, the anchor points are the same per-residue representative points for every
K/R/H in the pocket residue list, capped at the **10** nearest the pocket centroid.

**Fit.** For every (template, query) pair, every injective map from a 4-subset of the
template anchors to the query anchors is scored by Kabsch superposition of those 4 points.
A placement is **admissible** when no ligand heavy atom lands within **2.2 Å** of a query
protein heavy atom. The pair's score is the **lowest anchor RMSD over admissible
placements**; if none is admissible, or none reaches RMSD ≤ **4.0 Å**, the pair is censored
at 4.0 Å.

**The query's score** is the minimum over all eligible templates. Lower is better;
throughout, the reported statistic is the negated RMSD so that larger means better and the
direction of every interval is unambiguous.

**Template leakage, excluded three ways.** A template is ineligible for a query when
(a) it shares a UniProt accession with the query, (b) its accession lies in the query's
MMseqs2 30 % cluster, or (c) the query is itself a PDB-derived binder in the template's
strict homology group. Counts of templates removed by each rule are reported.

**Secondaries, fixed here.** The same procedure at 5 matched anchors, and the **buried
fraction** of the placed ligand — the fraction of its heavy atoms with at least 60 protein
heavy atoms within 8 Å — are computed and reported as descriptors. Neither can change a
decision.

## Decisions

Intervals are 95 % bootstrap intervals over **2,000 resamples of whole MMseqs2 30 %
clusters**, seed **20261006**. Per-copy and per-group estimands are both reported. **A
stratum with fewer than 5 clusters yields no evidence and is reported as not evaluable.**
Holm correction runs across N1 and N2, this study's two primary tests.

**N1 — calibration guard, and it runs first.** AUC of the fit score separating annotated
binders from their matched controls.

* The guard **passes** if the 95 % lower bound of the AUC exceeds 0.5.
* If it does not, the method cannot recognise a site it is shown, so it cannot recognise
  one it is not. **The candidate arm is then not run at all**: the workflow gates the
  candidate job on this verdict, so no candidate ranking exists to be reported
  selectively. The study reports the guard's failure and stops.

**N2 — primary, only if N1 passes.** Paired difference in mean fit score, each candidate
minus its own matched control, cluster bootstrap.

* **better in candidates** if the lower bound > 0;
* **worse** if the upper bound < 0;
* **no difference** if the whole interval lies within ±**0.25 Å**;
* **inconclusive** otherwise.

**N3 — the per-candidate shortlist, reported only if N1 passes.** A candidate is
**template-compatible** when its fit score exceeds the **95th percentile of the control
arm's** scores. The shortlist is exactly those candidates, each reported with its best
template (PDB id and copy), the matched anchor count, the anchor RMSD, and the placed
ligand's buried fraction. N3 is descriptive: it is a filter applied to an already-published
shortlist, not a new hypothesis test, and it carries no p-value.

## What this study cannot conclude

* Nothing about affinity, kinetics, or whether the protein ever meets IP6 in a cell.
* Nothing about specificity between inositol phosphates, or between IP6 and other
  polyanions — study C already showed the descriptors cannot make that call, and geometry
  alone cannot either.
* A candidate that passes N3 is a protein whose pocket **can** host IP6, which is a
  hypothesis worth an experiment and is not a finding of binding. The reporting language
  is fixed here: such proteins are called **template-compatible**, never *predicted
  binders*.
