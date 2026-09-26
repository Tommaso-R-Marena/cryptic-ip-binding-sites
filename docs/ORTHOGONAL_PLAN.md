# Orthogonal sequence evidence on the candidates: plan

Pre-registered on 2026-09-26, before metapredict or SHARK was run on any candidate.

**What was known when this was written.** Studies C and H are complete and read, so this
is a further **filter on a list already in hand**, like study H, not a test of the screen.

- 75 candidates; 24 have a conserved basic pocket, 30 are explained by another ligand,
  2 are not conserved and 19 could not be evaluated for want of UniRef50 homologues.
- Study H's conservation is **enriched** against matched controls (+0.646 [0.485, 0.808]),
  and its filter keeps 0.714 [0.524, 0.905] of annotated binders.
- Every *Dictyostelium* candidate is explained or unevaluable, so study H's shortlist is
  human and yeast only.

## Why these two tools, and what each is for

Study H's evidence has two weaknesses that sequence tools can address independently of
anything the screen or the docking used.

1. **A pocket predicted in a disordered region is an artefact.** AlphaFold pLDDT is a
   structural proxy for this, and study H's controls are pLDDT-matched, but pLDDT and
   sequence-based disorder are not the same measurement. **metapredict** (v3) gives a
   per-residue disorder score from sequence alone, so it can flag a pocket that pLDDT
   scored confidently but that sits in a region predicted to be disordered.
2. **MAFFT conservation assumes the alignment columns are right.** For divergent
   homologues that assumption is weak, and it is why 19 candidates are unevaluable.
   **SHARK** (`bio-shark` 2.0.6, Toth-Petroczy lab) is alignment-free:
   - `shark-capture` finds conserved short motifs across a divergent set without a global
     alignment, giving a second opinion on whether a pocket's basic residues sit in
     something conserved;
   - `shark-dive` detects remote homology by k-mer similarity, so it can test whether the
     unevaluable candidates have homologues that UniRef50's 50 % identity clustering
     missed.

## Data

- **Proteins and roles.** Exactly study H's three arms, so the comparisons are
  commensurable: 75 candidates, their pLDDT- and depth-matched controls, and the
  annotated IP binders as positive controls.
- **Pocket residues.** `top_pocket_residues` as in studies C and H.
- **Sequences and homologues.** The same UniProt sequences and UniRef50 member sets that
  study H used, fetched through the same code.
- **SHARK-dive target set.** A seeded sample of **3,000** of the screened human, yeast and
  *Dictyostelium* proteins (seed 20261002), not the whole proteomes: SHARK-dive is
  quadratic in sequence count and the full set does not fit the job budget. A remote
  homologue found there is evidence the protein is not orphan; it is **not** a substitute
  for orthologues across many species, and I4 is reported as such, against the sample.

## Metrics (fixed here)

- **Pocket disorder**: the mean metapredict v3 disorder score over the pocket's residues.
  A residue is disordered at a score > 0.5, and a pocket is **disordered** when its mean
  exceeds 0.5.
- **Motif support**: with `shark-capture` run on the protein's homologue set
  (`--k_min 3 --k_max 8`, defaults otherwise, and at most **40** homologues per protein —
  a runtime cap fixed here, since capture is superlinear in set size), the fraction of the
  pocket's basic positions that fall inside a captured consensus motif's match in the
  target sequence. A match span is SHARK's own 1-based inclusive `start`/`end`.
- **Remote homologues**: the number of target sequences `shark-dive` returns for a query,
  counted at SHARK's own reported similarity ranking, top 100 retained.

## Decisions

- **I1, primary (the only hypothesis test).** The paired difference in mean pocket
  disorder, candidates − matched controls, resampled together over MMseqs2 clusters
  (2,000 resamples, seed 20261002).
  - **candidates more ordered** if the 95 % upper bound < 0.
  - **no difference** if the interval lies inside ±0.05.
  - **inconclusive** otherwise. **Not evaluable** below 5 clusters.
  - Because the controls are already pLDDT-matched, this asks specifically whether
    sequence disorder carries information *beyond* AlphaFold's confidence. A null here is
    a real and reportable answer, not a failure.
- **I2, disorder QC and its calibration (descriptive).** The share of candidates whose
  pocket is disordered, and the same share among annotated binders. Any candidate with a
  disordered pocket is **demoted from the shortlist**. The filter is declared informative
  only if at most **0.25** of annotated binders are flagged; otherwise nothing is demoted
  on it — the same guard study H used for conservation.
- **I3, motif support (descriptive).** Motif support for candidates against matched
  controls, reported with a group-bootstrap interval. Orthogonal to study H's MAFFT
  conservation, and reported beside it rather than combined into a score.
- **I4, remote-homology rescue (descriptive).** For each of the 19 candidates study H
  could not evaluate, the number of remote homologues `shark-dive` finds. Reported as a
  count only: rescuing a protein from "no homologues" does **not** re-run the
  conservation criterion here, because a within-proteome remote hit is not an orthologue
  set, and doing so would be a post-hoc loosening of study H's rule.

## What this can and cannot support

- It can remove candidates whose pockets are predicted-disordered, and add an
  alignment-free second opinion on conservation.
- It cannot say anything about ligand identity. Study C measured that and found the
  descriptors too weak to change a ranking; no sequence tool changes that.
- It cannot make a candidate a demonstrated IP6 binder. Only the assays already named can.

## Tools explicitly considered and not used

- **BindCraft** designs *de novo protein binders* (diffusion plus ProteinMPNN plus an
  AlphaFold filter, needing a GPU and PyRosetta). IP6 is a small molecule, so the tool
  does not apply to the question "does this protein bind IP6". The one legitimate adjacent
  use — designing a protein binder as a laboratory tool reagent against a validated
  pocket — is downstream of an assay that has not happened, and is not attempted.
- **More AutoDock Vina at a larger budget** is not run: study G established that the
  ceiling is search-saturated, so repeating it would re-derive a finished negative. The
  untested docking variant is a **flexible receptor**, which is pre-registered separately
  in `docs/FLEXIBLE_PLAN.md`.

## Execution and outputs

- **Workflow.** `.github/workflows/orthogonal.yml`: a disorder job, a sharded
  `shark-capture` matrix, a `shark-dive` job, then a report.
- **Code.** `scripts/orthogonal.py`. SHARK steps are best-effort per protein and record
  their error without failing the run; the metapredict step is required.
- **Results.** Printed between `BEGIN_ORTHOGONAL_JSON` and `END_ORTHOGONAL_JSON`,
  extracted into `results/orthogonal/`.
