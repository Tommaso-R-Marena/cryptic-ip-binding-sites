# Study Q8: the within-protein control, built rather than wished for

Dated 2026-10-20, **before any structure has been fetched for this study and
before any sibling pocket has been passed to `describe_architecture`**.

## The hole this fills

Two amendments to the architecture plan name the same missing control, and
`docs/DETECTOR.md` repeats it:

> The real-lining comparison in `results/architecture/real/REAL.md` sets 75
> candidate pocket linings against random position sets matched on size and
> span. A random position set is not a protein, so **any** real pocket lining
> would beat that null. The control that would mean something is *other pockets
> from the same proteins*, and no committed table carries residue lists for
> them.

That last clause was checked, not assumed:
`results/ml_training/validation_pocket_features.csv` carries 12,190 non-selected
pockets across 136 structures but descriptors only, and the proteome-screen
shards that did carry `pocket_residues` are not in this repository.

**The control does not exist as a file. It can be made.** The same detector that
produced the candidates will re-detect every pocket in the same 75 proteins and
keep the residue list for each. Nothing new is invented: `ProteinAnalyzer`
already returns pockets and their lining residues, and the AlphaFold client
already fetches the models.

## The question, fixed before any of it runs

> Within one protein, is the **candidate** pocket's lining architecture
> different from that protein's **other** pockets of the same lining size?

Same protein, same fold, same model, same detector, same lining size. What
varies is which pocket. That is the comparison the decoy arm could not make.

## Why size matching is the whole design

The candidate pocket was selected by a scorer built on burial, basicity, volume
and conservation. Architecture is not among its features, so the comparison is
not circular by construction — but lining **size** correlates with nearly every
pocket property, and both halves of the architecture decomposition depend on it
(the topological null holds n and span fixed precisely because they are not the
question). An unmatched comparison would therefore rediscover "the selected
pocket is bigger", dressed as architecture.

So the primary comparison uses siblings matched on lining size, and the
unmatched comparison is reported beside it as the secondary that shows what
size alone buys.

## What will be computed

1. **Pocket detection.** For each of the 75 proteins in
   `results/specificity/candidates.csv`, fetch the current AlphaFold model, run
   `ProteinAnalyzer.detect_pockets()` and record, for every pocket, its lining
   residue numbers. The model version actually used is recorded per protein.
2. **Re-identifying the candidate.** Pocket indices do not survive a re-run: a
   newer model release or a different fpocket build renumbers them. The
   candidate is therefore the pocket with the **highest Jaccard overlap**
   against the committed `top_pocket_residues`, and that overlap is recorded. A
   protein whose best overlap is below `MIN_JACCARD = 0.5` is **excluded with a
   counted reason**, because the pocket being compared would not be the pocket
   that was selected.
3. **Matched siblings.** Every other pocket of the same protein whose lining
   size is within `SIZE_TOLERANCE = 2` residues of the candidate's, capped at
   `MAX_SIBLINGS = 40` by a seeded random sample when more qualify. A protein
   with fewer than `MIN_SIBLINGS = 8` matched siblings is excluded with a
   counted reason.
4. **Per protein**, the candidate's **mid-p quantile** among its matched
   siblings, separately for the topology and register halves, with an
   **unresolved half scored as unresolved and never as a p of 1** — the error
   that produced a spurious result in the paired within-lining study and is now
   pinned by its tests.
5. **Aggregated** with a group bootstrap over `cluster` (MMseqs2 30 %), the
   standing five-group floor, and Holm across the two primary tests.

## Decision rules, fixed now

* The null is a uniform quantile, so the test is whether the mean departs from
  0.5. **Below 0.5 means the candidate pocket is more structured than its own
  protein's pockets of the same size.**
* **Primary family:** two tests, topology and register, Holm-corrected.
* **Secondary, reported whatever the primary shows:** the same two quantiles
  against *unmatched* siblings; the fraction of candidates and of siblings
  called `helix_face` or `strand_face`; the distribution of matched-sibling
  counts; and the re-identification overlaps.
* **A difference that tracks lining size is not a finding.** Size is matched by
  construction, the residual size difference is reported, and a primary result
  accompanied by a residual size difference in the same direction is reported
  as confounded rather than as a finding.
* **Power is stated, not assumed.** The minimum detectable difference is
  computed from the same bootstrap, and the transfer ceiling stands over all of
  it: at the φ/ψ scatter real structures show, the method detects about a third
  of what it detects on ideal geometry, so a null here is weak evidence of
  absence and is to be reported with that number beside it.

## What this still cannot conclude

* These are candidates, not annotated sites. A result describes where the
  *selected* pocket sits relative to its protein's other pockets.
* It is not evidence of binding, and no inositol phosphate is involved anywhere
  in the measurement.
* A sibling pocket is not a negative. Some of a protein's other pockets may be
  real sites for something else; the comparison is to the protein's own
  background, not to a labelled non-site.
* AlphaFold models are predictions. Pocket detection on a predicted structure
  inherits whatever the prediction got wrong, equally for candidate and sibling.
