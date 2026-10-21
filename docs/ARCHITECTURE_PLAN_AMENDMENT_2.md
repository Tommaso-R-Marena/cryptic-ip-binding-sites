# Architecture plan, amendment 2: do the basic residues sit on a face?

Dated 2026-10-13, **before any basic-residue subset has been passed to
`describe_architecture`**. Written after amendment 1 and before its successor
analysis is run.

## Why a second amendment

Amendment 1 reported real candidate linings as structured at 0.33 against 0.04
for size- and span-matched decoys, and said plainly that the comparison is
nearly tautological: a random position set is not a protein, so *any* real
pocket lining beats that null. It named the control that would fix it — other
pockets from the same proteins — and said no committed table carried them.

**That has now been verified rather than assumed, and it holds, for a specific
reason worth recording.** `results/ml_training/validation_pocket_features.csv`
does carry 12,190 non-selected pockets across 136 structures, a median of 83
pockets per structure, with a label marking the 5 ligand sites. It carries
**descriptors only and no residue lists**, so the within-protein control exists
as numbers and not as positions, and the architecture statistics need
positions. The control remains unavailable, now for a checked reason.

## The better design that is available

`results/triage/triaged_candidates.csv` carries, for each of the same 75
pockets, the **basic residues of that pocket's lining** as `basic_positions`.
Checked: it is a subset of `top_pocket_residues` in **75 of 75** rows, median 7
residues against a median 15-residue lining, and 55 pockets have at least 6.

That permits a **paired, within-lining** comparison which controls everything
the decoy comparison could not:

> For each pocket, are its basic residues more phase-coherent, or more
> face-spaced, than a **random subset of the same lining, of the same size**?

Same protein, same fold, same lining, same subset size. The only thing that
varies is *which* residues were chosen, so a difference cannot be protein
identity, pocket selection, lining length or span. It is also the question this
project cares about most: inositol phosphate binds through basic residues, so
whether those residues sit on one face of a secondary-structure element is a
statement about the binding mode rather than about pocket-ness.

## What will be computed, fixed before any result is read

1. **Per pocket**, the architecture of the basic subset, and of
   `WITHIN_DRAWS` random subsets of the same lining with the same size. The
   per-pocket statistic is the **quantile of the basic subset among its own
   matched subsets** — a number on [0, 1] that is uniform under the null that
   basic residues are placed at random within the lining.
2. **Aggregated** over pockets with a group bootstrap on `cluster`, with the
   standing five-group floor, reporting the mean quantile and its interval. A
   mean quantile of 0.5 is the null; below 0.5 means the basic residues are
   *more* structured than chance within their own lining.
3. **Both halves separately**, topology and register, never pooled, because
   they were shown to detect different pocket shapes.
4. **The same for `conserved_positions`**, as the second and last subset
   tested, with Holm across the two subsets times the two halves.

## Decision rules

* The null is a uniform quantile, so the test is whether the mean departs from
  0.5. Reported with a group-bootstrap interval and a permutation p.
* Pockets with fewer than `MIN_SUBSET` basic residues are excluded with a
  counted reason, not scored, because a register null is unresolved below 20
  gap orderings and that is a property of the data.
* **A difference that tracks subset size is not a finding.** Size is matched by
  construction here, and the matched draw is what makes that true, but the
  subset sizes are reported anyway.
* The transfer result stands over all of it: at realistic φ/ψ scatter the
  method detects about a third of what it detects on idealised geometry, so a
  null here is weak evidence of absence and is to be reported with that number
  beside it.

## What this still cannot conclude

* These are candidates, not annotated sites. A result here describes where
  basic residues sit in **candidate** pocket linings.
* It is not evidence of binding. A basic face is a necessary-looking feature,
  not a demonstrated interaction.
* It says nothing about whether the pockets are real pockets; the within-lining
  design deliberately conditions on the lining and therefore cannot.
