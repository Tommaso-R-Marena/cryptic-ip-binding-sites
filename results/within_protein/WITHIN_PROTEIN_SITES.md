# Is the candidate pocket different from its own protein's pockets?

**Run.** `scripts/within_protein.py analyse`, seed 20261020, on pockets re-detected by `within_protein.py pockets`.

**Pre-registered** in `docs/WITHIN_PROTEIN_PLAN.md`, dated before any structure was fetched for this study.

## The control that did not exist

`results/architecture/real/REAL.md` compares 75 candidate linings with random position sets matched on size and span, and says plainly that the comparison is nearly tautological: a random position set is not a protein, so any real lining beats that null. The control that would mean something is *other pockets from the same proteins*, and no committed table carries residue lists for them.

So the control was built. The same detector re-detects every pocket in the same proteins and keeps each pocket's lining residues. Same protein, same fold, same model, same detector, same lining size - what varies is **which pocket**.

The candidate is re-identified by Jaccard overlap against the committed lining (floor 0.5), because pocket indices do not survive a re-run. Median overlap 1.00, lowest kept 0.60.

Siblings match the candidate's lining size to within 2 residues, at least 8 of them and at most 40, sampled with the seed above. Median 18 siblings per protein.

5 of 28 proteins measured.

**Which proteins survive, and why it matters.** The measured candidate linings run 5 to 15 residues (median 8). The exclusions above are almost entirely a size effect: a large lining has few pockets of its own size to be compared with, so the proteins that clear the sibling floor are the smaller-lining end of the table. That is a property of the matched design rather than a result, it was not chosen after seeing anything, and loosening the tolerance afterwards to recover the rest would be exactly the tuning this plan fixed in advance to avoid.

Excluded, by cause:

- 22 — fewer than 8 siblings within 2 residues of the candidate
- 1 — selected pocket not re-identified (recall below 0.3)

## Result

The per-protein statistic is the candidate's **mid-p quantile** among its matched siblings. A smaller architecture p means more structure, so the quantile counts how many siblings are *more* structured than the candidate: **below 0.5 means the candidate is more structured than its protein's same-size pockets**. An unresolved half is unresolved, never a p of 1, and ties are split because permutation p-values are discrete.

| half | proteins | groups | ties | mean quantile | difference from 0.5 | 95% interval | p | Holm | finding |
|---|---|---|---|---|---|---|---|---|---|
| topology | 5 | 4 | 0.00 | 0.279 | -0.221 | [-0.444, +0.038] | 0.1404 | 0.2809 | **no evidence** — 4 groups |
| register | 3 | 3 | 0.00 | 0.328 | -0.172 | [-0.500, +0.038] | 0.4128 | 0.4128 | **no evidence** — 3 groups |

**This arm did not have the power to decide anything, and the table above is not a null result.** the topology half rests on 4 resampling groups and could only have detected a difference of 0.37; the register half rests on 3 resampling groups and could only have detected a difference of 0.47. The project's standing rule is that a stratum with fewer than five groups yields no evidence, and that rule is applied here to the arm as a whole rather than quietly to a footnote. Whatever direction the means point in, this arm neither supports nor refutes the hypothesis.

Holm across 2 primary tests, fixed by the plan before any of it was read.

**Secondary.** Candidates called a helical or strand face: 0.20; their matched siblings: 0.09.

**Residual size.** Matching is to within a tolerance, not exact, so the residual is reported: topology -0.10 residues; register +0.45 residues. A result accompanied by a residual difference in the same direction is confounded rather than a finding.

**By burial class**, declared secondary and never promoted: a stratum with fewer than five groups gives no evidence, which is the project's standing rule and not a judgement about this table.

| class | proteins | groups | topology | register | evidence |
|---|---|---|---|---|---|
| surface | 5 | 4 | 0.273 | 0.328 | no — too few groups |

## How to read a null here

The transfer measurement in `results/architecture/ARCHITECTURE.md` stands over all of this: at the φ/ψ scatter real structures show, the method detects about **a third** of what it detects on idealised geometry. A null here is therefore weak evidence of absence, and the minimum detectable difference in the table is the honest statement of what this design could have seen.

## What this cannot conclude

* These are candidates, not annotated sites.
* A sibling pocket is not a negative. The comparison is to the protein's own background, not to a labelled non-site.
* AlphaFold models are predictions, and pocket detection inherits whatever the prediction got wrong — equally for candidate and sibling, which is the point of a within-protein design.
* Nothing here involves an inositol phosphate.
