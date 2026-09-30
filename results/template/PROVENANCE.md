# Study N: the two primary tests pass and the shortlist is degenerate

**Run.** Template fit run 36743030815, commit `549d08e`, branch
`claude/ip-binding-studies-26p6u4`. Report job 109987757920, guard job 109986231214.
Templates fetched 77 RCSB entries; 5 shards per binder arm, 10 per candidate arm; 276
fit records. Plan: `docs/TEMPLATE_PLAN.md`, pre-registered before any template was
extracted.

**Extraction.** `scripts/extract_log_block.py read <log> TEMPLATE_JSON`, then a script
check that every member's `score` is exactly the negated `rmsd` and that the members are
ranked. No figure was retyped. The artifact itself is not reachable from the analysis
sandbox (blob storage returns 403 through the proxy).

## What passed

- **N1, the calibration guard: pass, weakly.** AUC **0.599 [0.546, 0.652]**, p = 0.001,
  62 annotated binders against 63 matched controls in 93 clusters. The lower bound clears
  0.5, which is the pre-registered bar, and 0.60 is weak discrimination. The guard says
  the method is better than chance at recognising a real IP6 site. It does not say it is
  good at it.
- **N2: better in candidates**, +**1.453 Å [1.056, 1.834]** per group over 75 pairs in 65
  clusters, Holm p = 0.0. Per copy, +1.420 [1.051, 1.793].

## What that does not mean, and the reader should see this before the shortlist

**N3 is degenerate and is not reported as a shortlist.** The pre-registered threshold is
the 95th percentile of the control arm's scores. That percentile came out at **−4.0**,
which *is* the censoring floor: at least 95 % of the 75 controls (≥ 71 of them) received
**no admissible placement at all**. So "above the 95th percentile of controls" reduces to
"was not censored", and the 35 candidates it selects are simply every candidate that got
any clash-free placement. That is not a discriminating cut, and calling those 35 proteins
template-compatible would be reporting an artifact of the threshold rule.

**The placed ligand is not buried.** Across those 35 placements the buried fraction has
**median 0.000**, maximum 0.278, and is **exactly zero for 22 of 35**. A rigid transplant
that must avoid clashing tends to find room on the *surface*. The study therefore does not
demonstrate that any candidate pocket can host IP6 in a buried site, which is the only
thing the project cares about.

**N2 most likely measures censoring, not fit quality.** Because controls are almost all
censored at −4.0 and candidates are often not, the +1.45 Å difference is largely the
difference between "a placement exists" and "none does", not between a good fit and a
worse one. There is a plain confound for that difference: the fit score requires at least
**4 basic residues** in the pocket, and the matched controls were matched on pLDDT and hull
depth only — never on basic-residue count — while the candidates were *selected by a screen
that rewards basic pockets*. A control with three lysines cannot be scored at all.

## Status

Study N's primary tests are reported as run. **No candidate shortlist is published from
it**, and the Tier I/II list is unchanged by this study. The honest one-line summary is:
a rigid IP6 transplant lands on candidate pockets more often than on matched controls,
almost entirely on the surface, by a margin that is probably explained by how many basic
residues each pocket has.

The diagnostics that would settle it are recorded in
`docs/TEMPLATE_PLAN_AMENDMENT_1.md` and are explicitly post-hoc: they were written after
these results were read and can support no hypothesis test.
