# Study N, amendment 1: diagnostics, written after the results were read

Dated 2026-09-30, **after** Template fit run 36743030815 was read. Everything below is
therefore **post-hoc**. None of it is a hypothesis test, none of it can change N1 or N2,
and no p-value or interval from it may be reported as evidence. It exists to explain a
result the pre-registered design could not interpret.

## What the run showed that the plan did not anticipate

N3's threshold — the 95th percentile of the control arm — landed exactly on the censoring
floor of −4.0 Å, because at least 71 of 75 matched controls admitted no placement at all.
The pre-registered rule therefore selects "not censored" rather than "fits well", and the
buried fraction of the surviving placements has median 0.000. Both facts are recorded in
`results/template/PROVENANCE.md`.

## The confound this exposes

The fit score is undefined for a pocket with fewer than 4 basic residues: there is nothing
to map the template's four anchors onto. Controls were matched to candidates on organism,
`combined` score below the organism median, mean pLDDT within 5 and hull depth within 3 Å —
**never on the number of basic residues in the pocket**. The candidates, meanwhile, were
selected by a screen whose score rewards buried basic pockets. So the arms differ in the
one quantity the score structurally requires, and N2's +1.45 Å may be no more than that.

## Diagnostics to compute, all descriptive

1. **Split the censoring.** Report, per arm, how many queries were censored for *too few
   basic residues* versus *no admissible placement within the ceiling*. These are already
   distinguished in each fit record's `reason` field and were simply not tabulated.
2. **Basic-residue counts per arm.** The distribution of `n_query_anchors` for candidates
   and for controls, which measures the confound directly.
3. **Restricted comparison.** N2 recomputed over only those pairs where *both* members have
   at least 4 basic pocket residues, so censoring cannot drive it. Reported as a
   descriptive difference with an interval, explicitly not as a test.
4. **Burial of the placement.** The buried-fraction distribution for both arms, and how
   many placements exceed any sensible burial threshold at all.

## What would make this a real test again

A control arm matched on basic-residue count as well as pLDDT and depth, and a burial
requirement on the *placed* ligand rather than on the template it came from. Both are
changes to the design, so they belong to a **new pre-registration** written before that
arm is drawn — not to this study, whose results are now read.
