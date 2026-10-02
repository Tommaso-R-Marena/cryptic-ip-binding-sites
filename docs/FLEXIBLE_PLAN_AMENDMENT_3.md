# Study J, amendment 3: the per-copy compute cap was choosing the analysed set

Dated 2026-10-01, **before any J1, J2 or J3 estimate has been computed or read**. The only
study J output read so far is the per-copy `ok` / error lines in four shard logs of run
36850024210; no success rate, no paired difference, no interval. This amendment changes a
compute budget, not a scientific constant: the estimands, the 2 Å success criterion, the
±0.05 margin, seed 20261003, `homology_group_strict` resampling, the 5-group floor and
amendments 1 and 2 are all unchanged.

## What the run showed

Run 36850024210 docked each copy in a child process with a 5400-second (90-minute) cap,
three copies to a shard. Across the four shards whose logs were read — 12 copies — the
outcomes were:

| outcome | copies |
|---|---|
| docked | 5 |
| **timed out after 5400 s** | **5** |
| declined by amendment 2 (untypable residue within 8.0 Å of the box) | 2 |

Shard 11 timed out on all three of its copies. The copies that did finish took **30 to 86
minutes**, so the 90-minute cap was not sitting beyond the tail of the cost distribution —
it was sitting inside its body.

## Why that is not acceptable, and not merely wasteful

A cap that cuts 5 of 12 copies would be a problem even if it cut at random. It does not cut
at random. Flexible-receptor docking cost grows with the size of the receptor and the number
of movable side chains, so the copies the cap removes are systematically the **larger
structures with more crowded pockets** — which is the direction of the buried, cryptic sites
this project exists to study. Estimating J1 on the survivors would answer a question about
small, fast receptors and report it as a question about all of them.

This is the same defect as study O's control arm, in a different place: a comparison whose
population is decided by which members happened to be scoreable.

## The change

One copy per shard (272 shards), with the per-copy cap raised from 5400 s to **18000 s**
(300 minutes), inside the shard's existing 19800-second budget and the job's 355-minute
timeout. A copy that still exceeds 300 minutes is recorded as timed out, as before.

The cap is an operational limit, not a pre-registered scientific parameter, and it is being
raised rather than lowered — it cannot manufacture a result, only admit copies that were
previously excluded for cost. It is being set **before** any estimate is computed, and it is
fixed now: it may not be moved again on the basis of what the estimates turn out to be.

## What is still censored, and how it will be reported

Three exclusion routes remain, and the report names each separately rather than folding them
together:

1. **Timed out** at 300 minutes — recorded per copy, with the count in `error_classes`.
2. **Declined by amendment 2** — an untypable residue or a PQR/template hydrogen discrepancy
   within 8.0 Å of the docking box.
3. **PDB2PQR failure** — study A's own class, which it hits on 19 of 272 copies.

If the timed-out count after this change is still a material fraction of the 272, that is a
finding about the affordability of flexible-receptor docking for this ligand and will be
reported as one, rather than being absorbed into J1's denominator. J1, J2 and J3 will state
how many copies they are computed on next to the 272 the plan selected.

## The cost, stated plainly

272 jobs of up to 5 hours, at 20 concurrent, is roughly a day of wall-clock time. Run
36850024210 is cancelled rather than left to spend eleven more hours producing a dataset
whose inclusion is decided by a timer.
