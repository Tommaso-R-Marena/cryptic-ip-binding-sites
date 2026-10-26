# Study Q10 — amendment 2

**Dated 2026-10-26**, after the function arm reported twice. Both readings have
been seen; this records a reproducibility defect found by comparing them and
changes no decision rule. Nothing here was written before the second reading,
and it is marked as such rather than presented as pre-registration.

## 1. The function arm's control pool is not reproducible across machines

Study Q10 ran twice with the SAME inputs: same census, same structures, same
`within_protein.py`, same fpocket build (`fpocket-4.2.3-h5fd1fdb_2` in both).
Only `cryptic_ip/detect/sequence.py` differed, by a renamed constant and some
comments, which cannot reach a pocket. The two runs disagree:

| | run 40275846268 (committed) | run 40307921444 |
|---|---|---|
| copies measured | 121 | 120 |
| deposited pockets | 9654 | 9352 |
| linings clearing the floor | 5088 | 4924 |
| fpocket exceeded 900 s | 8 | 9 |
| structure not downloaded | 8 | 8 |

The excluded sets are identical but for **`6ZWM_A_2602`**, which finished inside the budget on the first machine and not on the second.

| statistic | committed | second run | finding, both |
|---|---|---|---|
| `basic_fraction` | +3.994 [+3.239, +4.747] | +3.457 [+2.350, +4.460] | yes |
| `composition_divergence` | +2.097 [+1.401, +2.845] | +2.091 [+1.242, +2.975] | yes |
| `nonbasic_divergence` | +2.519 [+1.658, +3.358] | +2.388 [+1.303, +3.422] | yes |

The site side barely moved (+3.924 to +4.080). The CONTROL side moved from -0.071 to +0.623, and the minimum detectable difference with it, 1.137 to 1.615.

**The cause is a wall-clock budget deciding dataset membership.**
`within_protein.py sites` gives fpocket 900 seconds per structure and records a
structure that exceeds it as excluded. That is a cutoff on elapsed time, so
which copies are in the control pool is a property of **how fast the runner
was**, not of the data. The second run's pocket step took over four hours
against the first's 3 h 59 m, and one more structure fell over the line.

## 2. What this does and does not put in doubt

* The **gate is untouched**: 95 shells, 13 groups, mean z +4.28, identical in
  both runs. It compares shells against subsets of their own chains and never
  reads a detected pocket.
* The **formation and type arms are bit-identical** across the two runs. They
  read the candidate-model pocket table, which is built in a different job and
  did not move: 4267 linings both times.
* Every **verdict holds**: four findings, five nulls, in both runs. The arm's
  site side barely moves; what moves is its control, and in the direction that
  makes the finding smaller rather than larger.
* What is in doubt is the **precision** of the function arm, not its sign. An
  interval that shifts by half a z when one of 121 structures drops out is an
  interval whose width is partly an accident of scheduling.

So no claim in `results/composition/COMPOSITION.md` is withdrawn, and the
committed artefact stays as it is: it is the run with more data, 121 copies
against 120, and replacing it with the second run would hide this rather than
record it.

## 3. What changed in the code, and what deliberately did not

Changed: the artefact now carries `n_deposited_copies`, `deposited_copies` and
`deposited_copies_digest`, and the report states in words that the pool is
machine-dependent and that two artefacts are comparable only where the digest
agrees. This defect was invisible in the artefact and visible only by diffing
two job logs; it will not be invisible again. The committed artefact predates
the fields and therefore makes no claim about its own pool, which is correct.

**Not changed, because it is a design decision rather than a bug fix:** the
900-second budget itself, and the rule that exceeding it excludes a copy. The
two defensible repairs point different ways and the choice has scientific
consequences, so it is not being made here on a maintainer's own initiative:

1. **Pin the excluded set.** Commit the list of copies the budget excludes and
   make the driver refuse to proceed when the set it observes differs, so a
   slower machine produces a loud failure to reproduce instead of a quietly
   different control pool. Keeps the dataset fixed; makes CI sensitive to
   machine speed.
2. **Stop the budget deciding membership.** Raise or remove it and let the
   resumable table carry a structure across several attempts, so a copy is
   excluded only when it genuinely cannot be processed. Fixes the dataset
   properly; costs much more CI time and may never terminate on some
   structures.

Option 1 is the cheaper and more honest of the two in the short run, and option
2 is the right answer if the budget is ever to be defended as a scientific
choice. Until one is taken, the function arm's intervals should be read with
section 1 in view.

---

## 4. Postscript, added 2026-10-26: the cache is what currently makes it stable

Run 40348053052, on head `6a10380`, touched nothing the pocket table's digest
covers, so the table was restored from cache instead of detected afresh. Its
pocket step took **one second** against 3 h 59 m and 4 h 20 m cold, and
`--resume` skipped all 137 copies.

The restored table reproduces run 40307921444 **exactly**: 120 copies, 4924
linings clearing the floor, digest `c026eea964496d04`, every top-level field
identical and every arm's difference, interval and verdict identical to the
last 15 significant figures. The only differences are last-bit floating-point
noise in a few means and interval endpoints, which comes from reduction order
across machines and appears between any two runs on this project.

So the position is now precise, and it is two statements rather than one:

* **fresh pocket detection is not reproducible** — section 1;
* **the same table read twice is** — this section.

Which means the cache is not only an economy. It is, at present, the only thing
making the function arm's control pool stable from one run to the next, and
that is not a property anyone designed. **It also expires**: GitHub evicts a
cache that goes unread, so the stability has a shelf life, after which the next
cold run lands on whatever pool that machine's speed produces.

This does not change which artefact is committed. `composition.json` stays the
121-copy run, for the reason section 2 gives, and swapping in a 120-copy run
because it happens to match today's cache would trade one arbitrary pool for
another while implying the newer one is authoritative. Both readings are on
record above, which is the honest state of it.

It does sharpen the choice in section 3. Option 1, pinning the excluded set, is
close to what the cache achieves by accident — so taking it deliberately would
cost little and would survive eviction, which the accident does not.
