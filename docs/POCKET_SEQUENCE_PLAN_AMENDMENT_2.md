# Study Q9 — amendment 2

**Dated 2026-10-20**, after reading the run of head `923a100`
(`results/pocket_sequence/`). It records what that run said, explains the one
thing its own audit flagged, and tightens that audit. Nothing here changes a
design or adds an arm.

## 1. All three questions now have an answer, and the interesting one is negative

The hydropathy row amendment 1 section 1 disclosed — top-ranked pockets at
`composition_z` +1.02 against bottom-ranked at −0.30, difference +1.317 — **does
not survive power.** Every better-powered form of the same contrast collapses it:

| arm | windows | groups | difference | could have detected | status |
|---|---|---|---|---|---|
| `formation` (pre-registered) | 10 | 9 | **+1.310** | 1.417 | under-powered, as before |
| `formation_matched` | 400 | 70 | +0.098 | 0.310 | exploratory |
| `formation_dense` | 375 | 75 | +0.079 | 0.280 | exploratory |
| `formation_matched_dense` | 400 | 75 | +0.081 | 0.321 | exploratory |
| **`formation_heldout`** (confirmatory) | **400** | **63** | **+0.014** | 0.300 | Holm 1.0000 |

The confirmatory arm is the one that matters: it runs on the 63 resampling
groups the pre-registered rule selected nothing from, so on data no decision in
amendment 1 could have been tuned to, with forty times the windows and seven
times the groups. It gives **+0.014**, with an interval of [−0.198, +0.215] and a
minimum detectable difference of 0.300. The effect is not small; it is absent.

**The power floor did its job.** A difference of +1.31 that clears a corrected p
is what a ten-window arm produces when the estimate is noise, and the plan
refused to call it a finding before any of this ran. Had the floor not been
there, the first run would have reported a sequence predictor of pocket
formation that does not exist.

The other two arms:

* **Type.** `type_dense`, 400 windows across 56 groups: +0.032 (charge) and
  −0.064 (hydropathy), against minimum detectable differences of 0.356 and
  0.428. Null, and powered enough to say so. The architecture class a lining's
  *positions* earn is not anticipated by the arrangement of its letters.
* **Function.** `function_dense`, 129 shells across 16 groups: +0.061 (charge)
  and +0.191 (hydropathy), against 0.596 and 0.934. **Null but not resolved** —
  the observed differences are well under what this design could detect, so this
  arm has not answered its question, and section 5 says what it would take.

## 2. Why the run's own audit said the look did not reproduce

The report says, correctly by its own rule and wrongly as English, that "the
pocket table is therefore not the one that was looked at". The pocket table is
the same table. What changed is a line of amendment 1's implementation: a
window's null seed moved from **its row position in the pockets table** to **a
hash of its own letters**, so that a window's null is a property of the window
rather than of whichever arm happened to select it. The null is therefore a
different draw, on the same windows.

What reproduced exactly:

* the selection — **10 windows across 9 groups**, identical, as it must be,
  since the rule reads ranks and lengths and never a seed;
* the window set the held-out arm excludes — **9 groups**, the same 9.

What moved: the difference, from 1.317 to 1.309855, a gap of **0.00715**. A z
against 2000 shuffles carries Monte-Carlo error of about
`sqrt((1 + z²/2) / 2000)`, so a fresh draw of the whole ten-window arm is
expected to move the difference by about **0.0173**. The observed gap is **0.41
of that** — under half of one standard error, and in the only quantity that a
re-seeding can move.

**So the held-out arm is held out from exactly what amendment 1 says it is**,
and `LOOK_TOLERANCE = 0.001` was simply the wrong unit: a Monte-Carlo quantity
compared against a fixed constant.

## 3. The audit is tightened, and it cannot be self-serving here

The reproduction check now reports, and decides on, three things separately:

* **selection** — windows and groups must match exactly, with no tolerance;
* **the difference** — compared against the Monte-Carlo spread implied by the
  run's own shuffle count and window count, not a fixed constant, with three
  times that spread as the bound;
* **the gap and the spread themselves**, printed, so the reader checks the
  arithmetic rather than trusting the verdict.

It is worth being explicit that loosening a fence here gains nothing: the
confirmatory arm it guards found **+0.014 at Holm 1.0000**. No verdict in this
study turns on the audit's wording. A tolerance fixed after seeing that a
confirmatory arm had fired would be a different matter, and this is not that.

## 4. The pre-registered row moved from Holm 0.037 to 0.080 without its data changing

In the first run the type and function arms had no input, so the plan's primary
family had **two** runnable members instead of four, and Holm corrected over
two. Feeding the function arm made two more runnable, and the same difference —
+1.317 against +1.310 — now corrects to 0.0800.

That is the correction working as specified, not a result changing. It is
recorded because the first run's 0.037 was in part an artefact of two arms being
starved, and anyone reading amendment 1 section 1 beside this run would
otherwise have to work out the discrepancy themselves.

## 5. The plan's window cannot resolve the function or type arm at this size

The yield diagnostics amendment 1 section 6 added:

| source | reason | linings |
|---|---|---|
| pockets | span over the forty-residue cap | 5443 |
| pockets | under the eight-residue floor | 963 |
| contact shells | span over the cap | 114 |
| contact shells | under the floor | 6 |

The plan's `local` window admits **2001 of 8407** pockets and **9 of 129**
annotated contact shells. Nine shells fall in four homology groups, and the
pre-registered type contrast, after length matching, falls in **one**. Both are
below the plan's five-group floor, so both read "no evidence" — which is the
plan's correct verdict and not a result.

An inositol-phosphate contact shell is a roughly spherical shell around a
ligand, so its residues come from several helices and strands at once and its
sequence span is usually far more than forty residues. **That is a property of
the object, not of this dataset**, and it means the plan's contiguous-span
window will never carry the function arm, at any sample size.

Making the function question confirmatory would need a window definition
pre-registered for it — the `dense` window is the obvious candidate and is
exploratory here by declaration — on data disjoint from this run. That is a
separate study with its own pre-registration, and **it is not begun here**,
because the arm that would motivate it is itself unresolved: +0.191 against a
minimum detectable 0.934.
