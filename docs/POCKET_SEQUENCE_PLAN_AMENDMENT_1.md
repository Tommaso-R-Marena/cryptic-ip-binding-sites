# Study Q9 — amendment 1

**Dated 2026-10-20**, after the first run of `docs/POCKET_SEQUENCE_PLAN.md` and
**before any result of anything this amendment adds has been read.**

## 1. What has already been looked at

The pre-registered workflow ran once (head `b2d6b97`). Two of its three arms
produced nothing at all, and the third produced ten windows. In full, so that
nothing below can be read as having been chosen in ignorance:

| arm | statistic | windows | groups | mean z (top) | mean z (bottom) | difference | Holm | could have detected | finding |
|---|---|---|---|---|---|---|---|---|---|
| formation | `helix_charge_periodicity` | 10 | 9 | −0.34 | −0.29 | −0.058 | 0.900 | 1.285 | no |
| formation | `helix_hydropathy_periodicity` | 10 | 9 | +1.02 | −0.30 | **+1.317** | **0.037** | 1.406 | **no — under-powered** |

The hydropathy row clears Holm and its interval excludes zero, and it still is
not a finding, because the plan's power floor says an observed difference
smaller than the minimum detectable difference decides nothing. That rule was
fixed before any data and it is not being revisited here. **The number above is
a look, not a result**, and everything in section 3 is written in the knowledge
of it.

The type arm skipped all eight of its statistics with `no group has both arms`
and the function arm reported zero site windows. Neither was a property of the
data: the workflow never passed `--classes` or `--sites`, so no class was ever
computed and no annotated site was ever read. **Nothing about the type arm or
the function arm has been seen.** The pre-registration of those two arms is
therefore intact, and section 5 only records the inputs being supplied.

## 2. Why the formation arm had ten windows out of two thousand

Two thousand and one pocket windows survived, and the arm used ten of them.
The cause is a collinearity the plan did not foresee:

* fpocket ranks pockets by volume, so a top-ranked pocket has a large lining.
* The plan's `local` window is the **contiguous stretch the lining spans**,
  capped at forty residues, and refused outright when the span exceeds the cap.
* So a top-ranked pocket's window is usually refused, and the few that survive
  sit against the cap while the bottom-ranked survivors sit near the floor.
* The plan then matches the two arms on window length to within two residues —
  and in this dataset rank *is* length, so almost no pair matches.

**That is itself a result about the design space, independent of any signal:**
an absolute-rank contrast on span-limited windows has almost no support,
because the covariate the design must control for is the treatment. It is
reported as such.

## 3. The formation arm is re-specified, and the re-specification is exploratory

The pre-registered arm stays exactly as it is, reported unchanged under the
name `formation`. Three arms are added:

* **`formation_matched`** — within each protein, its surviving windows are
  ordered by fpocket rank and split at the median; the better-ranked half is
  length-matched against the worse-ranked half under the plan's own two-residue
  tolerance. Every protein with two surviving windows contributes, and rank is
  no longer collinear with length because the split is taken *within* the
  window population that exists. **Exploratory. In no Holm family. Never a
  finding, whatever it shows** — it was specified after the look in section 1.
* **`formation_heldout`** — the identical rule, restricted to the resampling
  groups that the *pre-registered* rule selected nothing from. The script
  computes that set by running the pre-registered rule first, so it is defined
  by the plan rather than by me, and it is disjoint from every window whose
  value section 1 reports. **This arm is confirmatory**, in its own Holm family
  over the plan's two primary statistics, with the plan's five decision rules
  applied unchanged. It is the only new test that may become a finding.
* **`formation_dense`** — the pre-registered rule on the window of section 4.
  Exploratory.

A protein contributes at most **twelve** matched pairs, chosen with the study's
seed, so one large protein with two hundred pockets cannot carry an arm. Any
arm that is not per-protein contributes at most **four hundred** pairs. Both
caps are cost ceilings, fixed here, before any of these arms has run.

## 4. A second window definition, declared exploratory

The `local` window refuses a lining whose span exceeds forty residues, which is
most of them. Added beside it, never replacing it:

* **`dense`** — the forty-residue contiguous stretch of the chain that contains
  the most of the lining's positions, ties going to the lowest start. It is
  complete and contiguous by construction, like `local`, so a period still
  means what it means; it reads only positions, so it cannot see the letters it
  will be scored on; and **every such window is exactly forty residues long**,
  so the length confound the plan worries about is zero by construction rather
  than by matching.

It is a different object from the plan's window and it is **exploratory on
every arm**. No `dense` row is ever a finding and no `dense` row enters a Holm
family. Where a `local` arm is unresolved for want of windows, the `dense` row
beside it says whether the shortage was the only thing in the way.

## 5. The type and function arms, as planned

Neither arm's design changes. The first run simply never fed them, so two data
steps are added:

* **classes** — `describe_architecture` on each pocket's lining positions, 800
  permutation draws, seed fixed at the study's seed, emitting
  `uniprot_id,pocket_id,element`. The class reads positions only, never
  letters, which is what makes Q9b's two sides independent.
* **sites** — each annotated inositol-phosphate copy's contact shell, computed
  from its deposited structure with the census filter the within-protein study
  uses, and **no pocket detection at all**, because Q9c needs the shell and the
  chain sequence and nothing else.

## 6. Window-yield diagnostics

Every refusal of a window is counted by reason — span over the cap, a gap in
the model, under the floor, a non-standard letter — and the counts are reported.
This is reporting, not a design change. Had it been there the first time,
section 2 would not have needed a diagnosis.

## 7. Performance, which is not a method change

`periodicity_power`, `band_periodicity` and `dominant_period` now build their
cos/sin bases once per window length and period grid instead of once per probe
per draw, and solve the two-column least-squares fit through the pseudo-inverse
of its Gram matrix, which is the same minimum-norm solution `lstsq` returns.
One window's eight statistics at two thousand shuffles fell from 35.3 s to
0.69 s. Agreement with the least-squares reference over four hundred random
series at nine periods each is 1.1e-15, pinned by a test. **No committed number
changes**, and the budget this frees is what pays for section 3.
