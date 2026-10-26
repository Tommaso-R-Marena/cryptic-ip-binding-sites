# Study Q10 — amendment 1

**Dated 2026-10-26**, after two runs of `docs/COMPOSITION_PLAN.md` and **before
the function arm has been read on the code that will report it.**

No decision rule in the plan changes here. This exists because the function arm
has already been seen once, on code that has since been superseded, and a
result that has been looked at must be declared as a look before it is read
again as a result.

## 1. What has already been looked at

The workflow ran twice. The first run (head `720ee40`, run 40231209197)
completed all four arms. Its artefact is **not** the study's result, for two
reasons found afterwards:

* It predates `1ed8f5f`, which replaced the positive control's p-value. That
  run's p is the label-permutation p, measured at the time to reject at 7.25%
  of a nominal 5% against a zero-variance control arm. Its artefact carries no
  `p_rule` field, which is how it was identified.
* It predates `d02c8db`, which split pocket detection out of the contact-shell
  step. In that run the two shared one job, so eight annotated shells were lost
  along with the structures fpocket failed on: **87 shells in 12 homology
  groups**, against **95 in 13** once the steps were independent. The gate's
  own number moved with them, from mean z +4.06 to +4.28.

The second run (head `d02c8db`, run 40250237023) is on current code and is
committed as `results/composition/composition.json`. Its `analyse` job finished at
01:09 and the pocket table it needed was uploaded at 05:00, so its function arm
is **not runnable** and says so.

**The look.** Run 40231209197's function arm — annotated shells against the
other pockets of the same deposited structure, 87 shells in 10 groups — gave:

| statistic | mean z (shells) | mean z (control) | difference | 95% interval | Holm |
|---|---|---|---|---|---|
| `basic_fraction` | +3.96 | −0.01 | **+3.967** | [+3.176, +4.780] | 0.0030 |
| `composition_divergence` | +2.22 | +0.13 | **+2.094** | [+1.392, +2.862] | 0.0030 |
| `nonbasic_divergence` | +2.60 | +0.13 | **+2.478** | [+1.611, +3.326] | 0.0030 |

All three cleared every floor on that run. **These numbers are a look, not a
result.** When the arm reports on current code it will be read as the result,
and this amendment is what makes it impossible to present that reading as a
first look.

## 2. What is not revisited

* The gate, its direction, its floors and the rule that a shut gate makes every
  arm of this study and of study Q9 *instrument not demonstrated*. The gate has
  fired on both runs, at +4.06 and +4.28, and the plan's reading of that —
  Q9's nulls are statements about order, not about sensitivity — stands.
* The Holm families. The committed artefact's primary family has six members
  rather than nine because the function arm did not run; when it runs the
  family returns to the nine the plan named. Nothing is added or dropped.
* Every statistic, the subset null, the window, the size matching, the caps and
  the seed, all unchanged.

## 3. Why the arm is expected to read the same

Nothing between the two heads touches what the function arm measures. The
p-value fix applies to the one-sample control arm only; the function arm is a
two-arm contrast, whose calibration was checked over 400 null pools and left
alone deliberately. The eight recovered shells enter the arm's site side, so
the arm's n should rise from 87 and its groups from 10.

If it reads materially differently, that is itself the finding to report, and
the two tables above are what it would be compared against.

## 4. Limits, unchanged

A compositional difference between a contact shell and its own chain is not a
discovery, on either run. The function arm is the one arm of this study whose
control is hard enough to mean something, and even there it says that an IP
contact shell is more basic and more compositionally divergent from its chain
than the other pockets of the same structure — which is textbook, and is the
instrument check this study was built to be. It is not a sequence-only
predictor of anything and this study does not attempt one.

---

## 5. Postscript, added 2026-10-26 after the arm reported

Nothing above was changed. This section records what the function arm did on
current code (run 40275846268, head `66b400a`), against the look section 1
declared, because a declared prediction that is never checked is not a
discipline.

| statistic | the look (87 shells, 10 groups) | the result (95 shells, 10 groups) |
|---|---|---|
| `basic_fraction` | +3.967 [+3.176, +4.780] | **+3.994** [+3.239, +4.747] |
| `composition_divergence` | +2.094 [+1.392, +2.862] | **+2.097** [+1.401, +2.845] |
| `nonbasic_divergence` | +2.478 [+1.611, +3.326] | **+2.519** [+1.658, +3.358] |

All three are findings, all at Holm 0.0030, as on the superseded run. Section 3
said the arm was expected to read the same and it does; the largest move is
0.041 in `nonbasic_divergence`, well inside every interval.

**One prediction in section 3 was wrong.** It said the arm's "n should rise
from 87 and its groups from 10". The n rose to 95; the groups did not move. The
eight recovered shells fell in homology groups the arm already had, so the
function arm still rests on **ten** groups while the gate rests on thirteen.
That matters more than the reproduction does: the arm with the only control in
this study worth anything is also the arm with the fewest groups, and ten is
twice the plan's floor of five rather than comfortably above it.

Every `formation` and `type` row is identical to the run committed before this
one, to the last digit, which is the property the workflow's two tests exist to
hold: adding the expensive table changes the function arm and nothing else.
