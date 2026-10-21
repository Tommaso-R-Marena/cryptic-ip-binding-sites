# Study O: the method improved, the comparison is still confounded

**Run.** Template transplant 2 run 36826844329, commit `2fa4d57`, branch
`claude/ip-binding-studies-26p6u4`. Report job 110263024339. 46 jobs: templates, pool,
10 count shards, select, 10 binder shards, guard, 20 candidate shards, report. 275 fit
records; one protein had no AlphaFold model.

**Extraction.** `scripts/extract_log_block.py read <log> TEMPLATE2_JSON`, then a script
check that every member's `score` is exactly the negated `deviation` and that members are
ranked. No figure retyped.

**Plan.** `docs/TEMPLATE2_PLAN.md`, pre-registered before any fit under this design.

## What improved, genuinely

| | study N | study O |
|---|---|---|
| guard AUC | 0.599 [0.546, 0.652] | **0.667 [0.576, 0.758]** |
| search | pruning unsound under RMSD | **lossless** by `pruning_complete` |
| burial | inherited from the template | **required of the placed ligand** |
| shortlist cut | 95th percentile of controls (degenerate) | **absolute**, from crystal burial classes |
| controls censored | ≥ 71/75 (≥ 94.7 %) | **59/75 (78.7 %)** |

- **O1: pass.** AUC **0.667 [0.576, 0.758]**, p = 0.0, 62 annotated binders against 62
  matched controls in 99 clusters. Better separation than study N's guard, and the gate
  opened on a firmer footing.
- **O2: better in candidates**, +**1.014 Å** [0.849, 1.181] per group over 75 pairs in 65
  clusters; per copy +0.986 [0.823, 1.151]. Holm p = 0.0 for both tests.
- **O3:** **3 of 75 candidates** meet all three absolute criteria, against **0 of 75**
  matched controls.

## Why the shortlist is still not a finding

**The basic-count matching did not do its job.** This is the defect study O was built to
remove, and the balance block says it survived:

| arm | mean basic, query | mean basic, control | gap | censored query | censored control |
|---|---|---|---|---|---|
| candidates | 7.31 | 4.95 | **2.36** | 7/75 | **59/75** |
| annotated | 7.19 | 4.45 | **2.74** | 29/62 | 45/62 |

A pool of five is not enough. Controls are drawn on pLDDT, hull depth and score first, and
within five such proteins there is often none with ~7 basic pocket residues — because the
candidates come from a screen that *rewards* basic pockets, so the two populations differ
in exactly the quantity the fit score requires. Censoring fell from ≥ 94.7 % to 78.7 %,
which is real progress and not a fix.

So **"3 versus 0" cannot be read as biology.** With 59 of 75 controls unscoreable for want
of four basic residues, zero controls passing is close to guaranteed by the imbalance, and
O2's +1.01 Å remains substantially a difference in censoring rather than in fit quality —
the same criticism that applied to study N, reduced but not retired.

**One of the three sits on the threshold.** Q55FR9's placed ligand has relative SASA
**0.23615** against a cut of **0.2363** — a margin of **0.00015**. That is not a robust
pass and should not be presented as one.

**None of the three is cryptic.** All three fail the stricter cryptic cut (≤ 0.1044), so
at best they are semi-cryptic. The project's target is the buried, cryptic case.

## The three, recorded without promotion

| accession | worst anchor (Å) | basic residues | template | placed relative SASA | cryptic |
|---|---|---|---|---|---|
| Q54BM7 | 1.693 | 17 | 9R4I:A:501 | 0.193 | no |
| Q03619 | 2.169 | 5 | 5HP2:A:801 | 0.183 | no |
| Q55FR9 | 2.203 | 7 | 9D5J:A:801 | 0.236 | no (on the cut) |

These are **template-compatible** in the plan's fixed sense: a real crystal IP6 site can be
rigidly transplanted onto the pocket, clash-free, within 2.5 Å on every matched anchor, with
the placed ligand buried. They are **not predicted binders**, the comparison against controls
is confounded as above, and Q54BM7's 17 basic pocket residues make it the kind of pocket a
transplant finds room in for reasons that have little to do with IP6 specifically.

## What would actually settle it

Not more compute — a better control arm. Three options, in increasing honesty and cost:

1. **A much larger pool** (25–50 rather than 5), so basic count can be matched within it.
2. **Match on basic count first**, then on pLDDT and depth within that stratum, inverting the
   current order.
3. **Restrict O2 to pairs where both members are scoreable**, which removes censoring from the
   estimand entirely and measures fit quality alone.

All three change the design, so any of them is a new pre-registration, not an amendment to
this study. None has been run, and no number above may be re-read after choosing one.
