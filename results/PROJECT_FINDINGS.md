# What this project establishes

This page is **generated** by `scripts/project_findings.py` from the committed artefacts below. It recomputes nothing, no number on it is retyped, and a claim its artefact stops supporting fails the render rather than standing.

It exists because `results/publication/RESULTS_SUMMARY.md` is generated from the control benchmark, the validation dataset, the ML benchmark and the yeast pilot alone, so it predates and does not cover most of the work. **Where the two differ, this page is the current one.**

| line | artefact | command |
|---|---|---|
| controls | `results/publication/validation/positive_controls.csv` | `scripts/run_publication_package.py` |
| negatives | `results/publication/validation/negative_controls.csv` | `scripts/run_publication_package.py` |
| dataset | `data/validation/ip_binding_validation_dataset.csv` | `scripts/build_validation_dataset.py` |
| specificity | `results/specificity/specificity.json` | `scripts/specificity.py run` |
| candidates | `results/specificity/candidates.csv` | `scripts/specificity.py run` |
| architecture | `results/architecture/real/real.json` | `scripts/real_architecture.py run` |
| within | `results/within_protein/within_protein.json` | `scripts/within_protein.py analyse --arm candidate` |
| sites | `results/within_protein/within_protein_sites.json` | `scripts/within_protein.py analyse --arm site` |
| order | `results/pocket_sequence/pocket_sequence.json` | `scripts/pocket_sequence.py analyse` |
| gate | `results/composition/validation.json` | `scripts/pocket_composition.py validate` |
| census | `results/redocking/census.csv` | `scripts/redocking.py census` |
| coevolution | `results/coevolution/coevolution.json` | `scripts/coevolution.py report` |
| revision | `manuscript/RESULTS_DRAFT.md` | `hand-written; its revision note is quoted, not re-encoded` |
| composition *(not yet reported)* | `results/composition/composition.json` | `scripts/pocket_composition.py analyse` |

Pre-registered plans are in `docs/`, each with its dated amendments. The lining-architecture line has its own detailed page at `results/architecture/FINDINGS.md`.

## 0. Which artefacts are superseded, and by what

`results/publication/` holds outputs from before several pipeline defects were found. The project recorded the corrections in the working draft rather than by regenerating those outputs, so **the stale numbers are still in the files** and a reader meeting them first will take them as current. Three of them are still live in committed artefacts today: `results/yeast_pilot/yeast_pilot_hits.csv` carries the withdrawn candidate, with the misread input visible beside the true one in the same row; the validation dataset's burial counts are the per-entry ones; and the control panel's separation was measured on holo structures.

This is the manuscript's own revision note, quoted and not re-encoded, from `manuscript/RESULTS_DRAFT.md`:

> **Revision note (reconciliation with PR #41).** This draft was written from
> outputs of the pipeline before several defects in it were corrected. Each
> claim below was re-checked against a measurement, and five were changed:
>
> 1. **The yeast candidate is withdrawn.** LEU1's composite of 0.955 depended on
>    a mislabelled input: the scorer's "depth" term was reading fpocket's *mean
>    local hydrophobic density* (32.1 for this pocket) as a distance, and awarding
>    full credit for anything over 15. Re-running the scorer this draft used,
>    changing only that input to the true burial depth (2.84 Å), gives **0.711** —
>    below the 0.75 threshold. Re-screening the identical 499 proteins under the
>    corrected scorer yields no candidates.
> 2. **Burial depth discriminates weakly, not cleanly.** The "~3× difference"
>    came from three structures. Across 129 deposited entries depth separates
>    buried from exposed sites with AUROC 0.72 (95 % CI 0.58–0.84); enclosure
>    does so with 0.99.
> 3. **The benchmark is not 133 / 2 / 1.** Those counts came from ligand SASA
>    summed over every copy in an entry, which scales with copy number. Measured
>    per copy it is 89 surface, 26 semi-cryptic, 14 cryptic, 6 crystal artefacts.
> 4. **The ROC-AUC ≈ 0.5 was an artefact of the labels**, not evidence that the
>    task is under-powered.
> 5. **The tier-1 separation of 0.581 was measured with the ligand bound.**
>    Scored as a proteome target is scored - ligand removed - it is 0.081 on the
>    crystal structures and 0.069 on the AlphaFold models (see below).
>
> One claim is independently confirmed: HDAC1's inositol phosphate is surface-
> exposed. Evidence for each change is in `docs/METHODS.md`,
> `scripts/burial_survey.py` and the `yeast-rescreen` workflow.

**Everything in sections 1 onward is read from artefacts that postdate those corrections**, except where a sentence says otherwise and points back here.

## The short answer

* **The detector's own screening filter rejects the case it was built around.** Scored the way a proteome target is scored — ligand removed — ADAR2's site separates from the PH-domain negatives only narrowly, and the plan's strict filter (score ≥ 0.75, lining SASA ≤ 10 Å²) **rejects ADAR2's site on both gates**. What does separate it is depth inside the convex hull. Section 0 has the supersession; the 1 of 3 pass rate in `results/publication/validation/positive_controls.csv` was measured with the ligand bound and is not the number to quote.
* **Buried examples are scarce, and the counts to use are the per-copy ones.** The control benchmark's validation set reads 1 cryptic of 136 when ligand accessibility is summed over every copy in an entry, which scales with copy number; the revision note in section 0 corrects that to a per-copy count. The redocking census is per copy, and holds **15 cryptic annotated copies across 4 homology groups** among its 269 selected copies (45 cryptic rows in 662), plus 54 semi-cryptic across 14 groups — and that is what bounds the ground-truth arms of every study below. Small, but not one.
* **A learnable inositol-phosphate-versus-other-polyanion signal** at 0.828 [0.756, 0.886] ROC-AUC against 34 annotated binders among 33084 unseen proteins — which **does not beat its own same-protein baseline** (-0.015 [-0.044, 0.015], Holm 0.321, *not supported*).
* **75 candidate proteins that are hypotheses, not findings.** The best precision at any rank is 0.25, and it is 0.00 at rank 1.
* **Two structural decompositions that separate real linings from random ones and not from other pockets.** The architecture line (positions) and the sequence line (letters) account for most of the project's statistical work. Against random position sets the architecture contrast is large; against the protein's own other pockets, matched on size, both lines are null. Sections 2 to 4.

Nothing on this page is evidence that any candidate binds an inositol phosphate. The deposited cryptic copies above are structures with a buried inositol phosphate, classified by this project's own burial metric; **one of them, ADAR2, is the case with published evidence that the protein cannot fold without the cofactor**, and it is the case the detector was built around.

## 1. Can the score tell an IP site from another polyanion site?

Pre-registered in `docs/SPECIFICITY_PLAN.md`. The gate opened on the `no_supergroup` variant at a held-out point of 0.812, after the all-data variant was ruled **not evaluable** because one homology group held more than 40% of the IP pockets — which is the kind of thing that makes an impressive number worthless, and was caught before the number was used.

| test | ROC-AUC [95%] | Holm | decision |
|---|---|---|---|
| L1, same proteins | 0.844 [0.771, 0.902] | — | — |
| L3a, combined | 0.828 [0.756, 0.886] | 0.000 | **supported** |
| L3b, combined − L1 | -0.015 [-0.044, 0.015] | 0.321 | **not supported** |

**The third row is the one that matters and the easiest to skip.** The combined score ranks annotated binders well above chance, and it does **not** rank them better than the same-protein baseline does. Everything the pipeline adds beyond that baseline is, on this evidence, unmeasurable.

Of the 75 candidates, 30 are explained by an existing annotation and 45 are not. The unexplained ones are the hypotheses; they are also the ones with no independent support of any kind.

## 2. The lining-architecture line: positions

Pre-registered in `docs/ARCHITECTURE_PLAN.md` with five dated amendments, and set out in full at `results/architecture/FINDINGS.md`. It factors a pocket lining into **topology** (the gap multiset) and **register** (the arrangement within it), and tests each against its own permutation null.

75 of 75 real candidate linings were measurable. The method properties are real and calibrated, and the detailed page carries them. What they amount to on real data is a per-lining verdict:

| verdict | linings |
|---|---|
| none | 50 |
| clear | 15 |
| weak | 10 |

Against size- and span-matched **decoy** position sets the contrast is large: 33.3% of real linings resolve on one half or the other against 4.3% of 375 decoys. That is the architecture line's strongest real-data number, and the line says in its own plan why it is weaker than it looks: **a random position set is not a protein**, so the comparison is close to tautological — a real lining is contiguous, secondary-structure-bound and solvent-excluded, and a random set of positions is none of those.

**That is exactly why section 3 exists**: it replaces the decoy with the protein's own other pockets, matched on lining size, so protein, fold, model, detector and size are all fixed and only the choice of pocket varies. On that comparison the line is null.

Rates are in the artefact (`results/architecture/real/real.json`) and on the detailed page. The transfer ceiling stated there — the method's agreement with ideal geometry falling as backbone jitter and bulges are introduced — bounds this whole line and is the reason a null here is weak evidence rather than evidence of absence.

## 3. The control that did not exist, built — and null in both arms

Pre-registered in `docs/WITHIN_PROTEIN_PLAN.md`. The architecture line's own comparison sets a candidate lining against random position sets, and says in print that a random position set is not a protein. The control that means something is **other pockets of the same protein, matched on lining size** — and no committed table carried their residues, so this study built it.

| arm | half | proteins | groups | mean quantile | difference | Holm | could have detected |
|---|---|---|---|---|---|---|---|
| candidate | topology | 36 | 33 | 0.434 | -0.066 | 0.457 | 0.152 |
| candidate | register | 36 | 32 | 0.538 | +0.038 | 0.481 | 0.150 |
| annotated site | topology | 27 | 8 | 0.336 | -0.164 | 0.207 | 0.245 |
| annotated site | register | 27 | 8 | 0.592 | +0.092 | 0.379 | 0.308 |

Below 0.5 means more structured than the protein's own same-size pockets. **Nothing is a finding, and in every row the difference observed is smaller than the difference the design could have detected** — so this is weak evidence, not evidence of absence.

The annotated-site arm decides how to read the other: it runs the identical statistic on deposited inositol-hexakisphosphate sites, so it could have fired, and a failure there would mean the candidate arm's null is about the instrument rather than about candidates. It does not separate real sites from their own proteins' pockets at the power available.

One thing it established for free: re-detection reproduced the committed candidate linings at a median Jaccard overlap of 1.00, so the candidate table is reproducible from current models with current fpocket — which nothing had checked.

## 4. The sequence line: letters, in order

Pre-registered in `docs/POCKET_SEQUENCE_PLAN.md` with two dated amendments. A window of sequence factors as **composition** (which letters) **⊕ order** (in what arrangement), the same way a lining factors into topology and register. This line takes the quotient by composition and asks what the ordering carries.

**The result worth carrying off this line is a correction, not a discovery.** An underpowered form of the formation contrast gave +1.31 on 10 windows at a corrected p the plan refused, because the observed difference was smaller than the minimum detectable one. The re-specified rule, on the resampling groups the pre-registered rule selected nothing from, gives **+0.014** on 400 windows across 63 groups, interval [-0.198, +0.215], Holm 1.0000. Not a small effect — an absent one.

Without the power floor, fixed before any sequence was read, this project would have reported a sequence predictor of pocket formation that does not exist.

Architectural type is null and powered enough to say so. Binding is **null but not resolved**: the plan's contiguous-span window admits 9 of 129 annotated contact shells, because a shell draws residues from several secondary-structure elements at once and so has no short contiguous span. That is a property of the object and no sample size fixes it.

## 5. The instrument check the sequence line did not have

Pre-registered in `docs/COMPOSITION_PLAN.md`. Sections 3 and 4 return nulls, and a design that returns nulls has either found a real absence or has no sensitivity. **Nothing in those studies distinguishes the two**: their only positive control was a synthetic planted signal.

So this study takes the other half of the decomposition — composition, which the order null holds fixed and cannot see — and puts a known-true signal through the same machinery. An inositol-phosphate contact shell being arginine-rich is textbook and is **not offered as a finding**; it is offered as the check.

The gate is calibrated before any real lining reaches it: a planted basic patch gives mean z **+5.99** (p 0.008319) and opens it, while 12 composition-matched noise pools average +0.0470 against an expected zero and **0** of them open it. All 5 checks pass.

**What the gate decides.** If the positive control fails on real shells, every arm of that study reads *instrument not demonstrated* and **the sequence line's nulls above cannot be read as null results either**. That is fixed in code, not promised in prose.

**The study itself has not reported yet.** When it does, this section gains its result and section 4's nulls either firm up into statements about order or become uninterpretable. Until then the question is open, and this page says so rather than leaving the line out.

## 6. Feature families, and which have actually been tried

The two decompositions above cover a lining's **positions** and a window's **letters**. It is easy to read their nulls as though nothing about a sequence or a structure could predict an IP site, so here is what has and has not been put to a test.

| family | status |
|---|---|
| lining topology and register | tested, null against the control that matters (section 3) |
| sequence order | tested, null; binding arm unresolved (section 4) |
| sequence composition | instrument check calibrated, study in flight (section 5) |
| cross-organism IP-concentration coevolution | tested, **inconclusive** — -0.017 pp [-0.341, +0.269] over 722 matched bins and 20080 clusters (`docs/COEVOLUTION_PLAN.md`) |
| a learned pocket score | exists, and is a component of the combined score; the protein-level classifier benchmarked at chance |
| orthogroup-level conservation of hits | computed in the comparative analysis |
| **per-residue conservation from an alignment** | **not tested** |
| **residue-residue coevolution as a site feature** | **not tested** |
| **a learned representation of the site itself** | **not tested** |

The last three are where a signal could still live, and the nulls above say nothing about them. Note that the coevolution row is a different question from the last two: it asks whether an organism richer in inositol phosphate buries more of them, not whether a particular site's residues covary.

## What is not established, and what would change it

* **That any candidate binds an inositol phosphate.** No candidate has independent support. The honest description is a ranked list of hypotheses.
* **That the pipeline generalises to cryptic sites.** One cryptic example in the validation set cannot support a performance estimate for cryptic sites, and the deposited archive cannot supply more, because a buried site with its ligand in place is hard to crystallise. This is the project's binding constraint and it is not a statistical one.
* **That the machinery adds anything beyond the same-protein baseline.** Section 1's third row says it does not, at the power available.
* **That structure or sequence carries a detectable signature of an IP site.** Sections 2 to 4 are null or unresolved, and section 5 is what decides whether those nulls mean anything.

What would change the picture, in the order it would be worth doing: an experimentally tested candidate; published functional evidence for one of the deposited cryptic copies other than ADAR2, which would turn the ground-truth arms from structural classification into established biology; and one of the three untested feature families in section 6.

Every claim above is read from the artefact beside it. Where a number here and a number in `results/publication/RESULTS_SUMMARY.md` disagree, that file is older.
