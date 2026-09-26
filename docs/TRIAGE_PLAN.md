# Triaging the screen's IP candidates by pocket conservation: plan

Pre-registered on 2026-09-26, before any homologue of any candidate was fetched or
aligned for this study.

**What was known when this was written.** The candidate list
(`results/specificity/candidates.csv`, study C) had been read, so this study is a
**filter applied to a list already in hand, not a test of the screen**. Specifically:

- 75 candidates, 25 per organism, all unseen and none annotated as IP binders.
- Two annotated IP binders sit inside the human top 27 and none inside the yeast or
  *Dictyostelium* top 25, against a background rate of 34 in 33,084 unseen proteins.
- Study B applied a conservation criterion to two α-arrestin leads. ARRDC2 failed it: its
  mapped site had no basic residue at all. That criterion is reused here unchanged.

Nothing in studies A–G changes. This study adds an independent line of evidence —
evolutionary rather than structural — to candidates the screen has already ranked.

## Why

The screen ranks a pocket by geometry and charge in a single AlphaFold model. Nothing in
it asks whether that pocket is conserved. A pocket whose basic residues are preserved
across a protein's orthologues has been under selection; one whose basic residues are
idiosyncratic to a single species is more likely a surface accident of the model. This is
cheap, it is independent of every feature the screen uses, and study B already showed it
can kill a lead.

It cannot confirm a binder. A conserved basic pocket is necessary-ish, not sufficient:
plenty of conserved basic pockets bind other polyanions, which is what study C measured.

## Data

- **Targets.** All 75 candidates, using the `top_pocket_residues` recorded for each in
  `results/specificity/candidates.csv` (UniProt numbering, the model version screened).
- **Positive controls.** The annotated IP binders among the same three proteomes that the
  screen also scored (`annotated` true in study C's ranking). These are proteins the
  filter *should* mostly pass. They are the calibration: a filter that fails known
  binders is too strict to use.
- **Negative controls.** For each candidate, one protein drawn from the same organism's
  unseen set, matched on mean pLDDT (within 5) and hull depth (within 3 Å) and ranked
  below the 50th percentile of the combined score, chosen by a fixed rule in `copy_key`
  order with seed 20261001. These are pockets with the same model quality and burial but
  no ranking signal.
- **Homologues.** The UniProtKB members of each target's UniRef50 cluster, as in
  `docs/ARRESTIN_PLAN.md`: deduplicated at 100 % identity, at most 500 sequences in
  accession order. Fewer than 10 sequences: *not evaluable* for that protein.
- **Alignment.** MAFFT 7.526 `--auto`, through the same code as study B.

## The criterion (fixed here, unchanged from study B)

At each pocket position that is K, R or H in the target, the fraction of homologues
(target excluded, gaps count as not basic) carrying K, R or H in that column. A pocket is
**conserved-basic** when at least 3 pocket positions are basic in the target and at least
3 of them have a basic fraction ≥ 0.80.

## Decisions

- **H1, primary (the only hypothesis test).** The difference in the conserved-basic rate
  between candidates and their matched negative controls, paired on the matching.
  - Intervals: 2,000 bootstrap resamples of MMseqs2 30 % clusters (seed 20261001),
    resampling candidate and its control together.
  - **enriched** if the 95 % lower bound > 0.
  - **not enriched** if the upper bound < 0.05.
  - **inconclusive** otherwise. **Not evaluable** with fewer than 5 clusters.
- **H2, filter calibration (descriptive).** The conserved-basic rate among the annotated
  positive controls, with its interval. Reported beside H1. If the positive controls fail
  the filter at a high rate, H1's result is reported but the filter is declared
  uninformative for triage, and no candidate is promoted or demoted on it.
- **H3, the triaged list (the deliverable, descriptive).** Every candidate with its
  decision, the number of homologues, the basic positions and their conservation
  fractions. Candidates are labelled:
  - **conserved basic pocket** — passes, and is not explained by another ligand under
    study C's rule;
  - **explained** — study C's rule attributes it to another ligand class, regardless of
    conservation;
  - **not conserved**, or **not evaluable** (too few homologues).

Only the first label is offered as a candidate worth an experiment, and even then as a
hypothesis: nothing here measures binding.

## What the result can and cannot support

- **enriched, with positive controls passing.** The screen's top pockets are conserved
  more often than matched pockets of the same model quality, so the ranking is picking up
  something under selection rather than model noise. It still says nothing about *which*
  polyanion binds.
- **not enriched, or positive controls failing.** Conservation adds no information at
  this depth, the candidate list stands unfiltered on study C's evidence alone, and the
  honest summary is that we cannot separate the 45 unexplained candidates.
- No outcome makes any candidate a demonstrated IP6 binder. The only things that would
  are the assays study B's dossiers already name: ITC or a fluorescence binding assay on
  the purified protein with IP6, an ATP control, and charge-reversal mutants of the
  pocket's basic residues.

## Execution and outputs

- **Workflow.** `.github/workflows/triage.yml`: one fetch-and-align job, then a report.
- **Code.** `scripts/triage.py`. Conservation reuses `scripts/arrestin.py`'s functions by
  import; that script is not edited, because editing a finished study's script would
  re-run its workflow.
- **Results.** Printed between `BEGIN_TRIAGE_JSON` and `END_TRIAGE_JSON`, extracted into
  `results/triage/`.
