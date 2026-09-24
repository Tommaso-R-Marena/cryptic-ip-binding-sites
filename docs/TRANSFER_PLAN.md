# Transfer plan: buried phosphate-dense ligand sites

Pre-registered on 2026-09-24, before any entry of this dataset was searched,
downloaded or measured. Changes are appended as amendments, dated, with
reasons, before the results they could affect exist.

## Why

The pre-registered inositol phosphate benchmark (`docs/ANALYSIS_PLAN.md`, run
35949588200) could not evaluate its hull-depth hypotheses. The PDB's 60 buried
inositol phosphate pockets come from about six homology families, and one holds
87 % of them. No re-analysis of that data can fix this.

The question behind H1 and H2 is physical and not specific to inositol: does
depth from the protein's exterior help recognise a **buried** site for a
**phosphate-dense anion**? Many other ligands share that chemistry: nucleotide
di- and triphosphates, PRPP, sugar bisphosphates, isoprenoid pyrophosphates and
pyrophosphate itself. They are bound by many unrelated folds, so a buried subset
should span far more independent families.

This plan asks the question on that class, and then tests whether what is learned
there transfers to inositol phosphate sites.

## Ligand class (the definition is the rule, not a list)

A chemical component belongs to the class when **all** of these hold for its
formula in the PDB chemical component dictionary:

- it has at least 2 phosphorus atoms;
- its phosphorus atoms per heavy atom is at least 0.07;
- it is **not** an inositol phosphate (`looks_like_inositol_phosphate` is false);
- it is not lipid-linked.

For scale:

| component | phosphorus / heavy atoms | in the class |
|---|---|---|
| IP6 | 0.17 | reference only; excluded by the inositol rule |
| ATP | 0.10 | yes |
| ADP | 0.074 | yes |
| PRPP | 0.14 | yes |
| NADP | 0.06 | no: mostly nucleoside, not phosphate |
| CoA | 0.06 | no |
| FAD | 0.04 | no |

Candidate identifiers are seeded in
`cryptic_ip/database/polyanion_ligands.py`, resolved against the chemical
component API at run time, and kept only if they pass the rule. The resolved
registry and every rejection are written to the dataset manifest.

## Entries

- X-ray structures at 2.5 Å or better that contain at least one class component.
- Any entry that contains any inositol phosphate component is excluded, so the
  two datasets share no entry.
- A seeded random sample of **1,200** entries, with seed 20260925. The size is
  set by compute: 1,200 X-ray entries extract in about 60 shards.
- As in the benchmark's Amendment 1, an entry whose polymer atoms exceed 99,999
  is excluded and counted.

## Labels, descriptors, groups and holdout

These are identical to the benchmark, reusing its code unchanged:

- **Labels.** Pockets overlapping a ligand copy are positives. Burial class comes from
  that copy's relative SASA: cryptic ≤ 0.12, and crystal-artefact copies are excluded.
- **Descriptors.** 39 descriptors on ligand-free structures.
- **Groups.** MMseqs2 sequence groups and strict (sequence + Foldseek) groups.
- **Holdout.** The latest 20 % of strict groups.

The benchmark's task names are kept in the table so its code runs unchanged:

- `ip_site` here means "a pocket on a class ligand";
- `cryptic_ip_site` means "a pocket on a buried class ligand";
- `burial` means buried versus surface among class-ligand pockets.

Homology groups are computed **jointly** over this dataset and the inositol
phosphate benchmark's 367 entries, so that one namespace covers both.

## Hypotheses

**T1.** In this dataset, adding hull depth to the descriptors improves the
recognition of buried class-ligand sites. This is task `cryptic_ip_site`: paired
ROC-AUC, full arm minus `no_hull_depth`.

**T2.** Adding hull depth improves telling buried class-ligand sites from surface
ones. This is task `burial`, with the same paired comparison.

Both are decided exactly as H1 and H2 are, by the benchmark's code and rules:

- both groupings are required;
- the equivalence margin is 0.01;
- the result is Holm-corrected across T1 and T2;
- the holdout must agree when it holds at least 10 positive groups;
- the permutation control must pass;
- if one group holds more than 40 % of a task's positives, that task is not evaluable.

Following diagnostic D1 of the benchmark, the permutation control is judged on
**10 permutations**, not one. It fails if their mean pooled ROC-AUC exceeds 0.52,
or if more than 1 of the 10 exceeds 0.60.

**T3 (external, descriptive).**

- **Training set.** A model is trained on this dataset's development rows
  (`cryptic_ip_site` task, full arm), after removing every entry whose strict
  group contains any inositol phosphate benchmark entry.
- **Model selection.** Family, hyperparameters, threshold and calibration are
  selected by grouped inner cross-validation on those rows alone.
- **Scoring.** The locked model then scores every pocket of the inositol phosphate
  benchmark. That table is fixed by its SHA-256 from run 35949588200.
- **What is reported.**
  - ROC-AUC and PR-AUC for `cryptic_ip_site` and `ip_site` on that table;
  - intervals from the inositol table's strict groups;
  - the rule-based score beside them.
- **Status.** Five strict groups hold all buried inositol phosphate positives, so T3
  is **descriptive**: it has no p-value and supports no claim of transfer. It can
  still show that transfer fails, if the point estimate sits at chance.

## What would count as a finding

- T1 or T2 **supported**: hull depth helps recognise buried phosphate-dense sites.
  That would be the first evaluable answer to the question the benchmark could not settle.
- T1 or T2 **refuted**, meaning the interval lies inside ±0.01: hull depth adds
  nothing, and the screen's hull-depth gate should be reconsidered.
- **Not evaluable** again means that buried phosphate-dense sites are rare across
  families as well. That is itself informative for the proteome screen's prior.

## Ledger

Every run of this plan's workflow is listed in the report, with its commit and
table SHA-256. No result is dropped.
