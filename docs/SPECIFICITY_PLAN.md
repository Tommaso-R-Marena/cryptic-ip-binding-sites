# Inositol phosphate versus other polyanion sites: plan

Pre-registered on 2026-09-24, before any model was fitted to tell inositol phosphate
(IP) sites from other phosphate-dense ligand sites, and before the two pocket tables
were joined. Changes are dated amendments in separate files, written before the results
they could affect are read.

## Why

The learned `ip_site` model recognises buried polyanion sites in general: its top
proteome candidates are dominated by PAPS, nucleotide and sugar-phosphate binders
(`results/learned_screen/LEARNED_SCREEN.md`). This study asks whether the 39 pocket
descriptors can tell an IP site from a site for another phosphate-dense ligand, and, if
they can, whether a ranking that uses that distinction finds IP binders better.

## Data

- **IP pockets (label 1):** rows of the IP benchmark table (run 35949588200,
  `benchmark-table/table.csv.gz`, SHA-256
  a21d92fd0f57306848c260d1f163017f0bf1e82a5dfcf1f3be69e0c20b698a1d) with
  `label_ip_site == 1`: pockets on any non-artefact IP copy.
- **Other-polyanion pockets (label 0):** rows of the transfer table (run 35972729955,
  `transfer-table/table.csv.gz`, SHA-256 recorded at run time) with
  `label_ip_site == 1`: pockets on a class ligand (nucleotide di- and triphosphates,
  PRPP, sugar bisphosphates, pyrophosphates; 35 components by formula rule). The two
  datasets share no entry (checked again).
- Pockets on neither ligand are not used: the question is conditional on a polyanion
  site.
- **Descriptors:** the benchmark's 39 (`protocol.ARMS["full"]`), computed by the same code
  on ligand-free structures in both tables.

## Grouping and holdout

- **Groups:** the joint homology groups of `transfer-dataset/joint_grouped.csv`
  (MMseqs2 `homology_group` and sequence-plus-Foldseek `homology_group_strict`, computed
  over both datasets together), mapped by entry. The benchmark table's own group columns
  are replaced; a row without a joint group stops the run.
- **Temporal holdout:** the latest 20 % of joint strict groups by earliest release date
  among the entries in this table (`protocol.temporal_holdout`), chosen from dates and
  groups alone.
- **Evaluability (the 40 % rule):** if one group holds more than 40 % of the
  development IP pockets under a grouping, that grouping cannot support the evaluation.
  The same share for the other-polyanion class is reported.

## Protocol

The benchmark's machinery unchanged (`scripts/benchmark.py run`, task column
`label_ip_site`, arm `full`): nested grouped CV, 5 outer × 3 inner folds, family and
hyperparameters chosen jointly in the inner loop (10 draws per family), threshold and
Platt calibration inside the training folds. Runs: sequence grouping repeats 0-2 (repeat
0 also fits the locked model and scores the holdout once), strict grouping repeat 0,
and 10 label permutations (sequence grouping, repeats 0-9). Intervals: 2,000
group-bootstrap resamples over the grouping of each run.

## Variants (all pre-declared; nothing else will be tried)

- **S1 (primary): all rows.**
- **S1b: without the nucleotide-fold super-group.** Every row of the one joint strict
  group holding the most other-polyanion development pockets is removed from
  development and holdout (in the transfer study this was G:10JT, 61 % of the buried
  positives). Its identity and size are reported.
- **S1r: comparable data.** IP rows restricted to X-ray entries at 2.5 Å or better, the
  transfer set's inclusion rule, so that a model cannot separate the classes by
  dataset (method, resolution) instead of chemistry.

## Decision rule (each variant)

- **Permutation control:** the ten permutations are judged by the transfer plan's rule:
  the control fails if their mean pooled ROC-AUC exceeds 0.52 or more than 1 of 10
  exceeds 0.60. A failed control makes the variant *not evaluable*.
- **Not evaluable** also if the 40 % rule fails under either grouping, or a grouping
  cannot be split.
- **Learnable** if the CV ROC-AUC lower bound exceeds 0.5 under both the sequence and
  the strict grouping, and, where the holdout holds at least 5 IP-positive groups and 5
  other-polyanion groups, the holdout ROC-AUC lower bound exceeds 0.5 as well.
- **Not learnable** if the upper bound is at most 0.60 under both groupings.
- **Inconclusive** otherwise.

**Gate for the re-ranking.** The re-ranking (below) runs if S1 is *learnable*, or if S1
is *not evaluable* only because of the 40 % rule and S1b is *learnable*; and in either
case S1r's point estimate exceeds 0.5 (a specificity that vanishes on comparable data
is a dataset effect). Otherwise it does not run, and that is reported.

## Re-ranking the screened proteomes (conditional)

- **Models.** P(polyanion site) is the locked `ip_site` model of the learned screen
  (`learned-screen/model/ip_site_locked.joblib`, run 35992676292). P(IP | polyanion site)
  is the specificity model locked on the whole table of the variant that opened the gate
  (development and holdout; strict groups for inner CV; `protocol.select_candidate`, 10
  draws, seed 0), calibrated.
- **Combined pocket score** = P(polyanion site) × P(IP | polyanion site). Protein score:
  the best confident pocket (mean pLDDT ≥ 70), as in the learned screen.
- **Screen data:** shards and catalogues of run 35935291031.
- **Unseen proteins:** no MMseqs2 hit (≥ 30 % identity over ≥ 50 % of the shorter
  sequence, E ≤ 1e-3) against the UniProt sequences of the IP benchmark's entries **or**
  of the transfer set's entries, because the specificity model has seen both.
- **L3a.** Among unseen proteins, the combined score ranks UniProt-annotated IP binders
  (`known_ip_annotation`) above other proteins: pooled ROC-AUC, cluster bootstrap (2,000,
  MMseqs2 clusters at 30 %/50 %). *Supported* if the lower bound exceeds 0.5.
- **L3b.** The combined ranking beats the learned `ip_site` ranking (L1) on the same
  unseen proteins: paired ROC-AUC difference, same resampling. *Supported* if the lower
  bound exceeds 0.
- **Multiplicity.** Holm across the study's primary tests: S1 (p-value of the sequence-
  grouping CV ROC-AUC against 0.5), L3a and L3b. Support requires the Holm-adjusted p
  below 0.05 as well.
- **Candidates.** The top 25 unseen, non-annotated proteins per proteome by combined
  score, each with the empirical precision at its rank (annotated binders at or above
  it, a lower bound) and a rule-based classification: **explained** if its UniProt
  keywords, name or function match
  `ATP-binding|GTP-binding|Nucleotide-binding|NAD|FAD|FMN|Coenzyme A|sulfotransferase|glycosyltransferase|UDP-|diphosphate|pyrophosphate|bisphosphate|phosphoglycer|mitochondrial carrier|Solute carrier family 25|ATPase|kinase|3'-phosphoadenos`
  (case-insensitive), and **unexplained** otherwise. The same rule is applied to the
  learned screen's candidates for comparison.

## Descriptor extension

If S1 is not learnable or inconclusive, a descriptor extension (for example the count
and spread of basic nitrogens over a 9-11 Å ring footprint, the Lys/Arg ratio, and the
absence of adenine-stacking aromatics and Mg-coordination motifs) may be pre-registered
in a separate plan before any extended descriptor is computed. Any such descriptor must
be computed identically on the IP table, the transfer table and the proteome screen.
It is not part of this plan.

## Execution and outputs

`.github/workflows/specificity.yml`. Results printed between `BEGIN_SPECIFICITY_JSON`
and `END_SPECIFICITY_JSON` and extracted into `results/specificity/`.
