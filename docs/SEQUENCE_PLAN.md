# Sequence-signature plan: pre-registration for study Q6

Dated 2026-10-12, **before any annotated sequence has been read by any function
in `cryptic_ip/detect/sequence.py`**. The module exists and is unit-tested; it
has been run on strings written by hand to probe its estimators and on nothing
else. Every number that appears below is a property of the method or of the
dataset's size, computed without looking at any outcome. Changes after today
are dated amendments in separate files, written before the results they could
affect are read.

This plan covers study Q6 only. It does not amend studies A–P or Q1–Q5; those
stand as reported.

## The question

Is there a signature *in the sequence* at an annotated inositol-phosphate site
— in its entropy, its compressibility, its periodicity, the spacing of its
coordinating residues — that says a pocket of a particular kind forms there?

Two halves, and they are pre-registered separately because they need different
controls:

* **Q6a — association.** Do sequence descriptors of a site's lining window
  differ from matched controls?
* **Q6b — type.** Among sites, do the descriptors differ between *kinds* of
  pocket (buried vs. solvent-exposed; cryptic vs. constitutively open, in the
  sense `detect/ensemble.py` measures)?

Q6b is the half the question was really asking — a signature that "leads to a
certain type of pocket" — and it is the half with the weaker design, because
stratifying a 29-group dataset by pocket type is exactly the operation that
drops strata below the evidence floor. The power calculation below is done per
stratum for that reason, and the expected answer for most strata is "no
evidence", stated in advance.

## Three ways to answer this wrongly

Written down first because the module is built around them and because each one
constrains what may be reported.

1. **Kolmogorov complexity is uncomputable.** No quantity in this study is
   named after it or claimed to estimate it. `lempel_ziv_complexity` is the
   exact LZ76 parse count; `compressed_bits` is a real compressor's output
   minus its measured empty-input header. Both bound K from above up to a
   constant. A result will not be worded as a statement about algorithmic
   complexity.
2. **Plug-in Shannon entropy is badly biased at the size of a pocket lining,
   and the bias depends on N.** Mean bias over 4000 uniform random windows per
   row, seed 20261012, against a uniform source over the 20 amino acids (true
   entropy 4.322 bits): plug-in is −1.408 bits low at N=10, −0.799 at N=20,
   −0.302 at N=50, −0.144 at N=100. A pocket lining is 10–20 residues.
   Comparing a 12-residue lining with a 20-residue one on plug-in entropy
   reports a difference that is pure sample size. Grassberger is the default
   estimator (−0.821, −0.316, −0.029, −0.002 at the same N) and window length
   is controlled explicitly as well — see *Windows*.
3. **Composition is the confounder that matters.** This project established
   long ago that these sites are arginine- and lysine-rich. Every entropy and
   complexity measure of a window is a function of its composition, so a
   "signature" can be nothing but "this window has a lot of arginine" wearing a
   different name.

   **Therefore the per-window quantity entering every primary test is not a raw
   measure. It is `composition_z`: the measure's z-score against 2000 shuffles
   of the same window.** That null holds the window's length and its letter
   counts fixed and destroys only the order, which makes the z-score free of
   the bias in (2) by construction and makes it answer the only interesting
   form of the question — whether the ordering carries information the
   composition does not.

   A difference in a *raw* measure between sites and controls is not a finding
   under this plan, however large or however small its p-value. It is reported
   in the secondary table with the composition z beside it, and if the z is
   null the raw difference is reported as composition-explained.

## Data

* **Sites (arm 1).** The IP benchmark table (run 35949588200,
  `benchmark-table/table.csv.gz`, SHA-256
  `a21d92fd0f57306848c260d1f163017f0bf1e82a5dfcf1f3be69e0c20b698a1d`), rows
  with `label_ip_site == 1`. The lining window comes from
  `pocket_residue_numbers` and the deposited sequence of the same chain.
* **Controls (arm 0).** Two control sets, both required:
  * **Matched non-site pockets.** Other fpocket pockets on the *same* chains,
    matched to each site pocket on lining length (±2 residues) and on basic
    residue count (±1), by the same matching code as `matched_controls`. Same
    chain means the same amino-acid usage and the same structure-prediction
    pipeline, so the comparison is within-protein.
  * **Composition-matched shuffles.** Each site window's own
    `composition_z` null. This is the control that makes a positive
    interpretable, and it is per-window, so it needs no matching at all.
* **Sequences.** The chain sequences already present in the cached structures
  used by the benchmark run. No new fetch: the sandbox cannot reach RCSB,
  UniProt or AlphaFold, every data step runs in GitHub Actions, and this study
  introduces no new download.
* **Excluded.** Windows shorter than 8 residues after dropping non-standard
  residues (the estimators need at least 4, but an 8-residue floor is set here
  so a `k=3` block entropy is not read off 6 points), and any chain where more
  than 10% of the lining residues are non-standard.

## Windows

Three window definitions, fixed now, all three reported:

| window | definition | what it can see |
|---|---|---|
| `lining` | the lining residues in sequence order, concatenated, gaps closed | composition and letter statistics only; **no periodicity** |
| `local` | the contiguous stretch spanned by the lining residues, capped at 40 | periodicity, spacing, block structure |
| `spacing` | the lining residue *numbers* alone, via `spacing_spectrum` | gap structure with no letters at all |

The `lining` window must not be used for any periodicity statistic: closing
the gaps destroys the positional information a period means. This is written
down because it is the easiest mistake to make with this data, and the
analysis code will refuse it.

Window length is reported as a covariate in every comparison, and a comparison
whose two arms differ in mean window length by more than 2 residues is
reported as length-confounded regardless of its p-value.

## Grouping and estimands

* **Groups.** `homology_group_strict` (sequence plus Foldseek), with the
  MMseqs2 30% clustering as the sensitivity analysis. Resampling is over whole
  groups, never over copies.
* **Per-copy estimand.** The mean descriptor over pocket copies. Reported for
  completeness; not the basis of any claim, because 29 groups hold 209 copies
  unevenly.
* **Per-group estimand.** Each group contributes its own mean; the comparison
  is between group means. This is the primary estimand.
* **Evidence floor.** A stratum with fewer than 5 groups on either side gives
  no evidence, whichever way it points. `FeatureComparison.is_evidence`
  enforces this and the floor is applied before any difference is read.

## Primary tests

Six, fixed now, no substitutions. Each is a difference in per-group mean of a
per-window `composition_z`, sites against matched non-site pockets, with a
group bootstrap interval (2000 resamples) and a group-label permutation
p-value. **Holm correction across all six.**

| # | window | statistic under `composition_z` | the hypothesis |
|---|---|---|---|
| 1 | `local` | `entropy_rate_k2` | the ordering is more predictable at a site than its letters imply |
| 2 | `local` | `helix_charge_periodicity` | a basic helix face lines the pocket |
| 3 | `local` | `strand_charge_periodicity` | a basic strand face lines it instead |
| 4 | `local` | `dominant_period_power` | *some* charge period dominates, whichever it is |
| 5 | `lining` | `lz_normalised` | the lining letters are more repetitive than chance |
| 6 | `spacing` | `fraction_helical` (no null; a proportion) | the coordinating residues are 3–4 apart |

Test 6 has no shuffle null because it has no letters; it is tested as a
proportion difference between arms with the same group bootstrap.

One note on test 5, written today while the module was being validated and
before any site sequence was read. `lz_normalised` divides the LZ76 parse count
by its asymptotic expectation `n / log_a(n)`, and that expectation is not
reached at the length of a pocket lining: measured on 300 random windows each,
the mean is 0.88 at N=20, 0.84 at N=36, 0.83 at N=100 and 0.88 at N=500 — below
1 and not monotone in N. The raw value is therefore not comparable across
windows of different length, which is a second reason test 5 is run on the
composition z rather than on the raw ratio, and the length-confound rule above
applies to it with no exception.

Everything else the module computes — the four reduced alphabets, the raw
measures, `excess_entropy`, `compressed_bits`, the per-alphabet entropies — is
**secondary and exploratory**. It is reported in full, uncorrected, labelled
exploratory, and no sentence in the results will call any of it a finding.

## Power, computed before any result

The design has 29 strict groups, split unevenly. With the group counts
available and `minimum_detectable_difference` at 80% power and two-sided
α = 0.05, the smallest detectable difference in mean `composition_z` is
**(z_α + z_β) × sqrt(s²₁/n₁ + s²₀/n₀) = 2.80 × sqrt(s²₁/n₁ + s²₀/n₀)**. The
spread terms are not known until the descriptors are computed, so the
numerical MDE is computed and **printed in the results file before the
difference itself is read**, as `FLEXIBLE.md` does. The analysis script emits
the MDE block first and the comparison block second, in that order, so the
ordering is in the artefact and not merely asserted.

### Q6b, fixed now

The comparison shape, written down today before any descriptor exists so it
cannot be chosen to suit a result. For each pocket type present among the
sites, that type's sites are compared with the sites of **all other types** —
one-against-rest, not pairwise, because pairwise over k types multiplies the
family by k(k−1)/2 and this dataset cannot afford it. The features are the same
six. The MDE is computed per stratum and printed before the differences, as
above.

**Holm for Q6b runs across features × types**, so a study with three pocket
types and six features corrects over eighteen tests. That is a harsh
correction and it is the honest one: "which of these features differs for
which kind of pocket" is eighteen questions, not six.

With exactly **two** types the correction is six, not twelve: one-against-rest
tests the same partition twice there — the second type against the rest is the
first with its sign flipped — so only the first type is run and the report says
so. From three types up the rows are genuinely different partitions and all are
run.

A stratum below the 5-group floor is reported as "no evidence" with its group
count and no difference at all. **With 29 strict groups this is the expected
outcome for most strata**, and saying so in advance is the point: a stratified
null here means the design could not see, not that there is nothing there, and
the MDE beside it is what distinguishes the two.

Pocket types are assigned from existing descriptors before the split, and the
assignment is recorded per pocket in the inputs file rather than computed
inside the analysis, so the split cannot be adjusted after a difference is
seen.

## Decision rules

Written before the data so they cannot be chosen afterwards.

* **A signature is claimed** only if a primary test's Holm-corrected p is below
  0.05, its bootstrap interval excludes zero, its observed difference exceeds
  the pre-computed MDE, both arms clear the 5-group floor, and the arms are not
  length-confounded. All five.
* **A raw-measure difference with a null composition z** is reported as
  composition-explained. This is the most likely outcome and is a real result:
  it would say the sequence signature at these sites *is* their composition,
  which the project already knows.
* **A null result** is reported with its MDE, so "we found nothing" is
  separated from "we could not have found anything". The project has made this
  mistake once already — a held-out AUC of 0.499 was read as a modelling
  failure when at five groups nothing below 0.87 was detectable.
* **No feature selection on the outcome.** The six primary tests are the six
  above. If a secondary feature looks interesting, that is a hypothesis for a
  new pre-registration on new data, not a result of this one.
* **Pocket type is defined before the split**, from the existing descriptors
  (`relative_sasa` for burial, `crypticity` from `detect/ensemble.py` for
  cryptic vs. open), with the thresholds those modules already publish. No new
  threshold is introduced here, and a continuum is reported as a continuum:
  `relative_sasa` is known to run smoothly from 0.070 to 0.896 and the 0.12
  cutoff is a calibration, not a boundary in nature.

## What this study cannot conclude

* Nothing about causation. "Leads to pocket formation" is not testable from a
  cross-section of annotated sites; the strongest available claim is
  association between a sequence descriptor and a pocket property, measured
  with homology groups resampled.
* Nothing about Kolmogorov complexity, per (1).
* Nothing from a stratum below the evidence floor, per the floor.
* Nothing about prediction. A descriptor that separates sites from matched
  controls in this dataset is not a predictor until it is tested on groups held
  out from its fitting, which this study does not do and does not claim.

## Execution

The analysis runs in GitHub Actions (`.github/workflows/sequence.yml`), prints
its results base64-wrapped between `BEGIN_SEQUENCE_MD_B64` /
`BEGIN_SEQUENCE_JSON_B64` markers and their matching ends, and the artefact is
extracted programmatically into `results/sequence/`. No number is retyped by
hand.

Two jobs, and they are not the same kind of thing:

* **Method validation** (`sequence_signature.py validate`) runs on every push
  to this branch and on demand. It plants a helical-face ordering signature in
  synthetic windows and checks three things: that the machinery finds it, that
  it reports nothing for the same letters shuffled, and that it reports nothing
  for a basic-rich pool with no ordering. The three pools are matched on
  composition — identical basic fraction, entropy within 0.1 bits — so a
  separation between them cannot be composition. The committed run is
  `results/sequence/VALIDATION.md`: planted mean z **+4.05** firing in 100% of
  windows, shuffled control **−0.24** at 0%, basic-rich control **−0.13** at
  2%. A failed check fails the run. This is a property of the method and is not
  a study result.
* **The study** (`describe` then `report`) runs only when handed an inputs
  artefact, and emits the power block before the differences.

**What is not built yet, stated here rather than discovered later.** The inputs
file needs, per pocket, a chain sequence and the lining residues as
`(chain, residue_number)` pairs. The repo's existing lining residue numbers are
chain-blind — `PocketAnalyzer.get_pocket_residues` says so in its own docstring
— so a chain-aware pass over the structures is required to build them, written
against the benchmark table's actual columns. That pass is the remaining piece
of study Q6. Writing it from guessed column names would produce a workflow that
goes green on an empty analysis, which this project has had happen once already
(a 272-shard matrix was silently dropped and the run passed with no docking in
it), so the study job refuses to run without an inputs artefact instead.
