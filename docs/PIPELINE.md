# The pipeline, end to end: what we ran, what we ran it with, and why

Tommaso R. Marena · `cryptic-ip-binding-sites` · compiled 2026-10-01

This document states the whole computational pipeline for detecting buried inositol-phosphate
(IP) binding sites: every stage, the exact parameter values in the code, the software and
version each stage depends on, and the justification for each choice. It is written so that a
reader who has never seen the repository can tell what was done, and a reader who disagrees
with a choice can find the line that encodes it.

Two pipelines run in sequence, and the order matters:

- **Part A — discovery.** Measure pockets across predicted proteomes and rank them against
  criteria calibrated on structures where the answer is known.
- **Part B — validation.** Ask whether the docking that Part A would lean on to confirm a
  candidate can solve a problem whose answer we already have.

Part B's verdict governs how Part A's output may be read, which is why both are here.

---

## 1. Software inventory

Versions are the pinned values in `requirements.txt` and in the `.github/workflows/*.yml`
install steps. Everything runs on **Python 3.11** in CI (the package supports ≥ 3.8; CI also
tests 3.10).

### Python packages (pinned, `requirements.txt`)

| Package | Version | Used for | Why this one |
|---|---|---|---|
| numpy | 2.2.6 | all array geometry, SASA accumulation, bootstrap resampling | — |
| scipy | 1.15.3 | spatial queries (`cKDTree`) for contacts and neighbour search | KD-trees make per-atom neighbour queries tractable at proteome scale |
| pandas | 2.3.3 | the census, per-copy tables, all result frames | — |
| scikit-learn | 1.7.2 | the ML comparison baseline, calibration, ROC/PR metrics | the standard implementation; we needed nested grouped CV, which it supports directly |
| statsmodels | 0.14.5 | supporting statistical routines | — |
| xgboost | 3.2.0 | gradient-boosted comparison model | a second model family, so "ML does not help" is not a statement about one algorithm |
| shap | 0.49.1 | feature attribution on the fitted models | we needed to see *what* a model keyed on before trusting or discarding it |
| biopython | 1.86 | structure parsing, **Shrake–Rupley SASA**, sequence alignment | SASA from BioPython rather than the `freesasa` CLI so burial is computed in-process with no external binary to drift |
| prody | 2.6.1 | structure handling in parts of the analysis layer | — |
| matplotlib / seaborn / matplotlib-venn | 3.10.8 / 0.13.2 / 1.1.2 | figures | — |
| goatools | 1.6.4 | GO term enrichment on screen output | — |
| requests / aiohttp | 2.32.5 / 3.13.3 | RCSB, UniProt, AlphaFold DB fetches | async fetching for proteome-scale downloads |
| joblib | 1.5.3 | parallel map over structures | — |
| tqdm, pyyaml, jsonschema, typing-extensions, tabulate | 4.67.3, 6.0.3, 4.26.0, 4.15.0, 0.9.0 | progress, config, schema validation, tables | — |
| streamlit / py3Dmol | 1.55.0 / 2.5.4 | the single-protein web UI | — |

### Docking stack (pinned in the workflow env `PY_DOCKING`)

| Package | Version | Used for | Why this one |
|---|---|---|---|
| **AutoDock Vina** | 1.2.7 | all docking; `vina`, `vinardo` and `ad4` scoring | the Python API exposes a flexible receptor (`set_receptor(rigid, flex)`), which study J needs; 1.2.x is the maintained line |
| **RDKit** | 2026.3.6 | ligand construction from SMILES, ETKDGv3 embedding, substructure matching for symmetric RMSD | the graph-automorphism machinery the RMSD correction depends on |
| **Meeko** | 0.8.0 | side-chain torsion trees for the flexible receptor (study J only) | writing `ROOT`/`BRANCH` trees by hand per residue type risks poses that look plausible and are wrong; Meeko makes them correct by construction |
| **PDB2PQR** | 3.6.1 | receptor protonation, AMBER force field, PDB + PQR output | the standard tool; emits both the protonated PDB and per-atom charges we need |
| **PROPKA** | 3.5.1 | pKa assignment driving PDB2PQR at pH 7.4 | titratable residues — histidine above all — must be assigned, not guessed |
| **gemmi** | 0.7.5 | fast mmCIF/PDB parsing, symmetry operations | handles both deposited formats and the symmetry-mate test |

### External binaries (conda-forge / bioconda)

| Tool | Version | Used for | Why this one |
|---|---|---|---|
| **fpocket** | 4.2.3 | cavity detection and alpha-sphere volumes | the de facto standard for geometric pocket detection; its alpha-sphere volume is what our volume window is calibrated against |
| **AutoGrid** | 4.2.9 | AD4 scoring maps | required by the AD4 arm |
| **MAFFT** | 7.526 | multiple sequence alignment of orthologues | fast and standard for the hundreds-of-sequences case |
| **US-align** | 20241201 | structural superposition (AlphaFold model onto crystal) | needed to place a predicted model in the crystal's frame before cross-docking |
| **MMseqs2** | *unpinned* (conda) | sequence homology groups | orders of magnitude faster than BLAST at this scale |
| **Foldseek** | *unpinned* (conda) | structural homology groups | joins homologues too remote for sequence search, which sequence-only grouping would wrongly treat as independent |
| **APBS** | *unpinned* (conda) | Poisson–Boltzmann electrostatics on the control panel | the rigorous calculation, affordable on a handful of structures but not at proteome scale |

> **Reproducibility limitation, stated rather than hidden.** MMseqs2, Foldseek and APBS are
> installed unpinned from conda and will drift. Everything else is pinned exactly. A rerun
> should pin these three to the versions recorded in the run log of the relevant workflow.

### Infrastructure

| Component | Detail |
|---|---|
| Execution | GitHub Actions, `ubuntu-latest`, one job per shard |
| Why | the analysis container cannot reach RCSB, UniProt or AlphaFold, so every data step runs in CI rather than locally — and CI runs are logged, dated and re-runnable |
| Result transport | printed between `BEGIN_<NAME>` / `END_<NAME>` markers and extracted by `scripts/extract_log_block.py` |
| Why | artifacts cannot be downloaded into the analysis container; marker extraction means **no committed number is ever retyped by hand** |

### Considered and rejected

| Tool | Decision | Reason |
|---|---|---|
| **Boltz-2** co-folding | not used | a 100-residue IP6 co-fold did not finish in 5400 s on 4 CPUs at the cheapest settings. Measured, not assumed; reserved for GPU |
| **Meeko receptor builder** (for the main benchmark) | not used | it rejects residues with non-ideal geometry, which would drop copies non-randomly. We wrote our own deterministic PDBQT typer instead, so every structure takes one path |
| **freesasa** CLI | not used | SASA is computed in-process with BioPython's Shrake–Rupley, removing an external binary from the measurement path |

---

## 2. Data sources

| Source | What we took | Provenance |
|---|---|---|
| **RCSB PDB** | every deposited structure containing an inositol phosphate; CCD chemical component definitions | census run `36020459246`; 367 entries cached |
| **AlphaFold DB** | predicted proteomes: human, *S. cerevisiae*, *D. discoideum* | screen run `35935291031` |
| **UniProt / UniRef50** | annotations, and up to 500 orthologues per candidate for the conservation filter | study H |

---

## 3. Part A — discovery

### A1. Identify the ligands from coordinates, not from a code list

**What.** A molecule counts as an inositol phosphate when: six carbons form a ring at C–C
bonding distance (≤ 1.75 Å); an oxygen lies within 1.65 Å of at least 5 of the 6 ring carbons;
and phosphorus lies within 1.90 Å of one of those oxygens. The series (InsP3…InsP6) follows
from the phosphorus count. `cryptic_ip/analysis/inositol_detection.py`

**Why.** A hand-written list of PDB chemical component identifiers fails in two directions: a
site whose ligand code is missing from the list becomes invisible rather than negative, and a
code can name the wrong chemistry — our previous list contained `INS` (*myo*-inositol), which
carries no phosphate at all. Deciding from atoms needs no vocabulary and no network, so a
regioisomer, a pyrophosphate or a deoxy analogue is recognised on its structure. **This test
resolved 99% of deposited entries (135 of 136), which is the evidence that it is not too
strict.**

### A2. Measure burial per ligand copy

**What.** `relative_sasa = SASA(copy inside the complex) / SASA(same copy in isolation)`,
computed with Shrake–Rupley. Four complementary measures are reported: relative SASA, relative
phosphate SASA, burial depth, and enclosure. Classes, applied to the larger of the whole-ligand
and phosphate ratios: `cryptic ≤ 0.12`, `semi_cryptic ≤ 0.25`, `surface > 0.25`, and
`crystal_artifact` for a copy with fewer than 8 protein heavy-atom contacts within 4.5 Å.
`cryptic_ip/validation/burial_metrics.py`

**Why per copy, and why a ratio.** Absolute, copy-summed SASA — what earlier versions used — is
neither size- nor copy-number-independent: a structure with six InsP6 copies reported roughly
six times the SASA of one with a single copy, so a **crystallisation artefact was determining
the burial class**. The ratio removes both.

**Why the larger of the two ratios.** Conservative by design: a ligand buried to the ring but
with its phosphates in solvent should not be counted as sequestered, because sequestering the
phosphates is what a polyanion site has to do.

### A3. The 0.12 boundary, and the honest finding about it

**What.** The cutoff was calibrated on five deposited controls — ADAR2 (1ZY7) **0.093**; Btk PH
(1BWN) 0.253; PLCδ1 PH (1MAI) 0.373; HDAC1 (5ICN) 0.436; Pds5B (5HDT) 0.466 — then tested
against all 135 measurable deposited entries. `scripts/calibrate_controls.py`,
`scripts/burial_survey.py`

**What the survey found.** Two boundary estimators **disagree**: the density minimum falls at
0.138, Otsu's at 0.463. The quantiles run smoothly from 0.070 (q01) to 0.896 (q99). **There is
no natural two-class structure; burial is a continuum.**

**Why we kept 0.12 anyway.** Not because it is a natural break, but because it is the least
sensitive place to put a line: moving the cut across 0.10–0.15 changes the positive count by 5
entries, while the same 0.05 shift at 0.25–0.30 changes it by 18. It also sits within one
histogram bin of the density minimum, having been calibrated on ADAR2 alone — independent
corroboration.

**The limits we accept by doing this.** The positive class is *defined* by this cutoff rather
than discovered, so every downstream AUROC and enrichment is conditioned on it. The buried side
of the control panel is a **single structure** (ADAR2); HDAC1 was previously counted as a second
until we found that measurement had been made on component `6A0`, which carries no phosphate.
The original literature-derived cutoff of 0.05 would have yielded **zero** positives on the
deposited set.

### A4. Describe and score each pocket

**What.** Cavities from fpocket; 40 physically interpretable descriptors per pocket
(`cryptic_ip/analysis/features.py`); a rule-based composite of six smooth logistic components
(`cryptic_ip/analysis/scorer.py`):

| Component | Weight | Midpoint / window | Slope |
|---|---|---|---|
| Lining-residue SASA | 0.25 | 20 Å² | 0.12 |
| Burial depth | 0.22 | 12 Å | 0.45 |
| Basic residues | 0.20 | 3.5 | 1.1 |
| Enclosure | 0.13 | 0.75 | 12 |
| Cavity volume | 0.10 | 300–1600 Å³ (tolerance 400 Å³) | — |
| Electrostatic potential | 0.10 | 3.0 kT/e | 0.5 |

A component whose input is unavailable scores a neutral 0.5.

**Why rules rather than a fitted model.** The score has to be explainable to a reader without
reference to a fitted object, and the ML comparison (A6) needs a meaningful baseline.

**Why smooth components.** The original components were step functions — a basic-residue count
of 4 scored 0.8 and 3 scored 0.4. Discontinuities that large make the score unstable under any
measurement noise and are not justified by the biology. Every component is now a monotone
function of its input, so a small change in a measurement produces a small change in the score.

**Why the volume window is 300–1600 Å³.** It is a **cavity** window, not a ligand window. The
original 300–800 Å³ was the space InsP3–InsP6 themselves occupy, but it was applied to fpocket's
alpha-sphere volume, which measures the cavity and is systematically larger. Under the old
window ADAR2 scored 0.16 on volume while the Btk surface negative scored 1.00 — **the component
actively penalised the paradigm positive and rewarded a negative.** Measured real sites span
491–1525 Å³.

**Why a screened-Coulomb surrogate for electrostatics.** An APBS grid per structure is not
affordable at proteome scale. APBS is run on the controls, where the surrogate can be checked
against it.

**Why depth means burial depth.** `depth` is the distance from the pocket centre to the nearest
solvent-exposed atom. The previous pipeline passed fpocket's *mean local hydrophobic density*
into this slot and treated it as a distance. That single defect produced our one earlier
"candidate": P07264 scored 0.955 on a density of 32.1 read as a depth; its true depth is 2.84 Å
and the corrected score is 0.711, below threshold. Hull depth was tested as an alternative and
made the score **worse** on the held-out benchmark (AUROC 0.890 vs 0.934), so burial depth
stays.

### A5. Screen the predicted proteomes

**What.** A protein is a hit when its best confident pocket clears every gate
(`cryptic_ip/analysis/proteome_stats.py`, `PLAN_CRITERIA`):

`composite score ≥ 0.75` · `lining SASA ≤ 10 Å²` · `≥ 4 basic residues` ·
`volume 300–1600 Å³` · `mean pLDDT ≥ 70`

**Why every pocket is recorded, not only the passes.** Hit calling can then be varied — and the
threshold swept — without re-running a single structure.

**Why a pLDDT floor.** A pocket whose lining is a low-confidence prediction is not evidence
about the protein; it is evidence about the predictor.

**What the known-binder check says, including against us.** ADAR2 (0.575) and ADAR1 (0.614)
pass. ADAT1, PDS5B and HDAC1 are blocked at the score gate and HDAC3 at the volume gate. **Our
gates miss known binders**, so the screen's output is a hypothesis list, not a detector.

### A6. Triage by evolutionary conservation (study H)

**What.** Up to 500 orthologues per candidate from UniRef50, aligned with MAFFT; a pocket is
conserved-basic when **at least 3** basic pocket positions hold a basic residue in **at least
80%** of orthologues. `scripts/triage.py`, criterion imported unchanged from `scripts/arrestin.py`

**Result.** Candidates 0.925 [0.825, 1.000] against matched controls 0.365 [0.216, 0.514];
enrichment 0.646 [0.485, 0.808] over 38 matched pairs in 33 clusters. Annotated binders score
0.714, so the filter keeps most true positives.

**Why conservation.** A pocket that matters to the organism tends to be preserved; a basic patch
present only in one lineage is more likely an accident of surface composition.

**Why the criterion was imported, not reimplemented.** Study B had already fixed it. Rewriting
it would have let it drift, and editing study B's script would have re-run its workflow.

### A7. Test whether machine learning beats the rules — it does not

**What.** Nested, **group-aware** (homology groups never split across folds), calibrated
cross-validation against the rule-based baseline. Held-out test: random forest ROC-AUC **0.499**,
PR-AUC 0.0012; threshold score ROC-AUC **0.556**, PR-AUC 0.0017.

**Why we report this.** Severe pocket-level class imbalance — a handful of known buried sites
against millions of ordinary pockets — is the cause. The honest conclusion is that threshold
scoring stays the deployment mode until more positive structures exist, and that is what the
code does.

---

## 4. Part B — validation: can the docking be trusted?

### B1. Census

662 IP copies found across 367 entries → 462 eligible → **272 selected**. Exclusions are
itemised: 167 copies whose configuration differs from the CCD template, 32 crystal artefacts, 1
valence failure. `scripts/redocking.py census`

**Why exclude configuration mismatches.** If the deposited coordinates describe a different
stereochemistry from the reference definition, a "failed" redocking would be measuring the
mismatch, not the method.

### B2. Receptor preparation

**What.** (1) Strip every non-polymer atom — waters, ions, cofactors, and **every other IP
copy** — keeping all polymer chains of the asymmetric unit; selenomethionine written as
methionine, other modified residues removed and counted. (2) Protonate with PDB2PQR (AMBER)
driven by PROPKA at **pH 7.4**. (3) Type atoms for AutoDock with our own writer: polar hydrogens
kept as `HD`, non-polar merged into the parent, aromatic ring carbons `A`, acceptor nitrogens
`NA`, oxygens `OA`, sulfur `SA`. `cryptic_ip/docking/receptor.py`

**Why strip everything.** Anything left in the site templates the pose. A second IP copy is the
worst case: it is the answer.

**Why keep all chains.** These sites sit at interfaces; dropping chains would destroy them.

**Why our own typer rather than Meeko's receptor builder.** Meeko's builder rejects residues
with non-ideal geometry. The structures it rejects are not a random sample, so using it would
have dropped copies non-randomly — a bias in the benchmark itself.

### B3. Ligand construction

**What.** The docked ligand is **never built from crystal coordinates**. It is built from the
component's CCD SMILES, embedded with RDKit **ETKDGv3** at a fixed seed, randomly rotated, and
checked to start **more than 2 Å** from the crystal pose. The crystal copy is used only as the
RMSD reference. Two protonation states: `primary` deprotonates every terminal phosphate oxygen
then returns `floor(n_P / 2)` protons (InsP6 at **−9**, close to the measured charge near pH
7.4) and `deprotonated` leaves all oxygens charged (InsP6 at −12). `cryptic_ip/docking/ligand.py`

**Why.** This is what makes the test honest. Handing the program the ligand in the conformation
and place the crystal found is an open-book exam; the ≥ 2 Å start check enforces it mechanically
rather than by assertion.

### B4. Search

**What.** AutoDock Vina with every parameter fixed by the plan before any run: **exhaustiveness
32**, **20 poses**, **energy range 5 kcal/mol**, **seeds 1, 2, 3**. The box is centred on the
site, sized to the ligand's extent plus **8 Å per side**, minimum **22 Å** per side.
`cryptic_ip/docking/engine.py`

Arms: `primary` (Vina × 3 seeds + crystal-pose control), `secondary` (Vinardo, AD4, the
deprotonated ligand, metals kept), `pockets` (fpocket top-3 site finding and a decoy pocket),
`alphafold` (the same docking into the US-align-superposed AlphaFold model).

**Why three seeds.** To distinguish a real result from a lucky one, and to measure run-to-run
variability rather than assume it away.

**Why the decoy-pocket arm.** It asks a different question from pose accuracy — whether the
score can tell the true site from a plausible wrong one — and that turned out to be the one
thing the scores do well.

### B5. Scoring the answer — symmetry-corrected RMSD, no superposition

**What.** A pose is compared with the crystal pose **in the receptor frame, with no alignment**.
The core (every heavy atom except terminal oxygens on phosphorus) is matched by **all graph
automorphisms** (12 for an inositol phosphate core) found by RDKit substructure matching on an
element-only graph; for each mapping, each phosphorus's terminal oxygens are assigned by
**exhaustive permutation**. An incomplete crystal copy is matched as a substructure and the RMSD
runs over the atoms it has. Success is a top-ranked pose within **2.0 Å**.
`cryptic_ip/docking/rmsd.py`

**Why no superposition.** Aligning the two ligands would remove exactly the placement error
docking is being judged on.

**Why the symmetry correction.** Inositol phosphates have an automorphic heavy-atom graph (the
ring and its six phosphates can be relabelled) and each phosphate's terminal oxygens are
equivalent by resonance. Without the correction, a correct pose written in a different atom
order scores > 0 Å — we would be punishing correct answers for bookkeeping. We do not rely on
RDKit's `symmetrizeConjugatedTerminalGroups` because whether it covers P–O depends on the
version; a test checks our implementation against brute-force enumeration.

### B6. Statistics

**What.** Entries are grouped so homologues never count as independent:

- **sequence** grouping — MMseqs2, ≥ 30% identity over ≥ 50% of the shorter chain, E ≤ 1e-3;
- **strict** grouping — those links plus Foldseek structural links, TM-score ≥ 0.5 over ≥ 50%;
  groups are the connected components. Every chain participates, not only those touching the
  ligand. `cryptic_ip/benchmark/homology.py`

Every outcome is reported under **two estimands** — per copy and per group — with 95% intervals
from a **2000-resample bootstrap that resamples whole groups, never copies**. A stratum with
**fewer than 5 groups gives no evidence** and is reported as not evaluable. Multiple primary
tests within a study are corrected by **Holm**. `cryptic_ip/docking/stats.py`

**Why group resampling.** The PDB is not a fair sample: popular proteins are deposited dozens of
times. Resampling copies would give intervals that assume an independence the data does not
have.

**Why both estimands.** They answer different questions — "how often does this work on a
structure" and "how often does this work on a protein family" — and reporting only the flattering
one would be a choice made after seeing the data.

**Why a 5-group floor.** Fixed in advance so that a stratum cannot become evidence by being
small and lucky. It bites on our own most favourable number: cryptic-burial success is 0.343,
and with 4 groups we report it as **not evidence**.

---

## 5. What the validation found, and what follows

| Question | Result (group estimand) | Verdict |
|---|---|---|
| R1 top-pose success, redocking into the copy's own crystal | 0.110 [0.046, 0.193] | **unreliable** |
| Near-native pose anywhere in the list (E32) | 0.385 [0.256, 0.516] | — |
| G1 ceiling, E128 − E32 | +0.019 [−0.001, 0.047] | **search-saturated** |
| G2 top pose, E128 − E32 | +0.002 [−0.002, 0.009] | **no gain** |
| R4 true site vs decoy pocket (ROC-AUC) | 0.744 [0.658, 0.873] | **discriminates** |
| R3 docking into AlphaFold models | 0.045 [0.000, 0.122] | **not trustworthy** |
| F1 electrostatic re-ranking | +0.054 [−0.001, 0.114] | **no detectable difference** |

**The structural conclusion.** A near-native pose is found roughly four times more often than it
is ranked first, and quadrupling the search budget moves neither. **The bottleneck is the scoring
function, not the search.** The same function nonetheless separates a true site from a decoy at
AUC 0.744, so the scores are worth more as a site filter than as a pose ranker.

**Why this is unsurprising, and the standing caveat on all of it.** Vina's scoring has no
explicit electrostatic term. Scoring a −9 polyanion with a function that ignores charge is weak
evidence whatever the search budget.

**Study J** (flexible receptor) was the last structural explanation before the scoring
function stands alone as the suspect, and it came back the other way: letting the pocket's
side chains move makes redocking **worse**. Top-pose success falls to 0.026 from
0.096 rigid, a paired difference of -0.078 [-0.140, -0.031] per homology
group (Holm p 0.001); the best-of-list ceiling falls to 0.179 from
0.367. Both on 209 copies in 29 strict groups. The extra torsional freedom
enlarges the search space faster than the scoring function can exploit it, which leaves the
scoring function itself as the remaining suspect. Full figures and censoring in
`results/flexible/PROVENANCE.md`.

---

## 6. How the process was kept honest

1. **Pre-registration.** Every study has a plan in `docs/` — question, success criterion,
   estimand, interval method, decision rule — committed **before any data was fetched**. The
   biggest risk in a project like this is not a bug; it is deciding what counts as success after
   seeing the numbers.
2. **Dated amendments.** Changes are written as dated amendments before the affected results are
   read, with the reason (for example `docs/FLEXIBLE_PLAN_AMENDMENT_1..3.md`).
3. **No retyped numbers.** Results are extracted from run logs by `scripts/extract_log_block.py`;
   each result directory carries a `PROVENANCE.md` naming run, commit and job ids, and every
   committed figure is verified against the extracted JSON by script.
4. **Failures are reported.** Strata under the evidence floor, failed copies, censored
   comparisons and known binders our gates miss are all in the reports.

---

## 7. Limitations

- **Redocking is the easy case.** A copy is docked into its own crystal structure, so 0.110 is an
  upper bound on prospective performance, not an estimate of it.
- **The burial class is a cut on a continuum**, and the buried side of the control panel rests on
  one structure.
- **The screen ran on predicted models**, and it rewards pockets with the physical character of a
  polyanion site — which is why it surfaces proteins known to bind other negatively charged
  ligands. That is the search behaving as specified, not a false positive mode.
- **Study O's shortlist is not a shortlist.** Its control arm failed to match on basic-residue
  count (59 of 75 controls censored against 7 of 75 candidates), so "3 of 75 versus 0 of 75"
  cannot be read as biology. Three designs that would fix it are pre-registered; none has run.
- **Three external tools are unpinned** (MMseqs2, Foldseek, APBS).
- **One chemistry dominates the ground truth** (134 InsP6, 1 InsP5), so the calibration speaks to
  InsP6 sequestration and not to InsP3/InsP4 sites.
- **No wet-lab validation has been performed.** Everything here is computational, and the
  candidate list is a set of hypotheses to test, not a set of findings.
