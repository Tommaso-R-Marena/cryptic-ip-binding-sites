# Redocking benchmark for inositol phosphate ligands: plan

Pre-registered on 2026-09-24. When this plan was committed, no structure had been
fetched, prepared or docked for this study, and no docking of any inositol phosphate
had been run in this project. Changes are dated amendments in separate files
(`docs/REDOCKING_PLAN_AMENDMENT_n.md`), each written before the results it could affect
are read. Exploratory analyses go to the append-only ledger
`results/redocking/ledger.jsonl`, never into a decision.

## Question

When each crystal inositol phosphate (IP) ligand is put back into its own site, how
often does docking recover the crystal pose within 2 Å, and how does that depend on
how buried the site is? A second question follows from the proteome screen, which
runs on AlphaFold models: does docking into a predicted, ligand-free model recover
the crystal pose at all?

## 1. Data

**Entries.** The 367 PDB entries of the pre-registered IP benchmark table (run
35949588200, `benchmark-table/table.csv.gz`, SHA-256
a21d92fd0f57306848c260d1f163017f0bf1e82a5dfcf1f3be69e0c20b698a1d): exactly the
entries with at least one row in that table. Metadata (resolution, method, release
date, UniProt accessions, `homology_group`, `homology_group_strict`) come from
`benchmark-dataset/entries_grouped.csv` of the same run. Coordinates are fetched again
from RCSB with `scripts/fetch_structures.py rcsb` (legacy PDB preferred, mmCIF
otherwise), as the benchmark did.

**Copies.** Every IP copy in each entry's first model, found from coordinates by
`find_ligand_instances` (all phosphorylated inositols, InsP1-InsP8), with its burial
class from `compute_ligand_burial` (256 SASA points, exactly the benchmark's setting):
cryptic (relative SASA ≤ 0.12), semi-cryptic (≤ 0.25), surface, or crystal artefact
(fewer than 8 protein heavy-atom contacts within 4.5 Å). Recorded for every copy: IP
species (phosphorus count and component id), resolution, method, burial class,
relative SASA, strict group, interface (contacts from ≥ 2 chains), mean occupancy.

**Flags** (recorded, never used to drop a copy except where stated):

- *partial occupancy*: mean ligand-atom occupancy < 1.0;
- *alternate locations*: any ligand atom with an altloc; the highest-occupancy
  altloc is used (the project's parser does this);
- *incomplete*: the copy has fewer heavy atoms than its CCD component;
- *metal*: an Mg, Zn, Ca, Mn or Fe ion within 3.0 Å of any ligand heavy atom;
- *symmetry contact*: a protein atom of a crystal-symmetry mate (gemmi neighbour
  search over the unit cell, images other than the identity) within 4.5 Å of any
  ligand heavy atom. Not applicable to cryo-EM entries (no cell).

**Exclusions**, each counted with its reason:

1. crystal artefacts;
2. copies whose CCD entry has no parseable SMILES, or whose SMILES leaves a
   stereocentre unassigned;
3. copies whose configuration differs from the CCD's (section 3);
4. copies in entries whose receptor cannot be prepared (section 2);
5. copies not reached within a shard's time budget (section 9).

**Selection (compute cap).** At most **one copy per (entry, burial class)**. Among the
eligible copies of an entry and class, complete copies are preferred; among those, the
first by (chain, residue number, insertion code). This is fixed before any docking and
depends on no outcome. Every unselected copy is counted.

**Primary analysis set.** Selected copies that are complete and whose configuration
matches the CCD. Incomplete copies are docked and reported as a separate flagged
stratum, outside every decision.

## 2. Receptor

- **Strip** every non-polymer atom (waters, ions, all ligands, all other IP copies,
  cofactors). Keep every polymer chain of the deposited asymmetric unit. Selenomethionine
  is written as methionine; other modified residues are removed and counted.
- **Protonate** with PDB2PQR 3.6.1 (AMBER force field) and PROPKA 3.5.1 at pH 7.4
  (`--titration-state-method=propka --with-ph=7.4`). His tautomers and charges (HID,
  HIE, HIP) come from PROPKA's pKa estimates and PDB2PQR's hydrogen-bond optimisation.
  The counts of HID/HIE/HIP are recorded per structure. A structure PDB2PQR cannot
  process is excluded, with the error.
- **Type** for AutoDock with `cryptic_ip.docking.receptor`: polar hydrogens kept (HD),
  non-polar hydrogens merged; aromatic ring carbons A; histidine and nucleobase ring
  nitrogens without a hydrogen NA; oxygens OA; sulfur SA; phosphorus P. Charges are
  PDB2PQR's AMBER charges with merged hydrogen charges added (only AD4 uses charges).
  Meeko's receptor builder is not used for the receptor: it rejects residues with
  non-ideal geometry, which would drop structures non-randomly (it failed on a
  synthetic test peptide during development). Meeko prepares every ligand.
- **Metals sensitivity arm** (pre-declared): for copies flagged *metal*, the same
  docking is repeated with those ions (only those within 3.0 Å of the copy) kept in
  the receptor, typed by element (Mg, Zn, Ca, Mn, Fe) with charge +2. Stripping a
  bridging metal can make the true pose unreachable; the metal stratum is reported on
  its own.

## 3. Ligand

- **Source.** The copy's CCD entry (`files.rcsb.org/ligands/download/<ID>.cif`), fetched
  through `async_fetch`. SMILES preference: CACTVS SMILES_CANONICAL, OpenEye
  SMILES_CANONICAL, CACTVS SMILES, OpenEye SMILES. Never crystal coordinates.
- **Configuration check.** The crystal copy is rebuilt as a molecule (distance-based
  connectivity, bond orders from the CCD template with `AssignBondOrdersFromTemplate`,
  stereo from its 3-D coordinates). Its neutral isomeric canonical SMILES must equal
  the template's. This checks myo against scyllo and every other inositol
  configuration. A mismatch excludes the copy (exclusion 3). Incomplete copies cannot
  be checked and are flagged.
- **Protonation.** Primary state: every terminal phosphate oxygen deprotonated, then
  floor(n_P / 2) protons put back, one per phosphorus, on the phosphorus atoms of lowest
  RDKit canonical rank, each on that phosphorus's single-bonded terminal oxygen of
  lowest rank. Net charges: InsP6 −9, InsP5 −8, InsP4 −6, InsP3 −5 (InsP6 at pH ~7.4
  carries about −9 to −10). Sensitivity state: fully deprotonated (InsP6 −12).
  Vina's scoring has no electrostatic term and ignores partial charges; the two states
  differ in Vina only by which oxygens carry a donor hydrogen. They matter more for AD4.
- **Start pose.** RDKit ETKDGv3 with the docking seed as random seed; hydrogens
  explicit; a random rigid rotation (from the same seed) about the centroid; centroid on
  the box centre. The start pose's symmetry-corrected RMSD to the crystal copy is
  recorded and must exceed 2.0 Å; otherwise the next seed is used (up to 20).
- **PDBQT** by Meeko 0.8.0 (`MoleculePreparation(rigid_macrocycles=True)`: rigid ring,
  rotatable C-O, O-P and P-OH bonds, Gasteiger charges).

## 4. Search

AutoDock Vina 1.2.7 (Python bindings), all runner cores per run.

- **Box.** Cubic, centred on the crystal copy's centroid. Side = max(22 Å, d_max + 16 Å),
  where d_max is the largest heavy-atom distance in the ligand's seed-1 ETKDG
  conformer: the ligand's extent plus 8 Å on each side, independent of orientation.
  The same rule is used for every box in this plan and in docs/ARRESTIN_PLAN.md.
- **Primary.** Vina scoring, exhaustiveness 32, 20 poses, energy range 5 kcal/mol,
  seeds 1, 2 and 3, primary protonation state.
- **Secondary scoring.** Vinardo (seed 1). AD4 (seed 1) with autogrid4 4.2.9 maps
  (spacing 0.375 Å, the same box). If autogrid4 does not install or a map run fails,
  the AD4 arm is reported as not run for the affected copies, with the reason; nothing
  replaces it.
- **Protonation sensitivity.** Fully deprotonated ligand, Vina, seed 1.
- **Tertiary (site finding plus docking).** Per entry, for its first selected primary
  copy: fpocket 4.2.3 on the ligand-free structure through the project's analyser (the
  benchmark's code path). A box (same size rule) on each of fpocket's top 3 pockets by
  fpocket's own rank; Vina, seed 1. The best-scoring pose over the three boxes "lands in
  the true site" when its heavy-atom centroid is within 4.0 Å of the centroid of any
  non-artefact IP copy of the entry. Also recorded: whether any of the three pockets is
  a positive pocket (the benchmark's labelling rule).

## 5. RMSD

Heavy atoms only, **no superposition** (the pose is judged in the receptor frame),
implemented in `cryptic_ip.docking.rmsd`:

- the core (all heavy atoms except terminal phosphate oxygens) is matched through every
  graph automorphism of an element-only graph (bond orders, charges and hydrogens
  removed, so protonation cannot break a match);
- for each core mapping, each phosphorus's terminal oxygens are assigned by exhaustive
  permutation (exact, because different phosphorus atoms' oxygens contribute
  independent squared deviations);
- an incomplete crystal copy is matched as a substructure; its RMSD runs over the atoms
  it has.

RDKit's `CalcRMS(symmetrizeConjugatedTerminalGroups=True)` is not relied on, because
its coverage of P-O depends on the version; tests check this implementation against
brute-force enumeration. Tests show: a relabelled identical pose gives 0 Å while naive
atom-order RMSD gives more; swapping terminal oxygens gives 0 Å; a rigid 3 Å translation
gives 3 Å; a pose from a different conformer gives the value computed by brute force.

**Secondary metric:** RMSD over phosphorus atoms only (robust to ring flips), under the
same automorphisms.

## 6. Metrics

- **Primary.** Per copy, the fraction of the three seeds whose **top-ranked pose** is
  within 2.0 Å (0, 1/3, 2/3 or 1).
- **Secondary.** Best of the 20 poses within 2.0 Å (per seed, averaged over seeds);
  success at 1, 2 and 3 Å; the Spearman correlation between pose score and RMSD within
  each run; phosphorus-only success at 2.0 Å.
- **Estimands.** Per copy (each copy weighted equally) and per group (each strict
  homology group weighted equally: the mean of within-group means).
- **Intervals.** 2,000 group-bootstrap resamples over `homology_group_strict`
  (`protocol.group_bootstrap_weights`, seed 20260926); percentile 95 % intervals.
- **Strata.** Burial class (cryptic, semi-cryptic, surface); metal or not; interface or
  single-chain; IP species (InsP3, InsP4, InsP5, InsP6, other); resolution (≤ 2.0,
  2.0-2.5, 2.5-3.0, > 3.0 Å); X-ray or cryo-EM.
- **Small strata.** A stratum with fewer than 5 strict groups is reported with its
  interval and the words "fewer than 5 independent groups: this interval is not
  evidence". The buried (cryptic) class is expected to hold about 6 families, one of them
  (ADAR2 and relatives) holding most copies; its interval can at best separate a success
  rate near 0 from one near 1, and cannot support a claim about buried sites in general.

## 7. Decisions

Four primary questions, Holm-corrected together (p-values from the group bootstrap,
two-sided against the stated null). Each gets one of the labels below.

- **R1, protocol reliability** (all primary copies, group estimand).
  *Reliable*: lower bound ≥ 0.5. *Unreliable*: upper bound < 0.5. Otherwise
  *inconclusive*. p-value against 0.5.
- **R2, burial dependence**: success(cryptic) − success(surface), group estimand,
  both strata resampled together. *Not evaluable* if either stratum has fewer than 5
  strict groups. Otherwise *buried harder* (upper bound < 0, Holm p < 0.05),
  *buried easier* (lower bound > 0, Holm p < 0.05), else *no difference detected*.
- **R3, AlphaFold cross-docking** (section 8): success in the model, group estimand.
  *Trustworthy*: lower bound ≥ 0.5. *Not trustworthy*: upper bound < 0.5. Otherwise
  *inconclusive*. Reported per burial class as well (descriptive).
- **R4, discrimination**: does the Vina score separate the true site from the decoy
  pocket? ROC-AUC with true-site top scores as positives and decoy top scores as
  negatives (lower score ranks higher), group bootstrap. *Discriminates*: lower bound >
  0.5 and Holm p < 0.05. *Does not*: upper bound < 0.6. Otherwise *inconclusive*.

## 8. Controls and decompositions

- **Scoring or sampling failure.** For each primary copy whose seed-1 top pose is
  > 2.0 Å: the crystal pose (the docked ligand's atoms placed on the crystal
  coordinates, hydrogens rebuilt) is scored and locally minimised in the prepared
  receptor. If the minimised crystal pose scores better (lower) than the top docked
  pose, the failure is **sampling**; otherwise **scoring**. Fractions are reported, with
  the minimised pose's RMSD from the crystal.
- **Discrimination decoy.** The highest-ranked fpocket pocket of the ligand-free
  structure whose label is negative under the benchmark's rule (≤ 5 % overlap with
  every IP copy). Vina, seed 1, same box size. Feeds R4.
- **Noise.** Seed-to-seed standard deviation of the top score; agreement of the three
  seeds' success.
- **Protonation, Vinardo, AD4, metals.** Paired against the primary seed-1 outcome on
  the same copies; differences with group-bootstrap intervals; descriptive.
- **AlphaFold cross-docking.** For each primary copy, the chain with most contacts to
  the copy is matched to an AlphaFold DB model among the entry's UniProt accessions:
  the model whose sequence aligns to the chain with the highest identity, at least 90 %
  over at least 50 % of the chain (Biopython global alignment). Binding-site residues:
  chain residues with any heavy atom within 6.0 Å of the copy that are aligned to model
  residues. The model is superposed on the crystal by those residues' Cα atoms (Kabsch),
  prepared as a receptor exactly as above (the model is a single chain; no metals), and
  docked with the same box, Vina, seed 1. RMSD is against the crystal ligand in the
  crystal frame. The site Cα RMSD after superposition is recorded. Copies without a
  usable model (no accession, no model, fragmented model of a protein longer than 2,700
  residues, identity below the threshold, fewer than 3 site Cα pairs) are counted.
  Interface sites lose their partner chains in a single-chain model; they are
  reported as a stratum.

## 9. Execution

- One workflow, `.github/workflows/redocking.yml`: a copy census (fetch, copies, flags,
  CCD, selection), then docking in at most 30 shards per job, each with a 5 h 20 min
  budget after which unstarted copies are recorded as "not reached". The primary arm
  and each group of secondary arms run as separate jobs so no job exceeds 6 h.
- If more than 5 % of selected copies are not reached, the not-reached copies are
  docked in a second dispatch of the same workflow with the same code; the decision
  uses all copies. This rule is about completeness and depends on no outcome.
- The result JSON is printed between `BEGIN_REDOCKING_JSON` and `END_REDOCKING_JSON`,
  the per-copy table (gzip, base64) between `BEGIN_REDOCKING_COPIES_B64` and
  `END_REDOCKING_COPIES_B64`. Results in `results/redocking/` are extracted from the log
  by `scripts/extract_log_block.py`, never retyped.

## 10. Outputs

`results/redocking/REPORT.md` and `redocking.json`, the per-copy table, a forest-plot page
in the style of `scripts/benchmark_report_page.py`, and full accounting: copies found,
excluded (by reason), unselected, attempted, failed (by reason) and docked.

## What each outcome would mean

- R1 *reliable*: redocking recovers IP poses; later docking results (study B) rest on a
  validated protocol. *Unreliable*: docking poses for these polyanions should not be
  trusted, and docking scores in study B are weak evidence at best.
- R2 is expected to be *not evaluable* (about 6 buried families). That outcome is
  reported as such, not as "no difference".
- R3 *not trustworthy*: docking into AlphaFold models cannot confirm screen candidates.
- R4 *does not discriminate*: a good docking score at a candidate site says little.
