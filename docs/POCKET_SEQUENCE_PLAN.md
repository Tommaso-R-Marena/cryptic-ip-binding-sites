# Study Q9: does the sequence predict that a pocket forms, what kind, and what it binds?

Dated 2026-10-20, **before any sequence has been read from any structure for
this study and before any window has been passed to `composition_z`**. The
estimators in `cryptic_ip/detect/sequence.py` are unit-tested and have been
validated on synthetic pools (`results/sequence/VALIDATION.md`): they fire on a
planted helical face in 100% of windows and stay quiet on composition-matched
controls in 98–100%. They have never been run on a real sequence.

## The question, in three parts

> Are there patterns in the sequence that predict **pocket formation**, pocket
> **type**, or pocket **function**?

* **Q9a — formation.** Does a window's sequence say whether a pocket of
  substance forms there?
* **Q9b — type.** Among pockets, does the sequence say what *kind* of lining
  architecture the positions will show?
* **Q9c — function.** Among pockets, does the sequence distinguish the one an
  inositol phosphate occupies in a crystal?

## The decomposition, which is the whole design

A window of sequence factorises the same way a lining of positions does:

> **composition** (which letters) **⊕ order** (in what arrangement),

and the canonical null for order is permutation of the window's own letters —
the quotient by composition, exactly as the gap-permutation null in
`architecture.py` is the quotient by topology. That null holds length and
letter counts fixed and destroys only the arrangement.

**Every primary test in this study is on `composition_z`, never on a raw
measure.** A raw difference between arms is already explained by composition —
these sites are known to be basic, and pocket linings are known to be
hydrophobic — so a raw difference is reported in the secondary table as
composition-explained and is not a finding under this plan, whatever its
p-value. A non-zero z says the **arrangement** carries something the letters
alone do not. That is the only form of the question this study calls open.

## What is added to the module, and why

The existing statistics track **charge**, which is the right track for an
inositol-phosphate site. Pocket formation in general is a **hydropathy**
question — the amphipathic periodicity that puts one face of a helix against a
cavity and the other against solvent. The hydropathy track already exists;
three statistics over it do not, and are added before any data is read:
`helix_hydropathy_periodicity`, `strand_hydropathy_periodicity` and
`dominant_hydropathy_period_power`, defined exactly as their charge
counterparts with the Kyte–Doolittle scale already in the module.

Nothing else about the estimators changes.

## Data

* **Candidate proteins (Q9a, Q9b).** The 75 AlphaFold models of
  `results/specificity/candidates.csv`, with every pocket re-detected as study
  Q8 does it, and the chain sequence read from the same model file, so residue
  numbers and letters cannot disagree.
* **Deposited structures (Q9c).** The 137 annotated IHP copies study Q8's
  amendment 1 selects, with the site taken as the ligand's 4.5 Å contact shell
  in the ligand-bearing chain and the sequence read from the same file.
* **No new download.** Both sets are already fetched and cached by the
  within-protein workflow. This study adds no source.

## Windows

The three window definitions of `docs/SEQUENCE_PLAN.md` are inherited
unchanged, including its refusal: the `lining` window, which closes the gaps,
**must not** carry a periodicity statistic, because closing gaps destroys what
a period means. Periodicity statistics run on the `local` window — the
contiguous stretch the lining spans, capped at 40 residues.

Window length is a reported covariate in every comparison, and a comparison
whose arms differ in mean length by more than 2 residues is reported as
length-confounded regardless of its p-value.

## The three designs

**Q9a — formation.** Within each protein, the `local` windows of its
**top-ranked** pockets against the windows of its **bottom-ranked** pockets,
where rank is fpocket's own ordering, which is geometric and reads no sequence.
Top is rank ≤ `TOP_RANK` (5) and bottom is the lowest `TOP_RANK` ranks of the
same protein, matched on window length to within 2 residues. A protein
contributes only if it has both arms after matching.

*Also reported, where the data allow it:* windows of residues lining **no**
pocket at all, against lining windows of the same length. fpocket finds many
small pockets, so a protein may have too few non-lining stretches; the count of
proteins contributing is reported and the comparison is secondary for that
reason.

**Q9b — type.** Among the same pockets, the `local` windows of pockets the
architecture module calls `helix_face` against those it calls `none`. This
tests whether the **letters** anticipate what the **positions** say, which is
the one place the two lines of this project meet. The architecture class is
computed from positions alone and reads no letters, so the two sides are
independent by construction.

**Q9c — function.** The `local` window of each annotated IHP contact shell
against length-matched windows drawn from the same chain, outside the shell,
with the same number of windows per copy.

## Statistics, grouping and multiplicity

* Per window: `composition_z` with 2000 shuffles, seed 20261020.
* Per arm: the mean z, aggregated with a group bootstrap over **MMseqs2 30%
  clusters** for the candidate proteins and **`homology_group_strict`** for the
  deposited copies, with the project's standing five-group floor.
* **Primary family, fixed now:** four tests — `helix_charge_periodicity` and
  `helix_hydropathy_periodicity`, each on Q9a and Q9c — with Holm across those
  four. Q9b is reported with its own Holm across its two statistics, in its own
  family, because it asks a different question of a different contrast.
* Everything else in the module is secondary, reported in full, never promoted.
* The minimum detectable difference is computed from the same bootstrap and
  printed beside every row, as the within-protein study does.

## Decision rules

* A finding requires: Holm-corrected p below 0.05, the bootstrap interval
  excluding zero, the arm clearing the five-group floor, **and** the observed
  difference exceeding the minimum detectable difference. Any row failing the
  last two reads "no evidence", not "no".
* **A raw-measure difference is never a finding.** Only `composition_z` enters
  a primary test.
* A result that tracks window length is reported as length-confounded.
* The direction is pre-stated: a pocket-forming or site window is expected to
  show **more** periodicity than its own letters imply, so a positive z.

## What this cannot conclude, whatever it shows

* fpocket's rank is a geometric judgement of a predicted structure, not a
  measurement of a pocket's existence in solution. Q9a predicts fpocket, not
  nature.
* A crystal contact shell is one deposited conformation. Q9c predicts where an
  inositol phosphate sat in that crystal, not affinity and not specificity.
* Every arm is within-protein, so none of it generalises to sequences without
  a structure, and none of it is a sequence-only predictor of anything.
* The transfer ceiling of `results/architecture/ARCHITECTURE.md` stands over
  the periodicity statistics here as it does there.
