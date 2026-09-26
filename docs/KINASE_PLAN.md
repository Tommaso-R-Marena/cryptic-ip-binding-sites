# The inositol-phosphate kinases, and the missing cosubstrate: plan

Pre-registered on 2026-09-26, before any cofactor-holo receptor was built or docked.

## What was known when this was written, and what was read first

This plan was written **after** reading the IPK-superfamily slice of study A's
`results/redocking/copies.csv`. That slice is therefore a **re-analysis of data already
in hand (K0, K3 below), not a test**, and is reported as such. What is new, and blind, is
the holo arm of K1: no cofactor-bearing receptor has been built or docked.

Read first, and the reason this plan exists:

- The census holds **18 copies in the inositol-phosphate kinase (IPK) superfamily** —
  *Entamoeba* IP6KA (N9UNA8, the IP6 kinase itself), human PPIP5K2 (O43314), IPMK
  (Q8NFU5), ITPKA (P23677) and ITPKC (Q96DU7).
- **Every one of the 18 fails**: top-pose success 0 of 18, against 0.110 [0.046, 0.193]
  overall. Every failure is classified `scoring`, not sampling, and several have a
  near-native pose somewhere in the list (best-of-20 up to 1.0).
- **Every one is classified `surface`**, which is the burial class with the lowest
  best-of-list ceiling in study A (0.279).
- The diphosphoinositols are in the census too and also fail: InsP7 0 of 4, InsP8 0 of 2.

## The hypothesis

These are **catalytic** sites. In the crystal they hold the inositol phosphate *and* a
nucleotide cosubstrate with its catalytic metals: the IP6 that IP6K phosphorylates sits
against an ATP γ-phosphate and Mg²⁺. The redocking protocol strips every heteroatom
except metals within 2.8 Å, so it asks Vina to place a −9 ligand in a site whose
positive counter-charge has been deleted. Under that receptor the true pose is not
merely hard to rank — it is being scored in the wrong physical system.

This is a *different* explanation from the rigid receptor of `docs/FLEXIBLE_PLAN.md`.
There the missing ingredient is side-chain freedom; here it is a **missing molecule**.
The two are independent and this study does not test the other.

## Data

- **Copies.** Every primary-set copy with an outcome in study A (250) whose deposited
  structure contains a **nucleotide cofactor** — a residue whose CCD id is one of ATP,
  ADP, AMP, ANP, ACP, AGS, GTP, GDP, GNP, or ADX — with any cofactor heavy atom within
  **6.0 Å** of any ligand heavy atom. Which copies qualify is **not known** when this is
  written; the count is an output, not an input.
- **Pre-specified stratum.** The IPK superfamily, by the UniProt accessions listed above
  plus any census accession annotated with Pfam PF03770 (IPK domain).
- **Receptor, ligand, box, starting pose, RMSD, protonation, exhaustiveness.** Exactly as
  `docs/REDOCKING_PLAN.md`, through the same code, at exhaustiveness 32, seeds 1–3.
- **Pose lists.** Up to 40 poses within 10 kcal/mol, as `docs/RERANK_PLAN.md`, so the arms
  are comparable to studies F and G.

## Arms

| arm | receptor | seeds | copies |
|---|---|---|---|
| apo | as studies A/F/G: polymer plus metals within 2.8 Å | 1, 2, 3 | the qualifying copies — **reused, not re-docked** |
| holo | the same, plus the nucleotide cofactor as rigid receptor atoms | 1, 2, 3 | the qualifying copies |

**How the cofactor enters the receptor.** Its heavy atoms are appended to the receptor
PDBQT as rigid HETATM records. `cryptic_ip/docking/receptor.py` is **not edited** — a
finished study's module — so the appending lives in this study's own code, with a test
pinning the no-cofactor case to a byte-identical copy of the apo receptor.

**Atom typing, fixed here.** C→C, N→NA, O→OA, P→P, S→SA, halogens and metals to their
own types. Nitrogen and oxygen are typed as hydrogen-bond acceptors, which is what they
overwhelmingly are in a nucleotide. Aromatic carbons of the base are typed C rather than
A, so the base's stacking contribution is slightly under-rewarded; this is a stated
approximation, not a tuned choice.

**Charges, fixed here.** Gasteiger charges computed from the cofactor's CCD template and
matched to the crystal atoms by name; 0.0 where no template is available. Vina ignores
charges, so K1 and K2 are unaffected either way, and K4 reports how many copies carried
real charges — with none, K4 is declared not evidence about electrostatics.

**What the holo arm changes, and the decomposition.** The apo arm reused from study F has
neither the cofactor nor the metals. The holo arm adds **both**, because that is the
complex the crystal actually holds. To keep the two contributions separate the report
quotes study A's existing **metals-only** arm (`metals_success`) on the same copies rather
than re-docking it.

**A stated limitation of the primary engine.** Vina's scoring function is typed and
distance-based; it does **not** read partial charges. So in the Vina arm the cofactor acts
as a shaped, typed occluder — it restores the *shape and hydrogen-bonding* of the real
site but not its electrostatics. The charge part of the hypothesis is therefore tested
separately, in K4, through study F's screened-Coulomb term, where the cofactor's
phosphate oxygens carry the same formal charges the IP ligands are given.

## Decisions

- **K1, primary.** Paired difference in top-pose success at 2 Å, **holo − apo**, over
  qualifying copies, resampling `homology_group_strict` (2,000 resamples, seed 20261004).
  - **better** if the 95 % lower bound > 0; **worse** if the upper bound < 0.
  - **no material change** if the interval lies inside ±0.05, study G's margin.
  - **inconclusive** otherwise. **Not evaluable** below 5 groups.
  - Because the apo arm is known to be 0 of 18 on the IPK stratum, on that stratum the
    difference is simply the holo success rate. This is stated rather than hidden.
- **K2, the ceiling (secondary).** Best-of-list rate in each arm and their paired
  difference. A cofactor that raises the ceiling has changed what is *sampled*; one that
  raises top-pose success without moving the ceiling has changed what is *ranked*. These
  are different claims and are reported separately.
- **K3, the diphosphoinositols (descriptive, read first).** Top-pose success and
  best-of-list for InsP7 and InsP8 pooled as one PP-IP stratum (6 copies, 6 groups), with
  intervals, beside the other species. Descriptive because it was read before this plan.
- **K4, electrostatics with the cofactor present (secondary).** Study F's screened-Coulomb
  re-ranking recomputed with the cofactor's atoms included as fixed point charges, against
  the same re-ranking without it. Holm across K1 and K2 only; K3 and K4 are not corrected.

## What this can and cannot support

- It can say whether deleting the cosubstrate is what breaks docking at IP-kinase active
  sites, which is the one explanation study A's own diagnostics (scoring failure, pose
  present in the list) point at and no study here has tested.
- It cannot generalise to the screen's candidate pockets, which are not catalytic sites
  and mostly have no cosubstrate to restore. A win here narrows where the protocol may be
  trusted; it does not rescue the proteome ranking.
- It cannot make IP6K a novel finding: IP6K's IP6 site is known. The family is used here
  as a **hard, well-characterised positive control**, which is exactly what it is good for.

## Execution and outputs

- **Workflow.** `.github/workflows/kinase.yml`: a cofactor census, a sharded docking
  matrix, then a report.
- **Code.** `scripts/kinase.py`.
- **Results.** Printed between `BEGIN_KINASE_JSON` and `END_KINASE_JSON`, extracted into
  `results/kinase/`.
