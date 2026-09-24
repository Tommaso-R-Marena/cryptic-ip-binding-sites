# The α-arrestin lead: plan

Pre-registered on 2026-09-24, before any α-arrestin structure, sequence or pocket was
fetched, aligned or docked for this study. What was known when it was written: the
learned ranking (run 35992676292, `docs/LEARNED_SCREEN_PLAN.md`) placed human ARRDC2
(Q8TBH0) at rank 5 and yeast ART5 (P53244) at rank 27 among unseen proteins, with
their top pocket's score, hull depth and pLDDT (`results/learned_screen/LEARNED_SCREEN.md`).
Nothing else about α-arrestins had been computed. Changes are dated amendments in
separate files written before the results they could affect are read.

## Why

β-arrestins and visual arrestin are established IP6 binders. α-arrestins share the
arrestin fold but no detectable sequence homology with them, so the learned model could
not have ranked them from memory. If their top pocket is the structural counterpart of
the β-arrestin IP6 site, and that site keeps its basic residues across orthologues, the
lead is worth an experiment. If not, it is one more high-scoring polyanion pocket.

## Definitions (fixed here, before looking)

- **Arrestin-fold protein**: a UniProt entry with a Pfam cross-reference to
  Arrestin_N (PF00339) or Arrestin_C (PF02752).
- **Classic arrestin** (visual and β-arrestins): an arrestin-fold protein whose UniProt
  recommended name contains "arrestin" and contains neither "domain-containing" nor
  "trafficking adapter" (human: SAG, ARR3, ARRB1, ARRB2).
- **α-arrestin**: every other arrestin-fold protein (human ARRDC1-5 and TXNIP; yeast
  ARTs; Dictyostelium arrestin-domain proteins).
- Annotations come from the UniProt REST API (`organism_id` 9606, 559292 and 44689,
  fields accession, gene, protein name, Pfam), fetched through `async_fetch` inside the
  workflow. The family lists are written to the output before any score is read.

## B1. Family-level test (existing scores)

Data: `learned-screen/results/learned_screen/proteins.csv.gz` (run 35992676292): every
scored protein, with its learned score, `seen` flag (MMseqs2 homologue among benchmark
proteins) and its MMseqs2 cluster (30 % identity, 50 % coverage).

- **Test.** Among unseen scored proteins of the three proteomes pooled, α-arrestins
  against all other unseen proteins. Statistic: ROC-AUC of the learned score for
  α-arrestin membership (the Mann-Whitney U divided by n₁n₀), with a 2,000-resample
  bootstrap over MMseqs2 clusters (seed 20260927). One-sided, as pre-registered: the
  family ranks **high**.
- **Decision.** *Supported* if the bootstrap 5th percentile exceeds 0.5 (one-sided
  bootstrap p < 0.05). *Not supported* otherwise. *Not evaluable* if the α-arrestins
  fall in fewer than 5 clusters (the interval is then not evidence). The asymptotic
  one-sided Mann-Whitney p-value is reported beside it, with the caveat that it treats
  paralogues as independent. B1 is the study's only hypothesis test (Holm is trivial).
- Reported per organism and per protein: learned score and rank percentile among
  unseen proteins.
- **Positive control (descriptive).** Classic arrestins' scores and rank percentiles
  among all scored proteins. They are expected to rank high but may be excluded from
  B1 as "seen" (homologous to benchmark proteins); their seen status is reported.

## B2. Structural mapping

- **References.** Benchmark entries (`benchmark-dataset/entries_grouped.csv`, run
  35949588200) whose UniProt accessions include a classic arrestin (by the definition
  above, applied to the benchmark's accessions through the UniProt API), and that hold
  a non-artefact IP copy. Each such copy is a **reference site**; its residues are the
  protein residues of the arrestin chain with a heavy atom within 4.0 Å of the copy.
- **Targets.** AlphaFold DB models (current version, `scripts/fetch_structures.py
  alphafold`) of every human and yeast α-arrestin, plus ARRDC2 and ART5 by name if the
  definition missed them (reported), plus the four human classic arrestins as mapping
  controls. Dictyostelium α-arrestins are listed in B1 only.
- **Alignment.** US-align 20241201 (`USalign target.pdb reference.pdb -outfmt 1`) of
  each target model against each reference chain. The reference with the highest
  TM-score normalised by the reference length is used. A reference-site residue maps to
  the target residue aligned to it with ":" (distance < 5 Å after superposition); other
  positions are unmapped.
- **Reported:** TM-scores (both normalisations), RMSD, the per-residue mapping of every
  reference site, and the fraction of site residues mapped.
- A target with TM-score < 0.5 against every reference has **no mapped site**.

## B3. Pocket correspondence

- **Top learned pocket.** `top_pocket_residues` of the target in `proteins.csv.gz`
  (UniProt numbering, the model version screened).
- **Overlap metric.** Jaccard index between the mapped site residues and the top
  pocket's residues. **Overlap** when Jaccard ≥ 0.20 for at least one mapped site; the
  site with the highest Jaccard is the **lead site**. The distance between the site
  and pocket residue centroids (Cα) is reported as a secondary measure.

## B4. Conservation

- **Homologue set.** The UniProtKB members of the target's UniRef50 cluster (UniRef
  search `uniprot_id:<acc> AND identity:0.5`, then UniProtKB `uniref_cluster_50:<id>`),
  deduplicated at 100 % identity, at most 500 sequences in accession order, fetched
  through `async_fetch`. Fewer than 10 sequences: conservation *not evaluable* for that
  target.
- **Alignment.** MAFFT 7.526, `--auto`.
- **Score.** At each lead-site position that is K, R or H in the target, the fraction
  of homologues (target excluded; gaps count as not basic) with K, R or H in that column.
- **Conserved** when at least 3 lead-site positions are basic in the target and at
  least 3 of them have a basic fraction ≥ 0.80.

## B5. Docking

The protocol of `docs/REDOCKING_PLAN.md` (receptor preparation, ligand from the CCD,
primary protonation, box rule, Vina 1.2.7, exhaustiveness 32, 20 poses, energy range
5, seeds 1-3). Ligand: IP6 (CCD IHP). Decoy polyanion: ATP (CCD ATP), with the same protonation rule
applied to its phosphates (−3 in the primary state).

- **Sites.** For ARRDC2 (Q8TBH0) and ART5 (P53244): the lead site (or, with no mapped
  site, every mapped site; with none, none), and the top learned pocket. Box centre: the
  centroid of the site's (or pocket's) residue Cα atoms in the AlphaFold model.
- **Positive controls.** (a) IP6 redocked into every arrestin-IP6 reference structure
  (crystal receptor, crystal site centroid, the redocking protocol); (b) IP6 docked
  into the AlphaFold models of human ARRB1, ARRB2 and SAG at their own mapped sites.
- **Negative controls.** For each of ARRDC2 and ART5, 5 random pockets of the same
  screened model (from the screen's pocket table, seed 20260928): confident (mean pLDDT
  ≥ 70), and with Jaccard < 0.10 against the lead site and the top pocket.
- **Specificity (descriptive).** ATP docked into the lead site and the top pocket.

**Protocol validity for B.** Criteria 3 and 4 below use docking. They are evaluable
only if the positive-control crystal redocks (a) recover the crystal pose: mean over
reference sites of the three-seed top-pose success at 2.0 Å ≥ 0.5. Otherwise both are
*not evaluable*.

## B6. Decision rule (per protein: ARRDC2, ART5)

The lead is **supported for experiment** only if all four hold:

1. **Overlap**: the mapped site overlaps the top learned pocket (B3).
2. **Conservation**: the lead site's basic residues are conserved (B4).
3. **Convergence**: at the lead site, the top 5 poses of each of the 3 seeds (15
   poses), clustered in score order by leader clustering at 2.0 Å symmetric RMSD, have
   a largest cluster holding ≥ 50 % of them.
4. **Scores**: the best IP6 score at the lead site (lowest over seeds) is at least as
   good as the weakest positive control (≤ the highest positive-control best score) and
   better than every negative-control pocket's best score for the same protein.

Otherwise the lead is **not supported**, with the failing criteria named (a criterion
that is not evaluable counts as failing, and is named as not evaluable). Docking scores
for a −9 polyanion are weak evidence under any outcome; the report says so.

## Dossiers

One page per protein (ARRDC2, ART5): mapped residues, conservation, the contact table
of the best docked IP6 pose (residues with a heavy atom within 4.0 Å), the controls, the
verdict with each criterion, and what a wet-lab test would be (ITC or fluorescence
binding with IP6; charge-reversal mutants of the mapped basic residues).

## Execution and outputs

`.github/workflows/arrestin.yml`; results printed between `BEGIN_ARRESTIN_JSON` and
`END_ARRESTIN_JSON` and extracted into `results/arrestin/`.
