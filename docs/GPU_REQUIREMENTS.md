# What this project needs a GPU for

Written 2026-10-01. Every figure marked **measured** comes from a run in this
repository; every figure marked **estimate** does not and needs a probe before anyone
commits time to it. The distinction matters because the one thing we have already
measured overturned a plan.

## The measurement that forces the question

**Measured.** Cofold probe run 36729986766: Boltz-2 v2.2.1, one 100-residue protein with
IP6 (`ccd: IHP`), on 4 CPUs and 15 GB RAM, at the cheapest settings the model offers —
`--diffusion_samples 1`, `--recycling_steps 0`, `msa: empty`. It **did not complete in
5,400 s**. The probe stops at the first length that fails, so 300 and 600 residues were
never attempted. Boltz-2 installed in 1 m 39 s, so this is compute, not packaging.

**Measured.** The candidates are not small. Taking each candidate's highest pocket residue
number as a lower bound on its length: median **436**, mean 518, maximum **2,352**. Of 75
candidates, **52** are ≥ 300, **19** are ≥ 600 and **5** are ≥ 1,000 residues.

A 100-residue protein is smaller than almost every candidate. CPU co-folding is not
slow here, it is impossible.

## Priority 1 — co-folding: the only method that answers the actual question

Studies A–L never put an IP6 molecule into a candidate pocket. That is the gap:

- **Docking cannot do it.** Study F: top-pose success **0.110** over 242 crystal copies,
  with **230 of 242** failures attributable to scoring rather than search. Study G: a
  **16×** search budget moved success by **+0.019**. A screened-Coulomb re-rank recovered
  **+0.054** (permutation p = 0.0099) and no more.
- **Geometry alone cannot do it.** Study N rigidly transplanted 131 crystal IP6 sites onto
  each pocket. The guard passed only weakly (AUC **0.599 [0.546, 0.652]**), and the
  placements were **not buried**: median buried fraction **0.000**, zero for 22 of 35.

Co-folding (Boltz-2, Chai-1, AlphaFold3) is the right tool class for "does this fixed
protein bind this fixed ligand", and it is the one we cannot run without a GPU.

**Workload.** Four arms, so that "has a basic pocket" cannot masquerade as "binds IP6":
75 candidates, 75 matched controls, 62 annotated binders, plus ligand decoys (ATP,
sulfate, IP3) against the same receptors. Roughly **280 protein–ligand systems**, times 2–3
seeds for run-to-run variance: **~600–850 predictions**.

**The binding constraint is VRAM, not throughput.** Memory in AF3-class models grows
steeply with token count, so the 19 candidates ≥ 600 and especially the 5 ≥ 1,000 residues
are the problem cases, not the bulk. Either 80 GB, or those are truncated to a domain
around the pocket — and truncation is a modelling decision that must be pre-registered,
not improvised when a job runs out of memory.

**MSAs are a separate cost.** The probe used `msa: empty`, which degrades accuracy
substantially. Real alignments need either a GPU-side MSA server or precomputed
alignments on CPU plus storage. Whether MSAs are used must be stated in the plan, because
it changes what the result means.

## Priority 2 — molecular dynamics: the biggest standing weakness in the whole project

Every candidate is **a single apo AlphaFold conformer**. No induced fit, no alternative
rotamers, no evidence the pocket exists as anything but a static prediction. This is
precisely where an informed reader will push, and it is the caveat that most limits what
the candidate list can claim.

The code already exists: `cryptic_ip/validation/md_validation.py` (OpenMM),
`scripts/run_md_pilot_validation.py`, and `scripts/hpc/slurm_submit.sh`.

**Workload.** Pocket stability and cryptic-pocket opening across the shortlist, with and
without IP6 placed. A defensible design is 3 independent replicates × 100 ns per system.
**Estimate:** a ~50,000-atom system runs on the order of 100–250 ns/day on a current
datacentre GPU — to be confirmed by a short probe, not assumed. This is where GPU-hours
actually accumulate: 20 systems × 3 × 100 ns is the large ask, far larger than co-folding.

## What does *not* need a GPU — so the ask stays credible

- **Study J, flexible-receptor redocking** (pre-registered in `docs/FLEXIBLE_PLAN.md`,
  never run). AutoDock Vina is CPU-only. This is a CPU-core request, not a GPU one.
- **Model training.** Extra trees and XGBoost on 78,920 pockets × 39 descriptors is a CPU
  job. Asking for a GPU for this would weaken the case.
- **fpocket, SASA, MMseqs2, MAFFT, metapredict, the group bootstraps.** All CPU.
- **Study N itself.** It ran end to end in **15 minutes** on free shared runners.

## What a GPU will *not* fix

Worth saying out loud, because compute is often offered as an answer to these:

1. **Specificity.** Study C: re-weighting the proteome ranking by P(IP | site) did **not**
   improve it — ΔAUC **−0.015 [−0.044, 0.015]**. A better structure predictor does not make
   the descriptors distinguish IP6 from PAPS, UDP-sugar, ATP or sugar-phosphate sites, and
   30 of the 75 candidates are already flagged as explained by another phosphate- or
   sulfate-rich ligand.
2. **The control-matching defect.** Study N's controls were matched on pLDDT and hull depth
   but never on basic-residue count, while the fit score structurally requires ≥ 4 basic
   residues. That is a design fix, not a compute problem.
3. **Out-of-sample validation for buried IP6 sites specifically.** Buried-site positives
   occupy only **5 strict homology groups**, one holding **87%** of them, and none landed in
   the temporal holdout. No amount of GPU adds independent PDB entries.
4. **Wet-lab validation.** Nothing computational substitutes for it.

## The concrete ask

| item | request | basis |
|---|---|---|
| GPU | 1× 80 GB (A100/H100 class); 40 GB workable only if the ≥1,000-residue tail is truncated | VRAM, from the measured length distribution |
| Co-folding **probe** | ~2 GPU-hours, before anything else | the CPU probe cost 90 minutes and killed a study that would have wasted far more |
| Co-folding production | order 50–200 GPU-hours (**estimate**, set by the probe) | ~600–850 predictions |
| MD | order 500–1,500 GPU-hours (**estimate**) | systems × replicates × 100 ns |
| Storage | modest for co-folding; MD trajectories dominate, ~1–5 GB per 100 ns at a sane write frequency (**estimate**) | |
| CPU alongside | MSA generation, and study J's flexible redocking | both CPU-bound |

**Ask for the probe first.** That is the lesson of the CPU probe: it converted an
assumption into a measurement for 90 minutes of compute, and the measurement said do not
build the study. The same two GPU-hours will tell us what the production number actually
is instead of defending an estimate.

## Pre-registration, non-negotiable

Study M gets pre-registered **before any candidate is scored**, with the **calibration
guard first**: if the co-folding score cannot separate annotated binders from matched
controls, it cannot classify candidates either, and the study reports that and stops.
Study N's workflow already implements this as a hard gate — the candidate arm does not run
unless the guard passes — and study M should reuse that shape.
