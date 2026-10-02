# Study J: a flexible receptor makes redocking worse

**Run.** Flexible receptor run 36888975469, Report job 110917523300 (run attempt 2),
commit `67192a2`, branch `claude/ip-binding-studies-26p6u4`. Census from redocking run 36020459246. Started
2026-10-01 16:02 UTC, finished 2026-10-02 16:02 UTC.

**Plan.** `docs/FLEXIBLE_PLAN.md`, as amended by `docs/FLEXIBLE_PLAN_AMENDMENT_1.md`
(both arms docked fresh through one Meeko receptor), `_2.md` (what to do with residues
Meeko refuses) and `_3.md` (the per-copy compute cap). Seed 20261003,
2000 resamples of strict homology groups, exhaustiveness 32,
flexible side chains within 4.0 A of the ligand and at most 8.

Every figure below is read from `flexible.json` by script. Nothing is transcribed by hand.

## What the study found

**The flexible receptor is worse than the rigid one, on both measures.**

| | flex | rigid | paired difference, per group | per copy | Holm p | decision |
|---|---|---|---|---|---|---|
| J1 top pose within 2 A | 0.026 | 0.096 | -0.078 [-0.140, -0.031] | -0.070 [-0.106, -0.036] | 0.001 | **worse** |
| J2 best-of-list ceiling | 0.179 | 0.367 | -0.170 [-0.243, -0.099] | -0.188 [-0.266, -0.115] | 0.0 | **worse** |

Both on 209 copies in 29 strict homology groups, of the
272 the plan selected. Both intervals exclude zero on both estimands.

This is the opposite of what the study was built to find. Letting the pocket's side chains
move does not recover study A's failures; it loses about 7 points of
top-pose success and about 19 points of ceiling. The extra freedom
enlarges the search space faster than the scoring function can exploit it.

## J3, by burial class (exploratory, not Holm-corrected)

| burial class | decision | difference, per group | copies | groups | evidence |
|---|---|---|---|---|---|
| cryptic | not evaluable | -0.156 [-0.500, +0.177] | 11 | 4 | **no** |
| semi_cryptic | worse | -0.137 [-0.226, -0.053] | 44 | 11 | yes |
| surface | worse | -0.022 [-0.044, -0.005] | 154 | 22 | yes |

The cryptic stratum has 4 strict groups against the
pre-registered floor of 5, so it is **not evidence** whichever way it points. The two strata
that are evidence both say worse, and semi-cryptic loses more than surface.

## The receptor pipeline, audited separately

The rigid arm scores **0.096** under this study's Meeko-prepared
receptor against study F's **0.110** under the project-prepared one, over
209 copies. Amendment 1 required this comparison so that a change of receptor
preparation could not masquerade as an effect of side-chain freedom. The gap is about
14 thousandths and in the
same direction as study F's own interval, so the receptor change does not explain J1: the
flexible arm is worse than a rigid arm built by the same tool from the same protonation.

## What was censored, and why

Every one of the 272 selected copies returned a record. 209 were
scored and 63 failed. Dock 117 and Dock 245 were lost on the first attempt when GitHub
reclaimed their runners (exit 143 for 117; 245's log was unavailable). Both were re-run as
attempt 2, both returned a record, and both are scored here - so no copy is missing for a
reason other than the protocol's own.

| cause | copies | what it is |
|---|---|---|
| `ReceptorError` | 19 | PDB2PQR itself fails or disagrees with its own output; study A hits this class too |
| timed out after 18000 s | 17 | the 300-minute per-copy cap of amendment 3 |
| `RuntimeError` | 12 | amendment 2 declining a copy: PQR hydrogens disagree with Meeko's template within 4.0 A of the box |
| `PolymerCreationError` | 9 | Meeko cannot template a residue within 4.0 A of the box |
| `AtomValenceException` | 2 | RDKit sanitisation fails inside Meeko's PQR reader |
| `TypeError` | 2 | Meeko writes an atom serial of 100000, which overruns its PDBQT column and Vina rejects the file |
| `ValueError` | 2 | the PQR carries no chain identifier, so Meeko sees several residues sharing one key |

**17 copies hit the 300-minute cap.** That is
6.2% of the records, against about 42%
(5 of 12 readable) under the 90-minute cap that amendment 3 replaced. It is a finding about
cost rather than a footnote: flexible-receptor docking of an inositol phosphate is expensive
enough that a tail of large receptors cannot be scored within five hours each.

Two failure classes are column-overflow bugs in the toolchain rather than properties of the
structures, and both are recorded rather than worked around: Meeko's PDBQT writer emits
`ATOM  100000` for receptors with at least 100,000 atoms, which eats the column separator and
makes Vina reject the file; and for two entries the PQR reaches Meeko with no chain
identifier, so residue keys collide. Neither is the PQR column bug fixed in `cf3941d` - that
one produced `invalid literal for int()`, and it appears **nowhere** in these 272 copies.

## Receptor residues deleted under amendment 2

96 residues across 32 copies, each more than 8.0 A
beyond the docking box, which is Vina's own interaction cutoff, so none can affect a pose
inside the box. Both arms of a copy delete the same residues or the copy fails. The per-copy
lists are in `flexible.json` under `residues_dropped.by_copy`.

## Limits

* **Redocking is the easy case.** Each copy is docked into its own crystal structure, so
  these rates are upper bounds on prospective performance.
* **63 of 272 copies could not be scored**, and the causes are not random
  across structure size: the timeouts and both column-overflow classes concentrate in large
  receptors. J1 therefore describes the copies that are affordable to dock.
* **Vina's scoring has no explicit electrostatic term.** Scoring a -9 polyanion with a
  function that ignores charge is weak evidence whatever the receptor model.
* Flexible side chains per copy: mean 7.00, max
  8, with 2 copies
  having none, so the flexible arm really did differ from the rigid one almost everywhere.
