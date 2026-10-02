# Study J, amendment 1: both arms are docked fresh, with one receptor pipeline

Dated 2026-10-01, **before any flexible-receptor pose exists**. It changes how the arms are
produced, not what J1 estimates. J1's estimand, 2 Å success criterion, ±0.05 margin, seed
20261003, `homology_group_strict` resampling and 5-group floor are all unchanged, as are J2
and J3.

## What the plan assumed

`docs/FLEXIBLE_PLAN.md` says the rigid arm is **reused, not re-docked**, taken from the
`rerank-arms-*` artifacts of rerank run 36055694750, "and this is what makes the study
affordable". It also says the receptor is prepared "exactly as in `docs/REDOCKING_PLAN.md`,
through the same code".

Those two requirements turn out to conflict, and the conflict only becomes visible once the
flexible arm is actually built.

## Why they conflict

A flexible-receptor Vina run needs a **second** PDBQT holding each movable side chain as a
torsion tree (`ROOT`, `BRANCH`, `ENDBRANCH`), with the same atoms removed from the rigid
file. The project's own receptor path (`cryptic_ip.docking.receptor.write_pdbqt`) writes a
rigid receptor only; it has no torsion-tree writer.

There are two ways to get one, and both break something:

1. **Hand-write the torsion trees** from the project's typed atoms. This keeps the rigid
   receptor byte-identical to studies A/F/G, so the rigid arm could be reused. But it means
   hand-encoding side-chain rotatable-bond topology for every residue type, including which
   ring systems are internally rigid. A subtle error there produces poses that look
   plausible and are wrong, which is the failure mode this project is least willing to
   accept.
2. **Let Meeko build them** (`Polymer.flexibilize_sidechain`, verified to emit a correct
   `ROOT`/`BRANCH` tree). The trees are then right by construction, but Meeko also retypes
   and recharges the rigid part, so the rigid receptor is **not** atom-identical to the one
   studies A/F/G docked into.

Under option 2, reusing study G's rigid arm would compare a Meeko-prepared flexible receptor
against a project-prepared rigid one. Any difference in top-pose success would then mix the
thing being tested (side-chain freedom) with a thing nobody is asking about (which tool
assigned the atom types and charges). That is a confound, not a saving.

## The decision

**Option 2, with both arms re-docked through the same Meeko receptor.** Correct torsion trees
matter more than the saving, and an unconfounded comparison matters more than either.

- The flexible arm docks into a Meeko rigid + flex pair.
- The rigid arm docks into **the same Meeko rigid file with no flexible residues**, which the
  writer emits as an empty flex string — so the two arms differ in exactly one thing.
- Protonation still comes from the project's own PDB2PQR path at pH 7.4, so the two
  pipelines share their protonation input and differ only in typing and torsion handling.
- Study G's rigid arm is **not** reused, and the studies A/F/G numbers are untouched and
  remain the reference for the project-prepared receptor.

**Cost.** 250 copies × 3 seeds × 2 arms = **1,500 docking runs** rather than 750. That is the
price of removing the confound and is accepted here.

**Reported alongside J1**, so the change is auditable rather than assumed harmless: the rigid
arm's top-pose success under the Meeko receptor, next to study F's 0.110 under the
project-prepared receptor. If those differ materially, the receptor pipeline matters on its
own and that is a finding about the pipeline, reported as such and not folded into J1.

A pinning test still holds the two paths together where the plan asked: with no flexible
residues selected, `scripts/flexible.py` must produce an empty flex file and take the same
single-receptor code path, so "rigid" means rigid.
