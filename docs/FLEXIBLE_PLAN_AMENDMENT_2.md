# Study J, amendment 2: what to do with receptor residues Meeko refuses

Dated 2026-10-01, **before any flexible-receptor pose exists**. No J1, J2 or J3 figure has
been computed on real structures; the only docking run so far scored two copies and was
cancelled. This amendment adds a rule for receptors Meeko will not build and says what is
recorded when it applies. J1's estimand, the 2 Å success criterion, the ±0.05 margin, seed
20261003, `homology_group_strict` resampling and the 5-group floor are unchanged, as are J2
and J3 and amendment 1.

## What forced it

Run 36828409167 docked 2 of the 10 copies in shard 1 and failed the other 8. Three of those
eight are Meeko refusing the receptor outright, in two shapes:

- `Template matching failed for: ['B:351', 'B:423']` — Meeko has no chemical template for
  those residues.
- `Residue A:174 matched with template 'None' has H discrepancy: 3 missing, 0 excess` —
  the hydrogens PDB2PQR placed disagree with the template Meeko matched, which is Meeko
  saying the PQR's charges are not applicable to the receptor it would build.

Both are about residues **elsewhere in the protein**, not about the pocket. Studies A, F and
G never saw either, because the project's own receptor path types and charges the receptor
itself; amendment 1's move to Meeko is what exposes them. Left alone they would cost the
study roughly a third of its copies, and non-randomly: the structures with unusual residues
are not a random third.

## The rule

A residue Meeko cannot type, or whose hydrogens disagree with its template, is **deleted
from the receptor** if and only if every one of its atoms lies more than **8.0 Å beyond the
docking box**, measured from the nearest box face. Otherwise the copy is **not docked** and
is recorded as a failure with the residue named.

8.0 Å is AutoDock Vina's own interaction cutoff, not a tuned number: a receptor atom more
than 8 Å outside the box cannot contribute to the score of any pose inside it, so deleting
it cannot change a pose or an RMSD. The cutoff was fixed before any copy was scored and may
not be adjusted on results.

Two further requirements, because this touches the thing being compared:

1. **Both arms see the same receptor.** The rigid and flexible receptors are built by the
   same code from the same PQR, and a copy whose two arms delete different residues is
   **failed**, not reported. Without that, a deletion could masquerade as side-chain
   freedom.
2. **Every deletion is published.** The report states how many residues were deleted, in
   how many copies, with the per-copy lists, next to the drop radius. A reader who distrusts
   the rule can subtract those copies.

## What this does not fix

The remaining failure classes stay failures and are reported as such:

- **PDB2PQR itself failing** on a structure. Study A hits this too (19 of 272 copies) and
  it is nothing to do with this study.
- **An untypable residue near the pocket.** By the rule above the copy is dropped, which is
  the honest outcome: the receptor Meeko would build there is not the receptor the plan
  asked for.

## A correction, not an amendment

Four of shard 1's eight failures were a bug of mine, now fixed, and are recorded here so the
run history reads correctly rather than because the protocol changed. PDB2PQR writes the PQR
on PDB columns, where a four-digit residue number fills its field and runs into the chain —
`GLY A2401`. Meeko reads the PQR by splitting on whitespace, so it read the chain as
`A2401` and then the x coordinate as the residue number, raising `ValueError: invalid
literal for int() with base 10`. The file is now re-spaced from its columns before Meeko
sees it, which is what the plan always intended: the same atoms, the same charges, read
correctly. Exactly the four copies whose residue numbers reach four digits failed this way,
and the failure and the fix are both reproduced in the test suite against real PDB2PQR and
real Meeko.

Had that `int()` happened to succeed, the damage would have been worse than a crash: the
chain would have been `A2401`, every monomer key wrong, no flexible residue found, and both
arms would have docked a rigid receptor under a `flex` label — a null result that looked
like a measurement. That is the failure mode this project is least willing to accept, and
it is the reason the fix is pinned by tests rather than merely applied.
