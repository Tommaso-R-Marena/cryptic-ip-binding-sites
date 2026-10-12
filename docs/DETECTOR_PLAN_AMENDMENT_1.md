# Detector plan, amendment 1: a method-validation benchmark, and two new tiers

Dated 2026-10-12, **before any benchmark figure was read**. The suites were
written and the fixtures built; this document was committed before the numbers
in `results/detector/BENCHMARK.md` were looked at. Q1-Q5 are unchanged, their
decision rules are unchanged, and none of them has run.

## What is added, and why it is not a study

`docs/DETECTOR_PLAN.md` asks five biological questions of real structures. This
amendment adds something different in kind: a **method-validation benchmark**
on synthetic fixtures, covering determinism, invariance, robustness to
coordinate noise, discrimination on constructed cases, cost and its scaling,
statistical recovery against known parameters, and a head-to-head against the
deployed rule-based scorer.

It is not a study and it answers none of Q1-Q5. The distinction matters enough
to put in writing because a benchmark that blurred it would be the most
dangerous artefact in this repository: a table of numbers that looks exactly
like validation while saying nothing about biology. So every suite declares
what its inputs can support, and the report repeats that it says nothing about
whether the detector finds real binding sites.

## Why the existing synthetic structures could not be used

`cryptic_ip.testing.synthetic` places atoms on a lattice and labels them by
priority. Its own header says NOT A REAL STRUCTURE, and it is right: in the
files it writes a CA-CB "bond" can be 3.3, 11.3 or 14.7 angstrom long.

That is fine for what it was built for, which needs only positions. It is not
fine for anything depending on covalent geometry. Rotating a side chain about a
14.7 angstrom "bond" flings its atoms through the protein; every conformer is
then rejected as a clash, and the output looks exactly like a tightly packed
pocket rather than like invalid input. The first attempt at benchmarking the
rotamer tier on those structures produced contacts of 0.17 angstrom - two atoms
superimposed - and reported zero accepted conformers with no explanation.

Two things follow, and both are in this amendment rather than in a commit
message:

1. `cryptic_ip/detect/geometry_fixtures.py` builds backbones by internal
   coordinates with ideal bond lengths and angles, and attaches side chains
   with real chi torsions. Measured CA-CB length 1.530 angstrom with zero
   deviation across residues.
2. `cryptic_ip/detect/rotamer.py` refuses a torsion axis that is not a
   covalent bond, naming the offending residue and the axis length, rather than
   silently producing nonsense.

## The two new tiers, pre-registered before their effect was measured

**Bulk-solvent burial.** The cavity detector counted directional enclosure,
which cannot separate a sealed cavity from a deep groove: a surface channel is
enclosed in five of seven directions along most of its length and reads as
buried. Each cavity now also reports whether solvent can reach it from
outside, by a flood fill through free space eroded by about a water radius, and
publishes `sealed_fraction`, `mouth_area` and `reaches_bulk`. The 0.5 boundary
on `sealed_fraction` that makes `is_buried` a label is stated, not fitted, and
callers wanting the continuum read the fraction.

**Side-chain rotamers.** `docs/DETECTOR.md` named the sharpest limit of the
conformer ensemble: elastic-network modes move atoms rigidly with their residue
and so cannot open a site occluded by one side chain, which is a large share of
real cryptic sites. The ensemble now has a second tier that rotates
pocket-lining side chains about chi1 and chi2 at fixed offsets, screened for
clashes **relative to the input conformer** rather than against an absolute
cutoff - a predicted model routinely carries marginal contacts, and against a
fixed bar every rotation of such a side chain is a clash including the one that
opens the pocket.

Neither tier subsumes the other and the report keeps them apart.

## What the benchmark may not be used for

* It may not be cited as evidence that the detector finds inositol-phosphate
  sites. The fixtures are idealised barrels with chosen sequences and no
  energy minimisation.
* Its ground-truth labels are claims about construction, not about biology,
  and `GROUND_TRUTH_REASON` publishes the reason for each so they can be
  argued with. The first set of labels was wrong - two poly-alanine channels
  were labelled positive because they are buried, when burial is not
  site-ness - and the benchmark is what caught it. Burial is tested on its own
  in B4, where it belongs.
* Its cost figures are from barrels, not candidates. The candidates run to
  2,352 residues and the real per-structure cost is expected to be
  substantially higher; study Q's first stage is a timing probe on real
  structures for that reason.
* A threshold or constant may not be tuned on it. Every parameter in
  `docs/DETECTOR_PLAN.md` stays where that document fixed it.

## Decision rule for the benchmark itself

The benchmark is a regression gate, not a hypothesis test. A failing
measurement means the code changed behaviour, and the response is to fix the
code or to amend this document saying why the new behaviour is correct - never
to loosen the bound quietly. The suites that assert nothing (robustness, cost)
report their curves without a pass or fail, because there is no principled
threshold for them and inventing one would be the same mistake as a tuned
cutoff.
