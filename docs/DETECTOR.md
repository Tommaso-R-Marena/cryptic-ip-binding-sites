# The next-generation pocket detector

Code in `cryptic_ip/detect/`, driven by `scripts/detect.py`, pre-registered in
`docs/DETECTOR_PLAN.md`, and cluster-ready as study Q in
`hpc/studies/detector.yml`.

**Nothing here has been run on real candidate structures.** Everything below is
either a property of the method, measured on synthetic input, or a design
decision with its reason. The five questions that would make any of it a
finding are Q1-Q5 in the plan, and none has run.

---

## 1. What problem each piece solves

The deployed screen's own results constrain the design. Each row is a finding
from studies A-P and what it rules out.

| finding | consequence | what replaces it |
|---|---|---|
| burial is a continuum; two boundary estimators disagree (0.138 vs 0.463) and the positive class is *defined* by the 0.12 cut | any classifier against that boundary is circular | a monotone score over continuous terms, with no class boundary |
| the five-gate cascade blocks ADAT1, PDS5B, HDAC1 and HDAC3 — all annotated binders | a conjunction of five gates has the union of five false-negative rates | one score, one physical veto, and abstention for missing input |
| held-out ML ROC-AUC 0.499 against a 0.556 baseline | not a capacity problem: unlabelled pockets were used as negatives, and positives sit in 5 homology groups | positive-unlabelled learning, group-held-out, with the minimum detectable effect computed first |
| ΔAUC −0.015 for IP-specific re-weighting; 30 of 75 candidates explained by another polyanion | an absolute "is this a polyanion pocket" has the answer yes everywhere | a comparative margin against a competitor panel, with `polyanion_generic` as a real outcome |
| every candidate is a single apo conformer | a cryptic site is closed in the apo state, so burial measures the wrong variable | persistence and crypticity across an elastic-network conformer ensemble |

## 2. The modules

### `cavity.py` — detection in process, and true burial

LIGSITE's protein-solvent-protein count on a 0.8 Å grid: mark the free points,
count the directions in which protein is met on both sides, keep points
enclosed in at least 5 of 7. Nothing novel, deliberately — the novelty is what
runs on top, and a new cavity detector underneath it would be a second thing to
validate.

Two consequences: `fpocket` stops being a hard requirement (the analyze path
currently raises without it), and hundreds of conformers per protein become
affordable. **It is not fpocket and its volumes are not interchangeable with
fpocket's** — alpha spheres on a Voronoi tessellation against enclosed grid
points. Q5 measures the agreement; no threshold calibrated on one may be
applied to the other.

Verified against geometries whose answer is known by construction: a closed
shell gives exactly one cavity at the origin with 7/7 enclosure, cutting a hole
in it drops the enclosure to 6.21/7, and a flat slab gives none.

**Burial is now solvent connectivity, not a direction count.** The enclosure
count cannot separate a sealed cavity from a deep groove: a surface channel is
enclosed in five of seven directions along most of its length and reads as
buried. So each cavity also reports whether solvent can reach it from outside,
by a flood fill through free space eroded by about a water radius, and carries
`sealed_fraction`, `mouth_area` and `reaches_bulk`.

Measured on a block with a slot cut into one face, where the slot's clearance
is set by the atom lattice:

| clearance | sealed fraction | mouth area | verdict |
|---|---|---|---|
| 1.8 Å | 1.00 | 0 Å² | buried — a 1.8 Å throat is impassable to water |
| 5.8 Å | 0.30 | 82 Å² | open groove |
| 9.8 Å | 0.22 | 110 Å² | open groove |

The first row is the point: a slot that looks open on paper is sealed if water
cannot get down it, and the directional count calls all three buried.

### `ensemble.py` — crypticity, the signal that was missing

Displace along the softest elastic-network modes (both directions, ±1 and ±2 Å),
detect cavities in every conformer, match them across conformers by position,
and report per site:

* **persistence** — the fraction of conformers in which it exists;
* **crypticity** — `1 − apo_volume / best_volume`, so a site absent in the apo
  and open elsewhere scores 1.

A site with high persistence and low crypticity is an ordinary open pocket; one
with middling persistence and high crypticity is what this project set out to
find. Verified: an explicit ensemble whose apo has no cavity and which opens to
484 Å³ gives persistence 0.60, crypticity 1.00, classification `cryptic`.

**The limit this had, and no longer has alone:** elastic-network modes move
atoms rigidly with their residue, so they cannot open a site occluded by one
side chain — a large share of the real ones. The ensemble now has a second
tier.

### `rotamer.py` — the second tier

Pocket-lining side chains are rotated about chi1 and chi2 at fixed offsets and
screened for clashes. Composed with the backbone modes, this opens sites that
no amount of low-frequency backbone motion reaches.

Measured on helix barrels whose channel is plugged by tryptophans, comparing the
maximum crypticity each tier finds:

| barrel radius | backbone only | both tiers | gain |
|---|---|---|---|
| 6.5 Å | 0.130 | 0.130 | +0.000 — too tight for the TRP to rotate at all |
| 7.0 Å | 0.038 | **1.000** | **+0.962** |
| 7.5 Å | 0.046 | 0.067 | +0.021 |
| 8.0 Å | 0.028 | 0.056 | +0.028 |
| 9.0 Å | 0.012 | 0.041 | +0.029 |

At 7.0 Å the rotamer tier found a site the backbone tier missed entirely. At
6.5 Å no rotation is sterically allowed and the gain is honestly zero.

Two things here are not obvious and both were found by running it:

* **Clash screening is relative to the input, not absolute.** A predicted model
  routinely carries marginal contacts, and against a fixed 2.9 Å bar *every*
  rotation of such a side chain is a clash — including the one that opens the
  pocket. The first version produced an empty ensemble with no explanation.
* **An unphysical torsion axis is refused.** The project's existing synthetic
  structures place atoms on a lattice, where a CA–CB "bond" can be 14.7 Å.
  Rotating about that flings atoms through the protein, every conformer is
  rejected, and the output looks exactly like a tightly packed pocket rather
  than like invalid input — which is what it was.

That last point is why `geometry_fixtures.py` exists: it builds backbones by
internal coordinates with ideal bond lengths and angles (measured CA–CB 1.530 Å,
zero deviation) so geometry-dependent methods can be benchmarked on geometry
they can trust.

### `openability.py` — cost to open, not whether open

An anisotropic network model Hessian on Cα atoms, its lowest non-trivial modes,
projected onto the pocket lining's outward directions. The **opening overlap**
is how much of a mode is the opening motion; weighted by 1/λ, because a mode's
thermal amplitude goes as the inverse of its eigenvalue.

Correctness rests on the null space: a free elastic network has exactly six
zero modes. The tests check that, plus that a rigid translation and a rigid
rotation both cost zero energy, and that a disconnected network is refused
rather than analysed.

**The cost is in `gamma * A^2`, not kT.** The spring constant has never been
calibrated against a measured opening free energy. The field is named
`opening_cost_gamma_a2` for that reason — it was `opening_cost_kt` until the
mismatch was noticed, and this pipeline has already read a hydrophobic density
as a depth in angstrom once.

### `field.py` — rotation-invariant electrostatics

The pocket's potential, from formal charges under Debye-Hückel screening at
physiological ionic strength, expanded on six concentric shells in real
spherical harmonics. The per-shell, per-degree **power** is kept, and that is
exactly invariant to rotation about the pocket centre — a rotation mixes the
harmonics within a degree orthogonally and leaves the sum of their squares
alone.

This is what removes the alignment that study N needed. Rigid transplantation of
131 crystal IP6 sites gave AUC 0.599 and placements that were not buried
(median buried fraction 0.000); a rotation-invariant descriptor has no
orientation to get wrong.

Measured properties, which is how the constants were chosen:

| quantity | value | at |
|---|---|---|
| rotation discrepancy | 2.4e-2 | 192 quadrature points |
| | 2.7e-3 | **768 (default)** |
| | 6.0e-4 | 3072 |
| screening-cutoff error | 3.0e-4 | 5 Debye lengths (default) |
| cost, 78,920 pockets, 1 core | 3.0 h | 768 points |

The l=4 harmonics were hand-derived at first and were wrong. Every lower degree
was exactly orthonormal, so the obvious checks passed; the only symptom was a
rotation invariance stuck at 130% error however fine the grid. **An invariance
that does not improve with more sampling is a basis error, not a quadrature
error** — that distinction found it, and a Gram-matrix test now pins it.

### `specificity.py` — comparative, not absolute

IP6 against a panel of ATP, PAPS, UDP-glucose, sulfate, phosphate, IP3 and a
heparin tetrasaccharide, reporting one of four outcomes with a pre-registered
margin of 0.15:

* `ip6_selective` — IP6 fits better than every competitor by the margin
* `polyanion_generic` — a competitor is within the margin; it is a polyanion
  pocket and this method cannot say which polyanion
* `competitor_preferred` — a competitor fits better
* `no_polyanion_site` — the well is unfavourable, so no comparison is meaningful

The middle outcome is the one the deployed screen cannot express and is probably
the right answer for a large share of the 75 candidates. Reporting it is not a
weaker claim than a candidate list; it is the claim the evidence supports. The
summary names *which* competitor explains how many pockets, because "30 of 75
are explained by another ligand" is weaker than knowing whether they are ATP
sites or sulfate sites.

A margin sweep is published next to the result, so a reader can see how much
the answer depends on the threshold — a check nobody can run on a tuned one.

### `monotone.py` — abstention and one veto

Every feature carries a direction declared in advance, a monotone ramp, a
non-negative weight, and a reason. Monotonicity is checked over 2,000 random
perturbations rather than trusted.

Three things it does that the gate cascade cannot:

1. **A missing input is an interval, not a failure.** A crystal structure has no
   pLDDT; a pipeline without APBS has no electrostatic term. The score is
   computed over the range the missing feature could take, and **abstains** when
   that interval straddles the threshold. A protein is never rejected for a
   measurement nobody made.
2. **Alternative measurements substitute.** `relative_sasa` (measured on a
   ligand copy) and `cavity_openness` (measured by the cavity detector) are the
   same property by different means, under separate names because they are
   separate measurements. Having one drops the other from the score rather than
   widening it.
3. **One physical veto.** `field_depth_kt >= 0`: a well that repels the ligand
   cannot hold it, however buried the pocket. This was added after an acidic
   synthetic decoy — perfectly enclosed, zero basic residues, −3.6 kT/e —
   scored 0.554 against the 0.55 threshold and was called a site. Burial and
   model confidence were outvoting the physics. An unmeasured well *abstains*
   rather than clearing the veto.

It also answers the question the cascade could not: **which measurement is
holding this pocket back**, which is what an experimentalist needs.

### `pu.py` — learning without negatives, and honest power

Labelled positives plus an unlabelled pool, with the label frequency `c`
estimated per homology group and held out whole. Then
`P(positive | x) = g(x) / c`, clipped and visibly saturated rather than
exceeding one.

Elkan and Noto's mean estimator of `c` is **biased low by 30% at every true
value**, measured with the estimator isolated from any classifier by handing it
the exact `g(x)`:

| true c | mean | q90 | q95 | q99 | max |
|---|---|---|---|---|---|
| 0.2 | 0.70 | 0.96 | 0.98 | 0.99 | 1.00 |
| 0.5 | 0.69 | 0.96 | 0.98 | 0.99 | 1.00 |
| 0.9 | 0.70 | 0.96 | 0.98 | 0.99 | 1.00 |

(each entry is the estimate divided by the truth). The bias is the class
overlap, and it matters in the dangerous direction — the posterior divides by
`c`, so a `c` that is too low inflates P(positive). The default is the 95th
percentile within each group; max and q99 rest on one or two rows, and with 134
entries in 5 groups there are not enough labelled positives to spend on an
extreme quantile.

Also here: split-conformal prediction sets (which can contain *both* labels,
the method saying it does not know), a leakage audit that runs on every fit,
and the **minimum detectable AUC**:

| positive groups | min detectable AUC |
|---|---|
| 2 | 1.00 |
| 3 | 0.98 |
| **5** | **0.87** |
| 10 | 0.77 |
| 40 | 0.65 |

This project's buried-site positives occupy 5 strict homology groups, one
holding 87% of them. **Nothing below AUC 0.87 is detectable at that size**, so
the previously reported 0.499 is a statement about power, not about biology.
`python scripts/detect.py power` prints this, and Q4's decision rule requires
it to be published before any AUC.

### `sequence.py` and `signature.py` — information in the sequence, and a refusal

Pre-registered separately as study Q6 in `docs/SEQUENCE_PLAN.md`. The question:
is there a signature in the sequence at an annotated site — in its entropy, its
compressibility, its periodicity, the spacing of its coordinating residues —
that says a pocket of a particular kind forms there?

Three traps, each measured before the design was fixed, because all three
manufacture a signature that is not there.

| trap | the measurement | what the code does |
|---|---|---|
| **Kolmogorov complexity is uncomputable** | — | nothing is named after it; `lempel_ziv_phrases` is the exact LZ76 parse and `compressed_bits` is a compressor's output minus its measured empty-input header. Both bound K from above up to a constant |
| **plug-in Shannon entropy is biased, and the bias depends on N** | mean bias over 4000 random windows per row, seed 20261012, against a uniform source whose true entropy is 4.322 bits: −1.408 bits at N=10, −0.799 at N=20, −0.302 at N=50, −0.144 at N=100. A lining is 10–20 residues | Grassberger by default (−0.821, −0.316, −0.029, −0.002), window length controlled explicitly, and the plug-in value reported beside it as `entropy_bias_visible` |
| **composition is the confounder** | these sites are Arg/Lys-rich, which this project established long ago | the headline statistic is `composition_z`: the measure's z against 2000 shuffles of the same window. Same length, same letters, order destroyed — so the z is bias-free by construction and answers only whether the *ordering* carries information |

Four faults were found in the machinery by measuring it, and each one is now a
test:

* **The LZ76 parser was the wrong one.** The widely circulated Kaspar–Schuster
  snippet drops the longest-match bookkeeping, so it advances one symbol per
  phrase. It agreed with every short textbook case — `aaaaaaaa` to 2,
  `abababab` to 3, `abcdefgh` to 8 — while scoring a random 500-residue window
  at **485 phrases against the correct 213**. The textbook cases cannot see
  that fault, so the parse is returned and tested against its own definition.
* **Periodicity exceeded 1 at the Nyquist period.** Projecting onto two
  quadratures and normalising by `variance · n/2` returned 2.0 for a pure
  alternating signal, because the sine vanishes at period 2 and only one
  quadrature contributes. Replaced by a least-squares explained variance, which
  normalises by the basis it actually fitted.
* **A single exact probe has a blind spot.** A strictly period-4 charge train is
  *exactly* orthogonal to a 3.6-residue probe over 36 points, scoring 0.000 for
  a sequence that plainly has a periodic basic face. So `HELIX_BAND` is a band,
  3.0 to 4.2 — the lower edge is 3.0 and not 3.2 because a 3₁₀ helix is exactly
  3.0 residues per turn and the first version scored a perfect period-3 pattern
  at 0.099.
* **Band power is confounded by harmonics.** A period-7 pattern has content at
  3.5 and scored 0.33 in the helical band. The fundamental still wins the scan
  (7.1, not 3.5), so `dominant_period` is reported too and `helix_face`
  requires the dominant period to be *in* the band, not merely power there.

`lz_normalised` carries its own warning: `n / log_a(n)` is asymptotic, and on
random windows the ratio is 0.88 at N=20, 0.84 at N=36, 0.83 at N=100 and 0.88
at N=500 — below 1 and **not monotone in N**. Like the entropy bias, that makes
the raw value unsafe to compare across window lengths, and the composition null
is the answer to both.

`signature.py` decides *which* string, which is where the remaining mistakes
live:

* A pocket's lining gives two very different windows — the lining residues with
  the gaps closed, and the contiguous stretch they span — and **only the second
  can carry a period**. Closing the gaps destroys exactly the positional
  information a periodicity measures. Every function in `sequence.py` accepts
  either, because both are strings, so `statistic_for_window` refuses a
  positional statistic on the lining window.
* Lining residue numbers in this repo are **chain-blind** by their own
  docstring, so they arrive as `(chain, number)` pairs. A pocket with less than
  80% of its lining on one chain is an interface pocket, has no single sequence
  window, and is excluded with that reason rather than measured.
* Exclusions are reported by cause, and the non-standard-residue fraction is
  tested before the window-length floor, because dropped residues are what
  *made* a window short.

And the refusal: `ESTABLISHED_SIGNATURES` is empty, so `sequence_contribution`
returns `None` — an abstention, which `monotone.py` treats as "unmeasured"
rather than as a zero that would push a candidate down. The descriptors are
computed, recorded and available, and **they move no score until study Q6
reports**. A name may be added only alongside a results file that cleared all
five of the plan's decision rules. The project has a held-out AUC of 0.499 on
record from the last time an untested feature was trusted.

The method validation is committed at `results/sequence/VALIDATION.md`: a
planted helical-face signature fires at mean z **+4.05** in 100% of windows,
the same letters shuffled give **−0.24** at 0%, and a basic-rich pool with no
ordering gives **−0.13** at 2% — with all three pools holding an identical
basic fraction of 0.270 and entropies within 0.1 bits, so the separation cannot
be composition.

## 3. Running it

```bash
python scripts/detect.py spec                      # the specification, auditable
python scripts/detect.py power                     # what this dataset could detect
python scripts/detect.py describe model.pdb        # one structure, full detector
python scripts/detect.py describe model.pdb --no-ensemble   # skip most of the cost
python scripts/detect.py shard --structures DIR --shard 0 --shards 272 --out part.json

python scripts/sequence_signature.py validate                  # study Q6 method validation
python scripts/sequence_signature.py describe --inputs in.jsonl --out sites.json
python scripts/sequence_signature.py report --sites sites.json --out sequence.json
```

On the cluster, study Q in `hpc/studies/detector.yml`. It is **not** plan-gated,
because `docs/DETECTOR_PLAN.md` exists — unlike study M and study P, which
`rockfish.py submit` still refuses.

Measured cost: the full detector with a 13-conformer ensemble took **4.0 s** on
a 1,202-atom synthetic structure at 384 quadrature points. The candidates are
much larger (median 436 residues, maximum 2,352) and the cavity grid grows with
the bounding box, so the real per-structure cost is expected to be one to two
orders of magnitude higher on the long tail. That is an expectation, not a
measurement, which is why study Q's first stage is a timing probe on real
structures rather than a production sweep — the same discipline that let a
90-minute CPU probe kill a co-folding study before it wasted far more.

## 4. Benchmarks

`scripts/detect_benchmark.py run` executes seven suites; the committed result is
`results/detector/BENCHMARK.md` with its provenance alongside. 48 measurements,
0 failing. Pre-registered as a method validation in
`docs/DETECTOR_PLAN_AMENDMENT_1.md`, which also states what it may not be used
for.

| quantity | measured |
|---|---|
| field spectrum rotation discrepancy | 2.0e-03 |
| field spectrum translation discrepancy | 6.1e-16 |
| field discrepancy at 0.1 Å coordinate noise | 0.014 |
| field discrepancy at 0.5 Å coordinate noise | 0.138 |
| cavity volume change at 0.5 Å noise | 0.215 |
| basic minus acidic well depth | 3.43 kT/e |
| cavity detection scaling exponent | 1.02 |
| **field panel scaling exponent** | **0.032** |
| conformal coverage at α = 0.1 | 0.903 |
| minimum detectable AUC at 5 homology groups | 0.87 |
| monotonicity violations in 2000 perturbations | 0 |
| deployed scorer separation margin | 0.052 |
| monotone score separation margin | 0.095 |

The field panel's scaling exponent of 0.032 is the one worth a sentence: **the
descriptor's cost does not grow with the protein**, measured over barrels from
360 to 1,260 atoms whose radius grows with their helix count so the structure
really does get bigger. That is what the 5-Debye-length screening cutoff was
for, and it means the 2,352-residue tail of the candidate set costs what a
small protein costs.

**It says nothing about whether the detector finds real binding sites.** The
fixtures are idealised barrels with chosen sequences and no energy
minimisation. Q1–Q5 ask the biological questions of real structures and none
has run.

The benchmark also found two faults in itself, both recorded in
`results/detector/PROVENANCE.md`: its ground-truth labels were wrong (two
poly-alanine channels marked positive for being buried, when burial is not
site-ness — and both scorers duly reported a negative separation margin), and
its cost suite was timing its own warm-up.

## 5. What it cannot do

* No output is calibrated against a measured affinity. None is a prediction.
* The opening cost is ordinal at best and is not in kT.
* The electrostatic field is analytic, in uniform bulk-water permittivity. A
  buried pocket's lower local permittivity is not modelled, so well depths are
  *underestimated* for the most buried sites — the direction that loses
  candidates rather than inventing them.
* The ligand templates are spherically averaged. An inositol phosphate is not.
* Elastic-network modes cannot open a side-chain cryptic site.
* One chemistry dominates the ground truth: 134 InsP6 against 1 InsP5. IP3 is
  in the competitor panel, not the target, for that reason.
* No wet-lab validation. The output is a set of hypotheses.
