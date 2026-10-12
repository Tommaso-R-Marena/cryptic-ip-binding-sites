# Architecture homology: can the metric tell two folds apart?

**Run.** `scripts/homology_probe.py run`, seed 20261012, on the session container. Every figure below is computed by that command and written in the same pass, so nothing is transcribed.

Two fold families with ideal covalent geometry — helix barrels, whose face residues sit 3–4 apart, and beta sheets, whose face residues sit 2 apart. Each lining comes from the 3-D construction and is then cut to a spatial slab, which is what a bound ligand selects.

## What is matched, so the controls have nothing to exploit

- lining size: **18 positions, identical for every lining**
- residues per secondary-structure element: 18
- linker between elements: 10
- barrel spans: 91–200 (mean 145) — the two ranges overlap
- sheet spans: 83–199 (mean 143) — the two ranges overlap

30 linings: 15 barrel, 15 sheet

## Retrieval

| metric | P@1 | separation | permutation p | cluster purity |
|---|---|---|---|---|
| `size_only` *(control)* | 0.50 | +0.00 | 0.4918 | 0.53 |
| `span_only` *(control)* | 0.53 | -0.09 | 0.4033 | 0.50 |
| `topology` | 0.97 | +0.67 | 0.0005 | 0.67 |
| `register` | 1.00 | +3.09 | 0.0005 | 1.00 |
| `architecture` | 0.97 | +0.78 | 0.0005 | 0.67 |

Both controls are dead: `size_only` P@1 0.50 at p = 0.492; `span_only` P@1 0.53 at p = 0.403. So the separation is not size and not span.

**`register` separates the two folds cleanly** — cluster purity 1.00, separation +3.09. Note which half of the decomposition wins: for *detecting a face within one lining* the topology test is primary and register is useless, and for *comparing two linings* it is the other way round. The two halves have opposite utility on the two tasks, which is a reason to keep reporting both and not to collapse them.

The combined `architecture` metric is **worse** than its better half here, because its stated scale of 4 residues per unit of register distance leaves topology dominating about ten to one. That constant is a documented convention and has deliberately not been retuned against this result; the halves are always available separately and a reader should pick by task.

## Graded relationships, and an honest negative

| family | metric | Spearman(distance, element-count difference) |
|---|---|---|
| barrel | `span_only` *(control)* | +0.954 |
| barrel | `topology` | +0.948 |
| barrel | `register` | +0.798 |
| barrel | `architecture` | +0.946 |
| sheet | `span_only` *(control)* | +0.928 |
| sheet | `topology` | +0.861 |
| sheet | `register` | +0.478 |
| sheet | `architecture` | +0.853 |

Distance does track how many elements a fold has — but so does the span-only control, at +0.954 (barrel), +0.928 (sheet), as high as or higher than any real metric. **The graded relationship is therefore explained by span and is not a finding.** Family separation survives every control; the graded ordering does not.

## What this cannot be used for

* It is **architecture** homology, not sequence or structural homology. It does not show two proteins are related and it does not replace Foldseek or MMseqs2.
* It must never be used for the project's statistical grouping. The independence that `homology_group_strict` provides comes from sequence and structure comparison; a metric on lining spacing provides none of it.
* Two idealised folds are not a fold census. The claim is that the metric can separate two architectures that differ in face periodicity, which is the easiest possible discrimination, and it is an upper bound on real performance.
