# Exploration plan: zero-parameter physics descriptors

Pre-registered before any value below was computed. Commit this file, then run
`scripts/explore_descriptors.py`. Every run appends to a ledger; nothing is
reported that is not in the ledger, and nothing in the ledger is dropped.

## Why this exists

The confirmatory benchmark (`docs/ANALYSIS_PLAN.md`) asks one narrow question:
does hull depth add to a learned classifier? This plan asks a broader one:
does **any single physically motivated descriptor, with no fitted parameters**,
separate inositol-phosphate sites from other pockets beyond chance, once
homologous proteins are not counted as independent evidence?

A descriptor with no fitted parameters cannot overfit. The remaining risk is
the garden of forking paths: testing many descriptors and reporting the best.
The rules below control that risk. They are:

1. every descriptor's direction is predicted from physics before looking;
2. every descriptor is tested and reported, with false-discovery control across all of them;
3. selection happens on the development partition only;
4. confirmation happens once, on the temporal holdout, for a set fixed by rule
   before the holdout is read.

## Data

The table is the benchmark's `benchmark-table` artifact, built by
`scripts/benchmark.py prepare`. It holds one row per pocket, with descriptors,
labels, homology groups, a holdout flag and the rule-based score.

- **Development** means rows with `holdout == False`. The exploration stage
  refuses to run if it is given holdout rows, and asserts that none are present.
- **Holdout** means rows with `holdout == True`. Only the confirmation stage
  reads them, and only for the descriptors named in the confirmatory set.

Missing values (for example, basic-nitrogen distances when a pocket has none)
rank as the least binding-like value.

## Tasks

| task | positives | role |
|---|---|---|
| `ip_site` | pockets that overlap a crystal IP ligand | **primary**; inferential |
| `cryptic_ip_site` | IP sites that are buried | descriptive only |

Prepare found that the positives of `cryptic_ip_site` come from about six
independent families, and one family holds 87% of them. No interval computed
on that data can be trusted, so it is reported per group with no p-values and
no claims.

## Descriptor directions (a priori)

The score for each descriptor is `sign × value`, where the sign comes from this
table. A descriptor "confirms its prior" only if its development AUROC is above
0.5 in the declared direction.

**Higher predicts binding (+1).** IP ligands are polyanions of up to −12, held
by arginine, lysine and histidine and by hydroxyls, in enclosed, polar pockets.
The descriptors are:

- **Basic residues:** `n_basic_residues`, `n_strong_basic_residues`, `basic_fraction`,
  `n_basic_residues_core`, `n_basic_nitrogens`.
- **Positive charge and potential:** `net_formal_charge`, `positive_charge_density`,
  `charge_balance`, `coulomb_potential_kt`.
- **Burial:** `enclosure`, `burial_depth`, `buried_residue_fraction`, `hull_depth`.
- **Polar and stacking contacts:** `n_hydroxyl_residues`, `hydroxyl_fraction`,
  `polar_fraction`, `n_aromatic_residues`, `aromatic_fraction`.
- **Cavity packing:** `alpha_sphere_density`.

**Lower predicts binding (−1).** Acidic residues repel a polyanion, and
coordinating nitrogens must be close and clustered. Buried sites have low
solvent exposure, and IP sites are polar. The descriptors are:

- **Acidic residues:** `n_acidic_residues`, `acidic_fraction`, `n_acidic_residues_core`.
- **Basic-nitrogen geometry:** `basic_nitrogen_min_distance`,
  `basic_nitrogen_mean_distance`, `basic_nitrogen_dispersion`.
- **Solvent exposure:** `mean_relative_sasa`, `sasa_mean`, `sasa_median`, `sasa_min`, `sasa_max`.
- **Hydropathy:** `hydropathy_mean`.

**No prior (two-sided; 0).** Size and shape have no clear expected
direction: IP6 needs roughly 300–800 Å³, but fpocket over-segments large
cavities and under-segments small ones. These are scored with a sign of +1 and
reported, but they can never "confirm a prior". The descriptors are:
`pocket_volume`, `hull_volume`, `n_alpha_spheres`, `radius_of_gyration`,
`asphericity`, `max_extent`, `sasa_total`, `n_residues`.

### Pre-declared composites

These are fixed now and have no weights:

- `electropositive_enclosure`: the within-structure percentile rank of
  `coulomb_potential_kt`, plus that of `enclosure`.
- `basic_cluster`: the within-structure percentile rank of `n_basic_nitrogens`,
  plus that of `−basic_nitrogen_dispersion`.
- `rule_score`: the repository's existing hand-built score, reported as a reference.

Earlier work in this repository already evaluated `rule_score`, `burial_depth`
and `hull_depth` on these kinds of data. Their priors are therefore **not blind**.
They are reported, but excluded from the confirmatory set.

## Metrics

Two metrics are computed for each descriptor and composite on development rows.

1. **Pooled AUROC** over pockets, with a 95% interval and a two-sided p-value
   against 0.5 from a group bootstrap over `group_strict`. There are 2,000
   resamples, and the seed is 20260924.
2. **Within-structure recovery.** For each structure with at least one
   positive pocket, check whether the top-ranked pocket is a positive (top-1),
   and whether any of the top 3 is (top-3). A structure's chance rate is the
   expected hit rate under a random ordering of its pockets. The statistic is
   observed minus chance, averaged over structures, with each strict group
   weighted equally. It has a group-bootstrap interval and p-value against 0.

Top-1 recovery is the practical question: given a protein, can the descriptor
point at the site? Pooled AUROC mixes proteins with different pocket counts.

## Multiplicity and decisions

- The family of tests on development is 39 descriptors plus 2 composites, for
  both metrics: 82 tests. Benjamini–Hochberg is applied across all 82 at q = 0.05.
- **Confirmatory set.** Take the descriptors and composites that meet all four conditions:
  - their pooled-AUROC test passes BH;
  - their direction matches the prior (sign ≠ 0);
  - they are not in the non-blind list;
  - their top-1 recovery excess has a positive point estimate.

  Rank them by the lower bound of the pooled-AUROC interval and take at most 3.
  The set is written to the ledger and committed before confirmation runs.
- **Confirmation** runs on holdout rows, for the confirmatory set only, and runs once. A
  descriptor is **confirmed** if its holdout pooled-AUROC interval lies above
  0.5 and its Holm-adjusted p-value (across the set) is below 0.05. The holdout
  has positives in only about 8 strict groups, so a non-confirmation is
  reported as "not confirmed (underpowered)", not as a refutation.
- A discriminating descriptor is not a discovery. `basic_fraction` separating
  IP sites from random pockets is the expected textbook result, and it is reported as such.
  A result counts as **new** only if it beats `rule_score` in a paired
  comparison on development (group-bootstrap interval of the AUROC difference
  above 0), and that comparison is itself one more ledger entry.

## Ledger

`scripts/explore_descriptors.py` appends one JSON line per stage to
`results/exploration/ledger.jsonl`. Each line records:

- the table's SHA-256;
- the git commit;
- the stage;
- every statistic computed.

The ledger is committed as it grows. Any additional analysis requires an
amendment to this file that is committed first, and it adds to the multiplicity family.
