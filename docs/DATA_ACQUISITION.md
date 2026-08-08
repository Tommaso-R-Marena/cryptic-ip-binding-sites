# Data acquisition

How the ground-truth dataset is collected, what it contains, and what to do when
the structural databases are unreachable.

## What gets collected

| Stage | Source | Output |
|---|---|---|
| Ligand vocabulary | RCSB chemical component dictionary | Every validated inositol phosphate component |
| Ligand-bearing entries | RCSB search API | Every PDB entry containing one, all experimental methods |
| Decoy entries | RCSB search API | High-resolution entries containing **none** |
| Coordinates | `files.rcsb.org` | mmCIF (preferred) or PDB, SHA-256 verified |
| Annotations | `data.rcsb.org` GraphQL | UniProt, organism, resolution, method, R-free, release date |
| Measurements | local | Per-ligand-copy burial, contacts, coordination |

## Ligand vocabulary is discovered, not hard-coded

Earlier versions hard-coded eight component identifiers. That approach fails in
two ways that are invisible at run time: the PDB contains many more inositol
phosphate species than any hand-written list (regioisomers, inositol
pyrophosphates, non-*myo* stereoisomers), and a mistyped identifier contributes
nothing without raising an error.

`cryptic_ip.database.ip_ligands` instead:

1. Full-text searches the chemical component dictionary for inositol species.
2. Resolves **every** candidate — discovered or seeded — against the component
   API, so unknown identifiers are dropped with a recorded reason.
3. Derives the series (InsP1 … InsP8) by counting phosphorus atoms in the
   *reported formula*, never by parsing the identifier.
4. Excludes lipid-linked species (phosphatidylinositol phosphates) by default,
   since they bind at membranes rather than forming buried cofactor sites.

The seed list only widens the search; nothing in it is trusted. Rejected
candidates and their reasons are written to the manifest.

## All experimental methods, not just X-ray

The previous query filtered on `exptl.method == "X-RAY DIFFRACTION"`. Inositol
phosphates act as structural cofactors most often in large assemblies, which are
increasingly solved by cryo-EM, so that filter discarded exactly the structures
the project is about. The default now accepts every method; `--max-resolution`
and `--experimental-methods` remain available when a restriction is wanted.

## Pagination and batching

- The search API caps a response at 10 000 rows. The client reads `total_count`
  and pages until the set is exhausted, so large queries are not silently
  truncated.
- Entry metadata is fetched through batched GraphQL, 100 entries per request.
  Per-entry REST calls would cost thousands of requests for a proteome-scale
  query.
- mmCIF is preferred over legacy PDB because entries with more than 62 chains or
  99 999 atoms have **no** PDB-format file at all.

## Provenance and resumability

Every API response is cached on disk under a hash of the request, so re-running
is free and interrupted runs resume. The manifest records:

- the ligand registry and the queries that produced it, including rejections;
- per-file SHA-256 checksums, byte counts and retrieval timestamps;
- API request, retry and cache-hit counts;
- per-structure measurement failures;
- software versions and the exact command line.

## Usage

```bash
# Full collection with protein-level negatives
python scripts/build_ip_validation_dataset.py \
  --n-decoys 500 \
  --jobs 8 \
  --download-workers 4

# Fast smoke run
python scripts/build_ip_validation_dataset.py --max-entries 20 --jobs 2
```

Outputs:

| Path | Contents |
|---|---|
| `data/validation/ip_ligand_instances.csv` | One row per ligand copy |
| `data/validation/ip_binding_validation_dataset.csv` | One row per entry, most buried copy |
| `data/validation/dataset_manifest.json` | Full provenance |
| `data/validation/raw/structures/` | Coordinate files |
| `data/validation/raw/api_cache/` | Cached API responses |

## Network requirements

The build needs outbound HTTPS to:

- `search.rcsb.org`
- `data.rcsb.org`
- `files.rcsb.org`

If they are unreachable the script exits with status 2 and an explanatory
message rather than writing an empty dataset. Restricted environments — CI
sandboxes, air-gapped clusters, and network policies that allowlist only package
registries — are common, so:

- **Tests skip** rather than fail when the databases are unreachable. Set
  `CRYPTIC_IP_REQUIRE_NETWORK=1` to turn that back into a failure for a release
  check.
- **`cryptic-ip self-check`** verifies the whole measurement path offline against
  synthetic structures with known ground truth.
- To move a dataset between environments, copy `data/validation/` across; the
  manifest checksums let the recipient verify every file.

## Feature extraction

```bash
python scripts/extract_pocket_features.py --jobs 8
```

Per-structure results are cached as JSON, so this is resumable too. The run
prints pocket-detector recall on known ligand sites; a value below 80 % is
flagged, because every missed site is unreachable by the classifier no matter how
good it is.

See also: [METHODS.md](METHODS.md), [IP_VALIDATION_DATASET.md](IP_VALIDATION_DATASET.md).
