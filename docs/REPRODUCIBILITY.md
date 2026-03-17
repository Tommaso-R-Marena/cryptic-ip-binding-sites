# Reproducibility Guide

This document defines the end-to-end reproducibility workflow for all cryptic IP-binding site analyses.

## 1) Configuration management

### Files
- Default run configuration: `config/defaults/pipeline.yaml`
- Schema for validation: `config/schemas/pipeline.schema.yaml`

### Validate config before every run
```bash
python scripts/reproducibility_workflow.py validate-config \
  --config config/defaults/pipeline.yaml \
  --schema config/schemas/pipeline.schema.yaml
```

### Default parameters and when to change them
- `pipeline.score_threshold: 0.60` — increase (e.g., 0.7) to reduce false positives, decrease (e.g., 0.5) for exploratory sweeps.
- `pipeline.max_structures: null` — set to an integer for pilot runs or CI smoke tests.
- `pipeline.seed: 42` — only change when intentionally running sensitivity analyses.
- `pipeline.sort_columns` — keep as `composite_score, uniprot_id, pocket_id` for deterministic ranking.
- `pipeline.float_tolerance: 1e-6` — use for expected-value comparisons across hardware.
- `pipeline.output_precision_decimals: 6` — increase to 8+ for near-threshold score studies.
- `data_sources.alphafold_release_date` and `data_sources.pdb_fetch_date` — update at each data refresh.

### Version controlling config
Commit all config edits with analysis code. Suggested commit message format:
```text
repro(config): update score threshold to 0.65 for comparative run
```

## 2) Provenance tracking

Generate a JSON-LD manifest for every analysis:
```bash
python scripts/reproducibility_workflow.py manifest \
  --config config/defaults/pipeline.yaml \
  --inputs data/structures/AF-Q9Y2R2-F1-model_v4.pdb \
  --outputs results/screening_results.csv \
  --out results/provenance_manifest.jsonld
```

Each manifest captures:
- Input/output checksums (SHA256)
- Exact parameter values
- Software/runtime versions
- Pipeline git commit hash (`pipelineVersion`)
- Data source versions (`alphafold_release_date`, `pdb_fetch_date`)

## 3) Result reproducibility requirements

### Randomness control
CLI entry points now support deterministic seeds:
- `cryptic-ip analyze --seed 42 ...`
- `cryptic-ip screen --seed 42 ...`

### Deterministic ordering
Proteome screening results are stably sorted by:
1. `composite_score` (descending)
2. `uniprot_id` (ascending)
3. `pocket_id` (ascending)

### Floating-point notes
Small floating-point differences can occur across CPUs/BLAS libraries. Compare score outputs using tolerance (`abs(a-b) <= 1e-6`) rather than strict equality.

### Output checksums
Checksums are generated into:
- `results/repro_bundle/checksums.json` (bundle export)
- `results/provenance_manifest.jsonld` (per-run provenance)

## 4) Environment pinning

### Pinned dependencies
`requirements.txt` is fully pinned with `==` versions.

### Runtime verification
```bash
python scripts/reproducibility_workflow.py check-versions --requirements requirements.txt
```
Fail the run if mismatches are reported.

### Lockfile
`environment.lock.yml` captures the tested conda-style environment.

### Tested software matrix
- Python version and OS are written into the provenance manifest under `softwareVersions`.

## 5) Archival and sharing

### Generate methods section text automatically
```bash
python scripts/reproducibility_workflow.py methods \
  --config config/defaults/pipeline.yaml \
  --manifest results/provenance_manifest.jsonld \
  --out results/METHODS_AUTO.md
```

### Export complete analysis bundle
```bash
python scripts/reproducibility_workflow.py bundle \
  --config config/defaults/pipeline.yaml \
  --manifest results/provenance_manifest.jsonld \
  --methods results/METHODS_AUTO.md
```

Bundle contents include code/config snapshots, methods text, manifest, and checksums.

### Zenodo DOI workflow
1. Create a GitHub release with the corresponding commit/tag.
2. Connect repository to Zenodo and enable release archiving.
3. Upload `results/repro_bundle/` as release asset (if needed).
4. Record DOI in `CITATION.cff` and manuscript metadata.

### Research compendium checklist
- `README.md` + `docs/REPRODUCIBILITY.md`
- Raw data acquisition scripts in `scripts/`
- Config + schema in `config/`
- Provenance manifests in `results/`
- Validation scripts and expected baselines in `config/defaults/expected_validation.json`

## 6) Validation framework

### Reproduce published baseline from raw workflow
```bash
python scripts/reproducibility_workflow.py reproduce \
  --expected-json config/defaults/expected_validation.json
```

### Continuous validation against known structures
Run this command in CI after dependency install:
```bash
python scripts/reproducibility_workflow.py reproduce
```

If observed values diverge from expected baseline, the command exits with non-zero status and should alert CI maintainers.

## Recommended full reproducibility runbook

```bash
# 1) Validate config
python scripts/reproducibility_workflow.py validate-config

# 2) Run analysis
cryptic-ip screen data/structures --output results/screening_results.csv --seed 42

# 3) Create provenance
python scripts/reproducibility_workflow.py manifest \
  --inputs data/structures/AF-Q9Y2R2-F1-model_v4.pdb \
  --outputs results/screening_results.csv

# 4) Check environment
python scripts/reproducibility_workflow.py check-versions

# 5) Generate methods and export bundle
python scripts/reproducibility_workflow.py methods
python scripts/reproducibility_workflow.py bundle

# 6) Reproduce baseline validation
python scripts/reproducibility_workflow.py reproduce
```
