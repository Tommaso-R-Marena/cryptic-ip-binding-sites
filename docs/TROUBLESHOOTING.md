# Troubleshooting Guide

## Analysis module (`cryptic_ip.analysis`)

### `fpocket` executable not found
**Symptom:** `FileNotFoundError: fpocket` during `detect_pockets()`.
**Fix:** Install fpocket and ensure it is on `PATH`.
```bash
conda install -c conda-forge fpocket
which fpocket
```

### APBS electrostatics failures
**Symptom:** APBS wrapper exits non-zero, missing `.dx` output.
**Fix:** Install `apbs` + `pdb2pqr`, then test manually.
```bash
conda install -c conda-forge apbs pdb2pqr
apbs --version
```

### Empty candidate set
**Symptom:** no pockets exceed threshold.
**Fix:** lower score threshold (e.g. 0.60 → 0.50) and inspect individual score components before changing defaults globally.

## Database module (`cryptic_ip.database`)

### AlphaFold download throttling / HTTP 429
Use lower request rate in batch mode:
```python
AlphaFoldBatchDownloader(output_dir="data", requests_per_second=1.5)
```

### Corrupted cache/database
Run integrity checks:
```bash
cryptic-ip validate --database results/cache.sqlite
```

## Validation module (`cryptic_ip.validation`)

### OpenMM unavailable
**Symptom:** Import error for `openmm` while running `md-validate`.
**Fix:** install OpenMM from conda-forge.
```bash
conda install -c conda-forge openmm
python -c "import openmm; print(openmm.__version__)"
```

### MD runs too slowly
Reduce system size, simulation length, and top-n candidates; run on GPU where available.

## Dependency installation issues

### `fpocket` build issues on Linux/macOS
Prefer conda binaries over source builds.

### APBS package conflicts
Create a fresh env from `environment.yml` and then install APBS stack in one command to avoid ABI mismatch.

### OpenMM + CUDA mismatch
Pin OpenMM to a version compatible with your CUDA toolkit and verify with `openmm.testInstallation`.

## Performance optimization tips

- Use `--use-ml-model` to prioritize realistic hits earlier.
- Use `screen --max-structures` for calibration runs.
- Increase parallel workers only until I/O saturation.
- Use `ParallelProcessor` checkpoints for long proteome jobs.

## Docker-specific problems

### Container cannot find mounted data
Ensure absolute bind mounts and matching in-container paths in CLI flags.

### Slow I/O in bind-mounted directories
Write temporary outputs inside container filesystem and copy final CSV artifacts out.

### Missing external executables in custom image
Start from repo Dockerfile and keep `fpocket`, `freesasa`, and `apbs` in image build steps.

## Streamlit deployment issues

### App boots locally but fails in cloud
- Ensure `requirements.txt` includes `streamlit`, `cryptic_ip` dependencies, and plotting libraries.
- Add secrets via provider dashboard (never hardcode tokens).

### Memory limits on cloud runtime
- Avoid full-proteome jobs from UI; queue heavy jobs offline.
- Limit uploaded structure size and max candidates rendered in 3D.
