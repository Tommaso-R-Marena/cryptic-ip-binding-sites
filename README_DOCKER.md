# Docker Guide for `cryptic-ip-binding-sites`

This guide provides a reproducible, production-ready container workflow for the cryptic IP-binding pipeline, including:

- Python 3.10 runtime (compatible with Python 3.8+ requirement)
- fpocket, FreeSASA, APBS, PDB2PQR, OpenMM installed from conda-forge
- CLI-first container entrypoint (`cryptic-ip`)
- Optional Jupyter Lab service for interactive analysis

## Quick Start

### 1) Build image

```bash
./scripts/docker/build.sh
```

### 2) Run CLI command

```bash
docker compose run --rm analysis "cryptic-ip --version"
```

### 3) Analyze a structure

```bash
docker compose run --rm analysis \
  "cryptic-ip analyze data/structures/adar2.pdb --output results/adar2_scores.csv"
```

### 4) Start Jupyter (optional)

```bash
docker compose --profile notebook up jupyter
```

Then open: `http://localhost:8888` and use token from `JUPYTER_TOKEN` (default: `cryptic-ip`).

## Files

- `Dockerfile`: Multi-stage image definition
- `docker-compose.yml`: analysis + optional notebook services
- `.dockerignore`: context slimming for faster builds
- `scripts/docker/build.sh`: local and multi-arch builds
- `scripts/docker/test.sh`: in-container validation script

## Volume Mounting

The compose setup mounts:

- `./data -> /workspace/data`
- `./results -> /workspace/results`
- `./scripts -> /workspace/scripts` (read-only for analysis service)

Recommended host layout:

```text
project/
├── data/
│   ├── structures/
│   └── reference/
├── results/
└── scripts/
```

## Environment Variables

### Analysis service

- `ANALYSIS_MEM_LIMIT` (default: `8g`)
- `ANALYSIS_CPUS` (default: `4`)
- `OMP_NUM_THREADS` (default: `4`)
- `OPENMM_CPU_THREADS` (default: `4`)

### Notebook service

- `JUPYTER_PORT` (default: `8888`)
- `JUPYTER_TOKEN` (default: `cryptic-ip`)
- `JUPYTER_MEM_LIMIT` (default: `8g`)
- `JUPYTER_CPUS` (default: `4`)

Example:

```bash
ANALYSIS_MEM_LIMIT=16g ANALYSIS_CPUS=8 docker compose run --rm analysis "cryptic-ip --help"
```

## Build & Test Workflow

### Standard local build

```bash
./scripts/docker/build.sh
```

### Multi-architecture build (amd64 + arm64)

```bash
./scripts/docker/build.sh --multi-arch
```

> Requires Docker Buildx and a configured builder capable of multi-platform builds.

### Run in-container tests

```bash
./scripts/docker/test.sh
```

This script performs:

1. Container image build for `analysis`
2. Targeted pytest run inside the container
3. ADAR2 validation run (`scripts/phase1_validate_adar2.py`)

To skip ADAR2 validation (for offline CI):

```bash
RUN_ADAR2=0 ./scripts/docker/test.sh
```

## Common Workflows

### Validate dependencies inside container

```bash
docker compose run --rm analysis "cryptic-ip check-dependencies"
```

### Proteome screening

```bash
docker compose run --rm analysis \
  "cryptic-ip screen data/structures --output results/screening_results.csv --max-structures 100"
```

### MD validation

```bash
docker compose run --rm analysis \
  "cryptic-ip md-validate results/screening_results.csv --output-dir results/md_validation --top-n 20"
```

## Optimization Notes

- **Image size control**:
  - Multi-stage build keeps runtime focused on the conda env + code.
  - `.dockerignore` excludes large data/results and git metadata.
- **Cache-friendly layering**:
  - Dependency files copied before full source for better layer reuse.
- **Conda cleanup**:
  - `micromamba clean --all --yes` removes package caches.
- **Reproducibility**:
  - Bioinformatics dependencies are sourced from conda-forge.

## Troubleshooting

### Build fails with memory errors

- Increase Docker Desktop memory allocation.
- Reduce parallel load and rebuild.
- Try: `ANALYSIS_MEM_LIMIT=12g` (or higher) for heavy workflows.

### Multi-arch build fails

- Ensure Buildx is installed and active:

```bash
docker buildx ls
```

- Create/use a compatible builder:

```bash
docker buildx create --use --name cryptic-ip-builder
```

### ADAR2 validation fails in CI/offline

- ADAR2 script may depend on networked resources.
- Run with `RUN_ADAR2=0` to isolate container/functionality checks.

### Notebook not reachable

- Verify port mapping and token:

```bash
docker compose --profile notebook ps
```

- Check logs:

```bash
docker compose --profile notebook logs -f jupyter
```
