# Cryptic IP Binding Sites

[![Validation Pipeline](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/validate_pipeline.yml/badge.svg)](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/validate_pipeline.yml)
[![Run Tests](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/run_tests.yml/badge.svg)](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/run_tests.yml)
[![Docker Publish](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/docker-publish.yml/badge.svg)](https://github.com/OWNER/cryptic-ip-binding-sites/actions/workflows/docker-publish.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Checked with mypy](https://www.mypy-lang.org/static/mypy_badge.svg)](https://mypy-lang.org/)
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit)](https://github.com/pre-commit/pre-commit)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/OWNER/cryptic-ip-binding-sites/blob/main/notebooks/CrypticIP_FullPipeline.ipynb)

## Scientific Background

This project builds a computational pipeline to screen entire proteomes for proteins that bury inositol phosphates (IPs) in their structural core — not for signaling, but because they cannot fold without them. 

ADAR2 (UniProt P78563, PDB 1ZY7) is the gold-standard example: an IP6 molecule sits completely buried in its enzyme core (SASA ≈ 0 Å²), coordinated by key basic residues (K376, K519, R522, R651, K672, W687), and is absolutely required for catalytic activity. This project attempts to find similar *cryptic structural* IP-binding sites across three eukaryotic proteomes (*S. cerevisiae*, *H. sapiens*, and *D. discoideum*).

We distinguish between two paradigms:
1. **Surface IP Binding** (e.g., PH domains): Reversible signaling, high solvent accessibility, shallow pockets. We use these as *negative controls*.
2. **Buried IP Binding** (e.g., ADAR2, HDACs): Structural cofactors, deep pockets, extremely low solvent accessibility. These are our *positive controls*.

## Pipeline Overview

The pipeline analyzes each AlphaFold `.pdb` structure via the following 8 steps:

```text
AlphaFold Structure -> [1] Clean & prep (extract pLDDT) -> [2] Detect Pockets (fpocket)
-> [3] SASA Filter (FreeSASA) -> [4] Electrostatics (pdb2pqr+APBS) -> [5] Basic Residue Census
-> [6] Composite Scoring -> [7] pLDDT Quality Filter -> [8] Output JSON & CSV
```

## Quick Start

### Validation (Local, Fast)
Run the pipeline locally to process only the positive and negative controls (~10 min on CPU):
```bash
python -m pipeline --mode=validation
```

### Full Run (GitHub Actions)
The most efficient way to run this pipeline is using GitHub Actions (found in `.github/workflows/`):
1. Navigate to the **Actions** tab on GitHub.
2. Manually trigger `Fetch Yeast Proteome` (or Human / Dictyostelium), which downloads and validates the data.
3. Once fetched, manually trigger `Screen Yeast Proteome` to process the entire dataset in parallel.
4. Results are automatically published as GitHub Releases.

### Colab
You can test the core pipeline logic directly in Google Colab (no installation required) by clicking the badge above.

## Installation

### Dependencies
Required standard tools (Ubuntu):
```bash
sudo apt-get install -y fpocket apbs pdb2pqr
```

Python dependencies:
```bash
pip install -r requirements.txt
```

Alternatively, you can use the provided install script for an automated Ubuntu setup:
```bash
bash scripts/install_dependencies.sh
```

Or run via Docker:
```bash
docker build -t cryptic-ip .
docker run --rm cryptic-ip --help
```

## Configuration
All parameters are exposed in `config/config.yaml`.
This includes volume limits for `fpocket`, basic residue cutoffs, and the weights for the composite scoring formula:
`S = w1*(normalized_depth) + w2*(1 - SASA) + w3*(normalized_electrostatic) + w4*(normalized_basic_count)`

## Output Format
Master results are output as a CSV containing:
- `UniProtID`, `Organism`, `ProteinName`
- `PocketID`, `CompositeScore`, `PocketVolume`, `MeanSASA`, `ElectrostaticPotential`, `BasicResidueCount`, `MeanPLDDT`
- `PassAllFilters`, `ManualReviewFlag` (candidates > 0.7 composite score)

## Validation Targets

The pipeline executes a firm GO/NO-GO checkpoint. It will fail unless:
- ADAR2 IP6 site ranks in the top 3 detected pockets
- 4 Negative Controls (PH domains) score below 0.4
- No overlap exists between positive and negative control score distributions

## Contributing
This codebase is intended to scale to massive datasets automatically using CI/CD. Please ensure all modifications pass the `validate_pipeline.yml` tests.
