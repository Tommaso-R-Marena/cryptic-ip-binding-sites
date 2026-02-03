# Getting Started

## Installation

### Prerequisites

- Python 3.8+
- fpocket (pocket detection)
- APBS (electrostatics, optional)
- FreeSASA (solvent accessibility, handled by ProDy)

### Installing fpocket

**Linux/Mac:**
```bash
# Install fpocket from GitHub
git clone https://github.com/Discngine/fpocket.git
cd fpocket
make
sudo make install
```

**Or using conda:**
```bash
conda install -c conda-forge fpocket
```

### Installing the Pipeline

```bash
# Clone repository
git clone https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites.git
cd cryptic-ip-binding-sites

# Create environment
conda create -n cryptic-ip python=3.8
conda activate cryptic-ip

# Install
pip install -r requirements.txt
pip install -e .
```

## Quick Test

Validate installation by testing on ADAR2:

```bash
cryptic-ip validate
```

This will:
1. Download ADAR2 AlphaFold structure
2. Detect pockets
3. Score and rank them
4. Report if IP6 site was found

## Basic Usage

### Analyze a Single Structure

```bash
cryptic-ip analyze path/to/structure.pdb
```

### Download Proteomes

```bash
# Download yeast proteome (~15 GB)
cryptic-ip download yeast

# Download all proteomes
cryptic-ip download yeast human dictyostelium
```

### Screen Proteome

```bash
cryptic-ip screen data/structures/yeast --output yeast_results.csv
```

## Python API

```python
from cryptic_ip.analysis import ProteinAnalyzer

# Analyze structure
analyzer = ProteinAnalyzer('structure.pdb')
pockets = analyzer.detect_pockets()
scored = analyzer.score_all_pockets()

# Get top candidates
candidates = scored[scored['composite_score'] >= 0.70]
print(candidates)
```

## Next Steps

- See [notebooks/](../notebooks/) for detailed examples
- Read [PIPELINE.md](PIPELINE.md) for technical details
- Check [FAQ.md](FAQ.md) for common issues
