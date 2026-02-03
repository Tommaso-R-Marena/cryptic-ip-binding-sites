# Installation Guide

## System Requirements

### Operating System
- Linux (tested on Ubuntu 20.04+)
- macOS (tested on 10.15+)
- Windows (via WSL2)

### Minimum Hardware
- **CPU**: 4 cores recommended
- **RAM**: 8 GB minimum, 16 GB recommended for proteome screening
- **Storage**: 150 GB for complete proteome databases
  - Yeast: ~20 GB
  - Human: ~65 GB
  - Dictyostelium: ~40 GB
  - Plus working space for results

### Software Dependencies
- **Python**: 3.8 or higher
- **fpocket**: 3.0+ (pocket detection)
- **FreeSASA**: 2.0+ (solvent accessibility)
- **APBS**: 3.0+ (electrostatics)
- **PyMOL**: (optional, for visualization)

---

## Quick Install

For users who just want to get started:

```bash
# Clone repository
git clone https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites.git
cd cryptic-ip-binding-sites

# Install Python package
pip install -e .

# Install external tools (Ubuntu/Debian)
sudo apt-get install fpocket freesasa apbs

# Verify installation
cryptic-ip --version
cryptic-ip check-dependencies
```

---

## Detailed Installation

### Step 1: Clone Repository

```bash
git clone https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites.git
cd cryptic-ip-binding-sites
```

### Step 2: Create Virtual Environment (Recommended)

```bash
# Using venv
python3 -m venv cryptic_ip_env
source cryptic_ip_env/bin/activate  # On Windows: cryptic_ip_env\Scripts\activate

# OR using conda
conda create -n cryptic_ip python=3.9
conda activate cryptic_ip
```

### Step 3: Install Python Dependencies

```bash
# Basic installation
pip install -e .

# Development installation (includes testing tools)
pip install -e ".[dev]"

# Full installation (includes Jupyter, visualization tools)
pip install -e ".[full]"
```

### Step 4: Install External Tools

#### fpocket (Pocket Detection)

**Ubuntu/Debian:**
```bash
sudo apt-get install fpocket
```

**macOS:**
```bash
brew install brewsci/bio/fpocket
```

**From source:**
```bash
wget https://github.com/Discngine/fpocket/archive/refs/tags/4.1.tar.gz
tar -xzf 4.1.tar.gz
cd fpocket-4.1
make
sudo make install
```

**Verify:**
```bash
fpocket --version
```

#### FreeSASA (Solvent Accessibility)

**Ubuntu/Debian:**
```bash
sudo apt-get install freesasa
```

**macOS:**
```bash
brew install freesasa
```

**From source:**
```bash
wget https://github.com/mittinatten/freesasa/releases/download/2.1.2/freesasa-2.1.2.tar.gz
tar -xzf freesasa-2.1.2.tar.gz
cd freesasa-2.1.2
./configure
make
sudo make install
```

**Verify:**
```bash
freesasa --version
```

#### APBS (Electrostatics)

**Ubuntu/Debian:**
```bash
sudo apt-get install apbs
```

**macOS:**
```bash
brew install apbs
```

**From conda:**
```bash
conda install -c conda-forge apbs
```

**Verify:**
```bash
apbs --version
```

#### pdb2pqr (PDB to PQR conversion for APBS)

```bash
pip install pdb2pqr

# Verify
pdb2pqr --version
```

### Step 5: Verify Installation

```bash
# Check that all tools are accessible
cryptic-ip check-dependencies

# Should output:
# ✓ Python: 3.9.x
# ✓ fpocket: 4.1
# ✓ FreeSASA: 2.1.2
# ✓ APBS: 3.0.0
# ✓ pdb2pqr: 3.1.0
# All dependencies satisfied!
```

---

## Optional Tools

### PyMOL (Structure Visualization)

**Open-source version:**
```bash
# Ubuntu/Debian
sudo apt-get install pymol

# macOS
brew install pymol

# From conda
conda install -c conda-forge pymol-open-source
```

**Commercial version** (better performance):
Download from [pymol.org](https://pymol.org)

### Jupyter (Interactive Notebooks)

```bash
pip install jupyter notebook ipywidgets
jupyter notebook
```

---

## Database Setup

### Download Validation Structures

```bash
# Create directory structure
mkdir -p structures/validation
cd structures/validation

# ADAR2 (positive control)
wget https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb
wget https://files.rcsb.org/download/1ZY7.pdb

# PH domains (negative controls)
wget https://files.rcsb.org/download/1MAI.pdb  # PLCδ1
wget https://files.rcsb.org/download/1BTK.pdb  # Btk

# Additional validation examples
wget https://files.rcsb.org/download/5HDT.pdb  # Pds5B
wget https://files.rcsb.org/download/5ICN.pdb  # HDAC1
```

### Download Proteome Databases (Optional)

Only needed for large-scale screening:

```bash
cd structures

# Yeast proteome (~6,000 proteins, ~15 GB)
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002311_559292_YEAST_v4.tar
tar -xf UP000002311_559292_YEAST_v4.tar -C yeast/

# Human proteome (~23,000 proteins, ~50 GB)
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000005640_9606_HUMAN_v4.tar
tar -xf UP000005640_9606_HUMAN_v4.tar -C human/

# Dictyostelium proteome (~12,600 proteins, ~30 GB)
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000002195_44689_DICDI_v4.tar
tar -xf UP000002195_44689_DICDI_v4.tar -C dictyostelium/
```

---

## Troubleshooting

### Common Issues

#### "fpocket: command not found"

**Solution**: Make sure fpocket is installed and in your PATH:
```bash
which fpocket
# If not found, install using method above

# Add to PATH if needed (add to ~/.bashrc or ~/.zshrc):
export PATH="/path/to/fpocket/bin:$PATH"
```

#### "ModuleNotFoundError: No module named 'cryptic_ip'"

**Solution**: Install the package in editable mode:
```bash
cd /path/to/cryptic-ip-binding-sites
pip install -e .
```

#### APBS calculation failures

**Symptom**: APBS runs but produces errors or warnings  
**Common causes**:
- Grid size issues: protein too large for default grid
- Charge assignment problems: unusual residues or ligands

**Solution**: Check APBS logs in output directory and adjust parameters:
```python
from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

calc = ElectrostaticsCalculator(grid_spacing=1.0)  # Adjust grid
```

#### Memory errors during proteome screening

**Symptom**: Process killed or "MemoryError"  
**Solution**: Process in smaller batches:
```bash
cryptic-ip screen --proteome yeast --batch-size 100
```

#### Slow performance

**Tips**:
- Use parallel processing: `--n-jobs 4`
- Use SSD for database storage
- Filter by protein size: `--max-length 1000`
- Skip low-confidence structures: `--min-plddt 70`

---

## Cluster/HPC Installation

For running on high-performance computing systems:

### Module-based Systems

```bash
# Load modules
module load python/3.9
module load fpocket
module load apbs

# Install in user space
pip install --user -e .
```

### Containerized Deployment

Dockerfile provided in repository:

```bash
# Build container
docker build -t cryptic-ip .

# Run analysis
docker run -v $(pwd)/data:/data cryptic-ip \
  cryptic-ip analyze /data/protein.pdb --output /data/results
```

Singularity image:
```bash
# Convert Docker image
singularity build cryptic-ip.sif docker://cryptic-ip:latest

# Run on cluster
singularity exec cryptic-ip.sif cryptic-ip screen --proteome yeast
```

---

## Testing Installation

```bash
# Run test suite
pytest tests/

# Quick validation test
cryptic-ip validate --structure structures/validation/AF-P78563-F1-model_v4.pdb

# Expected output:
# Running validation on ADAR2...
# ✓ Pocket detection: PASSED
# ✓ Scoring: PASSED (score=0.87)
# ✓ Separation test: PASSED
# Validation successful!
```

---

## Next Steps

After installation:
1. Read the [Tutorial](TUTORIAL.md)
2. Try the [Quick Start notebook](../notebooks/01_Quick_Start.ipynb)
3. Review the [API documentation](API.md)

## Getting Help

- **Documentation**: See `docs/` directory
- **Issues**: [GitHub Issues](https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites/issues)
- **Email**: Your contact info here
