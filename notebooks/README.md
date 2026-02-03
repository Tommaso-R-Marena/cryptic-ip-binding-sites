# Jupyter Notebooks

This directory contains interactive Jupyter notebooks for analyzing cryptic IP binding sites.

## Notebooks

### 01_Quick_Start.ipynb
**Purpose**: Get started quickly with the pipeline  
**Content**: Basic usage examples, simple workflow demonstration  
**Time**: 10-15 minutes  
**Prerequisites**: None

### 02_ADAR2_Analysis.ipynb
**Purpose**: Deep dive into ADAR2 validation  
**Content**: Detailed analysis of the ADAR2 positive control, pocket detection, scoring  
**Time**: 30 minutes  
**Prerequisites**: Understanding of protein structure basics

### 03_Proteome_Screening.ipynb
**Purpose**: Run proteome-wide screening  
**Content**: Batch processing, parallelization, progress monitoring  
**Time**: Variable (hours to days depending on proteome size)  
**Prerequisites**: Completed validation in notebook 02

### 04_Validation_Analysis.ipynb
**Purpose**: Comprehensive validation workflow  
**Content**:
- ADAR2 positive control analysis
- PH domain negative control testing
- Score separation visualization
- AlphaFold vs crystal structure comparison
- Parameter optimization

**Time**: 45 minutes  
**Prerequisites**: Validation structures downloaded

### 05_Results_Analysis.ipynb
**Purpose**: Analyze screening results across proteomes  
**Content**:
- Load and filter results from yeast, human, Dictyostelium
- Identify high-confidence candidates
- Compare hit rates across organisms
- Test evolutionary co-evolution hypothesis
- Functional enrichment analysis
- Export candidate list for experimental validation

**Time**: 1 hour  
**Prerequisites**: Completed proteome screening (notebook 03)

## Usage

### Starting Jupyter

```bash
# From repository root
jupyter notebook notebooks/
```

### Recommended Order

1. **Quick Start** (01) → Get familiar with the tools
2. **ADAR2 Analysis** (02) → Understand the validation target
3. **Validation Analysis** (04) → Validate the pipeline thoroughly
4. **Proteome Screening** (03) → Run large-scale screening
5. **Results Analysis** (05) → Analyze and interpret findings

### Data Requirements

Some notebooks require structure files that must be downloaded separately:

```bash
# Download validation structures
mkdir -p structures/validation
cd structures/validation

# ADAR2 AlphaFold
wget https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb

# ADAR2 crystal structure
wget https://files.rcsb.org/download/1ZY7.pdb

# PH domain negative control
wget https://files.rcsb.org/download/1MAI.pdb
```

### Output Directories

Notebooks will create output in:
- `results/validation/` - Validation analyses
- `results/analysis/` - Screening result analyses
- `results/figures/` - Generated plots and visualizations

## Tips

- **Save frequently**: Notebooks can take time to run
- **Check outputs**: Each cell should produce informative output
- **Restart kernel if needed**: If results seem inconsistent, restart the kernel and re-run
- **Modify parameters**: Feel free to experiment with thresholds and settings

## Troubleshooting

**Problem**: `ModuleNotFoundError: No module named 'cryptic_ip'`  
**Solution**: Install the package: `pip install -e .` from repository root

**Problem**: Structure files not found  
**Solution**: Check file paths and download structures as shown above

**Problem**: fpocket or APBS not found  
**Solution**: Install dependencies: see main README.md installation section

**Problem**: Jupyter kernel crashes during large screening  
**Solution**: Use the command-line interface instead: `cryptic-ip screen` (see CLI documentation)

## Questions?

See the main repository README or open an issue on GitHub.
