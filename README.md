# Cryptic IP Binding Sites

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Computational Pipeline for Identifying Buried Inositol Phosphate Binding Sites

### Overview

This project develops a computational pipeline to systematically screen protein structures for **cryptic inositol phosphate (IP) binding sites**—buried pockets where IPs function as structural folding cofactors rather than signaling molecules.

### The Scientific Question

Some proteins don't just bind inositol phosphates for signaling—they need them to fold. ADAR2 is the paradigm: an IP6 molecule sits completely buried in the enzyme core, invisible from the surface, and the protein has zero catalytic activity without it ([Macbeth et al., 2005, *Science*](https://doi.org/10.1126/science.1115248)).

**How many other proteins do this?** Nobody has looked systematically.

### Approach

We combine AlphaFold structures with structural bioinformatics tools:
- **fpocket**: Pocket detection and characterization
- **FreeSASA**: Solvent accessibility analysis
- **APBS**: Electrostatic potential mapping

The pipeline identifies buried, positively-charged cavities that match the chemical signature of cryptic IP sites.

### Three-Organism Comparative Screen

1. **Saccharomyces cerevisiae** (~6,000 proteins): Small proteome, fast genetics, immediate experimental validation
2. **Homo sapiens** (~23,000 proteins): Clinical relevance, known examples (ADAR2, HDACs)
3. **Dictyostelium discoideum** (~12,600 proteins): IP6 concentration ~520 µM—10× higher than mammals. Evolutionary test: does high IP abundance correlate with more buried sites?

### Project Phases

#### Phase 1: Tool Validation ✅
- Validate pipeline on ADAR2 (known buried IP6 site)
- Test against negative controls (PH domains—surface binding)
- Tune scoring parameters for buried vs. surface discrimination

#### Phase 2: Structural Databases 🔄
- Download and organize AlphaFold proteomes
- Quality control and formatting
- pLDDT confidence filtering

#### Phase 3: Screening at Scale 📊
- Run validated pipeline across three proteomes
- Rank and filter candidates
- Generate manuscript-ready datasets

### Repository Structure

```
cryptic-ip-binding-sites/
├── cryptic_ip/           # Main Python package
│   ├── analysis/         # Pocket detection and scoring
│   ├── database/         # Proteome management
│   ├── validation/       # Known site validation
│   └── visualization/    # PyMOL scripts and plotting
├── data/                 # Data directory (gitignored)
│   ├── structures/       # AlphaFold and PDB structures
│   ├── results/          # Analysis outputs
│   └── validation/       # Validation test cases
├── scripts/              # Standalone analysis scripts
├── notebooks/            # Jupyter notebooks
├── tests/                # Unit tests
├── docs/                 # Documentation
└── examples/             # Usage examples
```

### Installation

```bash
# Clone repository
git clone https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites.git
cd cryptic-ip-binding-sites

# Create conda environment
conda create -n cryptic-ip python=3.8
conda activate cryptic-ip

# Install dependencies
pip install -r requirements.txt
pip install -e .

# Install external tools
# fpocket: https://github.com/Discngine/fpocket
# APBS: http://www.poissonboltzmann.org/
# FreeSASA: https://freesasa.github.io/
```

### Quick Start

```python
from cryptic_ip.analysis import ProteinAnalyzer
from cryptic_ip.validation import validate_adar2

# Validate pipeline on ADAR2
results = validate_adar2()
print(f"IP6 binding site rank: {results['ip6_pocket_rank']}")
print(f"SASA at IP6 site: {results['sasa']:.2f} Ų")
print(f"Electrostatic potential: {results['potential']:.2f} kT/e")

# Analyze a protein structure
analyzer = ProteinAnalyzer("path/to/structure.pdb")
pockets = analyzer.detect_pockets()
scores = analyzer.score_pockets(pockets)
candidates = analyzer.filter_candidates(scores, threshold=0.7)
```

### Pipeline Criteria

Cryptic IP site candidates must satisfy:
1. **Pocket depth** >15 Å from protein surface
2. **Low solvent accessibility** (SASA <5 ų)
3. **Strong positive electrostatics** (>5 kT/e)
4. **Basic residue cluster** ≥4 Arg/Lys/His within 5 Å
5. **Appropriate volume** 300–800 ų (IP3–IP6 range)
6. **High structure confidence** (pLDDT ≥70 for pocket-lining residues)

### Validation Targets

**Positive Controls** (should detect):
- ADAR2 (PDB: 1ZY7) - buried IP6
- Pds5B (PDB: 5HDT) - buried IP6
- HDAC1/3 (PDB: 5ICN/4A69) - interface IP4

**Negative Controls** (should reject):
- PLCδ1 PH domain (PDB: 1MAI) - surface IP3
- Btk PH domain (PDB: 1BTK) - surface IP4

### Expected Outputs

- Ranked candidate list for each proteome
- Comparative hit rate analysis
- Structural visualizations
- Conservation analysis
- Manuscript-ready figures

### Citation

If you use this pipeline, please cite:

```bibtex
@software{marena2026cryptic,
  author = {Marena, Tommaso R.},
  title = {Cryptic IP Binding Sites: Computational Pipeline},
  year = {2026},
  url = {https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites}
}
```

Key references:
- Macbeth et al. (2005). *Science* 309:1534-1539. [ADAR2-IP6 structure]
- Dick et al. (2018). *Nature* 560:509-512. [HIV-1 capsid-IP6]

### Contributing

This is a research project under active development. Contributions welcome:
- Bug reports and feature requests via Issues
- Pull requests for improvements
- Validation on additional known IP-binding proteins

### License

MIT License - see LICENSE file

### Contact

Tommaso R. Marena  
The Catholic University of America  
GitHub: [@Tommaso-R-Marena](https://github.com/Tommaso-R-Marena)

---

**Status**: Phase 1 Implementation (Validation Pipeline)

**Last Updated**: February 2026
