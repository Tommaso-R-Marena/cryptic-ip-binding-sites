# Tutorial: Detecting Cryptic IP Binding Sites

This tutorial walks through the complete workflow from structure preparation to candidate identification.

## Table of Contents

1. [Quick Start](#quick-start)
2. [Understanding the Pipeline](#understanding-the-pipeline)
3. [Step-by-Step Workflow](#step-by-step-workflow)
4. [Interpreting Results](#interpreting-results)
5. [Advanced Usage](#advanced-usage)

---

## Quick Start

### Analyze a Single Protein

```bash
# Download ADAR2 structure
wget https://alphafold.ebi.ac.uk/files/AF-P78563-F1-model_v4.pdb

# Run analysis
cryptic-ip analyze AF-P78563-F1-model_v4.pdb --output results/

# View results
cat results/analysis_report.txt
```

**Expected output:**
```
Protein: ADAR2 (P78563)
Top-scoring pocket:
  Score: 0.87
  Volume: 642 ų
  Depth: 18.3 Å
  SASA: 2.1 ų
  Basic residues: 6
  ✓ HIGH CONFIDENCE buried IP binding site candidate
```

---

## Understanding the Pipeline

### What the Pipeline Does

The pipeline identifies **buried** inositol phosphate binding sites in protein structures through five steps:

1. **Pocket Detection** (fpocket): Find all cavities in the protein
2. **Solvent Accessibility** (FreeSASA): Calculate how buried each pocket is
3. **Electrostatics** (APBS): Map positive charge distribution
4. **Residue Analysis**: Count basic residues (Arg, Lys, His) near pockets
5. **Composite Scoring**: Combine metrics to identify cryptic IP sites

### What Makes a Good Candidate?

Cryptic IP binding sites have:
- **Deep burial**: >15 Å from protein surface
- **Low accessibility**: SASA <5 ų
- **Strong positive potential**: >5 kT/e
- **Multiple basic residues**: ≥4 within 5 Å
- **Appropriate volume**: 300-800 ų (fits IP3-IP6)

**Contrast with surface IP binding** (PH domains):
- Shallow pockets <8 Å deep
- High accessibility SASA >50 ų
- Exposed to solvent
- Reversible binding

---

## Step-by-Step Workflow

### Step 1: Prepare Your Structure

```python
from cryptic_ip.database.downloader import StructureDownloader

# Download from AlphaFold
downloader = StructureDownloader()
structure = downloader.download_alphafold('P78563')  # ADAR2 UniProt ID
print(f"Downloaded: {structure}")

# Or use your own PDB file
structure = 'my_protein.pdb'
```

### Step 2: Run Pocket Detection

```python
from cryptic_ip.analysis.pocket_detection import PocketDetector

detector = PocketDetector()
pockets = detector.detect_pockets(
    structure,
    output_dir='results/pockets'
)

print(f"Detected {len(pockets)} pockets")
for i, pocket in enumerate(pockets[:3]):
    print(f"Pocket {i+1}: Volume={pocket['volume']:.1f} ų, Depth={pocket['depth']:.1f} Å")
```

**Output:**
```
Detected 12 pockets
Pocket 1: Volume=642.3 ų, Depth=18.3 Å
Pocket 2: Volume=324.1 ų, Depth=12.7 Å
Pocket 3: Volume=215.8 ų, Depth=9.2 Å
```

### Step 3: Calculate Solvent Accessibility

```python
from cryptic_ip.analysis.sasa import SASACalculator

sasa_calc = SASACalculator()

# Calculate for each pocket
for pocket in pockets:
    pocket_sasa = sasa_calc.calculate_pocket_sasa(
        structure,
        pocket['residues']
    )
    pocket['sasa'] = pocket_sasa
    print(f"Pocket {pocket['id']}: SASA = {pocket_sasa:.2f} ų")
```

### Step 4: Compute Electrostatics

```python
from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

elec_calc = ElectrostaticsCalculator()

# Calculate potential at pocket centers
for pocket in pockets:
    potential = elec_calc.calculate_potential_at_point(
        structure,
        pocket['center']
    )
    pocket['potential'] = potential
    print(f"Pocket {pocket['id']}: Potential = {potential:.2f} kT/e")
```

### Step 5: Score and Rank

```python
from cryptic_ip.analysis.scoring import CompositeScorer

scorer = CompositeScorer()

# Score all pockets
for pocket in pockets:
    score = scorer.calculate_composite_score(pocket)
    pocket['score'] = score

# Sort by score
pockets.sort(key=lambda x: x['score'], reverse=True)

print("\nTop 3 candidates:")
for i, pocket in enumerate(pockets[:3]):
    print(f"{i+1}. Score: {pocket['score']:.3f} | "
          f"Volume: {pocket['volume']:.0f} ų | "
          f"Depth: {pocket['depth']:.1f} Å | "
          f"SASA: {pocket['sasa']:.1f} ų")
```

**Output:**
```
Top 3 candidates:
1. Score: 0.872 | Volume: 642 ų | Depth: 18.3 Å | SASA: 2.1 ų
2. Score: 0.543 | Volume: 324 ų | Depth: 12.7 Å | SASA: 8.4 ų
3. Score: 0.312 | Volume: 215 ų | Depth: 9.2 Å | SASA: 15.2 ų
```

### Step 6: Visualize Results

```python
import matplotlib.pyplot as plt
import pandas as pd

# Convert to DataFrame
df = pd.DataFrame(pockets)

# Plot score vs. key metrics
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

axes[0].scatter(df['depth'], df['score'])
axes[0].set_xlabel('Depth (Å)')
axes[0].set_ylabel('Composite Score')
axes[0].axhline(0.7, color='red', linestyle='--', label='Threshold')
axes[0].legend()

axes[1].scatter(df['sasa'], df['score'])
axes[1].set_xlabel('SASA (ų)')
axes[1].set_ylabel('Composite Score')
axes[1].axhline(0.7, color='red', linestyle='--')

axes[2].scatter(df['potential'], df['score'])
axes[2].set_xlabel('Potential (kT/e)')
axes[2].set_ylabel('Composite Score')
axes[2].axhline(0.7, color='red', linestyle='--')

plt.tight_layout()
plt.savefig('results/score_analysis.png', dpi=300)
plt.show()
```

---

## Interpreting Results

### Score Interpretation

| Score Range | Interpretation | Action |
|-------------|----------------|--------|
| **>0.7** | High confidence | Priority for experimental validation |
| **0.5-0.7** | Moderate confidence | Manual inspection recommended |
| **0.3-0.5** | Low confidence | Likely false positive |
| **<0.3** | Very low confidence | Exclude |

### What to Look For

**High-confidence candidate checklist:**
- [ ] Score >0.7
- [ ] Pocket volume 300-800 ų (appropriate for IP3-IP6)
- [ ] Depth >15 Å (deeply buried)
- [ ] SASA <5 ų (minimal solvent exposure)
- [ ] ≥4 basic residues within 5 Å
- [ ] Strong positive electrostatic potential (>5 kT/e)
- [ ] High pLDDT confidence (>70) if using AlphaFold

### Manual Inspection

Even high-scoring pockets should be visually inspected:

```bash
# Open in PyMOL
pymol results/pockets/protein_pocket1.pdb

# In PyMOL console:
show surface
color blue, resn ARG+LYS+HIS
zoom pocket
```

**Look for:**
- Pocket is genuinely internal (not just a surface cleft)
- Basic residues point inward toward pocket center
- No obvious structural artifacts
- Pocket shape compatible with IP geometry

---

## Advanced Usage

### Custom Scoring Weights

Adjust parameter importance:

```python
custom_weights = {
    'sasa': 0.35,      # Emphasize burial
    'depth': 0.25,
    'potential': 0.20,
    'basic_residues': 0.15,
    'volume': 0.05
}

scorer = CompositeScorer(weights=custom_weights)
```

### Batch Processing

```python
from cryptic_ip.pipeline import ScreeningPipeline
from pathlib import Path

pipeline = ScreeningPipeline(n_jobs=4)

# Process multiple structures
structures = list(Path('structures/').glob('*.pdb'))
results = pipeline.screen_batch(structures, output_dir='results/batch')

# Export results
results.to_csv('batch_results.csv', index=False)
```

### Filtering by pLDDT

For AlphaFold structures, filter low-confidence regions:

```python
from cryptic_ip.validation.quality import pLDDTFilter

filter = pLDDTFilter(threshold=70.0)

# Only analyze high-confidence pockets
for pocket in pockets:
    if filter.passes(structure, pocket['residues']):
        # Proceed with analysis
        pass
```

### Conservation Analysis

Check if predicted binding residues are conserved:

```python
from cryptic_ip.analysis.conservation import ConservationAnalyzer

analyzer = ConservationAnalyzer()

# Requires multiple sequence alignment
conservation = analyzer.calculate_conservation(
    pocket_residues=[376, 519, 522, 651, 672],
    alignment='alignments/ADAR2_orthologs.aln'
)

print(f"Mean conservation score: {conservation.mean():.2f}")
```

---

## Example: Complete Analysis

Putting it all together:

```python
from cryptic_ip.pipeline import AnalysisPipeline

# Initialize pipeline
pipeline = AnalysisPipeline(
    score_threshold=0.7,
    plddt_threshold=70.0,
    output_dir='results/complete_analysis'
)

# Run full analysis
result = pipeline.analyze(
    structure='AF-P78563-F1-model_v4.pdb',
    protein_name='ADAR2',
    uniprot_id='P78563'
)

# Generate report
pipeline.generate_report(
    result,
    output_file='results/ADAR2_report.pdf',
    include_figures=True,
    include_pymol_session=True
)

print(f"Analysis complete. Report saved to results/ADAR2_report.pdf")
```

---

## Next Steps

- **Validation**: See [Validation Analysis notebook](../notebooks/04_Validation_Analysis.ipynb)
- **Proteome screening**: See [Proteome Screening notebook](../notebooks/03_Proteome_Screening.ipynb)
- **API reference**: See [API documentation](API.md)
- **Examples**: Check `examples/` directory

## Troubleshooting

Common issues and solutions in [INSTALLATION.md](INSTALLATION.md#troubleshooting)

## Citation

If you use this tool in your research, please cite:

```
[Citation to be added upon publication]
```
