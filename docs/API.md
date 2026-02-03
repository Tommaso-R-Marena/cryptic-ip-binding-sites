# API Reference

Complete API documentation for `cryptic-ip-binding-sites`.

## Table of Contents

1. [Pipeline Classes](#pipeline-classes)
2. [Analysis Modules](#analysis-modules)
3. [Database Management](#database-management)
4. [Validation Tools](#validation-tools)
5. [Utilities](#utilities)

---

## Pipeline Classes

### `AnalysisPipeline`

Main entry point for analyzing single structures.

```python
from cryptic_ip.pipeline import AnalysisPipeline

pipeline = AnalysisPipeline(
    score_threshold=0.7,
    plddt_threshold=70.0,
    output_dir='results/'
)
```

**Parameters:**
- `score_threshold` (float): Minimum composite score for candidates (default: 0.7)
- `plddt_threshold` (float): Minimum AlphaFold confidence (default: 70.0)
- `output_dir` (str | Path): Output directory for results
- `n_jobs` (int): Number of parallel processes (default: 1)

**Methods:**

#### `analyze(structure, **kwargs)`

Analyze a single protein structure.

```python
result = pipeline.analyze(
    structure='protein.pdb',
    protein_name='MyProtein',
    uniprot_id='P12345'
)
```

**Parameters:**
- `structure` (str | Path): Path to PDB file
- `protein_name` (str, optional): Protein name for reporting
- `uniprot_id` (str, optional): UniProt identifier

**Returns:**
- `AnalysisResult`: Object containing pockets, scores, and metadata

**AnalysisResult attributes:**
- `pockets` (list): List of detected pockets with scores
- `top_candidate` (dict): Highest-scoring pocket
- `metadata` (dict): Protein information
- `passed_threshold` (bool): Whether any pocket passed threshold

---

### `ScreeningPipeline`

Batch processing for proteome-wide screening.

```python
from cryptic_ip.pipeline import ScreeningPipeline

pipeline = ScreeningPipeline(
    n_jobs=8,
    batch_size=100,
    output_dir='results/screening'
)
```

**Parameters:**
- `n_jobs` (int): Number of parallel processes
- `batch_size` (int): Structures to process per batch (for memory management)
- `output_dir` (str | Path): Output directory
- `checkpoint_interval` (int): Save progress every N structures (default: 500)

**Methods:**

#### `screen_proteome(proteome, **kwargs)`

Screen entire proteome.

```python
results_df = pipeline.screen_proteome(
    proteome='yeast',  # or 'human', 'dictyostelium'
    structure_dir='structures/yeast/',
    resume_from_checkpoint=True
)
```

**Parameters:**
- `proteome` (str): Proteome identifier ('yeast', 'human', 'dictyostelium')
- `structure_dir` (str | Path): Directory containing PDB files
- `resume_from_checkpoint` (bool): Resume interrupted screening (default: False)
- `max_structures` (int, optional): Limit number of structures (for testing)

**Returns:**
- `pandas.DataFrame`: Results table with all scored proteins

---

## Analysis Modules

### Pocket Detection

#### `PocketDetector`

Wrapper for fpocket with additional analysis.

```python
from cryptic_ip.analysis.pocket_detection import PocketDetector

detector = PocketDetector(
    fpocket_path='/usr/local/bin/fpocket',
    min_volume=100.0
)
```

**Parameters:**
- `fpocket_path` (str, optional): Path to fpocket executable
- `min_volume` (float): Minimum pocket volume in ų (default: 100.0)
- `max_pockets` (int): Maximum pockets to return (default: 20)

**Methods:**

#### `detect_pockets(structure, output_dir)`

```python
pockets = detector.detect_pockets(
    'protein.pdb',
    output_dir='results/pockets'
)
```

**Returns:**
- `list[dict]`: List of pocket dictionaries with keys:
  - `id` (int): Pocket identifier (1-indexed)
  - `volume` (float): Pocket volume in ų
  - `depth` (float): Depth from surface in Å
  - `score` (float): fpocket druggability score
  - `center` (tuple): (x, y, z) coordinates of pocket center
  - `residues` (list): Residue IDs lining the pocket

---

### Solvent Accessibility

#### `SASACalculator`

Calculate solvent-accessible surface area.

```python
from cryptic_ip.analysis.sasa import SASACalculator

calc = SASACalculator(
    probe_radius=1.4,  # Water probe
    n_points=100       # Resolution
)
```

**Methods:**

#### `calculate_pocket_sasa(structure, residues)`

Calculate SASA for pocket-lining residues.

```python
sasa = calc.calculate_pocket_sasa(
    'protein.pdb',
    residues=[123, 145, 178, 201]
)
# Returns: 3.42 (very buried)
```

#### `calculate_per_residue_sasa(structure)`

Get SASA for all residues.

```python
residue_sasa = calc.calculate_per_residue_sasa('protein.pdb')
# Returns: dict {residue_id: sasa_value}
```

---

### Electrostatics

#### `ElectrostaticsCalculator`

Calculate electrostatic potential using APBS.

```python
from cryptic_ip.analysis.electrostatics import ElectrostaticsCalculator

calc = ElectrostaticsCalculator(
    apbs_path='/usr/local/bin/apbs',
    grid_spacing=0.5,
    temperature=298.15
)
```

**Parameters:**
- `apbs_path` (str, optional): Path to APBS executable
- `grid_spacing` (float): Grid spacing in Å (default: 0.5)
- `temperature` (float): Temperature in K (default: 298.15)
- `ion_concentration` (float): Salt concentration in M (default: 0.15)

**Methods:**

#### `calculate_potential_at_point(structure, point)`

Get potential at specific coordinates.

```python
potential = calc.calculate_potential_at_point(
    'protein.pdb',
    point=(12.3, 45.6, 78.9)
)
# Returns: 7.2 (in kT/e, strong positive)
```

#### `calculate_potential_map(structure)`

Generate full 3D potential map.

```python
map_data = calc.calculate_potential_map('protein.pdb')
# Returns: numpy array with potential values
```

---

### Scoring

#### `CompositeScorer`

Combine metrics into single score.

```python
from cryptic_ip.analysis.scoring import CompositeScorer

scorer = CompositeScorer(weights={
    'sasa': 0.30,
    'depth': 0.25,
    'potential': 0.25,
    'basic_residues': 0.15,
    'volume': 0.05
})
```

**Methods:**

#### `calculate_composite_score(pocket_data)`

```python
score = scorer.calculate_composite_score({
    'volume': 642.0,
    'depth': 18.3,
    'sasa': 2.1,
    'potential': 7.5,
    'basic_residues': 6
})
# Returns: 0.87 (high confidence)
```

**Score formula:**

\[
\text{score} = \sum_{i} w_i \cdot f_i(x_i)
\]

Where:
- \( w_i \) = weight for metric \( i \)
- \( f_i \) = normalization function (sigmoid or linear)
- \( x_i \) = raw metric value

---

## Database Management

### `ProteomeDatabase`

Manage structure databases and metadata.

```python
from cryptic_ip.database.manager import ProteomeDatabase

db = ProteomeDatabase(db_path='structures.db')
```

**Methods:**

#### `register_structure(uniprot_id, organism, file_path, **metadata)`

```python
db.register_structure(
    uniprot_id='P78563',
    organism='human',
    file_path='/path/to/AF-P78563-F1.pdb',
    avg_plddt=87.5,
    protein_length=702
)
```

#### `get_structures_by_organism(organism)`

```python
yeast_structures = db.get_structures_by_organism('yeast')
# Returns: list of structure records
```

#### `get_proteome_metadata(proteome)`

```python
meta = db.get_proteome_metadata('yeast')
print(meta)
# {'uniprot_id': 'UP000002311',
#  'organism': 'Saccharomyces cerevisiae',
#  'protein_count': 6049,
#  'ip6_concentration': 20.0}
```

---

## Validation Tools

### `ValidationSuite`

Validate pipeline performance.

```python
from cryptic_ip.validation.validation import ValidationSuite

validator = ValidationSuite()
```

**Methods:**

#### `validate_adar2(alphafold_pdb, crystal_pdb, output_dir)`

Run complete ADAR2 validation.

```python
results = validator.validate_adar2(
    'AF-P78563-F1-model_v4.pdb',
    '1ZY7.pdb',
    output_dir='validation_results/'
)

print(results)
# {'pocket_identified': True,
#  'score': 0.87,
#  'rmsd': 1.23,
#  'passed': True}
```

#### `validate_negative_control(structure, output_dir)`

Test that surface IP sites score lowly.

```python
results = validator.validate_negative_control(
    '1MAI.pdb',  # PH domain
    output_dir='validation_results/'
)

print(results['max_score'])
# 0.32 (correctly low)
```

#### `calculate_rmsd(structure1, structure2, selection=None)`

Align and calculate RMSD.

```python
rmsd = validator.calculate_rmsd(
    'alphafold.pdb',
    'crystal.pdb',
    selection='resid 376+519+522+651+672'  # Binding site only
)
# Returns: 1.23 (Å)
```

---

## Utilities

### Residue Analysis

```python
from cryptic_ip.utils.residues import count_basic_residues, get_residue_type

basic_count = count_basic_residues(
    'protein.pdb',
    pocket_center=(12.3, 45.6, 78.9),
    radius=5.0
)
# Returns: 6
```

### Structure I/O

```python
from cryptic_ip.utils.structure_io import load_structure, clean_structure

# Load structure
structure = load_structure('protein.pdb')

# Remove waters, ligands
cleaned = clean_structure(structure, remove_waters=True, remove_ligands=True)
```

### Coordinate Operations

```python
from cryptic_ip.utils.geometry import distance, center_of_mass

dist = distance(
    point1=(0, 0, 0),
    point2=(3, 4, 0)
)
# Returns: 5.0

com = center_of_mass(coordinates=[
    (0, 0, 0),
    (1, 1, 1),
    (2, 2, 2)
])
# Returns: (1.0, 1.0, 1.0)
```

---

## Command Line Interface

All functionality is also available via CLI:

```bash
# Single structure analysis
cryptic-ip analyze protein.pdb --output results/

# Proteome screening
cryptic-ip screen --proteome yeast --structures structures/yeast/ --jobs 8

# Validation
cryptic-ip validate --structure AF-P78563-F1-model_v4.pdb

# Database management
cryptic-ip db --init structures.db
cryptic-ip db --register P78563 human /path/to/structure.pdb
```

See `cryptic-ip --help` for complete CLI reference.

---

## Type Hints and IDE Support

All functions include complete type hints:

```python
from typing import List, Dict, Tuple, Optional
from pathlib import Path

def detect_pockets(
    structure: Path | str,
    output_dir: Path | str,
    min_volume: float = 100.0
) -> List[Dict[str, Any]]:
    ...
```

This enables full IDE autocomplete and type checking.

---

## Examples

More examples in `examples/` directory:
- `examples/basic_analysis.py`
- `examples/batch_screening.py`
- `examples/custom_scoring.py`
- `examples/conservation_analysis.py`

---

## Contributing

See [CONTRIBUTING.md](../CONTRIBUTING.md) for API development guidelines.
