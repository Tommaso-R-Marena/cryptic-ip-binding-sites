# Testing Guide

This document explains how to run tests and verify that the pipeline works with real data.

## Quick Test

Run all tests:

```bash
pytest tests/ -v
```

## Test Categories

### 1. Database Client Tests

Test real API integrations:

```bash
pytest tests/test_database_clients.py -v
```

These tests verify:
- ✅ AlphaFold API connectivity and structure downloads
- ✅ PDB API connectivity and crystal structure downloads
- ✅ UniProt API connectivity and protein annotations
- ✅ Complete workflow: UniProt → AlphaFold + PDB

**What gets downloaded:**
- ADAR2 AlphaFold structure (P78563, ~200 KB)
- ADAR2 crystal structure (1ZY7, ~100 KB)
- Protein metadata from all three databases

### 2. Structure Validation Tests

Test real structure parsing and analysis:

```bash
pytest tests/test_structure_validation.py -v
```

These tests verify:
- ✅ Downloading structures from AlphaFold and PDB
- ✅ Parsing PDB files with BioPython
- ✅ Extracting pLDDT confidence scores
- ✅ Finding basic residues (Arg, Lys, His)
- ✅ Verifying known IP6-binding residues exist
- ✅ API connectivity for all services

### 3. Network Tests

Tests that require internet connection are marked:

```bash
# Run only network tests
pytest tests/ -v -m requires_network

# Skip network tests
pytest tests/ -v -m "not requires_network"
```

## Notebook Tests

Verify notebooks execute correctly:

```bash
# Install nbconvert if needed
pip install nbconvert

# Execute all notebooks
jupyter nbconvert --execute --to notebook \
  --inplace notebooks/*.ipynb

# Or run specific notebook
jupyter nbconvert --execute --to notebook \
  --inplace notebooks/01_Quick_Start.ipynb
```

**What notebooks do:**

1. **01_Quick_Start.ipynb**
   - Downloads ADAR2 from AlphaFold (real data)
   - Loads and inspects structure with BioPython
   - Visualizes pLDDT confidence scores
   - Identifies basic residues for IP coordination
   - ✅ **Works in Colab and Binder**

2. **02_ADAR2_Analysis.ipynb**
   - Downloads both AlphaFold prediction and crystal structure
   - Compares structures
   - Finds IP6 molecule in crystal
   - Identifies known coordinating residues
   - Visualizes binding site region
   - ✅ **Works in Colab and Binder**

3. **03_Proteome_Screening.ipynb**
   - Demonstrates batch processing
   - Shows how to screen entire proteomes
   - Uses real proteome data sources

4. **04_Validation_Analysis.ipynb**
   - Tests on positive controls (ADAR2, IP3K, BTK)
   - Tests on negative controls
   - Validates scoring metrics

5. **05_Results_Visualization.ipynb**
   - Summarizes screening results
   - Creates publication-quality figures
   - Analyzes proteome-wide patterns

## Test Coverage

Generate coverage report:

```bash
pytest tests/ --cov=cryptic_ip --cov-report=html

# View report
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
start htmlcov/index.html  # Windows
```

## Continuous Integration

GitHub Actions automatically runs tests on:
- Push to any branch
- Pull requests
- Multiple Python versions (3.8, 3.9, 3.10, 3.11)
- Multiple operating systems (Ubuntu, macOS)

See `.github/workflows/ci.yml` for details.

## Manual Validation

Validate structures manually:

```bash
# Validate ADAR2 structures
python scripts/validate_structures.py --proteins P78563

# Check AlphaFold database versions
python scripts/check_alphafold_versions.py

# Validate specific PDB entries
python scripts/validate_pdb_entries.py --pdb-ids 1ZY7,5HDT,1MAI
```

## Expected Test Output

### Successful Test Run

```
$ pytest tests/ -v

======================== test session starts =========================
platform linux -- Python 3.9.0, pytest-7.4.0
cachedir: .pytest_cache
rootdir: /path/to/cryptic-ip-binding-sites
configfile: pytest.ini
testpaths: tests
collected 15 items

tests/test_database_clients.py::TestAlphaFoldClient::test_init PASSED
tests/test_database_clients.py::TestAlphaFoldClient::test_fetch_structure_adar2 PASSED
tests/test_database_clients.py::TestAlphaFoldClient::test_get_metadata PASSED
tests/test_database_clients.py::TestPDBClient::test_fetch_structure_1zy7 PASSED
tests/test_database_clients.py::TestPDBClient::test_get_entry_info PASSED
tests/test_database_clients.py::TestUniProtClient::test_get_protein_info PASSED
tests/test_database_clients.py::test_complete_data_retrieval_workflow PASSED
tests/test_structure_validation.py::TestStructureDownload::test_download_adar2_alphafold PASSED
tests/test_structure_validation.py::TestStructureDownload::test_download_adar2_crystal PASSED
tests/test_structure_validation.py::TestStructureParsing::test_structure_has_residues PASSED
tests/test_structure_validation.py::TestStructureParsing::test_structure_has_atoms PASSED
tests/test_structure_validation.py::TestStructureParsing::test_plddt_scores_exist PASSED
tests/test_structure_validation.py::TestStructureParsing::test_find_basic_residues PASSED
tests/test_structure_validation.py::test_api_connectivity PASSED

======================== 15 passed in 45.23s =========================
```

## Troubleshooting

### Network Errors

If tests fail with network errors:

```
requests.exceptions.ConnectionError: Failed to establish connection
```

**Solutions:**
1. Check internet connection
2. Try again (APIs may be temporarily down)
3. Check if APIs are blocked by firewall
4. Skip network tests: `pytest tests/ -m "not requires_network"`

### Missing Dependencies

If tests fail with import errors:

```
ModuleNotFoundError: No module named 'Bio'
```

**Solution:**
```bash
pip install -r requirements.txt
# or
pip install biopython requests pandas matplotlib
```

### Slow Tests

Tests download real data and may be slow:
- First run: ~1-2 minutes (downloads structures)
- Subsequent runs: ~10-20 seconds (uses cache)

Speed up with parallel execution:
```bash
pip install pytest-xdist
pytest tests/ -n auto  # Use all CPU cores
```

## Data Sources

All tests use **real data** from:

| Source | URL | Purpose |
|--------|-----|----------|
| **AlphaFold EBI** | `https://alphafold.ebi.ac.uk` | Predicted protein structures |
| **RCSB PDB** | `https://www.rcsb.org` | Experimental crystal structures |
| **UniProt** | `https://www.uniprot.org` | Protein annotations |

### Test Proteins

| Protein | UniProt | PDB | Description |
|---------|---------|-----|-------------|
| **ADAR2** | P78563 | 1ZY7 | Gold standard IP6-binding protein |
| **IP3K** | Q13572 | 1W2C | Positive control |
| **BTK** | Q06187 | 1MAI | Known IP3-binding kinase |

## CI/CD Pipeline

The GitHub Actions pipeline:

1. **Install dependencies** (fpocket, FreeSASA, APBS)
2. **Download test structures** (ADAR2, 1ZY7)
3. **Run linting** (black, flake8, isort, mypy)
4. **Run tests** (pytest with coverage)
5. **Execute notebooks** (verify they run)
6. **Security scan** (bandit, safety)
7. **Build package** (verify installation works)

See workflow status: [GitHub Actions](https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites/actions)

## Questions?

If tests fail unexpectedly:

1. Check [GitHub Issues](https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites/issues)
2. Verify API connectivity: `python tests/test_structure_validation.py::test_api_connectivity`
3. Open an issue with error logs

---

**All tests use real data - no mocks, no placeholders, no simplifications.**
