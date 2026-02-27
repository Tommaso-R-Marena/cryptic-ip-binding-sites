# Inositol Phosphate Validation Dataset Generation

Use `scripts/build_ip_validation_dataset.py` to build a reproducible validation set of
inositol-phosphate-bound crystal structures from RCSB PDB.

## What the script does

1. Queries RCSB for **X-ray diffraction** entries containing `IP3`, `IP4`, `IP5`, or `IP6`.
2. Downloads each `.pdb` file and saves ligand metadata JSON from the RCSB ChemComp API.
3. Computes ligand solvent-accessible surface area (SASA) with **FreeSASA**.
4. Classifies ligand burial:
   - `Cryptic`: SASA < 5 Å²
   - `Semi-cryptic`: SASA 5–20 Å²
   - `Surface`: SASA > 20 Å²
5. Writes a CSV with columns:
   `pdb_id, uniprot_id, ligand_type, sasa, classification, resolution, organism`.

## Reproducibility notes

- API calls are retried with exponential backoff and `429` rate-limit handling.
- Downloaded structures and ligand metadata are cached in a user-selected output directory.
- CSV rows are sorted by `(pdb_id, ligand_type)` for deterministic output order.
- A configurable minimum unique-protein check (`--min-proteins`, default `20`) is logged as a warning.

## Usage

```bash
python scripts/build_ip_validation_dataset.py \
  --output-csv data/validation/ip_binding_validation_dataset.csv \
  --download-dir data/validation/raw \
  --min-proteins 20 \
  --log-level INFO
```

## Dependencies

- `requests`
- `freesasa` Python bindings (and FreeSASA library)

If `freesasa` is unavailable, the script exits with a clear error when SASA computation is attempted.
