# Publication Package

Run the end-to-end manuscript preparation workflow:

```bash
# Requires fpocket (conda-forge), Python deps, and network access to RCSB/AlphaFold
python scripts/run_publication_package.py --output-dir results/publication
```

## What it produces

| Output | Description |
|--------|-------------|
| `data/validation/ip_binding_validation_dataset.csv` | RCSB-derived IP-binding structures with SASA-based burial labels |
| `results/publication/validation/control_benchmark.csv` | Positive/negative structural control scores |
| `results/publication/ml_training/` | Pocket features, ROC/PR curves, ML vs threshold comparison |
| `results/publication/figures/` | Figure 1–3 PDF/PNG panels + `figure_legends.txt` |
| `results/publication/RESULTS_SUMMARY.md` | Manuscript-ready results narrative |
| `results/publication/METHODS_AUTO.md` | Auto-generated methods text |
| `results/publication/provenance_manifest.jsonld` | Reproducibility manifest |

## Partial reruns

```bash
# Skip dataset rebuild if CSV already exists
python scripts/run_publication_package.py --skip-dataset-build

# Enable APBS electrostatics in controls and ML feature extraction
python scripts/run_publication_package.py --with-electrostatics

# Figures only (after CSV assets exist)
python scripts/generate_publication_figures.py \
  --config results/publication/figure_config.yaml --figures all
```

## Yeast pilot screen (500 proteins)

```bash
python scripts/run_yeast_pilot_screen.py --n-proteins 500 --workers 4
```

## MD pilot validation

```bash
python scripts/run_md_pilot_validation.py --top-n 5 --production-ns 1.0
```

## ML training (pocket-level labels, grouped splits)

```bash
python scripts/train_ml_classifier.py --skip-build-dataset --include-electrostatics
```

## External dependencies

Install structural biology tools via conda:

```bash
conda install -c conda-forge fpocket freesasa apbs pdb2pqr
```
