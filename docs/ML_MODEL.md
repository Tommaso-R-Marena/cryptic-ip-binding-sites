# ML Classifier Training and Usage

## Train the model

Run end-to-end training (dataset build + pocket feature extraction + model selection):

```bash
python scripts/train_ml_classifier.py
```

The training script performs:

- Validation-set refresh via `scripts/build_ip_validation_dataset.py` (unless `--skip-build-dataset`).
- Pocket-level feature extraction from all pockets detected in the validation structures.
- Hyperparameter grid search for **Random Forest** and **XGBoost**.
- Stratified 5-fold CV and confidence intervals for ROC/PR curves.
- Held-out test benchmark versus threshold scoring.
- Model serialization to `models/cryptic_ip_classifier_v1.pkl` (local artifact; not committed).
- SHAP feature-importance export (`models/cryptic_ip_classifier_v1_shap.png`, local artifact; not committed).
- Model card and metadata export under `models/`.

## Retraining workflow when dataset expands

```bash
python scripts/retrain_ml_classifier_if_needed.py \
  --dataset-csv data/validation/ip_binding_validation_dataset.csv
```

This command compares the dataset hash against `models/cryptic_ip_classifier_v1.metadata.json` and only retrains when needed (unless `--force`).

## Analyzer/CLI integration

Use ML scoring in CLI:

```bash
cryptic-ip analyze path/to/protein.pdb --use-ml-model
```

Or pass an explicit model path:

```bash
cryptic-ip analyze path/to/protein.pdb --use-ml-model --model-path models/cryptic_ip_classifier_v1.pkl
```

If the model is unavailable or fails to load, the analyzer falls back to threshold-based scoring automatically.

## Benchmark outputs for manuscript

After training, comparison outputs are written to `results/ml_training/`:

- `ml_vs_threshold_comparison.csv`
- `ml_vs_threshold_comparison.md`

These files summarize ROC AUC and PR AUC for ML vs threshold scoring on a held-out test split.
