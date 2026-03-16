# Model Registry

This directory stores versioned ML artifacts used by `ProteinAnalyzer` and CLI commands.

## Current versions

| Version | File | Status | Notes |
|---|---|---|---|
| v1 | `cryptic_ip_classifier_v1.pkl` | generated locally | Tree-based classifier selected from Random Forest/XGBoost grid search. Binary artifact is not committed. |

## Associated artifacts (v1)

- `cryptic_ip_classifier_v1.metadata.json`: training metadata (UTC timestamp, dataset reference/hash, hyperparameters, metrics).
- `cryptic_ip_classifier_v1_model_card.md`: concise model card for manuscript/reporting.
- `cryptic_ip_classifier_v1_shap.png`: SHAP global feature-importance bar chart (generated locally, not committed).

## Updating model versions

1. Refresh validation set: `python scripts/build_ip_validation_dataset.py ...`
2. Retrain when dataset changes:
   ```bash
   python scripts/retrain_ml_classifier_if_needed.py \
     --dataset-csv data/validation/ip_binding_validation_dataset.csv
   ```
3. If creating a new semantic model version (e.g., `v2`), write new files instead of overwriting `v1`, then append this table.


> Note: Binary artifacts (`.pkl`, `.png`) are intentionally excluded from version control. Generate them by running `python scripts/train_ml_classifier.py`.
