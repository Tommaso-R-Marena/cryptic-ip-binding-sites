# Cryptic IP Classifier v1

- Trained: 2026-02-28T20:21:59.841672+00:00
- Dataset: `data/validation/ip_binding_validation_dataset.csv`
- Structures: 24
- Pockets: 24
- Selected model: random_forest
- CV ROC AUC: 1.000
- CV PR AUC: 1.000

## Hyperparameters
```json
{
  "classifier__class_weight": null,
  "classifier__max_depth": null,
  "classifier__min_samples_leaf": 1,
  "classifier__min_samples_split": 2,
  "classifier__n_estimators": 200
}
```

## Held-out benchmark
| method | test_roc_auc | test_pr_auc |
|---|---|---|
| ML (random_forest) | 1.0 | 1.0 |
| Threshold scoring | 1.0 | 1.0 |
