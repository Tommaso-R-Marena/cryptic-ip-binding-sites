# Cryptic IP Classifier v1

- Trained: 2026-06-19T17:51:25.171931+00:00
- Dataset: `data/validation/ip_binding_validation_dataset.csv`
- Structures: 136
- Pockets: 12190
- Selected model: random_forest
- CV ROC AUC: 0.565
- CV PR AUC: 0.023

## Hyperparameters
```json
{
  "classifier__class_weight": null,
  "classifier__max_depth": null,
  "classifier__min_samples_leaf": 1,
  "classifier__min_samples_split": 4,
  "classifier__n_estimators": 200
}
```

## Held-out benchmark
| method | test_roc_auc | test_pr_auc |
|---|---|---|
| ML (random_forest) | 0.5814401814401815 | 0.026766709689335984 |
| Threshold scoring | 0.4470925470925471 | 0.013889916435016008 |
