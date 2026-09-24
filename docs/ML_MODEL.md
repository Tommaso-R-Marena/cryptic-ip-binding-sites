# ML Classifier Training and Usage

## Pipeline

Training has three separate stages, each resumable and independently inspectable:

```bash
# 1. Ground truth from the RCSB (requires network access)
python scripts/build_ip_validation_dataset.py --n-decoys 500 --jobs 8

# 2. Pocket descriptors and labels
python scripts/extract_pocket_features.py --jobs 8

# 3. Model comparison and selection
python scripts/train_ml_classifier.py
```

Separating them matters in practice: collection is slow and network-bound,
feature extraction is CPU-bound, and training is fast enough to iterate on. A
change to the model no longer forces a re-download of the PDB.

## What training does

- Loads the pocket table and **drops ambiguous pockets** (label `-1`) rather than
  relabelling them as negatives.
- Trains five model families on identical grouped splits: logistic regression
  (the interpretable floor), random forest, extremely randomised trees,
  histogram gradient boosting, and XGBoost when installed.
- Runs **nested** cross-validation: hyperparameters are searched in an inner
  grouped loop, and every reported number comes from outer folds the search
  never saw.
- Groups splits by protein (UniProt accession, falling back to PDB ID), so no
  protein appears on both sides of a fold.
- Calibrates probabilities on a group-disjoint slice of each training fold, and
  reports calibration quality (Brier score, expected calibration error).
- Selects the decision threshold on out-of-fold predictions only.
- Compares models with a paired DeLong test, and compares the winner against the
  rule-based scorer the same way.
- Writes a model card, full metadata, feature importances and diagnostic figures.

## Outputs

| Path | Contents |
|---|---|
| `models/<name>.pkl` | Serialised model with its feature schema and threshold |
| `models/<name>.metadata.json` | Seeds, versions, per-model nested-CV results, baseline comparison |
| `models/<name>_model_card.md` | Intended use, protocol, performance, limitations |
| `results/ml_training/model_comparison.csv` | Metrics with confidence intervals per model |
| `results/ml_training/model_delong_tests.csv` | Paired AUROC comparisons |
| `results/ml_training/ml_vs_threshold_comparison.csv` | Learned model vs rule-based scorer |
| `results/ml_training/feature_importance.csv` | Permutation importance |
| `results/ml_training/model_diagnostics.png` | ROC, precision-recall and calibration curves |

## Reported metrics

Discrimination is not enough on its own for a screening tool, so three families
are reported:

| Family | Metrics | Why |
|---|---|---|
| Discrimination | AUROC, average precision, MCC, F1, balanced accuracy | Ranking quality; average precision is the primary metric because positives are rare |
| Calibration | Brier score, expected calibration error | Downstream filters apply probability thresholds, which is meaningless if probabilities are not calibrated |
| Screening utility | Enrichment factor and precision at top 1 %, 5 %, 10 % | A proteome screen inspects the head of a ranked list |

Confidence intervals come from bootstrap resampling of **protein groups**, not
individual pockets: several pockets from one protein are not independent
observations, and resampling rows understates uncertainty.

## Using a trained model

```python
from cryptic_ip.analysis import ProteinAnalyzer

analyzer = ProteinAnalyzer(
    "structure.pdb", use_ml_model=True, model_path="models/cryptic_ip_classifier_v2.pkl"
)
hits = analyzer.run_pipeline().sort_values("composite_score", ascending=False)
```

`ProteinAnalyzer` reads the recorded validation AUROC before deploying a model
and **refuses to use one that performed at chance**, falling back to the
rule-based scorer with an explanation. Deploying a model with no demonstrated
skill is the worst available outcome: it produces confident-looking numbers with
nothing behind them.

## Interpreting a rule-based win

The rule-based scorer is included in every comparison, and the DeLong p-value for
learned-vs-rule is recorded. If the learned model does not beat the transparent
score significantly, the honest conclusion is that the extra complexity is not
buying anything — and the report should say so rather than presenting the
learned model as an advance.

## Backwards compatibility

The legacy six-feature schema is still supported for models serialised by earlier
versions:

```python
from cryptic_ip.analysis.ml_classifier import FEATURE_COLUMNS, CrypticSiteMLClassifier

clf = CrypticSiteMLClassifier(feature_names=FEATURE_COLUMNS)
```

One correction applies when scoring through the legacy schema: `pocket_depth` is
now filled from the genuine geometric burial depth. Earlier versions passed
fpocket's mean local hydrophobic density in that column, so a legacy model is now
fed the quantity its column name always claimed.

See also: [METHODS.md](METHODS.md), [DATA_ACQUISITION.md](DATA_ACQUISITION.md).
