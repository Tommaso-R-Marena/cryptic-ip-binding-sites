# Model Registry

Versioned ML artifacts used by `ProteinAnalyzer` and the CLI.

## Current versions

| Version | File | Status | Notes |
|---|---|---|---|
| v2 | `cryptic_ip_classifier_v2.pkl` | **produced by the current pipeline** | Nested, grouped, calibrated cross-validation over five model families on the 40-descriptor schema. Generate with `scripts/train_ml_classifier.py`. |
| v1 | `cryptic_ip_classifier_v1.pkl` | **superseded — do not use** | Trained under a defective labelling scheme; see below. |

## Why v1 must not be used

The v1 metadata records what the artifact actually achieved:

```
roc_auc      0.4964      (chance is 0.5)
pr_auc       0.000208
mcc          0.0
num_pockets  12190
num_positive 5
```

That is not underperformance, it is a broken pipeline, and the cause was in the
label definition rather than the model:

1. Burial class was assigned from ligand SASA **summed over every copy** of a
   ligand in an entry. Copy count is a crystallisation artefact, so nearly every
   entry was classified `Surface`.
2. Every pocket in a `Surface` entry was labelled negative — including the actual
   inositol phosphate sites. Most of the positive evidence in the dataset was
   therefore labelled as the negative class.
3. Positives were then limited to pockets within 8 Å of a ligand *centroid* in
   the few surviving entries, leaving 5 positives among 12 190 pockets.

`ProteinAnalyzer` now reads the recorded validation AUROC before deploying a
model and refuses to use one that performed at chance, falling back to the
rule-based scorer with an explanatory message. v1 is kept as a record of the
defect, not as a usable artifact.

## Associated artifacts

- `<name>.metadata.json` — training metadata: UTC timestamp, command line, seeds,
  software versions, per-model nested-CV results, chosen decision threshold,
  label counts, and the rule-based baseline comparison.
- `<name>_model_card.md` — intended use, out-of-scope uses, data description,
  evaluation protocol, performance with confidence intervals, feature importance
  and known limitations.

## Producing a model

```bash
python scripts/build_ip_validation_dataset.py --n-decoys 500   # needs RCSB access
python scripts/extract_pocket_features.py --jobs 8
python scripts/train_ml_classifier.py
```

Every stage is resumable. See [docs/METHODS.md](../docs/METHODS.md) for the
evaluation protocol and [docs/DATA_ACQUISITION.md](../docs/DATA_ACQUISITION.md)
for the collection design.

## Verifying without network access

```bash
cryptic-ip self-check
pytest tests/integration/test_synthetic_end_to_end.py
```

These exercise the full path on synthetic structures with known ground truth.

> Binary artifacts (`.pkl`, `.png`) are intentionally excluded from version
> control; regenerate them with the commands above.
