#!/usr/bin/env python3
"""Train, compare and select cryptic IP-site classifiers with nested grouped CV.

Protocol
--------
1. Load the pocket feature table produced by ``extract_pocket_features.py``.
2. Drop ambiguous pockets (label ``-1``) - they are excluded, not relabelled.
3. Train every candidate model under **nested, grouped, stratified**
   cross-validation on identical splits.
4. Report discrimination (AUROC, average precision), calibration (Brier, ECE) and
   screening utility (enrichment factor, precision at *k*), each with
   group-level bootstrap confidence intervals.
5. Compare models with a paired DeLong test, and compare the best model with the
   interpretable rule-based scorer the same way.
6. Persist the selected model, a model card, and a full metadata record.

The rule-based comparison is not decoration. If a transparent weighted score
matches the learned model, the learned model adds complexity without adding
information, and the honest report says so.
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.analysis.features import FEATURE_NAMES, feature_documentation  # noqa: E402
from cryptic_ip.analysis.ml_classifier import (  # noqa: E402
    CrypticSiteMLClassifier,
    classification_metrics,
    compare_models,
    default_model_specs,
    delong_roc_test,
    grouped_bootstrap_ci,
    model_comparison_table,
    pairwise_delong,
    select_threshold,
)
from cryptic_ip.analysis.scorer import PocketScorer  # noqa: E402

LOGGER = logging.getLogger("train_ml_classifier")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--features-csv",
        type=Path,
        default=Path("results/ml_training/pocket_features.csv"),
        help="Pocket feature table from extract_pocket_features.py",
    )
    parser.add_argument(
        "--work-dir",
        type=Path,
        default=Path("results/ml_training"),
        help="Directory for reports and curves",
    )
    parser.add_argument(
        "--model-dir", type=Path, default=Path("models"), help="Directory for model artefacts"
    )
    parser.add_argument(
        "--model-name",
        default="cryptic_ip_classifier_v2",
        help="Base name for the serialised model and its model card",
    )
    parser.add_argument(
        "--models",
        nargs="*",
        default=None,
        help="Candidate model names; defaults to the full comparison set",
    )
    parser.add_argument(
        "--group-column",
        default="group_key",
        help="Column identifying the grouping unit that splits must respect",
    )
    parser.add_argument("--label-column", default="label", help="Label column name")
    parser.add_argument("--n-splits", type=int, default=5, help="Outer CV folds")
    parser.add_argument("--inner-splits", type=int, default=3, help="Inner CV folds")
    parser.add_argument(
        "--n-search-iter", type=int, default=40, help="Randomised-search budget per inner loop"
    )
    parser.add_argument(
        "--n-bootstrap", type=int, default=2000, help="Group-bootstrap resamples for CIs"
    )
    parser.add_argument("--random-state", type=int, default=42, help="Master random seed")
    parser.add_argument(
        "--threshold-objective",
        default="mcc",
        choices=["mcc", "f1", "fbeta", "youden"],
        help="Objective used to pick the decision threshold on out-of-fold predictions",
    )
    parser.add_argument(
        "--no-figures", action="store_true", help="Skip ROC/PR/calibration figures"
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return parser.parse_args(argv)


def dataframe_to_markdown(frame: pd.DataFrame, *, float_format: str = "{:.4f}") -> str:
    """Render a DataFrame as a GitHub-flavoured Markdown table."""
    formatted = frame.copy()
    for column in formatted.columns:
        if pd.api.types.is_float_dtype(formatted[column]):
            formatted[column] = formatted[column].map(
                lambda v: "" if pd.isna(v) else float_format.format(v)
            )
    header = "| " + " | ".join(str(c) for c in formatted.columns) + " |"
    separator = "|" + "|".join(["---"] * len(formatted.columns)) + "|"
    rows = [
        "| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
        for row in formatted.astype(object).values
    ]
    return "\n".join([header, separator, *rows])


def load_training_frame(
    path: Path, *, label_column: str, group_column: str
) -> pd.DataFrame:
    """Load and validate the pocket feature table.

    Args:
        path: Feature CSV path.
        label_column: Label column name.
        group_column: Grouping column name.

    Returns:
        Rows with a trainable label (0 or 1).

    Raises:
        FileNotFoundError: If the table does not exist.
        ValueError: If required columns or both classes are missing.
    """
    if not path.exists():
        raise FileNotFoundError(
            f"Feature table not found: {path}. Run scripts/extract_pocket_features.py first."
        )
    frame = pd.read_csv(path)
    for column in (label_column, group_column):
        if column not in frame.columns:
            raise ValueError(f"Required column {column!r} missing from {path}")

    missing_features = [name for name in FEATURE_NAMES if name not in frame.columns]
    if missing_features:
        raise ValueError(
            f"Feature table is missing {len(missing_features)} descriptors: {missing_features[:8]}"
        )

    n_total = len(frame)
    ambiguous = int((frame[label_column] == -1).sum())
    frame = frame[frame[label_column].isin([0, 1])].copy()
    LOGGER.info(
        "Loaded %d pockets; %d ambiguous excluded; %d retained (%d positive)",
        n_total,
        ambiguous,
        len(frame),
        int((frame[label_column] == 1).sum()),
    )

    # A descriptor that is missing for every row carries no information and makes
    # the imputer warn on every fit. Report it once, clearly: an all-missing
    # column usually means an optional tool (APBS) was never run, and the reader
    # should know the model had no access to that evidence.
    all_missing = [
        name for name in FEATURE_NAMES if name in frame.columns and frame[name].isna().all()
    ]
    if all_missing:
        LOGGER.warning(
            "%d descriptor(s) are missing for every pocket and contribute nothing: %s. "
            "For 'electrostatic_potential' this means APBS was not run; the screened "
            "Coulomb surrogate 'coulomb_potential_kt' stands in for it.",
            len(all_missing),
            ", ".join(all_missing),
        )

    counts = frame[label_column].value_counts().to_dict()
    if len(counts) < 2:
        raise ValueError(
            f"Both classes are required for training; label counts are {counts}. "
            "Check that pocket labelling matched at least one ligand site."
        )
    return frame


def rule_baseline_scores(frame: pd.DataFrame) -> np.ndarray:
    """Score every pocket with the interpretable rule-based scorer."""
    return PocketScorer().score_frame(frame)


def plot_diagnostics(
    results: Dict[str, Any], work_dir: Path, *, best_name: str
) -> List[Path]:
    """Write ROC, precision-recall and calibration figures from out-of-fold data.

    Curves are drawn from out-of-fold predictions, so they describe generalisation
    rather than fit.

    Args:
        results: Training results keyed by model name.
        work_dir: Output directory.
        best_name: Name of the selected model.

    Returns:
        Paths of the written figures.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.calibration import calibration_curve
    from sklearn.metrics import precision_recall_curve, roc_curve

    work_dir.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    for name, result in results.items():
        if not result.oof_probabilities.size:
            continue
        y, p = result.oof_labels, result.oof_probabilities
        style = {"lw": 2.2 if name == best_name else 1.2, "alpha": 1.0 if name == best_name else 0.7}

        fpr, tpr, _ = roc_curve(y, p)
        axes[0].plot(fpr, tpr, label=f"{name} ({result.metrics.get('roc_auc', float('nan')):.3f})", **style)

        precision, recall, _ = precision_recall_curve(y, p)
        axes[1].plot(recall, precision, label=f"{name} ({result.metrics.get('pr_auc', float('nan')):.3f})", **style)

        try:
            observed, predicted = calibration_curve(y, p, n_bins=10, strategy="quantile")
            axes[2].plot(predicted, observed, marker="o", ms=3, label=name, **style)
        except ValueError:
            pass

    axes[0].plot([0, 1], [0, 1], "k--", lw=0.8)
    axes[0].set(xlabel="False positive rate", ylabel="True positive rate", title="ROC (out-of-fold)")
    axes[0].legend(fontsize=7, loc="lower right")

    if results:
        first = next(iter(results.values()))
        if first.oof_labels.size:
            axes[1].axhline(float(np.mean(first.oof_labels)), color="k", ls="--", lw=0.8)
    axes[1].set(xlabel="Recall", ylabel="Precision", title="Precision-recall (out-of-fold)")
    axes[1].legend(fontsize=7, loc="upper right")

    axes[2].plot([0, 1], [0, 1], "k--", lw=0.8)
    axes[2].set(
        xlabel="Predicted probability", ylabel="Observed frequency", title="Calibration"
    )
    axes[2].legend(fontsize=7, loc="upper left")

    figure.tight_layout()
    out_path = work_dir / "model_diagnostics.png"
    figure.savefig(out_path, dpi=200)
    plt.close(figure)
    return [out_path]


def write_model_card(
    path: Path,
    *,
    model_name: str,
    best_name: str,
    results: Dict[str, Any],
    comparison: pd.DataFrame,
    baseline_metrics: Dict[str, float],
    delong_table: pd.DataFrame,
    importance: Optional[pd.DataFrame],
    dataset_info: Dict[str, Any],
) -> None:
    """Write a model card documenting intended use, data, protocol and limits."""
    best = results[best_name]
    roc_point, roc_low, roc_high = best.metric_cis.get("roc_auc", (np.nan, np.nan, np.nan))
    pr_point, pr_low, pr_high = best.metric_cis.get("pr_auc", (np.nan, np.nan, np.nan))

    lines: List[str] = [
        f"# Model card: {model_name}",
        "",
        f"*Generated {datetime.now(timezone.utc).isoformat()}*",
        "",
        "## Intended use",
        "",
        "Ranks candidate pockets by their resemblance to buried inositol phosphate",
        "binding sites, to prioritise structures for manual inspection and",
        "experimental follow-up. It is a triage tool, not evidence of binding: a",
        "high-ranking pocket is a hypothesis to test, not a validated site.",
        "",
        "## Out of scope",
        "",
        "- Predicting binding affinity or occupancy.",
        "- Distinguishing inositol phosphates from other polyanions (nucleotides,",
        "  sulfate clusters, polyphosphates), which the descriptors do not separate.",
        "- Scoring structures whose fold is unreliable; check pLDDT for predicted models.",
        "",
        "## Training data",
        "",
        f"- Pockets: {best.n_positive + best.n_negative} "
        f"({best.n_positive} positive, {best.n_negative} negative)",
        f"- Groups (proteins): {best.n_groups}",
        f"- Ambiguous pockets excluded: {dataset_info.get('n_ambiguous', 'n/a')}",
        f"- Source table: `{dataset_info.get('features_csv', 'n/a')}`",
        "- Label rule: atom-level ligand overlap "
        "(positive at >= 30 %, negative at <= 5 %, ambiguous between)",
        "",
        "## Evaluation protocol",
        "",
        f"- Nested cross-validation: {len(best.folds)} outer folds, "
        f"{dataset_info.get('inner_splits', 'n/a')} inner folds",
        "- Splits are stratified **and grouped by protein**, so no protein appears",
        "  in both training and evaluation of the same fold.",
        "- Hyperparameters are selected in the inner loop only; every reported",
        "  number comes from outer folds the search never saw.",
        "- Probabilities are calibrated on a group-disjoint slice of each",
        "  training fold.",
        "- Confidence intervals come from bootstrap resampling of **groups**.",
        "",
        "## Performance (out-of-fold)",
        "",
        f"- AUROC: {roc_point:.3f} (95 % CI {roc_low:.3f}-{roc_high:.3f})",
        f"- Average precision: {pr_point:.3f} (95 % CI {pr_low:.3f}-{pr_high:.3f})",
        f"- Positive prevalence: {best.metrics.get('prevalence', float('nan')):.4f}",
        f"- MCC at the selected threshold ({best.selected_threshold:.3f}): "
        f"{best.metrics.get('mcc', float('nan')):.3f}",
        f"- Brier score: {best.metrics.get('brier', float('nan')):.4f}; "
        f"expected calibration error: {best.metrics.get('ece', float('nan')):.4f}",
        f"- Enrichment at 1 %: {best.metrics.get('enrichment_at_1pct', float('nan')):.2f}x",
        "",
        "### Model comparison",
        "",
        dataframe_to_markdown(comparison),
        "",
        "### Rule-based baseline",
        "",
        f"- AUROC {baseline_metrics.get('roc_auc', float('nan')):.3f}, "
        f"average precision {baseline_metrics.get('pr_auc', float('nan')):.3f}",
        "- The learned model must beat this to justify its complexity.",
        "",
    ]

    if not delong_table.empty:
        lines += ["### Paired DeLong tests", "", dataframe_to_markdown(delong_table), ""]

    if importance is not None and not importance.empty:
        lines += [
            "## Feature importance (permutation, average precision)",
            "",
            dataframe_to_markdown(importance.head(15)),
            "",
            "Correlated descriptors share credit, so these values rank contributions",
            "to the metric; they are not causal statements.",
            "",
        ]

    lines += [
        "## Known limitations",
        "",
        "- Pocket-detector recall bounds everything: a site fpocket never proposes",
        "  cannot be scored. See `labeling_summary.json` for the measured recall.",
        "- Positives come from solved structures, which are biased toward proteins",
        "  that crystallise well and toward ligands that were modelled.",
        "- The APBS electrostatic descriptor is often missing; the screened Coulomb",
        "  surrogate stands in and is a continuum approximation, not a",
        "  Poisson-Boltzmann solution.",
        "- pLDDT descriptors are meaningful only for predicted models; for",
        "  experimental structures the same column holds crystallographic B-factors.",
        "",
        "## Reproduction",
        "",
        "```bash",
        "python scripts/build_ip_validation_dataset.py",
        "python scripts/extract_pocket_features.py",
        "python scripts/train_ml_classifier.py",
        "```",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.model_dir.mkdir(parents=True, exist_ok=True)

    raw = pd.read_csv(args.features_csv) if args.features_csv.exists() else None
    n_ambiguous = int((raw[args.label_column] == -1).sum()) if raw is not None else 0

    frame = load_training_frame(
        args.features_csv, label_column=args.label_column, group_column=args.group_column
    )
    X = frame.loc[:, list(FEATURE_NAMES)]
    y = frame[args.label_column].astype(int).to_numpy()
    groups = frame[args.group_column].astype(str).tolist()

    model_names = list(args.models) if args.models else [spec.name for spec in default_model_specs()]
    LOGGER.info("Comparing models: %s", ", ".join(model_names))

    results = compare_models(
        X,
        y,
        groups,
        model_names=model_names,
        random_state=args.random_state,
        n_splits=args.n_splits,
        inner_splits=args.inner_splits,
        n_search_iter=args.n_search_iter,
        n_bootstrap=args.n_bootstrap,
    )
    if not results:
        LOGGER.error("No model trained successfully.")
        return 3

    comparison = model_comparison_table(results)
    comparison.to_csv(args.work_dir / "model_comparison.csv", index=False)
    (args.work_dir / "model_comparison.md").write_text(
        dataframe_to_markdown(comparison), encoding="utf-8"
    )

    best_name = str(comparison.iloc[0]["model"])
    best = results[best_name]
    LOGGER.info(
        "Selected %s: AUROC %.3f, AP %.3f",
        best_name,
        best.metrics.get("roc_auc", float("nan")),
        best.metrics.get("pr_auc", float("nan")),
    )

    # Rule-based baseline evaluated on the same out-of-fold rows.
    baseline_all = rule_baseline_scores(frame)
    baseline_threshold = select_threshold(y, baseline_all, objective=args.threshold_objective)
    baseline_metrics = classification_metrics(y, baseline_all, threshold=baseline_threshold)
    from sklearn.metrics import average_precision_score, roc_auc_score

    for metric_name, metric_fn in (
        ("roc_auc", roc_auc_score),
        ("pr_auc", average_precision_score),
    ):
        point, low, high = grouped_bootstrap_ci(
            y,
            baseline_all,
            groups,
            metric=lambda a, b, fn=metric_fn: float(fn(a, b)),
            n_bootstrap=args.n_bootstrap,
            random_state=args.random_state,
        )
        baseline_metrics[f"{metric_name}_ci_low"] = low
        baseline_metrics[f"{metric_name}_ci_high"] = high

    ml_vs_rule_difference, ml_vs_rule_p = delong_roc_test(
        best.oof_labels, best.oof_probabilities, baseline_all[: best.oof_labels.size]
    )
    comparison_rows = pd.DataFrame(
        [
            {
                "method": f"ML ({best_name})",
                "roc_auc": best.metrics.get("roc_auc", np.nan),
                "pr_auc": best.metrics.get("pr_auc", np.nan),
                "mcc": best.metrics.get("mcc", np.nan),
                "enrichment_at_1pct": best.metrics.get("enrichment_at_1pct", np.nan),
                "evaluation": "nested CV, out-of-fold",
            },
            {
                "method": "Rule-based scorer",
                "roc_auc": baseline_metrics.get("roc_auc", np.nan),
                "pr_auc": baseline_metrics.get("pr_auc", np.nan),
                "mcc": baseline_metrics.get("mcc", np.nan),
                "enrichment_at_1pct": baseline_metrics.get("enrichment_at_1pct", np.nan),
                "evaluation": "no fitting; applied directly",
            },
        ]
    )
    comparison_rows.to_csv(args.work_dir / "ml_vs_threshold_comparison.csv", index=False)
    (args.work_dir / "ml_vs_threshold_comparison.md").write_text(
        dataframe_to_markdown(comparison_rows)
        + f"\n\nDeLong AUROC difference {ml_vs_rule_difference:+.4f}, p = {ml_vs_rule_p:.3g}\n",
        encoding="utf-8",
    )

    delong_table = pairwise_delong(results)
    if not delong_table.empty:
        delong_table.to_csv(args.work_dir / "model_delong_tests.csv", index=False)

    estimator: CrypticSiteMLClassifier = getattr(best, "estimator")
    model_path = args.model_dir / f"{args.model_name}.pkl"
    estimator.save(str(model_path))
    LOGGER.info("Saved model to %s", model_path)

    importance: Optional[pd.DataFrame] = None
    try:
        importance = estimator.permutation_importance_(X, y, n_repeats=10)
        importance.to_csv(args.work_dir / "feature_importance.csv", index=False)
    except Exception as exc:  # noqa: BLE001 - importance is diagnostic, not essential
        LOGGER.warning("Permutation importance failed: %s", exc)

    if not args.no_figures:
        try:
            plot_diagnostics(results, args.work_dir, best_name=best_name)
        except Exception as exc:  # noqa: BLE001 - figures must not block training
            LOGGER.warning("Diagnostic figures failed: %s", exc)

    dataset_info = {
        "features_csv": str(args.features_csv),
        "n_ambiguous": n_ambiguous,
        "inner_splits": args.inner_splits,
    }
    write_model_card(
        args.model_dir / f"{args.model_name}_model_card.md",
        model_name=args.model_name,
        best_name=best_name,
        results=results,
        comparison=comparison,
        baseline_metrics=baseline_metrics,
        delong_table=delong_table,
        importance=importance,
        dataset_info=dataset_info,
    )

    metadata = {
        "model_version": args.model_name,
        "selected_model_type": best_name,
        "trained_at_utc": datetime.now(timezone.utc).isoformat(),
        "command": " ".join(sys.argv),
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "labeling": "atom-level ligand overlap; ambiguous pockets excluded",
        "split": f"nested StratifiedGroupKFold grouped by {args.group_column}",
        "random_state": args.random_state,
        "feature_columns": list(FEATURE_NAMES),
        "feature_documentation": feature_documentation(),
        "decision_threshold": best.selected_threshold,
        "n_pockets": int(len(frame)),
        "n_positive": int(best.n_positive),
        "n_negative": int(best.n_negative),
        "n_groups": int(best.n_groups),
        "n_ambiguous_excluded": n_ambiguous,
        "model_candidates": {name: result.to_dict() for name, result in results.items()},
        "rule_based_baseline": baseline_metrics,
        "ml_vs_rule_delong": {
            "auc_difference": ml_vs_rule_difference,
            "p_value": ml_vs_rule_p,
        },
    }
    try:
        from sklearn import __version__ as sklearn_version

        metadata["software"]["scikit-learn"] = sklearn_version
    except ImportError:
        pass

    (args.model_dir / f"{args.model_name}.metadata.json").write_text(
        json.dumps(metadata, indent=2, default=float), encoding="utf-8"
    )
    LOGGER.info("Wrote model card and metadata to %s", args.model_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
