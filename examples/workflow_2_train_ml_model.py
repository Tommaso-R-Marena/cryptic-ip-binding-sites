#!/usr/bin/env python3
"""Workflow 2: build validation dataset -> train classifier -> export model.

Expected outputs
- `ip_validation_dataset.csv`
- `cryptic_ip_rf.joblib`
- training metrics JSON

Typical runtime
- Dataset build: 5-20 minutes (network + FreeSASA)
- Model training: <2 minutes for small/medium datasets

Troubleshooting
- `freesasa is required`: install FreeSASA bindings.
- Sparse class labels: ensure both positive and negative classes exist.
- CV error about folds: provide at least `n_splits` samples per class.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cryptic_ip.analysis.ml_classifier import CrypticSiteMLClassifier
from scripts.build_ip_validation_dataset import build_dataset


CLASS_MAP = {"Cryptic": 1, "Semi-cryptic": 1, "Surface": 0}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Train cryptic IP ML model from validation dataset.")
    parser.add_argument("--output-dir", type=Path, default=Path("results/workflow2"))
    parser.add_argument("--build-dataset", action="store_true", help="Build dataset before training")
    parser.add_argument("--dataset-csv", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dataset_csv = args.dataset_csv or (args.output_dir / "ip_validation_dataset.csv")
    if args.build_dataset:
        build_dataset(
            output_csv=dataset_csv,
            download_dir=args.output_dir / "dataset_assets",
            min_proteins=20,
        )

    data = pd.read_csv(dataset_csv)
    features = pd.DataFrame({
        "pocket_depth": data.get("pocket_depth", 20.0 - data["sasa"].clip(upper=20.0)),
        "sasa": data["sasa"],
        "electrostatic_potential": data.get("electrostatic_potential", 8.0 - 0.1 * data["sasa"]),
        "n_basic_residues": data.get("n_basic_residues", (8.0 - 0.2 * data["sasa"]).clip(lower=1).round()),
        "pocket_volume": data.get("pocket_volume", 600.0 - 3.5 * data["sasa"]),
        "plddt_confidence": data.get("plddt_confidence", 85.0),
    })
    labels = data["classification"].map(CLASS_MAP).fillna(0).astype(int)

    classifier = CrypticSiteMLClassifier(model_type="random_forest", random_state=42)
    metrics = classifier.fit(features, labels)

    model_path = args.output_dir / "cryptic_ip_rf.joblib"
    classifier.save(str(model_path))

    with (args.output_dir / "training_metrics.json").open("w", encoding="utf-8") as handle:
        json.dump(metrics.__dict__, handle, indent=2, default=list)

    print(f"Saved model: {model_path}")


if __name__ == "__main__":
    main()
