#!/usr/bin/env python3
"""Retrain ML classifier when validation dataset has changed."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-csv", type=Path, default=Path("data/validation/ip_binding_validation_dataset.csv"))
    parser.add_argument("--model-dir", type=Path, default=Path("models"))
    parser.add_argument("--work-dir", type=Path, default=Path("results/ml_training"))
    parser.add_argument("--force", action="store_true", help="Retrain even if dataset hash is unchanged")
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    metadata_path = args.model_dir / "cryptic_ip_classifier_v1.metadata.json"

    if not args.dataset_csv.exists():
        raise FileNotFoundError(f"Validation dataset not found: {args.dataset_csv}")

    dataset_hash = sha256_file(args.dataset_csv)
    previous_hash = None
    if metadata_path.exists():
        previous_hash = json.loads(metadata_path.read_text(encoding="utf-8")).get("dataset_sha256")

    if not args.force and previous_hash == dataset_hash:
        print("Dataset hash unchanged. Skipping retraining.")
        return

    cmd = [
        "python",
        "scripts/train_ml_classifier.py",
        "--dataset-csv",
        str(args.dataset_csv),
        "--model-dir",
        str(args.model_dir),
        "--work-dir",
        str(args.work_dir),
        "--skip-build-dataset",
    ]
    subprocess.run(cmd, check=True)

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["dataset_sha256"] = dataset_hash
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print("Retraining complete; metadata hash updated.")


if __name__ == "__main__":
    main()
