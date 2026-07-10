#!/usr/bin/env python3
"""Export manuscript supplementary tables from the publication package."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


SUPPLEMENTARY_MANIFEST = [
    ("S1_validation_dataset.csv", "Validation structures with burial labels"),
    ("S2_control_benchmark.csv", "Tiered positive/negative control scores"),
    ("S3_ml_benchmark.csv", "ML vs threshold held-out comparison"),
    ("S4_comparative_hit_rates.csv", "Burial-class / organism hit-rate summary"),
    ("S5_gallery_candidates.csv", "Top structural control candidates for Figure 3"),
    ("S6_yeast_pilot_summary.json", "Yeast AlphaFold pilot screen summary (if available)"),
    ("S6_yeast_pilot_hits.csv", "Yeast pilot cryptic-site hits (if available)"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--publication-dir", type=Path, default=Path("results/publication"))
    parser.add_argument("--dataset-csv", type=Path, default=Path("data/validation/ip_binding_validation_dataset.csv"))
    parser.add_argument("--yeast-dir", type=Path, default=Path("results/yeast_pilot"))
    parser.add_argument("--output-dir", type=Path, default=None, help="Defaults to publication-dir/supplementary")
    return parser.parse_args()


def _copy_table(src: Path, dest: Path) -> bool:
    if not src.exists():
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return True


def _csv_row_count(path: Path) -> int:
    if not path.exists() or path.stat().st_size == 0:
        return 0
    try:
        return len(pd.read_csv(path))
    except pd.errors.EmptyDataError:
        return 0


def export_supplementary(
    publication_dir: Path,
    dataset_csv: Path,
    yeast_dir: Path,
    output_dir: Path | None = None,
) -> Path:
    out = output_dir or (publication_dir / "supplementary")
    out.mkdir(parents=True, exist_ok=True)

    exported: list[dict] = []

    mappings = [
        (dataset_csv, out / "S1_validation_dataset.csv"),
        (publication_dir / "validation" / "control_benchmark.csv", out / "S2_control_benchmark.csv"),
        (publication_dir / "ml_training" / "ml_vs_threshold_comparison.csv", out / "S3_ml_benchmark.csv"),
        (publication_dir / "comparative" / "hit_rates.csv", out / "S4_comparative_hit_rates.csv"),
        (publication_dir / "gallery" / "gallery_inputs.csv", out / "S5_gallery_candidates.csv"),
    ]
    for src, dest in mappings:
        if _copy_table(src, dest):
            exported.append({"table": dest.name, "source": str(src), "rows": _csv_row_count(dest)})

    if (yeast_dir / "yeast_pilot_summary.json").exists():
        shutil.copy2(yeast_dir / "yeast_pilot_summary.json", out / "S6_yeast_pilot_summary.json")
        exported.append({"table": "S6_yeast_pilot_summary.json", "source": str(yeast_dir)})
    if (yeast_dir / "yeast_pilot_hits.csv").exists():
        _copy_table(yeast_dir / "yeast_pilot_hits.csv", out / "S6_yeast_pilot_hits.csv")
        exported.append(
            {
                "table": "S6_yeast_pilot_hits.csv",
                "source": str(yeast_dir / "yeast_pilot_hits.csv"),
                "rows": _csv_row_count(out / "S6_yeast_pilot_hits.csv"),
            }
        )

    # Pocket-feature summary for ML supplement
    features = publication_dir / "ml_training" / "validation_pocket_features.csv"
    if features.exists():
        feat_df = pd.read_csv(features)
        summary = {
            "n_pockets": int(len(feat_df)),
            "n_positive_pockets": int(feat_df["label"].sum()) if "label" in feat_df.columns else None,
            "n_structures": int(feat_df["pdb_id"].nunique()) if "pdb_id" in feat_df.columns else None,
        }
        (out / "S3_ml_feature_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        exported.append({"table": "S3_ml_feature_summary.json", "summary": summary})

    index_lines = [
        "# Supplementary Tables",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "| Table | Description | Status |",
        "| --- | --- | --- |",
    ]
    for name, desc in SUPPLEMENTARY_MANIFEST:
        path = out / name
        status = "included" if path.exists() else "pending"
        index_lines.append(f"| `{name}` | {desc} | {status} |")

    (out / "SUPPLEMENTARY_INDEX.md").write_text("\n".join(index_lines) + "\n", encoding="utf-8")
    (out / "export_manifest.json").write_text(json.dumps(exported, indent=2), encoding="utf-8")
    return out


def main() -> int:
    args = parse_args()
    out = export_supplementary(args.publication_dir, args.dataset_csv, args.yeast_dir, args.output_dir)
    print(f"Supplementary tables written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
