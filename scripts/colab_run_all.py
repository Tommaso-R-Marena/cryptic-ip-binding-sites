#!/usr/bin/env python3
"""Run the complete cryptic IP pipeline (designed for Google Colab and local use).

Stages:
  install-check  Verify fpocket/APBS and Python imports
  structures     Download tier-1 validation PDBs (1ZY7, 1MAI)
  tier1          Phase 1 validation gate (ADAR2 vs PLCδ1)
  ml             Train/evaluate ML classifier on validation dataset
  yeast          Yeast AlphaFold pilot screen
  publication    Manuscript package (controls, figures, provenance)
  md             Optional short OpenMM MD pilot
  package        Zip results for download

Usage (Colab):
  python scripts/colab_run_all.py --preset quick
  python scripts/colab_run_all.py --preset pilot --output-dir /content/drive/MyDrive/cryptic_ip
  python scripts/colab_run_all.py --preset full --skip-md
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.colab_env import bootstrap_colab_runtime, reexec_if_needed

PRESETS = {
    "quick": {
        "n_proteins": 25,
        "workers": 2,
        "run_ml": True,
        "run_yeast": True,
        "run_publication": True,
        "run_figures": True,
        "run_md": False,
        "with_electrostatics": False,
    },
    "pilot": {
        "n_proteins": 500,
        "workers": 2,
        "run_ml": True,
        "run_yeast": True,
        "run_publication": True,
        "run_figures": True,
        "run_md": False,
        "with_electrostatics": False,
    },
    "full": {
        "n_proteins": 500,
        "workers": 2,
        "run_ml": True,
        "run_yeast": True,
        "run_publication": True,
        "run_figures": True,
        "run_md": True,
        "with_electrostatics": False,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--preset", choices=sorted(PRESETS), default="quick")
    parser.add_argument("--output-dir", type=Path, default=Path("results/colab_run"))
    parser.add_argument("--structures-dir", type=Path, default=Path("data/structures/yeast_pilot"))
    parser.add_argument("--n-proteins", type=int, default=None)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--score-threshold", type=float, default=0.75)
    parser.add_argument("--min-plddt", type=float, default=70.0)
    parser.add_argument("--with-electrostatics", action="store_true")
    parser.add_argument("--skip-electrostatics", action="store_true", default=True)
    parser.add_argument("--skip-ml", action="store_true")
    parser.add_argument("--skip-yeast", action="store_true")
    parser.add_argument("--skip-publication", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--skip-md", action="store_true")
    parser.add_argument("--skip-download", action="store_true", help="Reuse cached AlphaFold structures")
    parser.add_argument("--skip-package", action="store_true")
    return parser.parse_args()


def log(msg: str, log_path: Path) -> None:
    line = f"[{datetime.now(timezone.utc).isoformat()}] {msg}"
    print(line, flush=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


def run_cmd(cmd: Sequence[str], log_path: Path, *, cwd: Path = ROOT) -> None:
    log(f"CMD: {' '.join(cmd)}", log_path)
    env = os.environ.copy()
    bootstrap_colab_runtime()
    subprocess.run(cmd, check=True, cwd=cwd, env=env)


def stage_install_check(log_path: Path) -> None:
    log("STAGE install-check", log_path)
    import importlib
    import shutil

    bootstrap_colab_runtime()
    if shutil.which("fpocket") is None:
        raise RuntimeError("Missing fpocket. Run: bash scripts/colab_install.sh")

    for mod in ("numpy", "pandas", "Bio", "prody", "sklearn", "cryptic_ip"):
        importlib.import_module(mod)
    log(f"install-check OK (python={sys.executable})", log_path)


def stage_structures(log_path: Path) -> None:
    log("STAGE structures", log_path)
    validation_dir = ROOT / "data" / "validation"
    validation_dir.mkdir(parents=True, exist_ok=True)
    for pdb_id in ("1ZY7", "1MAI", "1BWN"):
        dest = validation_dir / f"{pdb_id}.pdb"
        if dest.exists():
            continue
        run_cmd(
            ["wget", "-q", "-O", str(dest), f"https://files.rcsb.org/download/{pdb_id}.pdb"],
            log_path,
        )
    log(f"Validation structures ready in {validation_dir}", log_path)


def stage_tier1(output_dir: Path, with_electrostatics: bool, log_path: Path) -> dict:
    log("STAGE tier1", log_path)
    from cryptic_ip.validation.validation_suite import ValidationSuite

    suite = ValidationSuite(data_dir=str(ROOT / "data" / "validation"), use_electrostatics=with_electrostatics)
    summary = suite.run_full_validation(output_dir=output_dir / "validation")
    sep = summary.get("separation_quality", {})
    if not sep.get("phase1_ready"):
        raise RuntimeError(
            f"Tier-1 gate FAILED: separation={sep.get('tier1_separation')} "
            f"(need >0.50 ADAR2 vs PLCδ1)"
        )
    log(f"tier1 PASSED separation={sep.get('tier1_separation'):.3f}", log_path)
    return summary


def has_current_feature_schema(path: Path) -> bool:
    """Whether a feature table matches the schema the trainer expects.

    Feature tables produced before the descriptor rework carry the legacy
    six-column schema and labels from the previous labelling scheme. Training on
    one would reproduce the defect it was built with, so such a table is
    regenerated rather than reused.

    Args:
        path: Candidate feature CSV.

    Returns:
        ``True`` when the table carries the current descriptor schema.
    """
    if path is None or not Path(path).exists():
        return False
    try:
        with Path(path).open("r", encoding="utf-8") as handle:
            header = handle.readline().strip().split(",")
    except OSError:
        return False
    required = {"group_key", "label", "burial_depth", "enclosure", "coulomb_potential_kt"}
    return required.issubset(set(header))


def stage_features(
    output_dir: Path,
    log_path: Path,
    *,
    structures_dir: Path | None = None,
    jobs: int = 2,
    sasa_points: int = 256,
) -> Path:
    """Extract pocket descriptors and labels from the control structures.

    Args:
        output_dir: Run output directory.
        log_path: Log file.
        structures_dir: Directory of structures; defaults to the validation set.
        jobs: Worker processes.
        sasa_points: SASA sample points per atom.

    Returns:
        Path to the feature table written.
    """
    log("STAGE features", log_path)
    python = bootstrap_colab_runtime()
    structures_dir = structures_dir or (ROOT / "data" / "validation")
    ml_dir = output_dir / "ml_training"
    features_csv = ml_dir / "pocket_features.csv"
    run_cmd(
        [
            python,
            "scripts/extract_pocket_features.py",
            "--structures-dir",
            str(structures_dir),
            "--entry-csv",
            str(ROOT / "data" / "validation" / "ip_binding_validation_dataset.csv"),
            "--output-csv",
            str(features_csv),
            "--cache-dir",
            str(ml_dir / "feature_cache"),
            "--summary-json",
            str(ml_dir / "labeling_summary.json"),
            "--jobs",
            str(jobs),
            "--sasa-points",
            str(sasa_points),
        ],
        log_path,
    )
    return features_csv


def stage_ml(
    output_dir: Path,
    with_electrostatics: bool,
    log_path: Path,
    *,
    features_csv: Path | None = None,
    structures_dir: Path | None = None,
    n_splits: int = 2,
    n_search_iter: int = 4,
    n_bootstrap: int = 100,
) -> None:
    """Train and compare classifiers on a current pocket feature table.

    Dataset building and feature extraction are separate, independently runnable
    stages now, so this consumes a feature table. When none is supplied - or the
    supplied one predates the descriptor rework - features are extracted first
    rather than training on a stale table.
    """
    if not has_current_feature_schema(features_csv):
        if features_csv is not None:
            log(
                f"features table {features_csv} predates the current schema; regenerating",
                log_path,
            )
        features_csv = stage_features(
            output_dir, log_path, structures_dir=structures_dir
        )

    log("STAGE ml", log_path)
    python = bootstrap_colab_runtime()
    # A deliberately light search budget. This stage demonstrates the pipeline
    # end to end on a small control set; the full default (five model families,
    # 5x3 nested folds, 40 search iterations, 2000 bootstrap resamples) takes
    # over ten minutes on a set this size and is meant for the real training run
    # via scripts/train_ml_classifier.py.
    cmd = [
        python,
        "scripts/train_ml_classifier.py",
        "--features-csv",
        str(features_csv),
        "--work-dir",
        str(output_dir / "ml_training"),
        "--model-dir",
        str(ROOT / "models"),
        "--model-name",
        "cryptic_ip_classifier_colab",
        "--models",
        "logistic_regression",
        "random_forest",
        "--n-splits",
        str(n_splits),
        "--inner-splits",
        "2",
        "--n-search-iter",
        str(n_search_iter),
        "--n-bootstrap",
        str(n_bootstrap),
    ]
    _ = with_electrostatics
    run_cmd(cmd, log_path)


def stage_yeast(
    output_dir: Path,
    structures_dir: Path,
    *,
    n_proteins: int,
    workers: int,
    score_threshold: float,
    min_plddt: float,
    with_electrostatics: bool,
    skip_download: bool,
    log_path: Path,
) -> None:
    log(f"STAGE yeast (n={n_proteins})", log_path)
    python = bootstrap_colab_runtime()
    cmd = [
        python,
        "scripts/run_yeast_pilot_screen.py",
        "--n-proteins",
        str(n_proteins),
        "--workers",
        str(workers),
        "--score-threshold",
        str(score_threshold),
        "--min-plddt",
        str(min_plddt),
        "--output-dir",
        str(output_dir / "yeast_pilot"),
        "--structures-dir",
        str(structures_dir),
    ]
    if skip_download:
        cmd.append("--skip-download")
    if with_electrostatics:
        cmd.append("--with-electrostatics")
    else:
        cmd.append("--skip-electrostatics")
    run_cmd(cmd, log_path)
    summary = json.loads((output_dir / "yeast_pilot" / "yeast_pilot_summary.json").read_text())
    log(f"yeast hit_rate={summary.get('hit_rate'):.4f} proteins_with_hits={summary.get('proteins_with_hits')}", log_path)


def stage_publication(
    output_dir: Path,
    *,
    with_electrostatics: bool,
    skip_figures: bool,
    log_path: Path,
    skip_ml_training: bool = False,
    skip_controls: bool = False,
) -> None:
    log("STAGE publication", log_path)
    python = bootstrap_colab_runtime()
    cmd = [
        python,
        "scripts/run_publication_package.py",
        "--output-dir",
        str(output_dir / "publication"),
        "--skip-dataset-build",
    ]
    if skip_figures:
        cmd.append("--skip-figures")
    if skip_ml_training:
        cmd.append("--skip-ml-training")
    if skip_controls:
        cmd.append("--skip-controls")
    if with_electrostatics:
        cmd.append("--with-electrostatics")
    else:
        cmd.append("--skip-electrostatics")
    run_cmd(cmd, log_path)


def stage_md(output_dir: Path, log_path: Path) -> None:
    log("STAGE md", log_path)
    python = bootstrap_colab_runtime()
    try:
        import openmm  # noqa: F401
    except ImportError:
        run_cmd([python, "-m", "pip", "install", "-q", "openmm", "mdtraj"], log_path)
    candidates = output_dir / "publication" / "gallery" / "gallery_inputs.csv"
    if not candidates.exists():
        candidates = ROOT / "results" / "publication" / "gallery" / "gallery_inputs.csv"
    run_cmd(
        [
            python,
            "scripts/run_md_pilot_validation.py",
            "--candidates-csv",
            str(candidates),
            "--output-dir",
            str(output_dir / "md_validation"),
            "--top-n",
            "3",
            "--production-ns",
            "0.5",
        ],
        log_path,
    )


def stage_package(output_dir: Path, log_path: Path) -> Path:
    log("STAGE package", log_path)
    archive = output_dir.parent / "colab_pipeline_results.zip"
    if archive.exists():
        archive.unlink()
    shutil.make_archive(str(archive.with_suffix("")), "zip", output_dir)
    log(f"Created {archive}", log_path)
    return archive


def resolve_config(args: argparse.Namespace) -> dict:
    cfg = dict(PRESETS[args.preset])
    if args.n_proteins is not None:
        cfg["n_proteins"] = args.n_proteins
    if args.workers is not None:
        cfg["workers"] = args.workers
    cfg["with_electrostatics"] = args.with_electrostatics and not args.skip_electrostatics
    if args.skip_ml:
        cfg["run_ml"] = False
    if args.skip_yeast:
        cfg["run_yeast"] = False
    if args.skip_publication:
        cfg["run_publication"] = False
    if args.skip_figures:
        cfg["run_figures"] = False
    if args.skip_md:
        cfg["run_md"] = False
    cfg["skip_download"] = args.skip_download
    cfg["score_threshold"] = args.score_threshold
    cfg["min_plddt"] = args.min_plddt
    cfg["skip_package"] = args.skip_package
    return cfg


def main() -> int:
    reexec_if_needed()
    args = parse_args()
    cfg = resolve_config(args)
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "colab_run.log"
    log_path.write_text("", encoding="utf-8")

    manifest = {"preset": args.preset, "started_at": datetime.now(timezone.utc).isoformat(), "config": cfg}
    log(f"Starting colab_run_all preset={args.preset} output={output_dir}", log_path)

    t0 = time.time()
    stage_install_check(log_path)
    stage_structures(log_path)
    tier1 = stage_tier1(output_dir, cfg["with_electrostatics"], log_path)

    if cfg["run_ml"]:
        stage_ml(output_dir, cfg["with_electrostatics"], log_path)

    if cfg["run_yeast"]:
        stage_yeast(
            output_dir,
            args.structures_dir,
            n_proteins=cfg["n_proteins"],
            workers=cfg["workers"],
            score_threshold=cfg["score_threshold"],
            min_plddt=cfg["min_plddt"],
            with_electrostatics=cfg["with_electrostatics"],
            skip_download=cfg["skip_download"],
            log_path=log_path,
        )

    if cfg["run_publication"]:
        stage_publication(
            output_dir,
            with_electrostatics=cfg["with_electrostatics"],
            skip_figures=not cfg["run_figures"],
            log_path=log_path,
        )

    if cfg["run_md"]:
        try:
            stage_md(output_dir, log_path)
        except Exception as exc:
            log(f"MD stage skipped/failed: {exc}", log_path)

    archive = None
    if not cfg["skip_package"]:
        archive = stage_package(output_dir, log_path)

    elapsed = time.time() - t0
    manifest.update(
        {
            "finished_at": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": elapsed,
            "tier1_separation": tier1.get("separation_quality", {}).get("tier1_separation"),
            "archive": str(archive) if archive else None,
        }
    )
    (output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    log(f"DONE in {elapsed/60:.1f} min — results in {output_dir}", log_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
