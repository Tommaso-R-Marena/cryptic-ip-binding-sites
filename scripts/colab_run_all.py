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
    subprocess.run(cmd, check=True, cwd=cwd)


def stage_install_check(log_path: Path) -> None:
    log("STAGE install-check", log_path)
    import importlib
    import shutil

    for tool in ("fpocket", "freesasa", "apbs", "pdb2pqr"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"Missing external tool: {tool}. Run: bash scripts/colab_install.sh")
    for mod in ("numpy", "pandas", "Bio", "prody", "sklearn"):
        importlib.import_module(mod)
    log("install-check OK", log_path)


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


def stage_ml(output_dir: Path, with_electrostatics: bool, log_path: Path) -> None:
    log("STAGE ml", log_path)
    cmd = [
        sys.executable,
        "scripts/train_ml_classifier.py",
        "--skip-build-dataset",
        "--work-dir",
        str(output_dir / "ml_training"),
        "--model-dir",
        str(ROOT / "models"),
    ]
    if with_electrostatics:
        cmd.append("--include-electrostatics")
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
    cmd = [
        sys.executable,
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
) -> None:
    log("STAGE publication", log_path)
    cmd = [
        sys.executable,
        "scripts/run_publication_package.py",
        "--output-dir",
        str(output_dir / "publication"),
        "--skip-dataset-build",
    ]
    if skip_figures:
        cmd.append("--skip-figures")
    if with_electrostatics:
        cmd.append("--with-electrostatics")
    else:
        cmd.append("--skip-electrostatics")
    run_cmd(cmd, log_path)


def stage_md(output_dir: Path, log_path: Path) -> None:
    log("STAGE md", log_path)
    try:
        import openmm  # noqa: F401
    except ImportError:
        run_cmd([sys.executable, "-m", "pip", "install", "-q", "openmm", "mdtraj"], log_path)
    candidates = output_dir / "publication" / "gallery" / "gallery_inputs.csv"
    if not candidates.exists():
        candidates = ROOT / "results" / "publication" / "gallery" / "gallery_inputs.csv"
    run_cmd(
        [
            sys.executable,
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
