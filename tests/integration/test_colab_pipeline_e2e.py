"""End-to-end integration test for the Colab orchestrator pipeline (minimal subset)."""

from __future__ import annotations

import json
import os
import shutil
from pathlib import Path

import pytest

from scripts.export_supplementary_tables import export_supplementary
from tests.integration.pipeline_helpers import (
    ensure_tier1_structures,
    resolve_ml_features_csv,
    seed_ml_structures,
    seed_yeast_structures,
)

pytestmark = [pytest.mark.integration, pytest.mark.slow]


@pytest.fixture(scope="module")
def fpocket_available() -> None:
    if shutil.which("fpocket") is not None:
        return
    micromamba_bin = Path("/tmp/mambaenv/bin")
    if (micromamba_bin / "fpocket").exists():
        os.environ["PATH"] = f"{micromamba_bin}{os.pathsep}{os.environ.get('PATH', '')}"
        return
    pytest.skip("fpocket not installed")


def test_colab_pipeline_minimal_e2e(tmp_path: Path, fpocket_available) -> None:
    """Run tier-1 -> ML -> yeast (n=2) -> publication -> supplements on a tiny subset."""
    from scripts.colab_run_all import (
        stage_install_check,
        stage_ml,
        stage_publication,
        stage_tier1,
        stage_yeast,
    )

    ensure_tier1_structures()
    features_csv = resolve_ml_features_csv(Path("data/validation/raw"))
    if features_csv is None and seed_ml_structures(Path("data/validation/raw")) < 2:
        pytest.skip("Need at least two cached validation PDBs or bundled ML features")

    output_dir = tmp_path / "colab_e2e"
    structures_dir = tmp_path / "yeast_structures"
    seed_yeast_structures(structures_dir, n_proteins=2)

    log_path = output_dir / "colab_run.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("", encoding="utf-8")

    stage_install_check(log_path)

    tier1 = stage_tier1(output_dir, with_electrostatics=False, log_path=log_path)
    sep = tier1.get("separation_quality", {})
    assert sep.get("phase1_ready"), f"Tier-1 gate failed: {sep}"

    stage_ml(
        output_dir,
        with_electrostatics=False,
        log_path=log_path,
        features_csv=features_csv,
    )
    ml_dir = output_dir / "ml_training"
    assert (ml_dir / "validation_pocket_features.csv").exists()
    assert (ml_dir / "ml_vs_threshold_comparison.csv").exists()
    assert (Path("models") / "cryptic_ip_classifier_v1.pkl").exists()

    stage_yeast(
        output_dir,
        structures_dir,
        n_proteins=2,
        workers=1,
        score_threshold=0.75,
        min_plddt=70.0,
        with_electrostatics=False,
        skip_download=True,
        log_path=log_path,
    )
    yeast_summary = json.loads((output_dir / "yeast_pilot" / "yeast_pilot_summary.json").read_text())
    assert yeast_summary["structures_screened"] == 2
    assert (output_dir / "yeast_pilot" / "yeast_pilot_hits.csv").exists()

    pub_ml = output_dir / "publication" / "ml_training"
    pub_ml.mkdir(parents=True, exist_ok=True)
    for artifact in (
        "validation_pocket_features.csv",
        "ml_vs_threshold_comparison.csv",
        "ml_vs_threshold_comparison.md",
    ):
        src = ml_dir / artifact
        if src.exists():
            shutil.copy2(src, pub_ml / artifact)

    pub_validation = output_dir / "publication" / "validation"
    pub_validation.mkdir(parents=True, exist_ok=True)
    for artifact in ("control_benchmark.csv", "positive_controls.csv", "negative_controls.csv"):
        src = output_dir / "validation" / artifact
        if src.exists():
            shutil.copy2(src, pub_validation / artifact)

    stage_publication(
        output_dir,
        with_electrostatics=False,
        skip_figures=True,
        log_path=log_path,
        skip_ml_training=True,
        skip_controls=True,
    )
    publication_dir = output_dir / "publication"
    assert (publication_dir / "RESULTS_SUMMARY.md").exists()
    assert (publication_dir / "validation" / "control_benchmark.csv").exists()

    supplementary_dir = export_supplementary(
        publication_dir,
        Path("data/validation/ip_binding_validation_dataset.csv"),
        output_dir / "yeast_pilot",
        publication_dir / "supplementary",
    )
    assert (supplementary_dir / "SUPPLEMENTARY_INDEX.md").exists()
    assert (supplementary_dir / "S2_control_benchmark.csv").exists()
    assert (supplementary_dir / "S6_yeast_pilot_summary.json").exists()

    log_text = log_path.read_text(encoding="utf-8")
    assert "STAGE tier1" in log_text
    assert "STAGE ml" in log_text
    assert "STAGE yeast" in log_text
    assert "STAGE publication" in log_text
