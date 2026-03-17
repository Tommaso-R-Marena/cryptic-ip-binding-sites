"""Reproducibility helpers for configuration, provenance, and archival workflows."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import random
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence

import numpy as np
import pandas as pd
import yaml
from jsonschema import Draft202012Validator


def set_global_seed(seed: Optional[int]) -> None:
    """Set deterministic seeds for Python and NumPy."""
    if seed is None:
        return
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def deterministic_sort_dataframe(
    frame: pd.DataFrame,
    sort_columns: Sequence[str],
    ascending: Sequence[bool] | bool = False,
) -> pd.DataFrame:
    """Sort with stable ordering to ensure deterministic output across runs."""
    return frame.sort_values(list(sort_columns), ascending=ascending, kind="mergesort").reset_index(
        drop=True
    )


def load_yaml(path: Path) -> Dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def validate_config(config: Mapping[str, Any], schema_path: Path) -> None:
    """Validate config dict against YAML or JSON schema."""
    schema = load_yaml(schema_path)
    validator = Draft202012Validator(schema)
    errors = sorted(validator.iter_errors(config), key=lambda err: list(err.path))
    if errors:
        messages = [f"{'/'.join(str(p) for p in err.path) or '<root>'}: {err.message}" for err in errors]
        raise ValueError("Configuration validation failed:\n" + "\n".join(messages))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8192), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_checksums(paths: Iterable[Path]) -> Dict[str, str]:
    return {str(path): sha256_file(path) for path in sorted(paths)}


def get_git_commit_hash(repo_dir: Path = Path(".")) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo_dir), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError:
        return "unknown"
    return result.stdout.strip()


def collect_software_versions(packages: Iterable[str]) -> Dict[str, str]:
    from importlib.metadata import PackageNotFoundError, version

    versions: Dict[str, str] = {
        "python": platform.python_version(),
        "platform": platform.platform(),
    }
    for package in packages:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def generate_provenance_manifest(
    *,
    config: Mapping[str, Any],
    inputs: Sequence[Path],
    outputs: Sequence[Path],
    parameters: Mapping[str, Any],
    data_sources: Mapping[str, Any],
    repo_dir: Path = Path("."),
) -> Dict[str, Any]:
    """Generate a JSON-LD provenance manifest for one analysis run."""
    packages = [
        "numpy",
        "pandas",
        "scipy",
        "biopython",
        "scikit-learn",
        "xgboost",
    ]
    manifest: Dict[str, Any] = {
        "@context": "https://w3id.org/ro/crate/1.1/context",
        "@type": "Dataset",
        "generatedAt": datetime.now(timezone.utc).isoformat(),
        "pipelineVersion": get_git_commit_hash(repo_dir),
        "config": config,
        "parameters": dict(parameters),
        "inputs": [{"path": str(path), "sha256": sha256_file(path)} for path in sorted(inputs)],
        "outputs": [{"path": str(path), "sha256": sha256_file(path)} for path in sorted(outputs)],
        "softwareVersions": collect_software_versions(packages),
        "dataSources": data_sources,
    }
    return manifest


def write_json(data: Mapping[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def parse_pinned_requirements(requirements_path: Path) -> Dict[str, str]:
    pinned = {}
    for raw in requirements_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if "==" not in line:
            continue
        name, version = line.split("==", 1)
        pinned[name.strip()] = version.strip()
    return pinned


def check_runtime_versions(requirements_path: Path) -> Dict[str, str]:
    pinned = parse_pinned_requirements(requirements_path)
    installed = collect_software_versions(pinned.keys())
    mismatches = {}
    for package, pinned_version in pinned.items():
        installed_version = installed.get(package)
        if installed_version != pinned_version:
            mismatches[package] = f"expected {pinned_version}, found {installed_version}"
    return mismatches


def export_analysis_bundle(
    *,
    bundle_dir: Path,
    include_paths: Sequence[Path],
    manifest: Mapping[str, Any],
    methods_text: str,
) -> Path:
    bundle_dir.mkdir(parents=True, exist_ok=True)
    copied: List[Path] = []
    for path in include_paths:
        if not path.exists():
            continue
        destination = bundle_dir / path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, destination)
        copied.append(destination)

    checksum_manifest = {str(path.relative_to(bundle_dir)): sha256_file(path) for path in sorted(copied)}
    write_json(dict(manifest), bundle_dir / "provenance_manifest.jsonld")
    write_json(checksum_manifest, bundle_dir / "checksums.json")
    (bundle_dir / "METHODS_AUTO.md").write_text(methods_text, encoding="utf-8")
    return bundle_dir


def generate_methods_text(config: Mapping[str, Any], manifest: Mapping[str, Any]) -> str:
    params = config.get("pipeline", {})
    data_sources = manifest.get("dataSources", {})
    return "\n".join(
        [
            "# Methods (auto-generated)",
            "",
            "We executed the cryptic IP-binding site pipeline with the following controlled parameters:",
            f"- Pocket score threshold: {params.get('score_threshold')}",
            f"- Maximum structures processed: {params.get('max_structures')}",
            f"- Random seed: {params.get('seed')}",
            f"- Deterministic sort columns: {', '.join(params.get('sort_columns', []))}",
            "",
            "Data sources:",
            f"- AlphaFold DB release date: {data_sources.get('alphafold_release_date', 'unknown')}",
            f"- PDB retrieval date: {data_sources.get('pdb_fetch_date', 'unknown')}",
            "",
            f"Pipeline commit: {manifest.get('pipelineVersion', 'unknown')}",
            "Floating point note: tiny score differences may occur across BLAS implementations;",
            "store outputs with >=6 decimal places and compare with tolerance (1e-6).",
        ]
    )
