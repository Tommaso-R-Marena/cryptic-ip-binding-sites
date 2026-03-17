#!/usr/bin/env python3
"""Utilities for reproducibility, provenance, and validation workflows."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.reproducibility import (
    check_runtime_versions,
    export_analysis_bundle,
    generate_methods_text,
    generate_provenance_manifest,
    load_yaml,
    validate_config,
    write_json,
)
from cryptic_ip.validation import validate_adar2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate_cmd = subparsers.add_parser("validate-config", help="Validate a YAML config")
    validate_cmd.add_argument("--config", type=Path, default=Path("config/defaults/pipeline.yaml"))
    validate_cmd.add_argument("--schema", type=Path, default=Path("config/schemas/pipeline.schema.yaml"))

    manifest_cmd = subparsers.add_parser("manifest", help="Generate JSON-LD provenance manifest")
    manifest_cmd.add_argument("--config", type=Path, default=Path("config/defaults/pipeline.yaml"))
    manifest_cmd.add_argument("--inputs", nargs="+", type=Path, required=True)
    manifest_cmd.add_argument("--outputs", nargs="+", type=Path, required=True)
    manifest_cmd.add_argument("--out", type=Path, default=Path("results/provenance_manifest.jsonld"))

    methods_cmd = subparsers.add_parser("methods", help="Generate methods text from config/manifest")
    methods_cmd.add_argument("--config", type=Path, default=Path("config/defaults/pipeline.yaml"))
    methods_cmd.add_argument("--manifest", type=Path, default=Path("results/provenance_manifest.jsonld"))
    methods_cmd.add_argument("--out", type=Path, default=Path("results/METHODS_AUTO.md"))

    bundle_cmd = subparsers.add_parser("bundle", help="Export full analysis bundle")
    bundle_cmd.add_argument("--config", type=Path, default=Path("config/defaults/pipeline.yaml"))
    bundle_cmd.add_argument("--manifest", type=Path, default=Path("results/provenance_manifest.jsonld"))
    bundle_cmd.add_argument("--methods", type=Path, default=Path("results/METHODS_AUTO.md"))

    version_cmd = subparsers.add_parser("check-versions", help="Validate pinned runtime versions")
    version_cmd.add_argument("--requirements", type=Path, default=Path("requirements.txt"))

    reproduce_cmd = subparsers.add_parser("reproduce", help="Reproduce baseline ADAR2 validation")
    reproduce_cmd.add_argument("--expected-json", type=Path, default=Path("config/defaults/expected_validation.json"))

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.command == "validate-config":
        config = load_yaml(args.config)
        validate_config(config, args.schema)
        print(f"Config is valid: {args.config}")
        return

    if args.command == "manifest":
        config = load_yaml(args.config)
        validate_config(config, Path("config/schemas/pipeline.schema.yaml"))
        manifest = generate_provenance_manifest(
            config=config,
            inputs=args.inputs,
            outputs=args.outputs,
            parameters=config["pipeline"],
            data_sources=config["data_sources"],
        )
        write_json(manifest, args.out)
        print(f"Wrote provenance manifest: {args.out}")
        return

    if args.command == "methods":
        config = load_yaml(args.config)
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        methods_text = generate_methods_text(config, manifest)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(methods_text, encoding="utf-8")
        print(f"Wrote methods text: {args.out}")
        return

    if args.command == "bundle":
        config = load_yaml(args.config)
        manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
        methods = args.methods.read_text(encoding="utf-8")
        include_paths = [Path(path) for path in config["archival"]["include_paths"]]
        out = export_analysis_bundle(
            bundle_dir=Path(config["archival"]["bundle_dir"]),
            include_paths=include_paths,
            manifest=manifest,
            methods_text=methods,
        )
        print(f"Exported analysis bundle: {out}")
        return

    if args.command == "check-versions":
        mismatches = check_runtime_versions(args.requirements)
        if mismatches:
            print("Version mismatch detected:")
            for package, mismatch in sorted(mismatches.items()):
                print(f"- {package}: {mismatch}")
            raise SystemExit(1)
        print("Runtime versions match pinned requirements.")
        return

    if args.command == "reproduce":
        expected = json.loads(args.expected_json.read_text(encoding="utf-8"))
        result = validate_adar2(use_alphafold=True)
        observed = {
            "validation_passed": bool(result.get("validation_passed")),
            "total_candidates": int(result.get("total_candidates", 0)),
        }
        if observed != expected:
            print(f"Expected {expected}, observed {observed}")
            raise SystemExit(1)
        print("Published baseline reproduced successfully.")


if __name__ == "__main__":
    main()
