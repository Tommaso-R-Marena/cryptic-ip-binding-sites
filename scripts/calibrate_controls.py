#!/usr/bin/env python3
"""Report burial and score measurements for the validation controls.

Why this exists
---------------
Several thresholds in this pipeline are only meaningful relative to the
measurements they are applied to: the burial-class boundary on relative SASA, the
control pass criteria on the composite score, and the tier-1 separation gate.
Choosing any of them without looking at what the paradigm structures actually
measure is guesswork, and adjusting one afterwards to make a test pass is worse
than guesswork.

This script prints, for each control structure, everything a threshold decision
needs:

* per-ligand-copy burial: relative SASA, relative phosphate SASA, depth,
  enclosure, contacts and coordination counts;
* the pocket that best overlaps the ligand, with the rule-based scorer's
  **component breakdown**, so a low composite can be attributed to a specific
  component rather than guessed at.

Run it whenever the scorer, its weights, or the burial definition changes. The
output is the evidence a threshold change should cite.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.analysis.analyzer import ProteinAnalyzer  # noqa: E402
from cryptic_ip.analysis.labeling import LigandSite, assign_pocket_labels  # noqa: E402
from cryptic_ip.analysis.scorer import PocketScorer  # noqa: E402
from cryptic_ip.analysis.structure_arrays import load_structure_arrays  # noqa: E402
from cryptic_ip.validation.burial_metrics import (  # noqa: E402
    compute_ligand_burial,
    find_ligand_instances,
)
from cryptic_ip.validation.structure_context import LIGAND_RESNAMES  # noqa: E402

LOGGER = logging.getLogger("calibrate_controls")

#: Control panel: identifier, role, and the expected qualitative burial state
#: from the primary literature.
CONTROLS: Sequence[Dict[str, str]] = (
    {"pdb_id": "1ZY7", "name": "ADAR2", "role": "positive", "expected": "cryptic"},
    {"pdb_id": "5HDT", "name": "Pds5B", "role": "positive", "expected": "surface_or_artifact"},
    {"pdb_id": "5ICN", "name": "HDAC1", "role": "positive", "expected": "semi_cryptic"},
    {"pdb_id": "1MAI", "name": "PLCd1_PH", "role": "negative", "expected": "surface"},
    {"pdb_id": "1BWN", "name": "Btk_PH", "role": "negative", "expected": "surface"},
)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--structures-dir",
        type=Path,
        default=Path("data/validation"),
        help="Directory containing the control structure files",
    )
    parser.add_argument(
        "--output-json",
        type=Path,
        default=Path("results/calibration/control_measurements.json"),
        help="Where to write the measurements",
    )
    parser.add_argument(
        "--sasa-points", type=int, default=512, help="SASA sample points per atom"
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )
    return parser.parse_args(argv)


def _resolve_structure(structures_dir: Path, pdb_id: str) -> Optional[Path]:
    """Find a control structure file, accepting either coordinate format."""
    for suffix in (".pdb", ".cif", ".ent"):
        candidate = structures_dir / f"{pdb_id}{suffix}"
        if candidate.exists():
            return candidate
    return None


def measure_control(
    path: Path, *, sasa_points: int = 512
) -> Dict[str, Any]:
    """Measure burial and scoring for one control structure.

    Args:
        path: Structure file.
        sasa_points: SASA sample points per atom.

    Returns:
        A record with burial measurements and the scorer's component breakdown
        for the pocket that best overlaps the ligand.
    """
    record: Dict[str, Any] = {"structure": str(path)}

    instances = compute_ligand_burial(path, n_points=sasa_points)
    record["n_ligand_instances"] = len(instances)
    record["instances"] = [
        {
            "instance_id": inst.instance_id,
            "comp_id": inst.comp_id,
            "relative_sasa": inst.relative_sasa,
            "relative_phosphate_sasa": inst.relative_phosphate_sasa,
            "sasa_complex": inst.sasa_complex,
            "sasa_isolated": inst.sasa_isolated,
            "burial_depth": inst.burial_depth,
            "enclosure": inst.enclosure,
            "n_protein_contacts": inst.n_protein_contacts,
            "n_basic_residues": inst.n_basic_residues,
            "n_basic_nitrogens": inst.n_basic_nitrogens,
            "burial_class": inst.burial_class,
            "is_probable_artifact": inst.is_probable_artifact,
        }
        for inst in instances
    ]
    if instances:
        record["most_buried"] = record["instances"][0]

    # Locate the pocket that actually holds the ligand, and break its score down.
    arrays = load_structure_arrays(path)
    found = find_ligand_instances(arrays, LIGAND_RESNAMES)
    if not found:
        record["error"] = "no inositol phosphate ligand parsed"
        return record

    sites = []
    for key, comp_id, atom_indices in found:
        suffix = f"{key[2]}{key[3]}".strip()
        sites.append(
            LigandSite(
                instance_id=f"{comp_id}_{key[1]}_{suffix}",
                comp_id=comp_id,
                coords=arrays.coords[atom_indices],
            )
        )

    analyzer = ProteinAnalyzer(str(path), skip_electrostatics=True)
    analyzer.detect_pockets()
    analyzer.calculate_sasa()
    record["n_pockets"] = int(len(analyzer.pockets))
    if record["n_pockets"] == 0:
        record["error"] = "fpocket detected no pockets"
        return record

    pocket_records = []
    geometry = []
    for pocket_id in analyzer.pockets["pocket_id"]:
        analysis = analyzer.analyze_pocket(int(pocket_id))
        pocket_records.append(analysis)
        centre = np.asarray(analysis["center"], dtype=float)
        spheres = analyzer._pocket_alpha_spheres(int(pocket_id))
        points = spheres if spheres is not None and len(spheres) else centre.reshape(1, 3)
        geometry.append((int(pocket_id), centre, points))

    assignments = assign_pocket_labels(geometry, sites)
    best = max(assignments, key=lambda a: a.overlap_fraction)
    record["site_pocket"] = {
        "pocket_id": best.pocket_id,
        "overlap_fraction": best.overlap_fraction,
        "label": best.label.name.lower(),
        "min_atom_distance": best.min_atom_distance,
    }

    analysis = next(r for r in pocket_records if int(r["pocket_id"]) == best.pocket_id)
    scorer = PocketScorer()
    components = scorer.component_scores(
        volume=analysis.get("pocket_volume"),
        depth=analysis.get("burial_depth"),
        sasa=analysis.get("sasa_mean"),
        basic_count=analysis.get("n_basic_residues"),
        potential=analysis.get("coulomb_potential_kt"),
        enclosure=analysis.get("enclosure"),
    )
    record["site_pocket"]["measurements"] = {
        "pocket_volume": analysis.get("pocket_volume"),
        "burial_depth": analysis.get("burial_depth"),
        "enclosure": analysis.get("enclosure"),
        "sasa_mean": analysis.get("sasa_mean"),
        "n_basic_residues": analysis.get("n_basic_residues"),
        "coulomb_potential_kt": analysis.get("coulomb_potential_kt"),
    }
    record["site_pocket"]["component_scores"] = components
    record["site_pocket"]["weights"] = scorer.weights
    record["site_pocket"]["composite_score"] = scorer.calculate_composite_score(
        volume=analysis.get("pocket_volume"),
        depth=analysis.get("burial_depth"),
        sasa=analysis.get("sasa_mean"),
        basic_count=analysis.get("n_basic_residues"),
        potential=analysis.get("coulomb_potential_kt"),
        enclosure=analysis.get("enclosure"),
    )

    # The highest-scoring pocket overall, for comparison with the true site.
    best_scoring = max(
        pocket_records,
        key=lambda r: scorer.calculate_composite_score(
            volume=r.get("pocket_volume"),
            depth=r.get("burial_depth"),
            sasa=r.get("sasa_mean"),
            basic_count=r.get("n_basic_residues"),
            potential=r.get("coulomb_potential_kt"),
            enclosure=r.get("enclosure"),
        ),
    )
    record["top_scoring_pocket_id"] = int(best_scoring["pocket_id"])
    record["site_is_top_scoring"] = int(best_scoring["pocket_id"]) == best.pocket_id

    analyzer.cleanup()
    return record


def _format_value(value: Any, width: int = 8, precision: int = 3) -> str:
    """Format a possibly-missing number for the summary table."""
    if value is None:
        return "n/a".rjust(width)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value).rjust(width)
    if not np.isfinite(number):
        return "n/a".rjust(width)
    return f"{number:{width}.{precision}f}"


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Entry point."""
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level), format="%(levelname)s | %(message)s"
    )

    results: List[Dict[str, Any]] = []
    for control in CONTROLS:
        path = _resolve_structure(args.structures_dir, control["pdb_id"])
        record: Dict[str, Any] = dict(control)
        if path is None:
            record["error"] = f"structure not found in {args.structures_dir}"
            LOGGER.warning("%s (%s): not available", control["name"], control["pdb_id"])
            results.append(record)
            continue
        try:
            record.update(measure_control(path, sasa_points=args.sasa_points))
        except Exception as exc:  # noqa: BLE001 - report and continue
            record["error"] = f"{type(exc).__name__}: {exc}"
            LOGGER.warning("%s failed: %s", control["name"], exc)
        results.append(record)

    print("\n" + "=" * 108)
    print("CONTROL BURIAL MEASUREMENTS (most buried ligand copy)")
    print("=" * 108)
    print(
        f"{'control':<12}{'pdb':<6}{'role':<10}{'expected':<20}{'ligand':<8}"
        f"{'relSASA':>9}{'relPSASA':>10}{'depth':>8}{'encl':>8}{'basic':>7}{'class':>16}"
    )
    for record in results:
        best = record.get("most_buried")
        if not best:
            print(
                f"{record['name']:<12}{record['pdb_id']:<6}{record['role']:<10}"
                f"{record['expected']:<20}  {record.get('error', 'no ligand')}"
            )
            continue
        print(
            f"{record['name']:<12}{record['pdb_id']:<6}{record['role']:<10}"
            f"{record['expected']:<20}{str(best.get('comp_id', '?')):<8}"
            f"{_format_value(best['relative_sasa'], 9)}"
            f"{_format_value(best['relative_phosphate_sasa'], 10)}"
            f"{_format_value(best['burial_depth'], 8, 2)}"
            f"{_format_value(best['enclosure'], 8)}"
            f"{_format_value(best['n_basic_residues'], 7, 0)}"
            f"{best['burial_class']:>16}"
        )

    print("\n" + "=" * 108)
    print("RULE-BASED SCORE AT THE TRUE LIGAND SITE (component breakdown)")
    print("=" * 108)
    for record in results:
        site = record.get("site_pocket")
        if not site:
            continue
        components = site["component_scores"]
        measurements = site["measurements"]
        print(
            f"\n{record['name']} ({record['pdb_id']}) "
            f"pocket {site['pocket_id']} of {record.get('n_pockets')} | "
            f"ligand overlap {site['overlap_fraction']:.0%} | "
            f"composite {site['composite_score']:.3f} | "
            f"site is top-scoring: {record.get('site_is_top_scoring')}"
        )
        for name, score in sorted(components.items()):
            weight = site["weights"].get(name, 0.0)
            raw = {
                "volume": measurements.get("pocket_volume"),
                "depth": measurements.get("burial_depth"),
                "sasa": measurements.get("sasa_mean"),
                "basic_residues": measurements.get("n_basic_residues"),
                "electrostatics": measurements.get("coulomb_potential_kt"),
                "enclosure": measurements.get("enclosure"),
            }.get(name)
            print(
                f"    {name:<16} score={score:5.3f}  weight={weight:5.3f}  "
                f"contribution={score * weight:5.3f}  measured={_format_value(raw, 9, 2)}"
            )

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(results, indent=2, default=float), encoding="utf-8")
    print(f"\nWrote {args.output_json}")

    measured = [r for r in results if r.get("most_buried")]
    if not measured:
        LOGGER.error("No control could be measured; check that structures are available.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
