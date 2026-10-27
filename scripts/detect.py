#!/usr/bin/env python3
"""The next-generation pocket detector: describe, score and rank one structure.

    python scripts/detect.py describe  <structure.pdb> [--out report.json]
    python scripts/detect.py shard     --structures DIR --shard 0 --shards 272 --out part.json
    python scripts/detect.py power     [--groups 5]
    python scripts/detect.py spec      [--json]

Every number this prints is computed from the structure given. Nothing is
calibrated against a measured affinity, and the pre-registered questions that
would make any of it a finding are in docs/DETECTOR_PLAN.md. Until those run,
treat the output as a hypothesis ordering and not as a prediction.

What it computes, and why each one is here
------------------------------------------
cavities      in-process LIGSITE-style detection, so a pocket can be found
              without fpocket and so hundreds of conformers are affordable
ensemble      persistence and crypticity across elastic-network conformers;
              this is the only part that can find a site that is closed in the
              apo model, which is what "cryptic" means
openability   whether a soft collective mode opens the pocket at all
field         the rotation-invariant electrostatic power spectrum, and its
              complementarity with an IP6 charge shell
specificity   IP6 against the competitor panel, with a pre-registered margin,
              reporting polyanion-generic rather than a candidate when that is
              what the evidence supports
score         a monotone additive score that abstains when an input is missing
              instead of failing the protein for a measurement nobody made
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis.structure_arrays import load_structure_arrays  # noqa: E402
from cryptic_ip.detect import cavity as CV  # noqa: E402
from cryptic_ip.detect import ensemble as EN  # noqa: E402
from cryptic_ip.detect import field as FI  # noqa: E402
from cryptic_ip.detect import monotone as MO  # noqa: E402
from cryptic_ip.detect import openability as OA  # noqa: E402
from cryptic_ip.detect import pu as PU  # noqa: E402
from cryptic_ip.detect import specificity as SP  # noqa: E402

#: The score above which a pocket is called a site, when nothing is missing.
#: Pre-registered in docs/DETECTOR_PLAN.md and not tuned on results.
DECISION_THRESHOLD = 0.55


def _echo(text: str = "") -> None:
    print(text)


def _polymer_atoms(arrays):
    """Polymer coordinates, elements, and the residue index of each atom."""
    mask = np.asarray(arrays.is_polymer, dtype=bool)
    if not np.any(mask):
        raise SystemExit("the structure has no polymer atoms")
    return (
        arrays.coords[mask].astype(float),
        [str(e) for e in np.asarray(arrays.elements, dtype=object)[mask]],
        np.asarray(arrays.residue_index, dtype=int)[mask],
        mask,
    )


def _node_map(arrays, mask):
    """Cα coordinates, and for each polymer atom the index of its own residue.

    The nearest-Cα shortcut is deliberately not used: inside a closed pocket
    the nearest alpha carbon to a side-chain atom is often a residue across the
    cavity, and displacing along a mode would then pull the two walls together
    rather than apart.
    """
    names = np.asarray(arrays.atom_names, dtype=object)
    polymer = np.asarray(arrays.is_polymer, dtype=bool)
    is_ca = polymer & (np.char.upper(names.astype(str)) == "CA")
    ca_residue = np.asarray(arrays.residue_index, dtype=int)[is_ca]
    node_coords = arrays.coords[is_ca].astype(float)
    position = {int(residue): index for index, residue in enumerate(ca_residue)}
    atom_residue = np.asarray(arrays.residue_index, dtype=int)[mask]
    missing = sorted({int(r) for r in atom_residue} - set(position))
    if missing:
        # A residue with no Cα cannot be moved by the network, so its atoms are
        # pinned to the nearest residue that can be. Reported rather than
        # silently done.
        nearest = {}
        ca_ids = np.array(sorted(position))
        for residue in missing:
            nearest[residue] = int(ca_ids[np.argmin(np.abs(ca_ids - residue))])
        for residue, substitute in nearest.items():
            position[residue] = position[substitute]
    return node_coords, np.array([position[int(r)] for r in atom_residue], dtype=int), len(missing)


def describe_structure(
    path: Path,
    *,
    spacing: float = 0.8,
    min_enclosure: int = 5,
    n_modes: int = 20,
    amplitudes: Sequence[float] = EN.DEFAULT_AMPLITUDES,
    which_modes: Sequence[int] = EN.DEFAULT_MODES,
    max_pockets: int = 5,
    field_points: int = FI.DEFAULT_POINTS_PER_SHELL,
    ensemble: bool = True,
) -> Dict[str, object]:
    """Run the whole detector on one structure."""
    arrays = load_structure_arrays(path)
    coords, elements, _, mask = _polymer_atoms(arrays)

    report: Dict[str, object] = {
        "structure": str(path),
        "n_atoms": int(len(coords)),
        "parameters": {
            "spacing": spacing, "min_enclosure": min_enclosure, "n_modes": n_modes,
            "amplitudes": list(amplitudes), "modes": list(which_modes),
            "field_points": field_points, "threshold": DECISION_THRESHOLD,
        },
    }

    apo = CV.detect_cavities(coords, elements, spacing=spacing, min_enclosure=min_enclosure)
    report["apo_cavities"] = apo.as_dict()

    tracks: List[EN.CavityTrack] = []
    if ensemble:
        try:
            node_coords, node_of_atom, pinned = _node_map(arrays, mask)
            modes = OA.normal_modes(node_coords, n_modes=n_modes)
            scan = EN.scan_with_nodes(
                coords, elements, node_coords, node_of_atom, modes,
                amplitudes=amplitudes, which_modes=which_modes,
                spacing=spacing, min_enclosure=min_enclosure,
            )
            report["ensemble"] = scan.as_dict()
            report["ensemble"]["residues_pinned_without_ca"] = pinned
            tracks = list(scan.tracks)
        except (OA.OpenabilityError, EN.EnsembleError, CV.CavityError) as exc:
            # An ensemble that cannot be built is reported, not silently
            # skipped: without it the crypticity column is absent rather than
            # zero, and a reader must be able to tell those apart.
            report["ensemble"] = {"unavailable": str(exc)}
            tracks = []

    sites: List[Dict[str, object]] = []
    candidates = tracks if tracks else [None] * min(len(apo.cavities), max_pockets)
    scorer = MO.MonotoneScorer()

    for index in range(min(max_pockets, max(len(candidates), 0))):
        track = tracks[index] if tracks else None
        cavity = apo.cavities[index] if (not tracks and index < len(apo.cavities)) else None
        centre = track.centre if track is not None else cavity.centre
        volume = track.max_volume if track is not None else cavity.volume
        crypticity = track.crypticity if track is not None else None
        persistence = track.persistence if track is not None else None

        field_report = FI.describe_pocket(arrays, centre, points_per_shell=field_points)
        ip6 = field_report.per_ligand[SP.TARGET]
        verdict = SP.verdict(field_report.per_ligand)

        try:
            open_report = OA.describe_pocket(arrays, centre, n_modes=min(n_modes, 20))
            openability = open_report.openability
            opening_cost = open_report.opening_cost_gamma_a2
        except OA.OpenabilityError as exc:
            openability, opening_cost = None, None
            open_report = None
            report.setdefault("warnings", []).append(f"openability at {centre}: {exc}")

        # Burial from this package's own detector, under its own name. The
        # ligand-copy relative SASA is a different measurement and abstains
        # unless a caller supplies it.
        enclosure = (track.mean_enclosure if track is not None
                     else (cavity.mean_enclosure if cavity is not None else None))
        openness = None if enclosure is None else float(
            np.clip(1.0 - enclosure / len(CV.SCAN_DIRECTIONS), 0.0, 1.0)
        )
        values: Dict[str, Optional[float]] = {
            "relative_sasa": None,
            "cavity_openness": openness,
            "field_depth_kt": float(np.clip(ip6["depth_kt"], -20.0, 20.0)),
            "field_shape": float(ip6["shape"]),
            "openability": openability,
            "n_basic_residues": _count_basic_nearby(arrays, centre),
            "pocket_volume": float(min(volume, 5000.0)),
            "mean_plddt": _mean_plddt_nearby(arrays, centre),
        }
        interval = scorer.score(values)
        sites.append({
            "centre": list(centre),
            "volume": float(volume),
            "crypticity": crypticity,
            "persistence": persistence,
            "field": field_report.as_dict(),
            "specificity": verdict.as_dict(),
            "openability": openability,
            "opening_cost_gamma_a2": opening_cost,
            "score": interval.as_dict(),
            "decision": interval.decision(DECISION_THRESHOLD),
            "blocking_feature": scorer.blocking_feature(values, DECISION_THRESHOLD),
            "features": {k: (None if v is None else float(v)) for k, v in values.items()},
        })

    report["sites"] = sites
    report["n_sites"] = len(sites)
    report["n_candidate_sites"] = sum(
        1 for site in sites
        if site["decision"] == "site" and site["specificity"]["is_candidate"]
    )
    return report


def _count_basic_nearby(arrays, centre, radius: float = 8.0) -> float:
    """Arg and Lys side-chain nitrogens within ``radius`` of the pocket centre."""
    coords = arrays.coords
    resnames = np.asarray(arrays.resnames, dtype=object)
    names = np.asarray(arrays.atom_names, dtype=object)
    centre = np.asarray(centre, dtype=float)
    distance = np.linalg.norm(coords - centre[None, :], axis=1)
    near = distance <= radius
    keys = {("ARG", "NH1"), ("ARG", "NH2"), ("ARG", "NE"), ("LYS", "NZ"), ("HIS", "NE2")}
    residues = set()
    for index in np.nonzero(near)[0]:
        if (str(resnames[index]).upper(), str(names[index]).upper()) in keys:
            residues.add(int(np.asarray(arrays.residue_index)[index]))
    return float(len(residues))


def _mean_plddt_nearby(arrays, centre, radius: float = 8.0) -> Optional[float]:
    """Mean B-factor near the pocket, which is pLDDT for a predicted model.

    ``None`` when the structure carries no usable B-factors, which is what the
    scorer's abstention is for: a crystal structure has no pLDDT and must not
    be failed for it.
    """
    bfactors = np.asarray(arrays.bfactors, dtype=float)
    if bfactors.size == 0 or np.allclose(bfactors, 0.0):
        return None
    centre = np.asarray(centre, dtype=float)
    near = np.linalg.norm(arrays.coords - centre[None, :], axis=1) <= radius
    if not np.any(near):
        return None
    value = float(np.mean(bfactors[near]))
    # A crystal structure's B-factors are not pLDDT and can exceed 100.
    return value if 0.0 <= value <= 100.0 else None


# ---------------------------------------------------------------- commands

def cmd_describe(args: argparse.Namespace) -> int:
    report = describe_structure(
        Path(args.structure),
        spacing=args.spacing,
        min_enclosure=args.min_enclosure,
        n_modes=args.modes,
        max_pockets=args.max_pockets,
        field_points=args.field_points,
        ensemble=not args.no_ensemble,
    )
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        _echo(f"wrote {args.out}")
    scorer = MO.MonotoneScorer()
    _echo(f"{report['structure']}: {report['n_atoms']} polymer atoms")
    apo = report["apo_cavities"]
    _echo(f"  apo cavities: {apo['n_cavities']}, total {apo['total_volume']:.0f} A^3")
    ens = report.get("ensemble") or {}
    if "unavailable" in ens:
        _echo(f"  ensemble: unavailable ({ens['unavailable']})")
    elif ens:
        _echo(f"  ensemble: {ens['n_conformers']} conformers, {ens['n_tracks']} tracked site(s), "
              f"{ens['n_cryptic']} cryptic or semi-cryptic")
    for index, site in enumerate(report["sites"]):
        score = site["score"]
        band = (f"{score['low']:.3f}" if score["point"] is not None
                else f"{score['low']:.3f}-{score['high']:.3f}")
        crypt = "n/a" if site["crypticity"] is None else f"{site['crypticity']:.2f}"
        _echo(f"  site {index}: score {band} -> {site['decision']}; "
              f"crypticity {crypt}; {site['specificity']['call']}")
        _echo(f"      {site['specificity']['explanation']}")
        for name in score.get("vetoed", ()):
            feature = next(f for f in scorer.features if f.name == name)
            _echo(f"      ruled out by {name}: {feature.veto_why.splitlines()[0]}")
        for name in score.get("veto_pending", ()):
            _echo(f"      cannot be cleared without measuring {name}")
        if site["blocking_feature"] and not score.get("vetoed"):
            _echo(f"      most limiting measurement: {site['blocking_feature']}")
        if score["missing"]:
            _echo(f"      not measured: {', '.join(score['missing'])}")
    _echo()
    _echo(f"{report['n_candidate_sites']} site(s) both scored above threshold and IP6-selective.")
    return 0


def cmd_shard(args: argparse.Namespace) -> int:
    """One shard of a structure directory, for the cluster."""
    directory = Path(args.structures)
    paths = sorted(p for p in directory.iterdir() if p.suffix.lower() in {".pdb", ".cif", ".ent"})
    if not paths:
        _echo(f"no structures in {directory}")
        return 1
    mine = [p for index, p in enumerate(paths) if index % args.shards == args.shard]
    _echo(f"shard {args.shard} of {args.shards}: {len(mine)} of {len(paths)} structures")
    out: List[Dict[str, object]] = []
    for path in mine:
        try:
            out.append(describe_structure(
                path, spacing=args.spacing, min_enclosure=args.min_enclosure,
                n_modes=args.modes, max_pockets=args.max_pockets,
                field_points=args.field_points, ensemble=not args.no_ensemble,
            ))
            _echo(f"  {path.name}: ok")
        except Exception as exc:  # noqa: BLE001 - one bad structure must not take the shard
            out.append({"structure": str(path), "error": f"{type(exc).__name__}: {exc}"})
            _echo(f"  {path.name}: {type(exc).__name__}: {exc}")
        # Written after every structure, so a task killed partway leaves what
        # it finished rather than nothing - the lesson of study J's shards.
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    _echo(f"wrote {args.out} ({len(out)} record(s))")
    return 0


def cmd_power(args: argparse.Namespace) -> int:
    """The minimum detectable effect, before any data is looked at."""
    _echo("Minimum detectable AUC at power 0.8, two-sided 0.05, by group count.")
    _echo("Computed from group counts rather than row counts, because copies of")
    _echo("one protein are not independent.")
    _echo()
    _echo(f"{'positive groups':>16} {'min detectable AUC':>20}")
    for groups in (1, 2, 3, 5, 10, 20, 40, 100):
        value = PU.minimum_detectable_auc(groups, args.negative_groups)
        _echo(f"{groups:>16} {value:>20.2f}")
    _echo()
    here = PU.minimum_detectable_auc(args.groups, args.negative_groups)
    _echo(f"This project's buried-site positives occupy {args.groups} strict homology "
          f"group(s), one holding 87% of them.")
    _echo(f"Nothing below AUC {here:.2f} is detectable at that size, so the previously "
          "reported 0.499 is a statement about power, not about biology.")
    return 0


def cmd_spec(args: argparse.Namespace) -> int:
    """Print the scoring specification and the competitor panel."""
    scorer = MO.MonotoneScorer()
    if args.json:
        _echo(json.dumps({
            "threshold": DECISION_THRESHOLD,
            "scorer": scorer.as_dict(),
            "templates": {k: {
                "net_charge": v.net_charge, "shell_charge": list(v.shell_charge),
                "radius_of_gyration": v.radius_of_gyration, "note": v.note,
            } for k, v in FI.LIGAND_TEMPLATES.items()},
            "specificity_margin": SP.DEFAULT_MARGIN,
            "shells": list(FI.DEFAULT_SHELLS),
        }, indent=2, sort_keys=True))
        return 0
    _echo(f"Decision threshold {DECISION_THRESHOLD}, specificity margin {SP.DEFAULT_MARGIN}.")
    _echo()
    _echo(f"{'feature':<20} {'dir':>4} {'midpoint':>10} {'weight':>7}  why")
    for feature in scorer.features:
        arrow = "up" if feature.direction == MO.UP else "down"
        _echo(f"{feature.name:<20} {arrow:>4} {feature.midpoint:>10g} {feature.weight:>7g}  "
              f"{feature.why.splitlines()[0][:70]}")
    _echo()
    _echo("Competitor panel (shell charge is cumulative, inside each shell radius):")
    _echo(f"  shells: {', '.join(str(s) for s in FI.DEFAULT_SHELLS)} A")
    for name, template in sorted(FI.LIGAND_TEMPLATES.items()):
        _echo(f"  {name:<14} net {template.net_charge:>5.1f}  Rg {template.radius_of_gyration:>4.1f}  "
              f"{template.note[:60]}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    def add_common(p):
        p.add_argument("--spacing", type=float, default=0.8, help="cavity grid spacing (A)")
        p.add_argument("--min-enclosure", type=int, default=5, help="of 7 directions")
        p.add_argument("--modes", type=int, default=20, help="elastic-network modes to compute")
        p.add_argument("--max-pockets", type=int, default=5)
        p.add_argument("--field-points", type=int, default=FI.DEFAULT_POINTS_PER_SHELL)
        p.add_argument("--no-ensemble", action="store_true",
                       help="skip the conformer ensemble, which is most of the cost")

    p = sub.add_parser("describe", help="run the detector on one structure")
    p.add_argument("structure")
    p.add_argument("--out")
    add_common(p)
    p.set_defaults(func=cmd_describe)

    p = sub.add_parser("shard", help="one shard of a structure directory")
    p.add_argument("--structures", required=True)
    p.add_argument("--shard", type=int, required=True)
    p.add_argument("--shards", type=int, required=True)
    p.add_argument("--out", required=True)
    add_common(p)
    p.set_defaults(func=cmd_shard)

    p = sub.add_parser("power", help="the minimum detectable effect at this dataset's size")
    p.add_argument("--groups", type=int, default=5, help="positive homology groups")
    p.add_argument("--negative-groups", type=int, default=100)
    p.set_defaults(func=cmd_power)

    p = sub.add_parser("spec", help="print the scoring specification and competitor panel")
    p.add_argument("--json", action="store_true")
    p.set_defaults(func=cmd_spec)

    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
