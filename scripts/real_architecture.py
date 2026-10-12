#!/usr/bin/env python3
"""Study Q7 on real pocket linings, against a decoy arm matched to each one.

Pre-registered in `docs/ARCHITECTURE_PLAN.md` as amended by
`docs/ARCHITECTURE_PLAN_AMENDMENT_1.md`, which fixes what is computed here
before any of it was read.

The linings are real: `results/specificity/candidates.csv` carries the lining
positions of 75 proteome-screen candidate pockets with a homology cluster each.
They are **candidates, not annotated sites** - `annotated` is False for every
row - so nothing here is evidence that annotated IP sites have a sequence-space
architecture. It is a description of what real candidate linings look like in
the invariant space, beside a null.

The null is the part that makes it readable. Every real lining is paired with
decoys drawn to the **same number of positions and the same span**, which is
the comparison that controls for the two things a detection rate is most easily
confounded by.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import architecture as ar  # noqa: E402
from cryptic_ip.detect import sequence as sq  # noqa: E402

SEED = 20261012
DECOYS_PER_LINING = 5


def load(path: Path) -> List[Dict[str, object]]:
    records = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line:
            records.append(json.loads(line))
    if not records:
        raise SystemExit(f"{path} holds no records")
    return records


def positions_of(record: Dict[str, object]) -> List[int]:
    return sorted({int(p) for _chain, p in record["residues"]})  # type: ignore[union-attr]


def measure(positions: Sequence[int], *, n_draws: int, seed: int) -> Dict[str, object]:
    """One lining's architecture, reduced to what the amendment asks for."""
    try:
        arch = ar.describe_architecture(positions, n_draws=n_draws, seed=seed)
    except ar.ArchitectureError as error:
        return {"measurable": False, "reason": str(error)}
    topology = arch.topology_p if arch.topology_p is not None else 1.0
    register = arch.register_p if arch.register_p is not None else 1.0
    return {
        "measurable": True,
        "n_positions": arch.n_positions,
        "span": arch.span,
        "segments": arch.segments,
        "occupancy": arch.occupancy,
        "gap_entropy": arch.gap_entropy,
        "dominant_period": float(arch.dominant.get("period", float("nan"))),
        "topology_p": float(topology),
        "register_p": float(register),
        "either_p": float(min(topology, register)),
        "verdict": arch.signature_verdict,
        "architecture_class": arch.architecture_class,
        "element": arch.element,
    }


def decoys_for(record: Dict[str, object], *, n_draws: int, rng: np.random.Generator,
               count: int = DECOYS_PER_LINING) -> List[Dict[str, object]]:
    """Random linings with this one's exact size and span."""
    positions = positions_of(record)
    n, span = len(positions), positions[-1] - positions[0]
    out = []
    for index in range(count):
        interior = rng.choice(np.arange(1, span), size=n - 2, replace=False)
        drawn = sorted(np.concatenate([[0], [span], interior]).astype(int).tolist())
        out.append(measure(drawn, n_draws=n_draws, seed=SEED + index))
    return out


def rate(rows: Sequence[Dict[str, object]], key: str = "either_p",
         alpha: float = 0.05) -> float:
    usable = [r for r in rows if r.get("measurable")]
    if not usable:
        return float("nan")
    return float(np.mean([float(r[key]) <= alpha for r in usable]))


def stratified(rows: Sequence[Dict[str, object]], records: Sequence[Dict[str, object]],
               n_bootstrap: int) -> Dict[str, object]:
    """The one pre-registered split: pockets explained by another polyanion.

    Group means over homology clusters, the project's standing five-group
    floor, Holm across the features tried.
    """
    features = ("either_p", "occupancy", "gap_entropy", "segments")
    arms: Dict[str, Dict[str, Dict[str, List[float]]]] = {}
    for row, record in zip(rows, records):
        if not row.get("measurable"):
            continue
        label = str(record.get("label", "unknown"))
        group = str(record["group"])
        for feature in features:
            arms.setdefault(feature, {}).setdefault(label, {}).setdefault(group, [])
            arms[feature][label][group].append(float(row[feature]))

    comparisons, p_values = [], {}
    for feature, by_label in arms.items():
        labels = sorted(by_label)
        if len(labels) != 2:
            comparisons.append({"feature": feature, "runnable": False,
                                "reason": f"{len(labels)} label(s), need exactly 2"})
            continue
        first, second = labels
        left = {g: float(np.mean(v)) for g, v in by_label[first].items()}
        right = {g: float(np.mean(v)) for g, v in by_label[second].items()}
        if not left or not right:
            comparisons.append({"feature": feature, "runnable": False,
                                "reason": "an arm has no group"})
            continue
        result = sq.compare_feature(feature, left, right,
                                    n_bootstrap=n_bootstrap, seed=SEED)
        entry = result.as_dict()
        entry.update({"runnable": True, "arms": [first, second]})
        comparisons.append(entry)
        p_values[feature] = result.p_value
    corrected = sq.holm(p_values) if p_values else {}
    for entry in comparisons:
        if entry.get("runnable"):
            entry["p_holm"] = corrected.get(str(entry["feature"]))
            entry["is_finding"] = bool(
                entry["p_holm"] is not None and entry["p_holm"] < 0.05
                and entry["is_evidence"] and entry["exceeds_power"]
                and not (entry["interval"][0] <= 0.0 <= entry["interval"][1]))
    return {"features": comparisons, "holm_family_size": len(p_values)}


def run(inputs: Path, *, n_draws: int, n_bootstrap: int) -> Dict[str, object]:
    records = load(inputs)
    rng = np.random.default_rng(SEED)
    real = [measure(positions_of(r), n_draws=n_draws, seed=SEED + i)
            for i, r in enumerate(records)]
    decoy: List[Dict[str, object]] = []
    for record in records:
        decoy.extend(decoys_for(record, n_draws=n_draws, rng=rng))

    usable = [r for r in real if r.get("measurable")]
    sizes = [int(r["n_positions"]) for r in usable]
    return {
        "seed": SEED,
        "inputs": str(inputs),
        "n_linings": len(records),
        "n_measurable": len(usable),
        "lining_size": {
            "min": int(min(sizes)), "median": float(np.median(sizes)),
            "max": int(max(sizes)),
        },
        "classes": dict(Counter(r.get("architecture_class", "unmeasurable") for r in real)),
        "verdicts": dict(Counter(r.get("verdict", "unmeasurable") for r in real)),
        "elements": dict(Counter(r.get("element", "unmeasurable") for r in real)),
        "rates": {
            "real_either": rate(real), "real_topology": rate(real, "topology_p"),
            "real_register": rate(real, "register_p"),
            "decoy_either": rate(decoy), "decoy_topology": rate(decoy, "topology_p"),
            "decoy_register": rate(decoy, "register_p"),
            "n_decoys": len(decoy),
        },
        "decoy_classes": dict(Counter(r.get("architecture_class", "unmeasurable") for r in decoy)),
        "stratified": stratified(real, records, n_bootstrap),
        "per_lining": [
            {"id": record["id"], "group": record["group"],
             "label": record.get("label"), **row}
            for record, row in zip(records, real)
        ],
    }


def markdown(report: Dict[str, object]) -> str:
    rates = report["rates"]
    size = report["lining_size"]
    lines = [
        "# Study Q7 on real pocket linings",
        "",
        f"**Run.** `scripts/real_architecture.py run`, seed {report['seed']}, on "
        f"`{report['inputs']}`. Every figure is computed by that command and written in "
        "the same pass.",
        "",
        "**Pre-registered** in `docs/ARCHITECTURE_PLAN.md` as amended by "
        "`docs/ARCHITECTURE_PLAN_AMENDMENT_1.md`, which fixed what is computed here "
        "before any of it was read.",
        "",
        "**These are candidates, not annotated sites.** `annotated` is False for every "
        "row, so nothing below is evidence that annotated IP sites have a sequence-space "
        "architecture. It is a description of real candidate linings beside a null.",
        "",
        f"{report['n_measurable']} of {report['n_linings']} linings measurable; "
        f"size {size['min']}–{size['max']} positions, median {size['median']:.0f}.",
        "",
        "## Real linings against decoys matched on size and span",
        "",
        "| arm | topology | register | **either** |",
        "|---|---|---|---|",
        f"| real ({report['n_measurable']}) | {rates['real_topology']:.2f} | "
        f"{rates['real_register']:.2f} | **{rates['real_either']:.2f}** |",
        f"| decoy ({rates['n_decoys']}) | {rates['decoy_topology']:.2f} | "
        f"{rates['decoy_register']:.2f} | **{rates['decoy_either']:.2f}** |",
        "",
        "Each decoy has the same number of positions and the same span as the real "
        "lining it was drawn for, so the comparison is not confounded by either.",
        "",
        "## Architecture classes",
        "",
        "| class | real | decoy |",
        "|---|---|---|",
    ]
    classes = set(report["classes"]) | set(report["decoy_classes"])  # type: ignore[arg-type]
    for name in sorted(classes):
        real_n = report["classes"].get(name, 0)  # type: ignore[union-attr]
        decoy_n = report["decoy_classes"].get(name, 0)  # type: ignore[union-attr]
        lines.append(f"| `{name}` | {real_n} | {decoy_n} |")
    lines += [
        "",
        "Verdicts: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report["verdicts"].items()))  # type: ignore[union-attr]
        + ". `unresolved` is reported separately from `none`, as the plan requires: they "
        "are different claims.",
        "",
        "Secondary-structure elements called: "
        + ", ".join(f"{k} {v}" for k, v in sorted(report["elements"].items())),  # type: ignore[union-attr]
        "",
        "## The one pre-registered split",
        "",
        "Pockets a competitor polyanion already explains, against the rest. Group means "
        "over homology clusters, Holm across the features, five-group floor.",
        "",
        "| feature | arms | difference | 95% interval | p | Holm | evidence | finding |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for entry in report["stratified"]["features"]:  # type: ignore[index]
        if not entry.get("runnable"):
            lines.append(f"| `{entry['feature']}` | — | — | — | — | — | — | "
                         f"not runnable: {entry['reason']} |")
            continue
        low, high = entry["interval"]
        lines.append(
            f"| `{entry['feature']}` | {' vs '.join(entry['arms'])} | "
            f"{entry['difference']:+.3f} | [{low:+.3f}, {high:+.3f}] | "
            f"{entry['p_value']:.4f} | {entry['p_holm']:.4f} | "
            f"{'yes' if entry['is_evidence'] else 'NO'} | "
            f"{'**yes**' if entry['is_finding'] else 'no'} |"
        )
    lines += [
        "",
        "## What this shows, and the control it does not have",
        "",
        "The decoys control size and span by construction, so the gap between the arms "
        "is not either of those. What it is, though, is **nearly tautological**: a random "
        "set of positions in a span is not a protein, and *any* real pocket lining - IP "
        "site or not - comes from residues that fold together and will therefore look "
        "structured against that null. This result establishes that the instrument "
        "registers real linings and gives the rate at which it does. It does **not** "
        "establish anything specific to inositol-phosphate sites.",
        "",
        "The control that would is **other pockets from the same proteins**: the same "
        "folds, the same screen, the same models, differing only in being the pocket the "
        "screen did not pick. No committed table carries them, so that comparison is "
        "named here as missing rather than approximated.",
        "",
        "Three further limits:",
        "",
        "* Not evidence about annotated IP sites: every lining here is a candidate, and "
        "the screen that chose them already favoured buried, basic pockets, so the "
        "selection and the signal cannot be separated in this data.",
        "* Not a predictor. Separating classes here is not prediction until it is tested "
        "on groups held out from fitting, which this does not do.",
        "* Lining size is reported above because a difference that tracks it is not a "
        "finding, and the decoys control for it by construction.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--inputs", type=Path,
                        default=Path("results/architecture/real/architecture_inputs.jsonl"))
    runner.add_argument("--draws", type=int, default=1000)
    runner.add_argument("--bootstrap", type=int, default=2000)
    runner.add_argument("--out", type=Path, required=True)
    runner.add_argument("--markdown", type=Path, default=None)
    args = parser.parse_args(argv)

    report = run(args.inputs, n_draws=args.draws, n_bootstrap=args.bootstrap)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
