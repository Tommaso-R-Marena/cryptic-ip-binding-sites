#!/usr/bin/env python3
"""Do a pocket's basic residues sit on a face, more than chance within its own lining?

Pre-registered in `docs/ARCHITECTURE_PLAN_AMENDMENT_2.md`, dated before any
basic-residue subset was passed to `describe_architecture`.

The design is paired and within-lining, which is what makes it worth running.
Amendment 1's comparison - real linings against random position sets - is
nearly tautological, because a random position set is not a protein. Here each
pocket's basic residues are compared with **random subsets of that same
pocket's lining, of the same size**. Protein, fold, lining and subset size are
all held fixed by construction, so a difference can only be *which* residues
were chosen.

The per-pocket statistic is the quantile of the basic subset among its own
matched draws, which is uniform on [0, 1] under the null that basic residues
sit at random within the lining. Below 0.5 means more structured than chance.
"""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import architecture as ar  # noqa: E402
from cryptic_ip.detect import sequence as sq  # noqa: E402

SEED = 20261013

#: Matched draws per pocket. The quantile's resolution is 1/(draws + 1), so 50
#: draws resolve to 0.02 - finer than any effect this design could establish at
#: 75 pockets.
WITHIN_DRAWS = 50

#: Fewest residues in a subset worth testing. Below the architecture module's
#: own ordering floor a register null is unresolved, which is a property of the
#: data rather than a result.
MIN_SUBSET = 6

#: A quantile needs this many **resolved** matched draws to mean anything. A
#: half whose draws are mostly unresolved is reported as unresolved for that
#: pocket rather than scored against a handful of survivors.
MIN_RESOLVED_DRAWS = 10

#: The subsets tested, in order. Two and only two, fixed by the amendment.
SUBSETS = ("basic_positions", "conserved_positions")

#: The halves, reported separately and never pooled.
HALVES = ("topology", "register")


def parse_positions(value: object) -> List[int]:
    """A list cell (``[1, 2]``) or a delimited one (``1,2``)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        return sorted({int(x) for x in ast.literal_eval(text)})
    return sorted({int(x) for x in text.replace(";", ",").split(",") if x.strip()})


def halves_of(positions: Sequence[int], *, n_draws: int,
              seed: int) -> Dict[str, Optional[float]]:
    """Both halves' Holm-corrected p, with **None for a half that is unresolved**.

    An unresolved half is not a p of 1. `architecture.py` separates the two
    deliberately - a lining whose gap multiset has too few orderings cannot be
    tested, which is a different claim from "no structure" - and an earlier
    version of this script collapsed them by substituting 1.0. The effect was
    not subtle: half of all register quantiles came back at exactly 1.000 and
    none at 0.000, because an unresolved subset tied with its equally
    unresolved draws. The distinction is kept here.
    """
    try:
        arch = ar.describe_architecture(positions, n_draws=n_draws, seed=seed)
    except ar.ArchitectureError:
        return {"topology": None, "register": None}
    return {
        "topology": None if arch.topology_p is None else float(arch.topology_p),
        "register": None if arch.register_p is None else float(arch.register_p),
    }


def quantile_within_lining(
    lining: Sequence[int],
    subset: Sequence[int],
    *,
    n_draws: int,
    draws: int,
    rng: np.random.Generator,
) -> Optional[Dict[str, float]]:
    """Where the subset's p sits among random subsets of the same lining.

    A lower p is more structured, so the quantile is the share of matched draws
    whose p is **at or below** the subset's. Under the null that the subset is
    an arbitrary choice of residues from the lining it is uniform on [0, 1].
    """
    observed = halves_of(subset, n_draws=n_draws, seed=SEED)
    if all(value is None for value in observed.values()):
        return None
    pool = sorted({int(p) for p in lining})
    size = len({int(p) for p in subset})
    if size >= len(pool):
        return None
    out: Dict[str, float] = {}
    collected: Dict[str, List[float]] = {half: [] for half in HALVES}
    for index in range(draws):
        drawn = sorted(rng.choice(pool, size=size, replace=False).tolist())
        values = halves_of(drawn, n_draws=n_draws, seed=SEED + index)
        for half in HALVES:
            if values[half] is not None:
                collected[half].append(values[half])
    for half in HALVES:
        if observed[half] is None:
            continue
        values = np.asarray(collected[half], dtype=float)
        if values.size < MIN_RESOLVED_DRAWS:
            continue
        # Mid-p quantile. Permutation p-values are discrete, so ties are common
        # and counting them as "at or below" pushes every quantile upward;
        # splitting them is the standard correction and makes the statistic
        # uniform under the null rather than merely bounded by it.
        below = float(np.sum(values < observed[half] - 1e-12))
        tied = float(np.sum(np.abs(values - observed[half]) <= 1e-12))
        out[f"{half}_p"] = float(observed[half])
        out[f"{half}_quantile"] = float((below + 0.5 * tied) / values.size)
        out[f"{half}_tied_fraction"] = float(tied / values.size)
        out[f"{half}_draws"] = float(values.size)
    if not any(key.endswith("_quantile") for key in out):
        return None
    out["subset_size"] = float(size)
    out["lining_size"] = float(len(pool))
    return out


def run(table: Path, *, n_draws: int, draws: int, n_bootstrap: int) -> Dict[str, object]:
    frame = pd.read_csv(table)
    required = {"top_pocket_residues", "cluster", *SUBSETS}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise SystemExit(f"{table} is missing {missing}")
    rng = np.random.default_rng(SEED)

    per_pocket: List[Dict[str, object]] = []
    excluded: Dict[str, int] = {}

    def exclude(reason: str) -> None:
        excluded[reason] = excluded.get(reason, 0) + 1

    for _index, row in frame.iterrows():
        lining = parse_positions(row["top_pocket_residues"])
        for name in SUBSETS:
            subset = parse_positions(row[name])
            if len(subset) < MIN_SUBSET:
                exclude(f"{name}: fewer than {MIN_SUBSET} residues")
                continue
            if not set(subset) <= set(lining):
                exclude(f"{name}: not a subset of the lining")
                continue
            result = quantile_within_lining(
                lining, subset, n_draws=n_draws, draws=draws, rng=rng)
            if result is None:
                exclude(f"{name}: unmeasurable subset or lining")
                continue
            per_pocket.append({
                "id": str(row.get("uniprot_id", _index)),
                "group": str(row["cluster"]),
                "subset": name,
                **result,
            })

    tests, p_values = [], {}
    for name in SUBSETS:
        for half in HALVES:
            rows = [r for r in per_pocket if r["subset"] == name]
            if not rows:
                tests.append({"subset": name, "half": half, "runnable": False,
                              "reason": "no measurable pocket"})
                continue
            by_group: Dict[str, List[float]] = {}
            for row in rows:
                if f"{half}_quantile" not in row:
                    continue
                by_group.setdefault(str(row["group"]), []).append(
                    float(row[f"{half}_quantile"]))
            if not by_group:
                tests.append({"subset": name, "half": half, "runnable": False,
                              "reason": "every pocket is unresolved for this half"})
                continue
            observed = {g: float(np.mean(v)) for g, v in by_group.items()}
            # The null is a uniform quantile with mean 0.5, so the comparison
            # arm is a constant 0.5 per group rather than a second sample.
            null_arm = {g: 0.5 for g in observed}
            result = sq.compare_feature(f"{name}:{half}", observed, null_arm,
                                        n_bootstrap=n_bootstrap, seed=SEED)
            entry = result.as_dict()
            resolved = [r for r in rows if f"{half}_quantile" in r]
            entry.update({"subset": name, "half": half, "runnable": True,
                          "mean_quantile": float(np.mean(list(observed.values()))),
                          "n_pockets": len(resolved),
                          "mean_tied_fraction": float(np.mean(
                              [r.get(f"{half}_tied_fraction", 0.0) for r in resolved]))})
            tests.append(entry)
            p_values[f"{name}:{half}"] = result.p_value
    corrected = sq.holm(p_values) if p_values else {}
    for entry in tests:
        if entry.get("runnable"):
            entry["p_holm"] = corrected.get(f"{entry['subset']}:{entry['half']}")
            entry["is_finding"] = bool(
                entry["p_holm"] is not None and entry["p_holm"] < 0.05
                and entry["is_evidence"] and entry["exceeds_power"]
                and not (entry["interval"][0] <= 0.0 <= entry["interval"][1]))
    return {
        "seed": SEED,
        "table": str(table),
        "within_draws": draws,
        "min_subset": MIN_SUBSET,
        "n_pockets_in_table": int(len(frame)),
        "n_measured": len(per_pocket),
        "excluded": dict(sorted(excluded.items())),
        "subset_sizes": {
            name: {
                "median": float(np.median([r["subset_size"] for r in per_pocket
                                           if r["subset"] == name] or [float("nan")])),
                "median_lining": float(np.median([r["lining_size"] for r in per_pocket
                                                  if r["subset"] == name] or [float("nan")])),
            } for name in SUBSETS
        },
        "tests": tests,
        "holm_family_size": len(p_values),
        "per_pocket": per_pocket,
    }


def markdown(report: Dict[str, object]) -> str:
    lines = [
        "# Do a pocket's basic residues sit on a face?",
        "",
        f"**Run.** `scripts/basic_face.py run`, seed {report['seed']}, on "
        f"`{report['table']}`. Every figure is computed by that command in the same pass.",
        "",
        "**Pre-registered** in `docs/ARCHITECTURE_PLAN_AMENDMENT_2.md`, dated before any "
        "basic-residue subset was measured.",
        "",
        "## The design, and why it is better than amendment 1's",
        "",
        "Amendment 1 compared real linings with random position sets and said plainly "
        "that the comparison is nearly tautological. Here each pocket's basic residues "
        "are compared with **random subsets of that same pocket's lining, of the same "
        "size**. Protein, fold, lining and subset size are held fixed by construction, "
        "so a difference can only be which residues were chosen.",
        "",
        f"The per-pocket statistic is the quantile of the subset among "
        f"{report['within_draws']} matched draws, uniform on [0, 1] under the null that "
        "basic residues sit at random within the lining. **Below 0.5 is more structured "
        "than chance.**",
        "",
        f"{report['n_measured']} subset measurements from "
        f"{report['n_pockets_in_table']} pockets.",
    ]
    if report["excluded"]:
        lines += ["", "Excluded, by cause:", ""]
        for reason, count in report["excluded"].items():  # type: ignore[union-attr]
            lines.append(f"- {count} — {reason}")
    lines += [
        "",
        "## Result",
        "",
        "An **unresolved** half is not scored as a p of 1. The architecture module "
        "separates \"cannot be tested\" from \"no structure\", and an earlier version of "
        "this script collapsed them by substituting 1.0 — which put half of all register "
        "quantiles at exactly 1.000 and none at 0.000, because an unresolved subset tied "
        "with its equally unresolved draws. The quantile is also a **mid-p** quantile, "
        "splitting ties, because permutation p-values are discrete and counting ties as "
        "\"at or below\" pushes every quantile upward.",
        "",
        "| subset | half | pockets | ties | mean quantile | difference from 0.5 | 95% interval | p | Holm | finding |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for entry in report["tests"]:  # type: ignore[union-attr]
        if not entry.get("runnable"):
            lines.append(f"| `{entry['subset']}` | {entry['half']} | — | — | — | — | — | — | — | "
                         f"not runnable: {entry['reason']} |")
            continue
        low, high = entry["interval"]
        lines.append(
            f"| `{entry['subset']}` | {entry['half']} | {entry['n_pockets']} | "
            f"{entry.get('mean_tied_fraction', float('nan')):.2f} | "
            f"{entry['mean_quantile']:.3f} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | {entry['p_value']:.4f} | "
            f"{entry['p_holm']:.4f} | {'**yes**' if entry['is_finding'] else 'no'} |"
        )
    lines += [
        "",
        f"Holm across {report['holm_family_size']} tests: two subsets times two halves, "
        "fixed by the amendment before any of it was read.",
        "",
        "## How to read a null here",
        "",
        "The transfer measurement in `results/architecture/ARCHITECTURE.md` stands over "
        "all of this: at the φ/ψ scatter real structures show, the method detects about "
        "**a third** of what it detects on idealised geometry. A null in this table is "
        "therefore weak evidence of absence, and the honest statement is that the design "
        "could not see an effect of the size it could see on ideal geometry.",
        "",
        "## What this cannot conclude",
        "",
        "* These are candidates, not annotated sites.",
        "* Not evidence of binding: a basic face is a suggestive feature, not an "
        "interaction.",
        "* Nothing about whether these are real pockets — the design conditions on the "
        "lining and therefore cannot ask that.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    runner = sub.add_parser("run")
    runner.add_argument("--table", type=Path,
                        default=Path("results/triage/triaged_candidates.csv"))
    runner.add_argument("--draws", type=int, default=WITHIN_DRAWS)
    runner.add_argument("--null-draws", type=int, default=400)
    runner.add_argument("--bootstrap", type=int, default=2000)
    runner.add_argument("--out", type=Path, required=True)
    runner.add_argument("--markdown", type=Path, default=None)
    args = parser.parse_args(argv)

    report = run(args.table, n_draws=args.null_draws, draws=args.draws,
                 n_bootstrap=args.bootstrap)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
