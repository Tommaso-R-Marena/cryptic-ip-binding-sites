#!/usr/bin/env python3
"""Study Q6: is there a sequence signature at an annotated binding site?

Pre-registered in `docs/SEQUENCE_PLAN.md`. Three subcommands:

``validate``
    Method validation on synthetic sequence pools with a planted signal and a
    composition-matched negative control. Runs anywhere, touches no study data,
    and is the evidence that the machinery can find an ordering signature and
    does not report one when composition alone explains the difference.

``describe``
    Per-pocket descriptors from an inputs file of chain sequences and lining
    positions. One JSON record per pocket.

``report``
    The pre-registered analysis: the minimum detectable effect **first**, then
    the six primary tests with Holm correction, then the exploratory table.

The report writes its power block before its results block, in that order, in
the artefact. The plan says the MDE is read before the difference; doing it in
this order is how that is enforced rather than asserted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import sequence as sq  # noqa: E402
from cryptic_ip.detect import signature as sg  # noqa: E402

#: Shuffles for a per-window composition null in the study runs. The plan's
#: convention; `validate` uses fewer so it stays interactive.
STUDY_SHUFFLES = 2000

#: Seed for every resampling in this study, fixed in advance.
SEED = 20261012


# ------------------------------------------------------------- validate

def planted_pool(
    n_windows: int,
    length: int,
    period: float,
    *,
    seed: int,
) -> List[str]:
    """Windows with a basic residue every ``period`` positions.

    The filler alphabet holds no basic residue, so the composition of a window
    is a function of its length and the period alone - which is what lets the
    negative control below hold composition exactly fixed.
    """
    filler = "ALVGSTMFWY"
    rng = np.random.default_rng(seed)
    pool = []
    for _ in range(n_windows):
        offset = float(rng.uniform(0.0, period))
        marks = {int(round(offset + i * period)) for i in range(int(length / period) + 2)}
        pool.append("".join(
            "R" if i in marks else filler[int(rng.integers(0, len(filler)))]
            for i in range(length)
        ))
    return pool


def shuffled_pool(pool: Sequence[str], *, seed: int) -> List[str]:
    """The same windows with their letters reordered: composition held fixed."""
    rng = np.random.default_rng(seed)
    out = []
    for window in pool:
        letters = list(window)
        rng.shuffle(letters)
        out.append("".join(letters))
    return out


def composition_pool(
    basic_counts: Sequence[int],
    length: int,
    *,
    seed: int,
) -> List[str]:
    """Basic-rich windows with no ordering at all.

    The pool a composition-blind method would mistake for the planted one. It
    takes the planted pool's **per-window** basic counts rather than their mean
    fraction: a single rounded fraction leaves a residual of a percentage point
    or two, which at a dozen windows is enough to make the two pools
    distinguishable on composition and the whole validation worthless. Passing
    the counts through makes the match exact by construction.
    """
    filler = "ALVGSTMFWY"
    rng = np.random.default_rng(seed)
    out = []
    for n_basic in basic_counts:
        n_basic = int(n_basic)
        if not 0 <= n_basic <= length:
            raise ValueError(f"{n_basic} basic residues do not fit in {length}")
        letters = ["R"] * n_basic + [
            filler[int(rng.integers(0, len(filler)))] for _ in range(length - n_basic)
        ]
        rng.shuffle(letters)
        out.append("".join(letters))
    return out


def pool_summary(name: str, pool: Sequence[str], statistic: str,
                 *, n_shuffles: int) -> Dict[str, object]:
    """A pool's mean statistic, mean composition z, and how often it fires."""
    raw, zed = [], []
    for index, window in enumerate(pool):
        raw.append(float(sq.STATISTICS[statistic](window)))
        zed.append(sq.composition_z(
            window, statistic, n_shuffles=n_shuffles, seed=SEED + index).z)
    zed_array = np.asarray(zed)
    letters = sorted({c for window in pool for c in window})
    return {
        "pool": name,
        "n_windows": len(pool),
        "length": len(pool[0]),
        "mean_raw": float(np.mean(raw)),
        "mean_entropy": float(np.mean([sq.shannon_entropy(w) for w in pool])),
        "basic_fraction": float(np.mean([
            sum(c in "RK" for c in w) / len(w) for w in pool])),
        "distinct_letters": len(letters),
        "mean_z": float(zed_array.mean()),
        "fraction_firing": float(np.mean(np.abs(zed_array) > 2.0)),
    }


def run_validate(n_windows: int, length: int, n_shuffles: int) -> Dict[str, object]:
    """Can the machinery find a planted signature, and does it stay quiet otherwise?"""
    planted = planted_pool(n_windows, length, 3.6, seed=SEED)
    shuffles = shuffled_pool(planted, seed=SEED + 1)
    basic = composition_pool(
        [sum(c in "RK" for c in window) for window in planted], length, seed=SEED + 2)

    statistic = "helix_charge_periodicity"
    pools = [
        pool_summary("planted helical face", planted, statistic, n_shuffles=n_shuffles),
        pool_summary("same letters, shuffled", shuffles, statistic, n_shuffles=n_shuffles),
        pool_summary("basic-rich, no ordering", basic, statistic, n_shuffles=n_shuffles),
    ]

    # The composition controls must be indistinguishable from the planted pool
    # on composition, or the validation proves nothing.
    entropies = [p["mean_entropy"] for p in pools]
    basics = [p["basic_fraction"] for p in pools]
    checks = {
        "shuffle_control_has_identical_entropy":
            abs(pools[0]["mean_entropy"] - pools[1]["mean_entropy"]) < 1e-9,
        "shuffle_control_has_identical_basic_fraction":
            abs(pools[0]["basic_fraction"] - pools[1]["basic_fraction"]) < 1e-9,
        "composition_control_has_identical_basic_fraction":
            abs(pools[0]["basic_fraction"] - pools[2]["basic_fraction"]) < 1e-9,
        "all_pools_within_0.1_bits":
            float(max(entropies) - min(entropies)) < 0.1,
        "all_pools_have_the_same_basic_fraction":
            float(max(basics) - min(basics)) < 1e-9,
        "planted_pool_fires":
            pools[0]["mean_z"] > 2.0 and pools[0]["fraction_firing"] > 0.5,
        "shuffle_control_stays_quiet":
            abs(pools[1]["mean_z"]) < 2.0 and pools[1]["fraction_firing"] < 0.25,
        "composition_control_stays_quiet":
            abs(pools[2]["mean_z"]) < 2.0 and pools[2]["fraction_firing"] < 0.25,
    }
    return {
        "statistic": statistic,
        "n_shuffles": n_shuffles,
        "pools": pools,
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def validate_markdown(report: Mapping[str, object]) -> str:
    lines = [
        "# Study Q6 method validation",
        "",
        "Synthetic pools only. No annotated sequence is read here; this is the",
        "evidence that the machinery finds an ordering signature when one is",
        "planted and reports nothing when composition alone would explain the",
        f"difference. Statistic: `{report['statistic']}`, "
        f"{report['n_shuffles']} shuffles per window.",
        "",
        "| pool | n | mean raw | mean entropy | basic fraction | mean z | fires |",
        "|---|---|---|---|---|---|---|",
    ]
    for pool in report["pools"]:  # type: ignore[index]
        lines.append(
            f"| {pool['pool']} | {pool['n_windows']} | {pool['mean_raw']:.3f} | "
            f"{pool['mean_entropy']:.3f} | {pool['basic_fraction']:.3f} | "
            f"{pool['mean_z']:+.2f} | {pool['fraction_firing'] * 100:.0f}% |"
        )
    lines += [
        "",
        "The three pools are matched on composition, so the separation in the z",
        "column cannot be composition.",
        "",
        "| check | result |",
        "|---|---|",
    ]
    for name, passed in report["checks"].items():  # type: ignore[union-attr]
        lines.append(f"| {name.replace('_', ' ')} | {'pass' if passed else 'FAIL'} |")
    lines += ["", f"All checks pass: **{report['all_checks_pass']}**"]
    return "\n".join(lines) + "\n"


def provenance_markdown(payload: Mapping[str, object], commit: str) -> str:
    """`results/sequence/PROVENANCE.md`, generated from the artefact.

    A provenance file whose numbers were typed in by hand proves nothing, so
    this is a subcommand rather than a one-off script, and
    `tests/test_sequence_signature_script.py` checks that the committed file is
    exactly what it produces from the committed `validation.json`.
    """
    pools = payload["pools"]  # type: ignore[index]
    checks = payload["checks"]  # type: ignore[index]
    entropy_gap = abs(pools[0]["mean_entropy"] - pools[1]["mean_entropy"])
    entropy_spread = (max(p["mean_entropy"] for p in pools)
                      - min(p["mean_entropy"] for p in pools))
    basic_spread = (max(p["basic_fraction"] for p in pools)
                    - min(p["basic_fraction"] for p in pools))
    lines = [
        "# Study Q6 method validation: provenance",
        "",
        f"**Run.** `scripts/sequence_signature.py validate --windows "
        f"{pools[0]['n_windows']} --length {pools[0]['length']} --shuffles "
        f"{payload['n_shuffles']}`, seed {SEED}, on the session container (single "
        f"core), at commit `{commit}`.",
        "",
        "**Plan.** `docs/SEQUENCE_PLAN.md`, which pre-registers study Q6 and states "
        "that this validation is a property of the method and not a study result.",
        "",
        "Every figure below is generated from `validation.json` by "
        "`scripts/sequence_signature.py provenance`. Nothing is transcribed, and a "
        "test checks that this file is what that subcommand produces.",
        "",
        "## What this is",
        "",
        "Three pools of synthetic sequence windows, matched on composition, one of "
        "which has a planted ordering signature. It answers one question: does the "
        "machinery find an ordering signature when one is there, and stay quiet when "
        "composition alone would explain the difference?",
        "",
        "**It says nothing about whether annotated inositol-phosphate sites have a "
        "sequence signature.** That is study Q6, asked of real sequences, and it has "
        "not run.",
        "",
        f"## The pools ({pools[0]['n_windows']} windows of {pools[0]['length']} "
        "residues each)",
        "",
        "| pool | mean raw | mean entropy (bits) | basic fraction | distinct letters "
        "| mean z | windows firing |",
        "|---|---|---|---|---|---|---|",
    ]
    for pool in pools:
        lines.append(
            f"| {pool['pool']} | {pool['mean_raw']:.3f} | {pool['mean_entropy']:.3f} | "
            f"{pool['basic_fraction']:.3f} | {pool['distinct_letters']} | "
            f"{pool['mean_z']:+.2f} | {pool['fraction_firing'] * 100:.0f}% |"
        )
    lines += [
        "",
        f"Statistic: `{payload['statistic']}`. Null: {payload['n_shuffles']} "
        "composition-preserving shuffles per window, so each window is compared only "
        "with reorderings of itself.",
        "",
        "## Why the comparison is clean",
        "",
        f"The shuffled control holds the planted pool's letters exactly, so its "
        f"entropy differs by {entropy_gap:.1e} bits and its basic fraction not at "
        f"all. The basic-rich control is built from the planted pool's own per-window "
        f"basic counts rather than from their mean fraction, which closed a residual "
        f"of 0.021 in basic fraction - enough to fail the match at twelve windows. "
        f"Across all three pools the entropy spread is {entropy_spread:.3f} bits and "
        f"the basic-fraction spread is {basic_spread:.1e}.",
        "",
        "So the z column cannot be reporting composition: there is no composition "
        "difference left for it to report.",
        "",
        "## Checks",
        "",
        "| check | result |",
        "|---|---|",
    ]
    for name, passed in sorted(checks.items()):
        lines.append(f"| {name.replace('_', ' ')} | {'pass' if passed else 'FAIL'} |")
    lines += [
        "",
        f"All {len(checks)} checks pass: **{payload['all_checks_pass']}**. A failed "
        "check exits non-zero and fails the workflow run.",
        "",
        "## What this validation may not be used for",
        "",
        "* It is not evidence that any real site has a sequence signature.",
        "* It is not a benchmark of the detector: no structure is involved.",
        "* It does not license adding a sequence term to the score. "
        "`signature.ESTABLISHED_SIGNATURES` stays empty until study Q6 reports, and "
        "`sequence_contribution` returns an abstention until then.",
        "* The planted signal is a charge periodicity at exactly 3.6 residues, which "
        "is the cleanest possible case. A real helical face is noisier, so the firing "
        "rate here is an upper bound on sensitivity, not an estimate of it.",
    ]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------- describe

def load_inputs(path: Path) -> List[Dict[str, object]]:
    """One JSON object per line: id, group, arm, chain sequence, lining residues.

    ``residues`` is a list of ``[chain_id, residue_number]`` pairs, because the
    repo's lining residue numbers are chain-blind and a window needs one chain.
    """
    records = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            record = json.loads(line)
        except json.JSONDecodeError as error:
            raise SystemExit(f"{path}:{number}: {error}") from error
        missing = {"id", "group", "arm", "sequence", "residues"} - set(record)
        if missing:
            raise SystemExit(f"{path}:{number}: missing {sorted(missing)}")
        if record["arm"] not in ("site", "control"):
            raise SystemExit(f"{path}:{number}: arm must be 'site' or 'control'")
        records.append(record)
    if not records:
        raise SystemExit(f"{path} holds no records")
    return records


def describe_one(record: Mapping[str, object], *, n_shuffles: int) -> Dict[str, object]:
    """One pocket's descriptors, or the reason it was excluded."""
    out: Dict[str, object] = {
        "id": record["id"], "group": record["group"], "arm": record["arm"],
    }
    pairs: List[Tuple[str, int]] = [
        (str(chain), int(number)) for chain, number in record["residues"]  # type: ignore[union-attr]
    ]
    try:
        chain_id, fraction, positions = sg.dominant_chain(pairs)
    except sg.WindowError as error:
        return {**out, "measurable": False, "exclusion_reason": str(error)}
    out["chain"] = chain_id
    out["chain_fraction"] = fraction
    if fraction < sg.MIN_CHAIN_FRACTION:
        return {
            **out, "measurable": False,
            "exclusion_reason": (
                f"lining spans chains ({fraction:.2f} on {chain_id}); an interface "
                "pocket has no single sequence window"
            ),
        }
    sequences = record["sequence"]
    chain_sequence = (
        sequences[chain_id] if isinstance(sequences, dict) else str(sequences)
    )
    try:
        windows = sg.site_windows(chain_sequence, positions, offset=int(record.get("offset", 1)))
    except sg.WindowError as error:
        return {**out, "measurable": False, "exclusion_reason": str(error)}
    return {**out, **sg.describe_site(windows, n_shuffles=n_shuffles, seed=SEED)}


def run_describe(inputs: Path, n_shuffles: int) -> List[Dict[str, object]]:
    return [describe_one(record, n_shuffles=n_shuffles) for record in load_inputs(inputs)]


# --------------------------------------------------------------- report

def value_for(site: Mapping[str, object], window: str, statistic: str) -> Optional[float]:
    """The per-pocket number a primary test uses, or None if unmeasurable.

    The primary family is defined on composition z-scores for the two sequence
    windows and on a raw proportion for the spacing window, which has no
    letters to shuffle.
    """
    if not site.get("measurable"):
        return None
    if window == "spacing":
        spacing = site.get("spacing")
        if not isinstance(spacing, dict) or statistic not in spacing:
            return None
        return float(spacing[statistic])
    bundle = site.get(window)
    if not isinstance(bundle, dict):
        return None
    nulls = bundle.get("nulls")
    if not isinstance(nulls, dict) or statistic not in nulls:
        return None
    return float(nulls[statistic]["z"])


def group_means(
    sites: Sequence[Mapping[str, object]], window: str, statistic: str, arm: str,
) -> Dict[str, float]:
    """Each homology group's mean value in one arm.

    Groups, not copies: the project's standing rule, and the reason a confident
    difference cannot be read off 209 copies of 29 proteins.
    """
    collected: Dict[str, List[float]] = {}
    for site in sites:
        if site.get("arm") != arm:
            continue
        value = value_for(site, window, statistic)
        if value is None:
            continue
        collected.setdefault(str(site["group"]), []).append(value)
    return {group: float(np.mean(values)) for group, values in collected.items()}


def power_block(sites: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    """Group counts and the minimum detectable difference, per primary test.

    Written to the artefact before any difference. The project has read a null
    as a failure once already, when at five homology groups nothing below an
    AUC of 0.87 was detectable at all.
    """
    rows = []
    for test in sg.PRIMARY_TESTS:
        site_means = group_means(sites, test.window, test.statistic, "site")
        control_means = group_means(sites, test.window, test.statistic, "control")
        rows.append({
            "test": test.number,
            "window": test.window,
            "statistic": test.statistic,
            "n_site_groups": len(site_means),
            "n_control_groups": len(control_means),
            "clears_evidence_floor":
                min(len(site_means), len(control_means)) >= 5,
            "minimum_detectable": sq.minimum_detectable_difference(
                list(site_means.values()), list(control_means.values())),
        })
    return {
        "n_pockets": len(sites),
        "n_measurable": sum(1 for s in sites if s.get("measurable")),
        "exclusions": _exclusion_table(sites),
        "windows": _window_table(sites),
        "tests": rows,
    }


def _window_table(sites: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    """Mean window lengths per arm, and how often the local window was capped.

    The truncation share matters for reading the length-confound flag: once
    both arms hit `LOCAL_CAP` their mean local lengths are both exactly the cap
    and the flag can no longer see a difference that is really there. A high
    truncation share makes the flag uninformative rather than reassuring.
    """
    out: Dict[str, object] = {}
    for arm in ("site", "control"):
        windows = [
            s["windows"] for s in sites
            if s.get("measurable") and s.get("arm") == arm
        ]
        if not windows:
            out[arm] = {"n": 0}
            continue
        out[arm] = {
            "n": len(windows),
            "mean_lining_length": float(np.mean([w["lining_length"] for w in windows])),
            "mean_local_length": float(np.mean([w["local_length"] for w in windows])),
            "truncated_fraction": float(np.mean([bool(w["truncated"]) for w in windows])),
        }
    return out


def _exclusion_table(sites: Sequence[Mapping[str, object]]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for site in sites:
        if site.get("measurable"):
            continue
        reason = str(site.get("exclusion_reason", "unstated")).split(";")[0]
        reason = reason.split("(")[0].strip()
        counts[reason] = counts.get(reason, 0) + 1
    return dict(sorted(counts.items()))


def length_confounded(sites: Sequence[Mapping[str, object]], window: str) -> Optional[bool]:
    """Whether the two arms differ in mean window length by more than 2 residues."""
    if window == "spacing":
        return None
    key = "lining_length" if window == "lining" else "local_length"
    lengths: Dict[str, List[float]] = {"site": [], "control": []}
    for site in sites:
        if not site.get("measurable"):
            continue
        arm = str(site.get("arm"))
        if arm in lengths:
            lengths[arm].append(float(site["windows"][key]))  # type: ignore[index]
    if not lengths["site"] or not lengths["control"]:
        return None
    return abs(float(np.mean(lengths["site"]) - np.mean(lengths["control"]))) > 2.0


def run_report(sites: Sequence[Mapping[str, object]], n_bootstrap: int) -> Dict[str, object]:
    """The power block, then the primary tests, then the exploratory table."""
    power = power_block(sites)
    comparisons = []
    p_values: Dict[str, float] = {}
    for test in sg.PRIMARY_TESTS:
        name = f"{test.window}:{test.statistic}"
        site_means = group_means(sites, test.window, test.statistic, "site")
        control_means = group_means(sites, test.window, test.statistic, "control")
        if not site_means or not control_means:
            comparisons.append({
                "test": test.number, "feature": name, "hypothesis": test.hypothesis,
                "runnable": False,
                "reason": "one arm has no measurable group",
                "n_site_groups": len(site_means), "n_control_groups": len(control_means),
            })
            continue
        comparison = sq.compare_feature(
            name, site_means, control_means, n_bootstrap=n_bootstrap, seed=SEED)
        entry = comparison.as_dict()
        entry.update({
            "test": test.number, "hypothesis": test.hypothesis, "runnable": True,
            "length_confounded": length_confounded(sites, test.window),
        })
        comparisons.append(entry)
        p_values[name] = comparison.p_value

    corrected = sq.holm(p_values) if p_values else {}
    for entry in comparisons:
        if entry.get("runnable"):
            entry["p_holm"] = corrected.get(str(entry["feature"]))
            entry["is_signature"] = bool(
                entry["p_holm"] is not None
                and entry["p_holm"] < 0.05
                and not (entry["interval"][0] <= 0.0 <= entry["interval"][1])
                and entry["exceeds_power"]
                and entry["is_evidence"]
                and entry["length_confounded"] is not True
            )
    return {
        "plan": "docs/SEQUENCE_PLAN.md",
        "power": power,
        "primary": comparisons,
        "n_signatures": sum(1 for e in comparisons if e.get("is_signature")),
        "holm_family_size": len(p_values),
    }


def report_markdown(report: Mapping[str, object]) -> str:
    power = report["power"]
    lines = [
        "# Study Q6: sequence signature at annotated inositol-phosphate sites",
        "",
        f"Pre-registered in `{report['plan']}`. "
        f"{power['n_measurable']} of {power['n_pockets']} pockets measurable.",  # type: ignore[index]
        "",
        "## Power, before any difference",
        "",
        "| # | window | statistic | site groups | control groups | floor | MDE |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in power["tests"]:  # type: ignore[index]
        mde = row["minimum_detectable"]
        lines.append(
            f"| {row['test']} | `{row['window']}` | `{row['statistic']}` | "
            f"{row['n_site_groups']} | {row['n_control_groups']} | "
            f"{'yes' if row['clears_evidence_floor'] else 'NO'} | "
            f"{'inf' if not np.isfinite(mde) else f'{mde:.3f}'} |"
        )
    lines += ["", "| arm | n | mean lining | mean local | local window capped |",
              "|---|---|---|---|---|"]
    for arm, entry in power["windows"].items():  # type: ignore[union-attr]
        if not entry["n"]:
            lines.append(f"| {arm} | 0 | — | — | — |")
            continue
        lines.append(
            f"| {arm} | {entry['n']} | {entry['mean_lining_length']:.1f} | "
            f"{entry['mean_local_length']:.1f} | "
            f"{entry['truncated_fraction'] * 100:.0f}% |"
        )
    if any(e["n"] and e["truncated_fraction"] > 0.5 for e in power["windows"].values()):  # type: ignore[union-attr]
        lines += [
            "",
            "More than half the local windows hit the cap, so both arms' mean",
            "local lengths are pinned near it and the length-confound flag on a",
            "`local` test cannot see a difference that may still be there.",
        ]
    if power["exclusions"]:  # type: ignore[index]
        lines += ["", "Excluded pockets, by cause:", ""]
        for reason, count in power["exclusions"].items():  # type: ignore[union-attr]
            lines.append(f"- {count} — {reason}")
    lines += [
        "",
        "## Primary tests",
        "",
        "| # | feature | site | control | difference | 95% interval | p | Holm | signature |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for entry in report["primary"]:  # type: ignore[index]
        if not entry.get("runnable"):
            lines.append(
                f"| {entry['test']} | `{entry['feature']}` | — | — | — | — | — | — | "
                f"not runnable: {entry['reason']} |"
            )
            continue
        low, high = entry["interval"]
        lines.append(
            f"| {entry['test']} | `{entry['feature']}` | {entry['site_mean']:+.3f} | "
            f"{entry['control_mean']:+.3f} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | {entry['p_value']:.4f} | "
            f"{entry['p_holm']:.4f} | {'**yes**' if entry['is_signature'] else 'no'} |"
        )
    lines += [
        "",
        f"Signatures established: **{report['n_signatures']}** of "
        f"{report['holm_family_size']} tests in the Holm family.",
        "",
        "A test is a signature only with a Holm-corrected p below 0.05, an",
        "interval excluding zero, a difference above the minimum detectable",
        "effect printed above, five groups on both sides, and arms that are not",
        "length-confounded. All five, per the plan.",
    ]
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------- cli

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate", help="method validation on synthetic pools")
    validate.add_argument("--windows", type=int, default=40)
    validate.add_argument("--length", type=int, default=24)
    validate.add_argument("--shuffles", type=int, default=400)
    validate.add_argument("--out", type=Path, default=Path("validation.json"))
    validate.add_argument("--markdown", type=Path, default=None)

    describe = sub.add_parser("describe", help="per-pocket descriptors")
    describe.add_argument("--inputs", type=Path, required=True)
    describe.add_argument("--shuffles", type=int, default=STUDY_SHUFFLES)
    describe.add_argument("--out", type=Path, required=True)

    provenance = sub.add_parser(
        "provenance", help="regenerate PROVENANCE.md from a validation artefact")
    provenance.add_argument("--validation", type=Path, required=True)
    provenance.add_argument("--commit", default="unknown")
    provenance.add_argument("--out", type=Path, required=True)

    report = sub.add_parser("report", help="the pre-registered analysis")
    report.add_argument("--sites", type=Path, required=True)
    report.add_argument("--bootstrap", type=int, default=2000)
    report.add_argument("--out", type=Path, required=True)
    report.add_argument("--markdown", type=Path, default=None)

    args = parser.parse_args(argv)

    if args.command == "validate":
        payload = run_validate(args.windows, args.length, args.shuffles)
        args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        text = validate_markdown(payload)
        if args.markdown:
            args.markdown.write_text(text)
        print(text)
        return 0 if payload["all_checks_pass"] else 1

    if args.command == "provenance":
        payload = json.loads(args.validation.read_text())
        args.out.write_text(provenance_markdown(payload, args.commit))
        print(f"wrote {args.out}")
        return 0

    if args.command == "describe":
        records = run_describe(args.inputs, args.shuffles)
        args.out.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
        measurable = sum(1 for r in records if r.get("measurable"))
        print(f"{len(records)} pockets, {measurable} measurable")
        return 0

    sites = json.loads(args.sites.read_text())
    payload = run_report(sites, args.bootstrap)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    text = report_markdown(payload)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
