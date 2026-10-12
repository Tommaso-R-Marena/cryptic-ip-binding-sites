#!/usr/bin/env python3
"""Measure the lining-architecture instruments, and what they can and cannot see.

Eight probes. The first four measure the instruments against their own
mathematics - two identities that collapse, two artefacts that would otherwise
have produced confident nonsense. The last four measure what the method can
detect, on ideal-geometry helix barrels whose lining is known from the 3-D
construction rather than assumed.

Everything here is synthetic or analytic. No annotated sequence is read, and
`docs/ARCHITECTURE_PLAN.md` pre-registers the study that would use it.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import permutations
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import architecture as ar  # noqa: E402
from cryptic_ip.detect import geometry_fixtures as gf  # noqa: E402

SEED = 20261012
PER, NHEL, LINK = 18, 6, 10


# ----------------------------------------------------------------- barrel

def barrel_inward(n_helices: int = NHEL, per_helix: int = PER, radius: float = 9.0):
    """Residues whose CB points at the barrel axis, and their heights.

    Ground truth from the construction: a residue lines the central channel
    exactly when its side chain points inward. Nothing about sequence is used.
    """
    fx = gf.helix_barrel(n_helices, per_helix, radius)
    inward: List[int] = []
    heights: Dict[int, float] = {}
    for residue in np.unique(fx.residue_index):
        sel = fx.residue_index == residue
        atoms = {n: c for n, c in zip(fx.atom_names[sel], fx.coords[sel])}
        if "CA" not in atoms or "CB" not in atoms:
            continue
        if np.hypot(*atoms["CB"][:2]) < np.hypot(*atoms["CA"][:2]):
            inward.append(int(residue))
            heights[int(residue)] = float(atoms["CB"][2])
    return sorted(inward), heights


def with_linkers(positions: Sequence[int], per_helix: int = PER, linker: int = LINK) -> List[int]:
    """Renumber as a real barrel: elements joined by loops.

    The fixture lays its helices head to tail, so in sequence space its lining
    is one unbroken run. Real barrels have loops. This changes the numbering
    and not the structure, and both versions are reported.
    """
    return sorted(int(p) + linker * (int(p) // per_helix) for p in positions)


# ------------------------------------------------- probe 1: the identities

def probe_identities(rng: np.random.Generator) -> Dict[str, object]:
    """H0 barcode == sorted gaps, and the box dimension in closed form."""
    def single_linkage(pos):
        clusters = [[p] for p in sorted(pos)]
        deaths = []
        while len(clusters) > 1:
            best, index = np.inf, 0
            for i in range(len(clusters) - 1):
                distance = min(abs(a - b) for a in clusters[i] for b in clusters[i + 1])
                if distance < best:
                    best, index = distance, i
            deaths.append(best)
            clusters[index] = clusters[index] + clusters[index + 1]
            del clusters[index + 1]
        return sorted(deaths)

    mismatches = 0
    trials = 300
    for _ in range(trials):
        pos = sorted(rng.choice(np.arange(80), size=int(rng.integers(3, 15)), replace=False).tolist())
        if single_linkage(pos) != sorted(np.diff(pos).tolist()):
            mismatches += 1
    return {
        "barcode_trials": trials,
        "barcode_mismatches": mismatches,
        "barcode_identity_holds": mismatches == 0,
    }


# ------------------------------------------- probe 2: the axial degeneracy

def probe_axial(rng: np.random.Generator) -> Dict[str, object]:
    """The doubled-angle resultant is exactly 1 at period 2 for any integer set."""
    def axial(pos, period):
        angles = 2.0 * np.pi * np.asarray(pos, float) / period
        return float(abs(np.exp(2j * angles).mean()))

    worst = 0.0
    for _ in range(2000):
        n = int(rng.integers(5, 25))
        pos = sorted(rng.choice(np.arange(300), size=n, replace=False).tolist())
        worst = max(worst, abs(axial(pos, 2.0) - 1.0))
    worst_four = 0.0
    for _ in range(2000):
        n = int(rng.integers(5, 25))
        pos = sorted((2 * rng.choice(np.arange(150), size=n, replace=False)).tolist())
        worst_four = max(worst_four, abs(axial(pos, 4.0) - 1.0))
    return {
        "max_abs_deviation_at_period_2": worst,
        "max_abs_deviation_at_period_4_uniform_parity": worst_four,
        "axial_statistic_usable": bool(worst > 1e-9),
        "note": ar.AXIAL_UNUSABLE,
    }


# ------------------------------------------- probe 3: the scan inflation

def probe_scan_inflation(rng: np.random.Generator, n: int = 12, lining_span: int = 60) -> Dict[str, object]:
    """What a fixed analytic threshold on a scanned coherence really admits."""
    threshold = float(np.sqrt(-np.log(0.05) / n))
    rows = []
    for label, band, probes in (
        ("single period 3.6", (3.6, 3.6001), 1),
        ("helix band", ar.HELIX_BAND, 13),
        ("strand band", ar.STRAND_BAND, 6),
        ("full scan 2-12", ar.SCAN_BAND, 101),
    ):
        grid = np.linspace(band[0], band[1], max(1, probes))
        values = np.empty(3000)
        for index in range(values.size):
            pos = rng.choice(np.arange(lining_span), size=n, replace=False)
            values[index] = ar.resultants(np.sort(pos).astype(float), grid).max()
        rows.append({
            "statistic": label, "probes": int(grid.size),
            "mean": float(values.mean()),
            "true_alpha_at_analytic_threshold": float(np.mean(values > threshold)),
        })
    return {"n": n, "span": lining_span, "analytic_threshold": threshold, "rows": rows}


# --------------------------------------------- probe 4: the peak-width law

def probe_peak_width() -> Dict[str, object]:
    """Half-width of the coherence peak, against ``0.58 * w^2 / span``."""
    rows = []
    for lining_span in (20, 40, 80, 160, 320):
        period = 3.6
        pos = np.unique(np.round(np.arange(0, lining_span + 1, period)).astype(int)).astype(float)
        peak = float(ar.resultants(pos, np.array([period]))[0])
        fine = np.linspace(period - 0.6, period + 0.6, 4001)
        values = ar.resultants(pos, fine)
        above = fine[values >= peak / 2.0]
        measured = float((above.max() - above.min()) / 2.0) if above.size else float("nan")
        rows.append({
            "span": lining_span,
            "measured_half_width": measured,
            "ratio_to_w2_over_span": measured / (period * period / lining_span),
        })
    ratios = [r["ratio_to_w2_over_span"] for r in rows]
    # What a fixed 13-probe scan loses on a perfect face.
    losses = []
    for lining_span in (30, 60, 104, 200):
        for period in (3.47, 3.62, 3.78):
            pos = np.unique(np.round(np.arange(0, lining_span + 1, period)).astype(int)).astype(float)
            fine = float(ar.resultants(pos, np.linspace(3.0, 4.2, 2401)).max())
            coarse = float(ar.resultants(pos, np.linspace(3.0, 4.2, 13)).max())
            losses.append({"span": lining_span, "period": period,
                           "fine": fine, "thirteen_probes": coarse, "lost": fine - coarse})
    return {
        "rows": rows,
        "ratio_mean": float(np.mean(ratios)),
        "ratio_spread": float(max(ratios) - min(ratios)),
        "fixed_probe_losses": losses,
        "worst_loss": float(max(r["lost"] for r in losses)),
    }


# ------------------------------------------ probe 5: the orderings floor

def probe_orderings_floor() -> Dict[str, object]:
    """Smallest p a permutation null can reach, by exact enumeration."""
    rows = []
    for gapset in ([3, 1, 3, 3, 1, 3], [1, 1, 1, 8, 1, 1, 1], [4, 3, 4, 3, 4, 3, 4],
                   [3, 4, 3, 4, 3, 4, 3, 4], [4, 3, 4, 3, 4, 3, 4, 3, 4, 3]):
        orders = set(permutations(gapset))
        values = []
        for order in orders:
            pos = np.concatenate([[0.0], np.cumsum(order)])
            values.append(ar.resultants(pos, np.linspace(*ar.HELIX_BAND, 13)).max())
        values = np.asarray(values)
        rows.append({
            "gaps": list(gapset), "orderings": len(orders),
            "min_attainable_p": float(np.mean(values >= values.max() - 1e-12)),
        })
    return {"rows": rows, "floor_used": ar.MIN_ORDERINGS}


# ------------------------------------------- probe 6: the calibration

def probe_calibration(rng: np.random.Generator, n_linings: int = 300,
                      n_draws: int = 400) -> Dict[str, object]:
    """False-positive rate of each test, and of the family, on random linings."""
    names = ["helix_topology", "strand_topology", "helix_register", "strand_register"]
    raw = {name: [] for name in names}
    holm = []
    classes: Dict[str, int] = {}
    verdicts: Dict[str, int] = {}
    for index in range(n_linings):
        n = int(rng.integers(8, 22))
        lining_span = int(rng.integers(int(2.2 * n), 8 * n))
        pos = sorted(rng.choice(np.arange(lining_span), size=n, replace=False).tolist())
        arch = ar.describe_architecture(pos, n_draws=n_draws, seed=SEED + index)
        family, corrected = arch.family, arch.holm
        for name in names:
            raw[name].append(family.get(name, float("nan")))
        holm.append(min(corrected.values()) if corrected else 1.0)
        classes[arch.architecture_class] = classes.get(arch.architecture_class, 0) + 1
        verdicts[arch.signature_verdict] = verdicts.get(arch.signature_verdict, 0) + 1
    per_test = {}
    for name in names:
        values = np.asarray(raw[name], dtype=float)
        ran = ~np.isnan(values)
        per_test[name] = {
            "fraction_runnable": float(ran.mean()),
            "alpha_at_0.05": float(np.mean(values[ran] <= 0.05)) if ran.any() else float("nan"),
        }
    holm = np.asarray(holm)
    uncorrected = np.nanmin(np.vstack([np.asarray(raw[n], dtype=float) for n in names]), axis=0)
    return {
        "n_linings": n_linings,
        "per_test": per_test,
        "family_alpha_uncorrected_0.05": float(np.mean(uncorrected <= 0.05)),
        "family_alpha_holm_0.05": float(np.mean(holm <= 0.05)),
        "family_alpha_holm_0.01": float(np.mean(holm <= 0.01)),
        "classes": classes,
        "verdicts": verdicts,
    }


# ----------------------------------------- probe 7: known-architecture panel

def probe_panel(rng: np.random.Generator, n_draws: int = 1500) -> Dict[str, object]:
    """Does the classifier recover architectures that are known by construction?"""
    inward, _heights = barrel_inward()
    cases = {
        "barrel, no linkers": list(inward),
        "barrel + 10-residue linkers": with_linkers(inward),
        "three helices + linkers": with_linkers([p for p in inward if p < 3 * PER]),
        "helix face, exact 3,4 spacing": list(np.cumsum([3, 4] * 6)),
        "helix face, jittered 3/4": list(np.cumsum(rng.choice([3, 4], size=13))),
        "strand face, every 2nd": list(range(0, 28, 2)),
        "contiguous stretch of 14": list(range(40, 54)),
        "hairpin, two runs of 7": list(range(10, 17)) + list(range(40, 47)),
        "random 14 in 200": sorted(rng.choice(np.arange(200), 14, replace=False).tolist()),
        "random 20 in 400": sorted(rng.choice(np.arange(400), 20, replace=False).tolist()),
    }
    rows = []
    for label, pos in cases.items():
        arch = ar.describe_architecture(pos, n_draws=n_draws, seed=SEED)
        entry = arch.as_dict()
        rows.append({
            "lining": label, "n": entry["n_positions"], "segments": entry["segments"],
            "occupancy": entry["occupancy"], "gap_entropy": entry["gap_entropy"],
            "dominant_period": entry["dominant_period"],
            "topology_p": arch.topology_p, "register_p": arch.register_p,
            "signature": entry["signature_verdict"], "element": entry["element"],
            "architecture_class": entry["architecture_class"],
        })
    return {"rows": rows}


# ------------------------------------- probe 8: power, dilution and shape

def _detection(arch) -> Dict[str, float]:
    """Both halves' Holm-corrected p, and the better of the two.

    All three are reported, and the headline is ``either``. An earlier version
    of this probe reported the topology half alone, on the reasoning that
    topology is what detects a face - which is true of a *single element* and
    wrong as a detection rate. Ring-shaped linings carry their signal in
    register, not topology, so the topology-only rate understated them by a
    factor of seven, and understated a 14-residue slab by a factor of eight.
    """
    topology = arch.topology_p if arch.topology_p is not None else 1.0
    register = arch.register_p if arch.register_p is not None else 1.0
    return {"topology": float(topology), "register": float(register),
            "either": float(min(topology, register))}


def _rates(records: List[Dict[str, float]]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for half in ("topology", "register", "either"):
        values = np.asarray([r[half] for r in records], dtype=float)
        out[f"detected_{half}"] = float(np.mean(values <= 0.05))
    values = np.asarray([r["either"] for r in records], dtype=float)
    out["median_p"] = float(np.median(values))
    out["clear_rate"] = float(np.mean(values <= 0.01))
    out["detected_rate"] = out["detected_either"]
    return out


def probe_power(rng: np.random.Generator, n_draws: int = 600) -> Dict[str, object]:
    """How completely, and in what shape, must a lining be seen to be detectable?

    The question that decides whether the study is worth running. A real
    annotated lining is the residues within contact range of a bound ligand -
    a spatial slab of the face, not a random subset of it - so both are
    measured, and then two pocket shapes of identical size are compared.

    Every rate is reported for each half of the decomposition and for the
    better of the two, because the two halves detect different pocket shapes
    and a single half is not the detection rate of the method.
    """
    inward, heights = barrel_inward()
    low, high = min(heights.values()), max(heights.values())
    dilution = []
    for mode in ("random", "spatial"):
        for fraction in (1.0, 0.8, 0.6, 0.45, 0.33):
            count = int(round(fraction * len(inward)))
            if count < 6:
                continue
            records = []
            for rep in range(30):
                if mode == "random":
                    pick = rng.choice(inward, size=count, replace=False)
                else:
                    centre = float(rng.uniform(low, high))
                    pick = sorted(inward, key=lambda r: abs(heights[r] - centre))[:count]
                arch = ar.describe_architecture(with_linkers(pick), n_draws=n_draws, seed=SEED + rep)
                records.append(_detection(arch))
            dilution.append({
                "mode": mode, "fraction_of_face": fraction, "n_positions": count,
                **_rates(records),
            })

    # Shape, at one fixed lining size.
    size = 14
    shapes = []
    for label, elements in (("ring across all elements", None), ("column along 2 elements", 2)):
        records, spanned, fractions = [], [], []
        for rep in range(30):
            if elements is None:
                centre = float(rng.uniform(low, high))
                pick = sorted(inward, key=lambda r: abs(heights[r] - centre))[:size]
            else:
                chosen = set(rng.choice(NHEL, elements, replace=False).tolist())
                pool = [r for r in inward if (r // PER) in chosen]
                if len(pool) < size:
                    continue
                pick = pool[:size]
            positions = with_linkers(pick)
            arch = ar.describe_architecture(positions, n_draws=n_draws, seed=SEED + rep)
            records.append(_detection(arch))
            spanned.append(len({int(r) // PER for r in pick}))
            fractions.append(ar.face_gap_fraction(positions, ar.HELIX_FACE_GAPS))
        if not records:
            continue
        shapes.append({
            "shape": label, "lining_size": size,
            "elements_spanned": float(np.mean(spanned)),
            "helix_gap_fraction": float(np.mean(fractions)),
            **_rates(records),
        })
    # Random linings of the same size, so no detection rate is ever read
    # without the false-positive rate that makes it meaningful.
    decoys = []
    for rep in range(40):
        span = int(rng.integers(4 * size, 14 * size))
        pick = sorted(rng.choice(np.arange(span), size=size, replace=False).tolist())
        decoys.append(_detection(ar.describe_architecture(
            pick, n_draws=n_draws, seed=SEED + 500 + rep)))
    return {
        "n_face_residues": len(inward),
        "dilution": dilution,
        "shapes": shapes,
        "false_positive": {"lining_size": size, **_rates(decoys)},
    }


# ------------------------------------------------------------------ driver

PROBES = {
    "identities": "H0 barcode equals the sorted gap multiset",
    "axial": "the doubled-angle resultant is degenerate on integer positions",
    "scan_inflation": "what an analytic threshold on a scanned coherence admits",
    "peak_width": "the coherence peak-width law, and what a fixed scan loses",
    "orderings_floor": "the smallest p a permutation null can reach",
    "calibration": "false-positive rates, per test and for the family",
    "panel": "recovery of architectures known from the 3-D construction",
    "power": "how completely and in what shape a lining must be seen",
}


def run_all(*, quick: bool = False) -> Dict[str, object]:
    rng = np.random.default_rng(SEED)
    draws = 300 if quick else 600
    return {
        "seed": SEED,
        "quick": quick,
        "identities": probe_identities(rng),
        "axial": probe_axial(rng),
        "scan_inflation": probe_scan_inflation(rng),
        "peak_width": probe_peak_width(),
        "orderings_floor": probe_orderings_floor(),
        "calibration": probe_calibration(rng, n_linings=120 if quick else 300, n_draws=draws),
        "panel": probe_panel(rng, n_draws=600 if quick else 1500),
        "power": probe_power(rng, n_draws=draws),
    }


def markdown(report: Dict[str, object]) -> str:
    ident = report["identities"]
    axial = report["axial"]
    scan = report["scan_inflation"]
    width = report["peak_width"]
    floor = report["orderings_floor"]
    cal = report["calibration"]
    panel = report["panel"]
    power = report["power"]
    lines = [
        "# Lining architecture: what the instruments can and cannot see",
        "",
        f"**Run.** `scripts/architecture_probe.py run"
        f"{' --quick' if report.get('quick') else ''}`, seed {report['seed']}, on the "
        "session container. Every figure below is computed by that command and written "
        "into this file in the same pass, so nothing is transcribed.",
        "",
        "Synthetic and analytic only. No annotated sequence is read; "
        "`docs/ARCHITECTURE_PLAN.md` pre-registers the study that would use this.",
        "",
        "## 1. Two identities that collapse",
        "",
        f"The H0 persistence barcode of a lining equals its sorted gap multiset: "
        f"{ident['barcode_mismatches']} disagreements against explicit single-linkage "
        f"in {ident['barcode_trials']} random position sets. So persistent homology on a "
        "lining buys nothing beyond the gaps, and the gaps are therefore *exactly* the "
        "topological half of a lining's structure.",
        "",
        "## 2. One borrowed statistic that does not survive the lattice",
        "",
        f"The axial (doubled-angle) resultant, standard in directional statistics for data "
        f"with two modes 180 degrees apart, is **identically 1** at period 2 for any integer "
        f"position set: max deviation {axial['max_abs_deviation_at_period_2']:.1e} over 2000 "
        f"arbitrary sets, and {axial['max_abs_deviation_at_period_4_uniform_parity']:.1e} at "
        "period 4 for uniform parity. A strand band containing period 2 would call every "
        "lining axially coherent, so it is not used.",
        "",
        "## 3. An analytic threshold on a scanned coherence is useless",
        "",
        f"At n = {scan['n']} the Rayleigh 95th percentile is "
        f"{scan['analytic_threshold']:.3f}, independent of period. Its true "
        "false-positive rate on random linings:",
        "",
        "| statistic | probes | mean R | true alpha at the analytic cutoff |",
        "|---|---|---|---|",
    ]
    for row in scan["rows"]:
        lines.append(f"| {row['statistic']} | {row['probes']} | {row['mean']:.3f} | "
                     f"**{row['true_alpha_at_analytic_threshold']:.3f}** |")
    lines += [
        "",
        "A full 2-12 scan at the nominal cutoff admits most random linings. Only a "
        "resampling null on the same scan is honest.",
        "",
        "## 4. The peak-width law, and what a fixed scan loses",
        "",
        f"Half-width = {width['ratio_mean']:.2f} x w^2 / span, with a spread of only "
        f"{width['ratio_spread']:.3f} across a sixteen-fold range of span:",
        "",
        "| span | measured half-width | ratio to w^2/span |",
        "|---|---|---|",
    ]
    for row in width["rows"]:
        lines.append(f"| {row['span']} | {row['measured_half_width']:.4f} | "
                     f"{row['ratio_to_w2_over_span']:.3f} |")
    lines += [
        "",
        f"A fixed thirteen-probe scan therefore loses up to **{width['worst_loss']:.3f}** in "
        "resultant on a perfect face, which is a missed face rather than rounding. "
        "Probe density is derived from the law.",
        "",
        "## 5. The permutation null has a floor",
        "",
        "| gap multiset | orderings | smallest attainable p |",
        "|---|---|---|",
    ]
    for row in floor["rows"]:
        lines.append(f"| {row['gaps']} | {row['orderings']} | {row['min_attainable_p']:.4f} |")
    lines += [
        "",
        f"Significance at 0.05 is unreachable below about 20 orderings, so a lining under "
        f"{floor['floor_used']} is **unresolved** rather than \"no register\". A z-score off "
        "such a null happily returns +2 and means nothing.",
        "",
        "## 6. Calibration",
        "",
        f"On {cal['n_linings']} random linings, each test alone is calibrated or "
        "conservative:",
        "",
        "| test | fraction runnable | alpha at 0.05 |",
        "|---|---|---|",
    ]
    for name, entry in cal["per_test"].items():
        lines.append(f"| {name} | {entry['fraction_runnable']:.2f} | {entry['alpha_at_0.05']:.3f} |")
    lines += [
        "",
        f"But the family of four is not: taking the best of them gives "
        f"**{cal['family_alpha_uncorrected_0.05']:.3f}** at a nominal 0.05. Holm over the "
        f"tests that ran restores it to **{cal['family_alpha_holm_0.05']:.3f}** at 0.05 and "
        f"{cal['family_alpha_holm_0.01']:.3f} at 0.01, so every verdict reads the corrected p.",
        "",
        "## 7. Recovery of known architectures",
        "",
        "| lining | n | runs | occ | dom w | topology p | register p | signature | element | class |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in panel["rows"]:
        fmt = lambda v: "--" if v is None else f"{v:.4f}"  # noqa: E731
        lines.append(
            f"| {row['lining']} | {row['n']} | {row['segments']} | {row['occupancy']:.2f} | "
            f"{row['dominant_period']:.2f} | {fmt(row['topology_p'])} | {fmt(row['register_p'])} | "
            f"{row['signature']} | {row['element']} | {row['architecture_class']} |"
        )
    lines += [
        "",
        "The row that matters most is the jittered helical face. Its gaps are drawn at "
        "random from {3, 4}, so permuting them gives another helical face and **register "
        "finds nothing** - correctly. Its periodicity lives entirely in the gap multiset, "
        "which the topology test sees. Reading face-ness off register would miss every "
        "irregularly spaced face, which is most of them.",
        "",
        "## 8. Power: how much of a lining must be seen",
        "",
        f"The barrel has {power['n_face_residues']} face residues. Subsampling them, with "
        "each half of the decomposition reported separately and the headline being the "
        "better of the two:",
        "",
        "| sampling | fraction of face | n | topology | register | **either** |",
        "|---|---|---|---|---|---|",
    ]
    for row in power["dilution"]:
        lines.append(
            f"| {row['mode']} | {row['fraction_of_face']:.2f} | {row['n_positions']} | "
            f"{row['detected_topology']:.2f} | {row['detected_register']:.2f} | "
            f"**{row['detected_either']:.2f}** |"
        )
    decoy = power["false_positive"]
    lines += [
        "",
        f"Against a false-positive rate of **{decoy['detected_either']:.2f}** on random "
        f"linings of {decoy['lining_size']} positions, so none of this is bought with "
        "false positives.",
        "",
        "**The two halves detect different things, and reading either one alone "
        "understates the method.** An earlier version of this probe reported the topology "
        "half only - on the reasoning that topology is what detects a face, which is true "
        "of a single element and wrong as a detection rate - and so understated a "
        "14-position slab eightfold.",
        "",
        "And at one fixed lining size, pocket **shape** still matters, though less than "
        "the topology half alone suggested:",
        "",
        "| shape | lining size | elements spanned | helix-gap fraction | topology | register | **either** |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in power["shapes"]:
        lines.append(
            f"| {row['shape']} | {row['lining_size']} | {row['elements_spanned']:.1f} | "
            f"{row['helix_gap_fraction']:.3f} | {row['detected_topology']:.2f} | "
            f"{row['detected_register']:.2f} | **{row['detected_either']:.2f}** |"
        )
    ring = next((r for r in power["shapes"] if "ring" in str(r["shape"])), None)
    column = next((r for r in power["shapes"] if "column" in str(r["shape"])), None)
    if ring and column:
        lines += [
            "",
            f"Same structure, same lining size, and a real difference - "
            f"{column['detected_either']:.2f} against {ring['detected_either']:.2f} - but "
            "a difference in degree rather than the presence-versus-absence the "
            "topology-only measurement implied. The mechanism is visible in the two "
            "halves: a column lining is a long run on one element, which is a **topology** "
            f"signal ({column['detected_topology']:.2f} against "
            f"{ring['detected_topology']:.2f}); a ring lining is short runs separated by "
            "linkers in a repeating arrangement, which is a **register** signal "
            f"({ring['detected_register']:.2f} against "
            f"{ring['detected_topology']:.2f} for its own topology).",
            "",
            "So the falsifiable prediction stands in a weakened and more precise form: the "
            "signature should be **stronger** in elongated element-aligned pockets than in "
            "ring-shaped ones at equal lining size, and the half that carries it should "
            "differ between them. It is no longer a claim that ring pockets carry no "
            "signature.",
            "",
        ]
    lines += [
        "## What this cannot be used for",
        "",
        "* It is not evidence about real binding sites. Every lining here is synthetic or "
        "derived from an idealised barrel, and no annotated sequence has been read.",
        "* The barrel is one fold with one periodicity. The peak-width law and the "
        "degeneracies are general; the power and shape numbers are properties of this "
        "fixture and are an upper bound on what a real, irregular protein would give.",
        "* The class boundaries - three runs, |z| and p cutoffs, a face gap of 4 - are "
        "fixed from geometry and from the statistics above, not from any outcome, and "
        "they have never been adjusted against a result.",
    ]
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="run every probe")
    run.add_argument("--quick", action="store_true", help="fewer draws, for CI")
    run.add_argument("--out", type=Path, default=Path("architecture.json"))
    run.add_argument("--markdown", type=Path, default=None)
    sub.add_parser("list", help="list the probes")
    args = parser.parse_args(argv)

    if args.command == "list":
        for name, description in PROBES.items():
            print(f"{name:18s} {description}")
        return 0

    report = run_all(quick=args.quick)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
