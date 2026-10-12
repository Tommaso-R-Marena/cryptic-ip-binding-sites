#!/usr/bin/env python3
"""Figures for the lining-architecture results.

Every number is read from the committed probe artefacts - `architecture.json`
and `homology.json` - and never recomputed here, so a figure cannot disagree
with the report it illustrates. If a probe is re-run, the figures change with
it or the script fails loudly on a missing key.

Three figures:

* **F1 power** - detection rate against how much of a face a lining captures,
  for random and spatial subsampling, with the band of real IP lining sizes
  marked. This is the figure that says where the method dies.
* **F2 shape** - detection rate for a ring-shaped and a column-shaped pocket of
  identical lining size. The falsifiable prediction, in one panel.
* **F3 homology** - fold-family retrieval beside its two controls, and the
  distance matrix that produced it.

Rendered light-mode: these are manuscript figures, and a vector PDF for a paper
has one surface. The PNG beside it is for viewing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

#: Categorical slots 1-3 of the validated palette. The first three slots are
#: the ones that clear the all-pairs colour-vision floors in both modes, which
#: is why no figure here carries a fourth series.
SERIES = ("#2a78d6", "#eb6834", "#1baf7a")

#: Sequential ramp for the distance matrix: one hue, light to dark. Never a
#: rainbow, and never a hue at a diverging midpoint - this quantity has no
#: midpoint.
SEQUENTIAL = "Blues"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8982"
GRID = "#e6e5e1"

#: Lining sizes the project's annotated IP sites actually have, marked on the
#: power figure because that is what makes the curve a verdict rather than a
#: curiosity.
REAL_LINING_RANGE = (10, 20)


def style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.labelcolor": INK_SECONDARY,
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.8,
        "xtick.color": INK_SECONDARY,
        "ytick.color": INK_SECONDARY,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
    })


def recede(ax, *, grid_axis: str = "y") -> None:
    """Grid and axes recessive; data in front."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(True, axis=grid_axis, zorder=0)
    ax.set_axisbelow(True)


# ------------------------------------------------------------------- F1

def figure_power(report: Dict[str, object], path: Path) -> Path:
    """Detection rate against the fraction of a face a lining captures."""
    power = report["power"]
    n_face = int(power["n_face_residues"])  # type: ignore[index]
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    recede(ax)

    modes = [("spatial", "spatial slab\n(what a ligand selects)", SERIES[0], "o", "-"),
             ("random", "random subset", SERIES[1], "s", "--")]
    rows_by_mode = {
        key: sorted((r for r in power["dilution"] if r["mode"] == key),  # type: ignore[union-attr]
                    key=lambda r: r["fraction_of_face"])
        for key, *_ in modes
    }

    def other(table, key, at_x):
        """The other series' value at the same x, for finding the widest gap."""
        name = "random" if key == "spatial" else "spatial"
        for row in table[name]:
            if row["n_positions"] == at_x:
                return row["detected_rate"]
        return 0.0

    for key, label, colour, marker, dash in modes:
        rows = rows_by_mode[key]
        x = [r["n_positions"] for r in rows]
        y = [r["detected_rate"] for r in rows]
        ax.plot(x, y, color=colour, marker=marker, linestyle=dash, label=label,
                markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=3,
                clip_on=False)
        # Direct-label where the curves are furthest apart, never at the shared
        # endpoint: both series finish at 100% and labels there collide.
        gaps = [abs(y[i] - other(rows_by_mode, key, x[i])) for i in range(len(x))]
        spot = int(np.argmax(gaps))
        ax.annotate(label.split("\n")[0], xy=(x[spot], y[spot]),
                    xytext=(9, -13),
                    textcoords="offset points", color=colour, fontsize=8,
                    fontweight="bold", va="center", ha="left")

    low, high = REAL_LINING_RANGE
    ax.axvspan(low, high, color=INK_MUTED, alpha=0.10, zorder=1, lw=0)
    ax.annotate("annotated IP linings in this project", xy=((low + high) / 2, 1.13),
                ha="center", va="bottom", fontsize=7.5, color=INK_SECONDARY)

    ax.set_xlabel(f"lining positions observed  (of {n_face} on the face)")
    ax.set_ylabel("detected at p ≤ 0.05")
    ax.set_title("The signature needs most of the face")
    ax.set_ylim(-0.03, 1.26)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.02))
    fig.tight_layout()
    return save(fig, path)


# ------------------------------------------------------------------- F2

def figure_shape(report: Dict[str, object], path: Path) -> Path:
    """The falsifiable prediction: same size, opposite answers."""
    shapes = report["power"]["shapes"]  # type: ignore[index]
    rows = sorted(shapes, key=lambda r: r["detected_rate"])  # type: ignore[call-overload]
    size = int(rows[0]["lining_size"])

    fig, ax = plt.subplots(figsize=(5.4, 2.5))
    recede(ax, grid_axis="x")
    positions = np.arange(len(rows))
    colours = [SERIES[1], SERIES[0]]
    bars = ax.barh(positions, [r["detected_rate"] for r in rows], height=0.52,
                   color=colours[: len(rows)], zorder=3)
    for bar, row in zip(bars, rows):
        bar.set_linewidth(2.0)
        bar.set_edgecolor(SURFACE)
        # The value sits outside the bar: a 7% bar is far too short to hold text.
        ax.annotate(f"{row['detected_rate'] * 100:.0f}%",
                    xy=(row["detected_rate"], bar.get_y() + bar.get_height() / 2),
                    xytext=(7, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=11, fontweight="bold", color=INK)

    ax.set_yticks(positions)
    ax.set_yticklabels(
        [f"{str(r['shape']).split(':')[0]}\n"
         f"{r['elements_spanned']:.0f} elements · {r['helix_gap_fraction']:.2f} on-face"
         for r in rows], fontsize=8, color=INK)
    for tick in ax.get_yticklabels():
        tick.set_linespacing(1.5)
    ax.set_xlim(0, 1.16)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.set_xlabel("detected at p ≤ 0.05")
    ax.set_title(f"Shape decides, not size — both linings are {size} residues")
    fig.tight_layout()
    return save(fig, path)


# ------------------------------------------------------------------- F3

def figure_homology(report: Dict[str, object], path: Path) -> Path:
    """Fold retrieval beside its controls, and the matrix behind it."""
    rows = list(report["retrieval"])  # type: ignore[arg-type]
    labels: List[str] = list(report["labels"])  # type: ignore[arg-type]
    matrix = np.asarray(report["distance_matrix"], dtype=float)
    best = str(report["best_metric"])

    fig, (left, right) = plt.subplots(
        1, 2, figsize=(8.4, 3.5), gridspec_kw={"width_ratios": [1.05, 1.0]})

    # --- left: precision at one, controls hatched
    recede(left, grid_axis="x")
    order = list(rows)
    positions = np.arange(len(order))[::-1]
    for position, row in zip(positions, order):
        control = bool(row["is_control"])
        colour = INK_MUTED if control else SERIES[0]
        bar = left.barh([position], [row["precision_at_one"]], height=0.58,
                        color=colour, zorder=3,
                        hatch="////" if control else None,
                        edgecolor=SURFACE, linewidth=2.0)
        left.annotate(f"{row['precision_at_one']:.2f}",
                      xy=(row["precision_at_one"], position), xytext=(6, 0),
                      textcoords="offset points", va="center", ha="left",
                      fontsize=9, fontweight="bold",
                      color=INK_SECONDARY if control else INK)
        del bar
    left.axvline(0.5, color=INK, linestyle=":", linewidth=1.2, zorder=4)
    left.set_ylim(float(positions[-1]) - 1.0, float(positions[0]) + 0.6)
    left.annotate("chance", xy=(0.5, float(positions[-1]) - 0.62), ha="center", va="center",
                  fontsize=7.5, color=INK_SECONDARY,
                  bbox=dict(boxstyle="round,pad=0.18", facecolor=SURFACE, edgecolor="none"))
    left.set_yticks(positions)
    left.set_yticklabels([f"{r['metric']}" + ("  (control)" if r["is_control"] else "")
                          for r in order], fontsize=8, color=INK)
    left.set_xlim(0, 1.18)
    left.set_xticks([0, 0.5, 1.0])
    left.set_xticklabels(["0", "0.5", "1.0"])
    left.set_xlabel("nearest neighbour is the same fold  (P@1)", labelpad=8)
    left.set_title("Both controls are dead")

    # --- right: the distance matrix, ordered by family
    index = sorted(range(len(labels)), key=lambda i: (labels[i], i))
    ordered = matrix[np.ix_(index, index)]
    image = right.imshow(ordered, cmap=SEQUENTIAL, origin="upper",
                         interpolation="nearest")
    boundary = sum(1 for i in index if labels[i] == sorted(set(labels))[0])
    for edge in (boundary - 0.5,):
        right.axhline(edge, color=SURFACE, linewidth=2.5)
        right.axvline(edge, color=SURFACE, linewidth=2.5)
    families = sorted(set(labels))
    centres = []
    start = 0
    for family in families:
        count = labels.count(family)
        centres.append(start + count / 2 - 0.5)
        start += count
    right.set_xticks(centres)
    right.set_xticklabels(families, fontsize=8, color=INK)
    right.set_yticks(centres)
    right.set_yticklabels(families, fontsize=8, color=INK, rotation=90, va="center")
    right.set_title(f"`{best}` distance, ordered by fold")
    for side in ("top", "right", "bottom", "left"):
        right.spines[side].set_visible(False)
    right.tick_params(length=0)
    bar = fig.colorbar(image, ax=right, fraction=0.046, pad=0.04)
    bar.set_label("distance", fontsize=8, color=INK_SECONDARY)
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=7, length=0, colors=INK_SECONDARY)

    handles = [Line2D([], [], marker="s", linestyle="none", color=SERIES[0],
                      markersize=7, label="architecture metric"),
               Line2D([], [], marker="s", linestyle="none", color=INK_MUTED,
                      markersize=7, label="control (knows one thing only)")]
    left.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.22),
                ncol=2, columnspacing=1.2, handletextpad=0.4)
    fig.tight_layout()
    return save(fig, path)


# ----------------------------------------------------------------- output

def save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path.with_suffix(".png")


FIGURES = {
    "power": "detection rate against the fraction of a face observed",
    "shape": "ring versus column at identical lining size",
    "homology": "fold retrieval beside its controls, and the distance matrix",
}


def build(architecture: Path, homology: Path, out_dir: Path) -> List[Path]:
    style()
    architecture_report = json.loads(architecture.read_text())
    homology_report = json.loads(homology.read_text())
    return [
        figure_power(architecture_report, out_dir / "F1_power"),
        figure_shape(architecture_report, out_dir / "F2_shape"),
        figure_homology(homology_report, out_dir / "F3_homology"),
    ]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--architecture", type=Path,
                        default=Path("results/architecture/architecture.json"))
    parser.add_argument("--homology", type=Path,
                        default=Path("results/architecture/homology.json"))
    parser.add_argument("--out-dir", type=Path,
                        default=Path("results/architecture/figures"))
    args = parser.parse_args(argv)
    for path in build(args.architecture, args.homology, args.out_dir):
        print(f"wrote {path} and {path.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
