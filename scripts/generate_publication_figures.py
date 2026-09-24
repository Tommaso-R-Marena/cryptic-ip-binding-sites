#!/usr/bin/env python3
"""Automated publication-ready figure generation for cryptic IP-binding site discovery.

Usage:
    python scripts/generate_publication_figures.py --config path/to/figure_config.yaml

The module renders a cohesive, journal-grade figure suite:

* Figure 1 - discovery overview (ADAR2 panel, annotated pipeline schematic, ROC with AUC).
* Figure 2 - comparative proteomics panels (hit rate, IP6 correlation, enrichment, conservation).
* Figure 3 - structural gallery of top candidates.
* Figure 4 - validation benchmark dashboard built from the real RCSB dataset + control run.

Figures share a single palette, typography, panel-label convention, and unit-annotated axes.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from textwrap import dedent
from typing import Dict, Optional

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patches

try:
    from PIL import Image
except ImportError:  # optional dependency at runtime
    Image = None

try:  # unit-annotated axis labels shared with the rest of the pipeline
    from cryptic_ip.analysis import units as _units

    ANGSTROM = _units.ANGSTROM
    ANGSTROM_SQ = _units.ANGSTROM_SQ
    ANGSTROM_CU = _units.ANGSTROM_CU
except Exception:  # pragma: no cover - keep the script runnable standalone
    ANGSTROM, ANGSTROM_SQ, ANGSTROM_CU = "\u00c5", "\u00c5\u00b2", "\u00c5\u00b3"


# --- Shared visual identity -------------------------------------------------

ACCENT = "#1f4e79"
ACCENT_LIGHT = "#d6eaf8"
CLASS_COLORS = {
    "Cryptic": "#c0392b",
    "Semi-cryptic": "#e67e22",
    "Surface": "#2e86c1",
}
CONTROL_COLORS = {"positive": "#c0392b", "negative": "#2e86c1"}
MODEL_COLORS = {"ML classifier": "#c0392b", "Threshold scoring": ACCENT}


DEFAULT_LEGENDS = dedent("""
    Figure 1 | Overview of cryptic IP-binding site discovery workflow.
    (A) ADAR2 structure with buried IP6 shown in stick representation and pocket surface shading.
    (B) End-to-end computational pipeline: AlphaFold structural input, fpocket cavity detection,
    FreeSASA solvent accessibility filtering, APBS electrostatic features, and final scoring.
    (C) Receiver operating characteristic (ROC) curves benchmarking candidate scoring performance.

    Figure 2 | Comparative proteomics of predicted IP-binding sites.
    (A) Hit-rate comparison across organisms with statistical significance annotations.
    (B) Correlation between IP6 concentration and predicted hit rate.
    (C) GO term enrichment heatmap for high-confidence candidates.
    (D) Phylogenetic conservation profile across top predicted hits.

    Figure 3 | Structural gallery of top predicted cryptic IP-binding candidates.
    Top 10 candidates shown with matched structural and electrostatic views,
    unified color scales, and consistent annotation/scale bars.

    Figure 4 | Validation benchmark on the curated RCSB IP-binding dataset.
    (A) Ligand solvent accessibility by burial class with cryptic/surface thresholds.
    (B) Dataset burial-class composition. (C) Positive vs negative control score separation.
    (D) Ligand accessibility versus crystallographic resolution, colored by burial class.
    """).strip()


def apply_publication_style(style_file: Path) -> None:
    """Apply journal-like plotting defaults."""
    plt.style.use(str(style_file))
    sns.set_theme(context="paper", style="whitegrid", font="DejaVu Sans")
    plt.rcParams.update(
        {
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
        }
    )


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_csv(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"Missing required CSV: {path}")
    return pd.read_csv(path)


def panel_label(ax: plt.Axes, letter: str) -> None:
    """Place a bold panel letter (A, B, ...) in the top-left corner of an axis."""
    ax.text(
        -0.12,
        1.08,
        letter,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def _annotate_bars(ax: plt.Axes, fmt: str = "{:.3g}", fontsize: int = 7) -> None:
    """Write the numeric value on top of each bar in an axis."""
    for patch in ax.patches:
        height = patch.get_height()
        if not np.isfinite(height):
            continue
        ax.annotate(
            fmt.format(height),
            (patch.get_x() + patch.get_width() / 2, height),
            ha="center",
            va="bottom",
            fontsize=fontsize,
            xytext=(0, 2),
            textcoords="offset points",
        )


def export_figure(fig: plt.Figure, basename: Path, dpi: int = 600) -> None:
    """Export vector PDF and high-resolution PNG; optionally write CMYK PNG copy."""
    fig.savefig(basename.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(basename.with_suffix(".png"), dpi=dpi, bbox_inches="tight")

    if Image is not None:
        png_path = basename.with_suffix(".png")
        cmyk_path = basename.with_name(f"{basename.stem}_CMYK.png")
        try:
            with Image.open(png_path) as img:
                img.convert("CMYK").save(cmyk_path)
        except OSError:
            pass


def run_pymol_render(structure_path: Path, output_png: Path, ip6_resn: str = "IP6") -> bool:
    """Render ADAR2 panel with PyMOL if available."""
    pymol_bin = shutil.which("pymol")
    if pymol_bin is None:
        return False

    pml_script = dedent(f"""
        reinitialize
        load {structure_path}
        hide everything
        show cartoon, all
        color gray80, all
        select ip6, resn {ip6_resn}
        show sticks, ip6
        color marine, ip6
        show surface, byres ip6 around 6
        set transparency, 0.35
        bg_color white
        ray 2000, 1600
        png {output_png}, dpi=600
        quit
        """)
    script_path = output_png.with_suffix(".pml")
    script_path.write_text(pml_script, encoding="utf-8")

    subprocess.run([pymol_bin, "-cq", str(script_path)], check=False)
    return output_png.exists()


def _adar2_placeholder(ax: plt.Axes) -> None:
    """Draw an informative ADAR2 fact card when a PyMOL render is unavailable."""
    ax.axis("off")
    card = patches.FancyBboxPatch(
        (0.04, 0.08),
        0.92,
        0.84,
        boxstyle="round,pad=0.02,rounding_size=0.03",
        linewidth=1.2,
        edgecolor=ACCENT,
        facecolor=ACCENT_LIGHT,
    )
    ax.add_patch(card)
    ax.text(0.5, 0.82, "ADAR2 (PDB 1ZY7)", ha="center", va="center", fontsize=11, fontweight="bold")
    ax.text(
        0.5,
        0.70,
        "Buried IP6 structural cofactor",
        ha="center",
        va="center",
        fontsize=9,
        style="italic",
    )
    facts = [
        f"Ligand SASA \u2248 0 {ANGSTROM_SQ} (fully buried)",
        f"Exterior window only 8.4 \u00d7 4.6 {ANGSTROM}",
        "6 basic coordinating residues",
        "No fold / catalysis without IP6",
    ]
    for i, fact in enumerate(facts):
        ax.text(0.12, 0.55 - i * 0.11, f"\u2022 {fact}", ha="left", va="center", fontsize=8.5)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def plot_pipeline_schematic(ax: plt.Axes) -> None:
    """Annotated end-to-end pipeline with the criterion each stage enforces."""
    steps = ["AlphaFold", "fpocket", "FreeSASA", "APBS", "Scoring"]
    criteria = [
        "pLDDT \u2265 70",
        f"depth > 15 {ANGSTROM}\nvol 300-800 {ANGSTROM_CU}",
        f"SASA < 5 {ANGSTROM_SQ}",
        "> +5 kT/e",
        "\u2265 4 basic",
    ]
    x_positions = np.linspace(0.1, 0.9, num=len(steps))

    for x, step, crit in zip(x_positions, steps, criteria):
        rect = patches.FancyBboxPatch(
            (x - 0.075, 0.55),
            0.15,
            0.20,
            boxstyle="round,pad=0.02",
            linewidth=1.2,
            edgecolor=ACCENT,
            facecolor=ACCENT_LIGHT,
        )
        ax.add_patch(rect)
        ax.text(x, 0.65, step, ha="center", va="center", fontsize=8.5, fontweight="bold")
        ax.text(x, 0.36, crit, ha="center", va="center", fontsize=6.8, color="0.25")

    for i in range(len(x_positions) - 1):
        ax.annotate(
            "",
            xy=(x_positions[i + 1] - 0.078, 0.65),
            xytext=(x_positions[i] + 0.078, 0.65),
            arrowprops=dict(arrowstyle="-|>", lw=1.5, color=ACCENT),
        )

    ax.text(0.5, 0.9, "Cryptic IP-site criteria", ha="center", va="center", fontsize=8, color="0.3")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def _auc(fpr: np.ndarray, tpr: np.ndarray) -> float:
    order = np.argsort(fpr)
    return float(np.trapz(np.asarray(tpr)[order], np.asarray(fpr)[order]))


def make_figure1(config: Dict, output_dir: Path) -> None:
    panel_a_path = output_dir / "figure1_panelA_adar2.png"
    structure_path = Path(config["figure1"]["adar2_structure"])
    ip6_resn = config["figure1"].get("ip6_resn", "IP6")
    rendered = run_pymol_render(structure_path, panel_a_path, ip6_resn=ip6_resn)

    roc_df = read_csv(Path(config["figure1"]["validation_roc_csv"]))
    if roc_df is None:
        raise ValueError("figure1.validation_roc_csv must be provided")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2), constrained_layout=True)

    if rendered:
        axes[0].imshow(mpimg.imread(panel_a_path))
        axes[0].axis("off")
    else:
        _adar2_placeholder(axes[0])
    axes[0].set_title("ADAR2 buried IP6 (validation gold standard)")
    panel_label(axes[0], "A")

    plot_pipeline_schematic(axes[1])
    axes[1].set_title("Discovery pipeline")
    panel_label(axes[1], "B")

    ax = axes[2]
    for model, grp in roc_df.groupby("model"):
        color = MODEL_COLORS.get(str(model))
        auc_val = _auc(grp["fpr"].to_numpy(), grp["tpr"].to_numpy())
        ax.plot(
            grp["fpr"], grp["tpr"], linewidth=2.2, color=color, label=f"{model} (AUC={auc_val:.2f})"
        )
        ax.fill_between(
            grp.sort_values("fpr")["fpr"], grp.sort_values("fpr")["tpr"], alpha=0.08, color=color
        )
    ax.plot([0, 1], [0, 1], "--", linewidth=1, color="0.5", label="Chance")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_title("Validation ROC")
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.legend(frameon=False, fontsize=7, loc="lower right")
    panel_label(ax, "C")

    export_figure(fig, output_dir / "Figure1_Overview")
    plt.close(fig)


def _add_sig_stars(ax: plt.Axes, x1: float, x2: float, y: float, stars: str) -> None:
    ax.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], lw=1.0, color="black")
    ax.text((x1 + x2) / 2, y + 0.015, stars, ha="center", va="bottom", fontsize=11)


def make_figure2(config: Dict, output_dir: Path) -> None:
    fig2_cfg = config["figure2"]
    hits_df = read_csv(Path(fig2_cfg["hit_rates_csv"]))
    corr_df = read_csv(Path(fig2_cfg["correlation_csv"]))
    go_df = read_csv(Path(fig2_cfg["go_heatmap_csv"]))
    phylo_df = read_csv(Path(fig2_cfg["phylo_conservation_csv"]))

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    ax = axes[0, 0]
    order = hits_df.sort_values("hit_rate", ascending=False)["organism"].tolist()
    sns.barplot(
        data=hits_df,
        x="organism",
        y="hit_rate",
        order=order,
        ax=ax,
        color=ACCENT,
        edgecolor="black",
        linewidth=0.6,
    )
    ax.set_title("Cryptic hit rate by group")
    ax.set_xlabel("Organism / ligand class")
    ax.set_ylabel("Cryptic fraction")
    if "n_structures" in hits_df.columns:
        rate_by = hits_df.set_index("organism")["n_structures"].to_dict()
        for patch, org in zip(ax.patches, order):
            ax.annotate(
                f"n={int(rate_by.get(org, 0))}",
                (patch.get_x() + patch.get_width() / 2, patch.get_height()),
                ha="center",
                va="bottom",
                fontsize=7,
                xytext=(0, 2),
                textcoords="offset points",
            )
    if {"group1", "group2", "stars"}.issubset(hits_df.columns):
        x_map = {k: i for i, k in enumerate(order)}
        y_top = hits_df["hit_rate"].max()
        for _, row in hits_df.dropna(subset=["group1", "group2", "stars"]).iterrows():
            _add_sig_stars(
                ax, x_map[row["group1"]], x_map[row["group2"]], y_top * 1.05, row["stars"]
            )
            y_top *= 1.08
    panel_label(ax, "A")

    ax = axes[0, 1]
    sns.regplot(
        data=corr_df,
        x="ip6_concentration",
        y="hit_rate",
        ax=ax,
        scatter_kws={"s": 55, "color": ACCENT, "edgecolor": "black", "linewidths": 0.5},
        line_kws={"color": "#c0392b"},
    )
    ax.set_title("IP6 concentration vs hit rate")
    ax.set_xlabel("IP6 concentration (\u00b5M)")
    ax.set_ylabel("Cryptic fraction")
    if len(corr_df["ip6_concentration"].unique()) >= 3:
        try:
            from scipy.stats import spearmanr

            rho, pval = spearmanr(corr_df["ip6_concentration"], corr_df["hit_rate"])
            if np.isfinite(rho):
                ax.annotate(
                    f"Spearman \u03c1={rho:.2f}\np={pval:.2g}",
                    (0.05, 0.9),
                    xycoords="axes fraction",
                    fontsize=8,
                    va="top",
                )
        except Exception:
            pass
    panel_label(ax, "B")

    ax = axes[1, 0]
    go_mat = go_df.set_index(go_df.columns[0])
    sns.heatmap(
        go_mat,
        cmap="rocket_r",
        linewidths=0.4,
        linecolor="white",
        annot=True,
        fmt=".3f",
        cbar_kws={"label": "Enrichment score"},
        ax=ax,
    )
    ax.set_title("Burial-class enrichment")
    ax.set_ylabel("")
    panel_label(ax, "C")

    ax = axes[1, 1]
    if {"species", "conservation_score"}.issubset(phylo_df.columns) and phylo_df[
        "conservation_score"
    ].notna().any():
        plot_df = phylo_df.dropna(subset=["species", "conservation_score"])
        if len(plot_df) > 1 and "candidate" in plot_df.columns:
            sns.pointplot(data=plot_df, x="species", y="conservation_score", hue="candidate", ax=ax)
        else:
            sns.barplot(data=plot_df, x="species", y="conservation_score", ax=ax, color=ACCENT)
        ax.set_ylabel("Conservation score")
    else:
        sns.heatmap(phylo_df.set_index(phylo_df.columns[0]), cmap="viridis", ax=ax)
    ax.set_title("Structural conservation")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=45)
    panel_label(ax, "D")

    export_figure(fig, output_dir / "Figure2_Comparative_Proteomics")
    plt.close(fig)


def make_figure3(config: Dict, output_dir: Path) -> None:
    gallery_df = read_csv(Path(config["figure3"]["top_candidates_csv"]))
    top_n = int(config["figure3"].get("top_n", 10))
    gallery_df = gallery_df.head(top_n)

    rows, cols = 5, 2
    fig, axes = plt.subplots(rows, cols, figsize=(10.5, 16), constrained_layout=True)
    axes = axes.flatten()

    for rank, (ax, (_, row)) in enumerate(zip(axes, gallery_df.iterrows()), start=1):
        struct_img = Path(row["structure_image"])
        elec_img = Path(row["electrostatic_image"])
        candidate = row["candidate"]

        ax.set_title(f"{rank}. {candidate}", fontsize=9, loc="left", fontweight="bold")
        ax.axis("off")

        if struct_img.exists():
            ax.imshow(mpimg.imread(struct_img), extent=[0, 0.48, 0, 1])
        else:
            ax.text(0.24, 0.5, "Missing\nstructure", ha="center", va="center", fontsize=8)

        if elec_img.exists():
            ax.imshow(mpimg.imread(elec_img), extent=[0.52, 1.0, 0, 1])
        else:
            ax.text(0.76, 0.5, "Missing\nelectrostatics", ha="center", va="center", fontsize=8)

        if "score" in row and pd.notna(row["score"]):
            ax.text(
                0.5,
                0.95,
                f"score = {float(row['score']):.3f}",
                ha="center",
                va="top",
                fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=ACCENT, lw=0.8),
            )

        ax.plot([0.82, 0.98], [0.06, 0.06], color="black", lw=2)
        ax.text(0.90, 0.08, f"5 {ANGSTROM}", ha="center", va="bottom", fontsize=7)
        ax.text(0.24, 0.02, "Structure", ha="center", fontsize=7)
        ax.text(0.76, 0.02, "Electrostatics", ha="center", fontsize=7)

    for ax in axes[len(gallery_df) :]:
        ax.axis("off")

    export_figure(fig, output_dir / "Figure3_Top_Candidates")
    plt.close(fig)


def make_figure4(config: Dict, output_dir: Path) -> None:
    """Validation benchmark dashboard grounded in the real RCSB dataset + control run."""
    fig4_cfg = config["figure4"]
    dataset = read_csv(Path(fig4_cfg["dataset_csv"]))
    controls_path = fig4_cfg.get("control_benchmark_csv")
    controls = (
        read_csv(Path(controls_path)) if controls_path and Path(controls_path).exists() else None
    )

    class_order = [
        c for c in ("Cryptic", "Semi-cryptic", "Surface") if c in set(dataset["classification"])
    ]
    palette = {c: CLASS_COLORS[c] for c in class_order}

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), constrained_layout=True)

    # A. Ligand SASA distribution by burial class, with cryptic/surface thresholds.
    ax = axes[0, 0]
    sns.boxplot(
        data=dataset,
        x="classification",
        y="sasa",
        order=class_order,
        palette=palette,
        showfliers=False,
        width=0.6,
        ax=ax,
    )
    sns.stripplot(
        data=dataset,
        x="classification",
        y="sasa",
        order=class_order,
        color="0.2",
        size=2.5,
        alpha=0.45,
        jitter=0.25,
        ax=ax,
    )
    # Symlog reveals the sparse low-SASA (cryptic) points despite the surface tail.
    ax.set_yscale("symlog", linthresh=20)
    ax.axhline(5, ls="--", lw=1, color="#c0392b")
    ax.axhline(50, ls="--", lw=1, color="#2e86c1")
    ax.text(
        0.02,
        5,
        f"cryptic \u2264 5 {ANGSTROM_SQ}",
        transform=ax.get_yaxis_transform(),
        va="bottom",
        ha="left",
        fontsize=6.5,
        color="#c0392b",
    )
    ax.text(
        0.02,
        50,
        f"surface \u2265 50 {ANGSTROM_SQ}",
        transform=ax.get_yaxis_transform(),
        va="bottom",
        ha="left",
        fontsize=6.5,
        color="#2e86c1",
    )
    ax.set_title("Ligand accessibility by burial class")
    ax.set_xlabel("")
    ax.set_ylabel(f"Ligand SASA ({ANGSTROM_SQ})")
    panel_label(ax, "A")

    # B. Dataset composition.
    ax = axes[0, 1]
    counts = dataset["classification"].value_counts().reindex(class_order).fillna(0)
    sns.barplot(
        x=counts.index.tolist(),
        y=counts.values,
        palette=palette,
        edgecolor="black",
        linewidth=0.6,
        ax=ax,
    )
    _annotate_bars(ax, fmt="{:.0f}")
    ax.set_title(f"Dataset composition (n={len(dataset)} structures)")
    ax.set_xlabel("")
    ax.set_ylabel("Structures")
    panel_label(ax, "B")

    # C. Positive vs negative control separation.
    ax = axes[1, 0]
    if controls is not None and {"control_type", "score"}.issubset(controls.columns):
        ctrl = controls.dropna(subset=["score"]).copy()
        type_order = [t for t in ("positive", "negative") if t in set(ctrl["control_type"])]
        sns.stripplot(
            data=ctrl,
            x="control_type",
            y="score",
            order=type_order,
            palette=CONTROL_COLORS,
            size=9,
            jitter=False,
            ax=ax,
        )
        for i, ctype in enumerate(type_order):
            mean_v = ctrl.loc[ctrl["control_type"] == ctype, "score"].mean()
            ax.plot([i - 0.25, i + 0.25], [mean_v, mean_v], color="black", lw=1.5)
        if {"positive", "negative"}.issubset(set(type_order)):
            pos_m = ctrl.loc[ctrl["control_type"] == "positive", "score"].mean()
            neg_m = ctrl.loc[ctrl["control_type"] == "negative", "score"].mean()
            ax.annotate(
                f"\u0394mean (pos\u2212neg) = {pos_m - neg_m:.3f}",
                (0.5, 0.94),
                xycoords="axes fraction",
                ha="center",
                fontsize=8,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=ACCENT, lw=0.8),
            )
        for _, r in ctrl.iterrows():
            ax.annotate(
                str(r.get("protein", "")),
                (r["control_type"], r["score"]),
                fontsize=6.5,
                xytext=(8, 0),
                textcoords="offset points",
                va="center",
            )
        ax.set_ylim(0, 1)
        ax.set_title("Control benchmark score separation")
        ax.set_xlabel("")
        ax.set_ylabel("Burial-aware validation score")
    else:
        ax.axis("off")
        ax.text(
            0.5, 0.5, "control_benchmark.csv\nnot available", ha="center", va="center", fontsize=9
        )
    panel_label(ax, "C")

    # D. Ligand SASA vs resolution, colored by class.
    ax = axes[1, 1]
    if "resolution" in dataset.columns:
        sub = dataset.dropna(subset=["resolution", "sasa"])
        for cls in class_order:
            pts = sub[sub["classification"] == cls]
            ax.scatter(
                pts["resolution"],
                pts["sasa"],
                s=22,
                alpha=0.7,
                color=CLASS_COLORS[cls],
                edgecolor="black",
                linewidth=0.3,
                label=cls,
            )
        ax.set_xlabel(f"Crystallographic resolution ({ANGSTROM})")
        ax.set_ylabel(f"Ligand SASA ({ANGSTROM_SQ})")
        ax.set_title("Accessibility vs resolution")
        ax.legend(frameon=False, fontsize=7, title="Burial class")
    else:
        ax.axis("off")
    panel_label(ax, "D")

    export_figure(fig, output_dir / "Figure4_Validation_Benchmark")
    plt.close(fig)


def write_legends(output_dir: Path, custom_text: Optional[str] = None) -> None:
    legends_path = output_dir / "figure_legends.txt"
    legends_path.write_text(
        custom_text.strip() if custom_text else DEFAULT_LEGENDS, encoding="utf-8"
    )


def load_config(path: Path) -> Dict:
    import yaml

    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate publication-ready figures from CSV inputs."
    )
    parser.add_argument(
        "--config", type=Path, required=True, help="YAML config file for figure inputs."
    )
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=["figure1", "figure2", "figure3", "figure4", "all"],
        default=["all"],
        help="Subset of figures to generate (default: all).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = load_config(args.config)

    output_dir = Path(config.get("output_dir", "figures/publication"))
    style_file = Path(config.get("style_file", "scripts/nature_publication.mplstyle"))
    ensure_dir(output_dir)
    apply_publication_style(style_file)

    targets = set(args.figures)
    if "all" in targets or "figure1" in targets:
        make_figure1(config, output_dir)
    if "all" in targets or "figure2" in targets:
        make_figure2(config, output_dir)
    if "all" in targets or "figure3" in targets:
        make_figure3(config, output_dir)
    if ("all" in targets or "figure4" in targets) and "figure4" in config:
        make_figure4(config, output_dir)

    write_legends(output_dir, config.get("figure_legends"))
    print(f"Figures written to: {output_dir}")


if __name__ == "__main__":
    main()
