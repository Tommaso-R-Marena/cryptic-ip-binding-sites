#!/usr/bin/env python3
"""Automated publication-ready figure generation for cryptic IP-binding site discovery.

Usage:
    python scripts/generate_publication_figures.py --config path/to/figure_config.yaml
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
from pathlib import Path
from textwrap import dedent
from typing import Dict, Iterable, Optional

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


DEFAULT_LEGENDS = dedent(
    """
    Figure 1 | Overview of cryptic IP-binding site discovery workflow.
    (A) ADAR2 structure with buried IP6 shown in stick representation and pocket surface shading.
    (B) End-to-end computational pipeline: AlphaFold structural input, fpocket cavity detection,
    FreeSASA solvent accessibility filtering, APBS electrostatic feature extraction, and final scoring.
    (C) Receiver operating characteristic (ROC) curves benchmarking candidate scoring performance.

    Figure 2 | Comparative proteomics of predicted IP-binding sites.
    (A) Hit-rate comparison across organisms with statistical significance annotations.
    (B) Correlation between IP6 concentration and predicted hit rate.
    (C) GO term enrichment heatmap for high-confidence candidates.
    (D) Phylogenetic conservation profile across top predicted hits.

    Figure 3 | Structural gallery of top predicted cryptic IP-binding candidates.
    Top 10 candidates shown with matched structural and electrostatic views,
    unified color scales, and consistent annotation/scale bars.
    """
).strip()


def apply_publication_style(style_file: Path) -> None:
    """Apply journal-like plotting defaults."""
    plt.style.use(str(style_file))
    sns.set_theme(context="paper", style="whitegrid", font="DejaVu Sans")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def read_csv(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None:
        return None
    if not path.exists():
        raise FileNotFoundError(f"Missing required CSV: {path}")
    return pd.read_csv(path)


def export_figure(fig: plt.Figure, basename: Path, dpi: int = 600) -> None:
    """Export vector PDF and high-resolution PNG; optionally write CMYK PNG copy."""
    fig.savefig(basename.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(basename.with_suffix(".png"), dpi=dpi, bbox_inches="tight")

    if Image is not None:
        png_path = basename.with_suffix(".png")
        cmyk_path = basename.with_name(f"{basename.stem}_CMYK.png")
        with Image.open(png_path) as img:
            img.convert("CMYK").save(cmyk_path)


def run_pymol_render(structure_path: Path, output_png: Path, ip6_resn: str = "IP6") -> bool:
    """Render ADAR2 panel with PyMOL if available."""
    pymol_bin = shutil.which("pymol")
    if pymol_bin is None:
        return False

    pml_script = dedent(
        f"""
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
        """
    )
    script_path = output_png.with_suffix(".pml")
    script_path.write_text(pml_script, encoding="utf-8")

    subprocess.run([pymol_bin, "-cq", str(script_path)], check=False)
    return output_png.exists()


def plot_pipeline_schematic(ax: plt.Axes) -> None:
    steps = ["AlphaFold", "fpocket", "FreeSASA", "APBS", "Scoring"]
    x_positions = np.linspace(0.08, 0.92, num=len(steps))

    for x, step in zip(x_positions, steps):
        rect = patches.FancyBboxPatch(
            (x - 0.08, 0.42),
            0.16,
            0.18,
            boxstyle="round,pad=0.02",
            linewidth=1.2,
            edgecolor="#1f4e79",
            facecolor="#d6eaf8",
        )
        ax.add_patch(rect)
        ax.text(x, 0.51, step, ha="center", va="center", fontsize=9, fontweight="bold")

    for i in range(len(x_positions) - 1):
        ax.annotate(
            "",
            xy=(x_positions[i + 1] - 0.09, 0.51),
            xytext=(x_positions[i] + 0.09, 0.51),
            arrowprops=dict(arrowstyle="->", lw=1.5, color="#1f4e79"),
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


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
        axes[0].set_title("A. ADAR2 buried IP6")
        axes[0].axis("off")
    else:
        axes[0].text(0.5, 0.5, "PyMOL unavailable\n(Add pre-rendered panel image)", ha="center", va="center")
        axes[0].set_title("A. ADAR2 buried IP6")
        axes[0].axis("off")

    plot_pipeline_schematic(axes[1])
    axes[1].set_title("B. Discovery pipeline")

    sns.lineplot(data=roc_df, x="fpr", y="tpr", hue="model", linewidth=2.0, ax=axes[2])
    axes[2].plot([0, 1], [0, 1], "k--", linewidth=1)
    axes[2].set_title("C. Validation ROC")
    axes[2].set_xlabel("False Positive Rate")
    axes[2].set_ylabel("True Positive Rate")
    axes[2].legend(title="Model", frameon=False, fontsize=8)

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

    sns.barplot(data=hits_df, x="organism", y="hit_rate", ax=axes[0, 0], palette="Blues_d")
    axes[0, 0].set_title("A. Hit rate by organism")
    axes[0, 0].set_xlabel("Organism")
    axes[0, 0].set_ylabel("Hit rate")
    if {"group1", "group2", "stars"}.issubset(hits_df.columns):
        x_map = {k: i for i, k in enumerate(hits_df["organism"].tolist())}
        y_top = hits_df["hit_rate"].max()
        for _, row in hits_df.dropna(subset=["group1", "group2", "stars"]).iterrows():
            _add_sig_stars(axes[0, 0], x_map[row["group1"]], x_map[row["group2"]], y_top * 1.05, row["stars"])
            y_top *= 1.08

    sns.regplot(data=corr_df, x="ip6_concentration", y="hit_rate", scatter_kws={"s": 45}, ax=axes[0, 1])
    axes[0, 1].set_title("B. IP6 concentration vs hit rate")
    axes[0, 1].set_xlabel("IP6 concentration")
    axes[0, 1].set_ylabel("Hit rate")

    go_mat = go_df.set_index(go_df.columns[0])
    sns.heatmap(go_mat, cmap="mako", linewidths=0.2, linecolor="white", ax=axes[1, 0])
    axes[1, 0].set_title("C. GO term enrichment")

    if {"species", "conservation_score"}.issubset(phylo_df.columns):
        sns.pointplot(data=phylo_df, x="species", y="conservation_score", hue="candidate", ax=axes[1, 1])
    else:
        sns.heatmap(phylo_df.set_index(phylo_df.columns[0]), cmap="viridis", ax=axes[1, 1])
    axes[1, 1].set_title("D. Phylogenetic conservation")
    axes[1, 1].tick_params(axis="x", rotation=45)

    export_figure(fig, output_dir / "Figure2_Comparative_Proteomics")
    plt.close(fig)


def make_figure3(config: Dict, output_dir: Path) -> None:
    gallery_df = read_csv(Path(config["figure3"]["top_candidates_csv"]))
    top_n = int(config["figure3"].get("top_n", 10))
    gallery_df = gallery_df.head(top_n)

    rows, cols = 5, 2
    fig, axes = plt.subplots(rows, cols, figsize=(10.5, 16), constrained_layout=True)
    axes = axes.flatten()

    for ax, (_, row) in zip(axes, gallery_df.iterrows()):
        struct_img = Path(row["structure_image"])
        elec_img = Path(row["electrostatic_image"])
        candidate = row["candidate"]

        ax.set_title(f"{candidate}", fontsize=9, loc="left")
        ax.axis("off")

        if struct_img.exists():
            left = mpimg.imread(struct_img)
            ax.imshow(left, extent=[0, 0.48, 0, 1])
        else:
            ax.text(0.24, 0.5, "Missing\nstructure", ha="center", va="center", fontsize=8)

        if elec_img.exists():
            right = mpimg.imread(elec_img)
            ax.imshow(right, extent=[0.52, 1.0, 0, 1])
        else:
            ax.text(0.76, 0.5, "Missing\nelectrostatics", ha="center", va="center", fontsize=8)

        ax.plot([0.82, 0.98], [0.06, 0.06], color="black", lw=2)
        ax.text(0.90, 0.08, "5 Å", ha="center", va="bottom", fontsize=7)
        ax.text(0.24, 0.02, "Structure", ha="center", fontsize=7)
        ax.text(0.76, 0.02, "Electrostatics", ha="center", fontsize=7)

    for ax in axes[len(gallery_df) :]:
        ax.axis("off")

    export_figure(fig, output_dir / "Figure3_Top_Candidates")
    plt.close(fig)


def write_legends(output_dir: Path, custom_text: Optional[str] = None) -> None:
    legends_path = output_dir / "figure_legends.txt"
    legends_path.write_text(custom_text.strip() if custom_text else DEFAULT_LEGENDS, encoding="utf-8")


def load_config(path: Path) -> Dict:
    import yaml

    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate publication-ready figures from CSV inputs.")
    parser.add_argument("--config", type=Path, required=True, help="YAML config file for figure inputs.")
    parser.add_argument(
        "--figures",
        nargs="+",
        choices=["figure1", "figure2", "figure3", "all"],
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

    write_legends(output_dir, config.get("figure_legends"))
    print(f"Figures written to: {output_dir}")


if __name__ == "__main__":
    main()
