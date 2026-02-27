"""Comparative multi-proteome analysis for cryptic IP-binding site prevalence."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from .statistical_validation import StatisticalValidation


@dataclass
class ComparativeResult:
    """Container for top-level comparative analysis tables."""

    hit_rate_table: pd.DataFrame
    correlation_table: pd.DataFrame
    pairwise_tests_table: pd.DataFrame
    ortholog_conservation_table: pd.DataFrame


class ComparativeIPAnalysis:
    """End-to-end comparative statistics and figure generation across proteomes."""

    REQUIRED_HIT_COLUMNS = {"organism", "protein_id", "is_hit"}
    REQUIRED_ORTHOLOG_COLUMNS = {"orthogroup", "organism", "protein_id"}

    def __init__(self, confidence_level: float = 0.95):
        self.confidence_level = confidence_level

    @staticmethod
    def _to_bool(series: pd.Series) -> pd.Series:
        if pd.api.types.is_bool_dtype(series):
            return series
        lowered = series.astype(str).str.lower()
        mapping = {
            "1": True,
            "0": False,
            "true": True,
            "false": False,
            "yes": True,
            "no": False,
        }
        if not lowered.isin(mapping.keys()).all():
            raise ValueError("is_hit column must contain boolean-like values.")
        return lowered.map(mapping)

    @staticmethod
    def _validate_columns(df: pd.DataFrame, required: Set[str], label: str) -> None:
        missing = required.difference(df.columns)
        if missing:
            missing_str = ", ".join(sorted(missing))
            raise ValueError(f"{label} is missing required columns: {missing_str}")

    def compute_hit_rates(self, hits_df: pd.DataFrame, ip6_map: Mapping[str, float]) -> pd.DataFrame:
        """Compute normalized hit rates and Wilson confidence intervals for each organism."""
        self._validate_columns(hits_df, self.REQUIRED_HIT_COLUMNS, "hits table")

        df = hits_df.copy()
        df["is_hit"] = self._to_bool(df["is_hit"])

        rows: List[Dict[str, float]] = []
        for organism, group in df.groupby("organism"):
            n_total = int(group["protein_id"].nunique())
            n_hits = int(group.loc[group["is_hit"], "protein_id"].nunique())
            p_hat = n_hits / n_total if n_total else np.nan

            ci_low, ci_high = stats.binomtest(n_hits, n_total).proportion_ci(
                confidence_level=self.confidence_level,
                method="wilson",
            )
            rows.append(
                {
                    "organism": organism,
                    "n_proteins": n_total,
                    "n_hits": n_hits,
                    "hit_rate": p_hat,
                    "hit_rate_percent": p_hat * 100,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "ci_low_percent": ci_low * 100,
                    "ci_high_percent": ci_high * 100,
                    "ip6_uM": ip6_map.get(organism, np.nan),
                }
            )

        return pd.DataFrame(rows).sort_values("organism").reset_index(drop=True)

    @staticmethod
    def spearman_ip6_correlation(hit_rates_df: pd.DataFrame) -> pd.DataFrame:
        """Test Spearman correlation between cellular IP6 concentration and hit rate."""
        valid = hit_rates_df.dropna(subset=["ip6_uM", "hit_rate"])
        rho, p_value = stats.spearmanr(valid["ip6_uM"], valid["hit_rate"])
        return pd.DataFrame(
            [
                {
                    "n_organisms": len(valid),
                    "spearman_rho": rho,
                    "spearman_p_value": p_value,
                }
            ]
        )

    @staticmethod
    def pairwise_organism_tests(hit_rates_df: pd.DataFrame) -> pd.DataFrame:
        """Run pairwise chi-square and Fisher exact tests between organisms."""
        rows = []
        by_org = hit_rates_df.set_index("organism")
        for org_a, org_b in combinations(by_org.index, 2):
            a = by_org.loc[org_a]
            b = by_org.loc[org_b]
            contingency = np.array(
                [
                    [int(a["n_hits"]), int(a["n_proteins"] - a["n_hits"])],
                    [int(b["n_hits"]), int(b["n_proteins"] - b["n_hits"])],
                ]
            )

            chi2, chi2_p, _, _ = stats.chi2_contingency(contingency)
            fisher_odds, fisher_p = stats.fisher_exact(contingency)
            rows.append(
                {
                    "organism_a": org_a,
                    "organism_b": org_b,
                    "chi2_stat": chi2,
                    "chi2_p_value": chi2_p,
                    "fisher_odds_ratio": fisher_odds,
                    "fisher_p_value": fisher_p,
                    "contingency_table": contingency.tolist(),
                }
            )

        out = pd.DataFrame(rows)
        if not out.empty:
            out["fisher_fdr_q"] = StatisticalValidation.benjamini_hochberg(
                out["fisher_p_value"].values
            )["fdr_q_value"].values
            out["chi2_fdr_q"] = StatisticalValidation.benjamini_hochberg(
                out["chi2_p_value"].values
            )["fdr_q_value"].values
        return out

    def ortholog_conservation(
        self,
        hits_df: pd.DataFrame,
        ortholog_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """Quantify orthogroup-level conservation of IP-binding hits across organisms."""
        self._validate_columns(hits_df, self.REQUIRED_HIT_COLUMNS, "hits table")
        self._validate_columns(ortholog_df, self.REQUIRED_ORTHOLOG_COLUMNS, "ortholog table")

        hits = hits_df.copy()
        hits["is_hit"] = self._to_bool(hits["is_hit"])

        merged = ortholog_df.merge(
            hits[["organism", "protein_id", "is_hit"]],
            on=["organism", "protein_id"],
            how="left",
        ).fillna({"is_hit": False})

        grouped = merged.groupby("orthogroup")
        rows = []
        for orthogroup, grp in grouped:
            org_count = grp["organism"].nunique()
            hit_orgs = grp.loc[grp["is_hit"], "organism"].nunique()
            rows.append(
                {
                    "orthogroup": orthogroup,
                    "organisms_present": org_count,
                    "organisms_with_hits": hit_orgs,
                    "hit_conserved_all_3": int(hit_orgs == 3),
                    "hit_in_any": int(hit_orgs >= 1),
                }
            )

        summary = pd.DataFrame(rows)
        if summary.empty:
            return summary

        conserved = int(summary["hit_conserved_all_3"].sum())
        partial = int(((summary["organisms_with_hits"] >= 1) & (summary["organisms_with_hits"] < 3)).sum())
        absent = int((summary["organisms_with_hits"] == 0).sum())

        contingency = np.array([[conserved, partial], [partial, absent]])
        _, fisher_p = stats.fisher_exact(contingency)

        return summary.assign(conservation_fisher_p_value=fisher_p)

    @staticmethod
    def _load_gene2go(path: str) -> Dict[str, Set[str]]:
        associations = pd.read_csv(path)
        required = {"protein_id", "go_id"}
        missing = required.difference(associations.columns)
        if missing:
            raise ValueError(f"GO association file {path} missing columns: {sorted(missing)}")

        mapping: Dict[str, Set[str]] = {}
        for _, row in associations.iterrows():
            mapping.setdefault(str(row["protein_id"]), set()).add(str(row["go_id"]))
        return mapping

    def run_go_enrichment(
        self,
        hits_df: pd.DataFrame,
        go_association_files: Mapping[str, str],
        obo_path: str,
    ) -> pd.DataFrame:
        """Run GO enrichment for each organism using GOATOOLS."""
        try:
            from goatools.obo_parser import GODag  # type: ignore
            from goatools.go_enrichment import GOEnrichmentStudy  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "GOATOOLS is required for GO enrichment. Install with `pip install goatools`."
            ) from exc

        self._validate_columns(hits_df, self.REQUIRED_HIT_COLUMNS, "hits table")
        hits = hits_df.copy()
        hits["is_hit"] = self._to_bool(hits["is_hit"])

        obodag = GODag(obo_path)

        all_rows: List[Dict[str, object]] = []
        for organism, assoc_path in go_association_files.items():
            assoc = self._load_gene2go(assoc_path)
            population = sorted(assoc.keys())
            study = sorted(hits[(hits["organism"] == organism) & hits["is_hit"]]["protein_id"].unique())

            enricher = GOEnrichmentStudy(
                population,
                assoc,
                obodag,
                methods=["fdr_bh"],
            )
            results = enricher.run_study(study)
            for rec in results:
                all_rows.append(
                    {
                        "organism": organism,
                        "go_id": rec.GO,
                        "go_term": rec.name,
                        "namespace": rec.NS,
                        "study_count": rec.study_count,
                        "population_count": rec.pop_count,
                        "p_uncorrected": rec.p_uncorrected,
                        "fdr_bh": rec.p_fdr_bh,
                        "enrichment": rec.enrichment,
                    }
                )

        out = pd.DataFrame(all_rows)
        if not out.empty:
            out = out.sort_values(["organism", "fdr_bh", "p_uncorrected"]).reset_index(drop=True)
        return out

    def generate_figures(
        self,
        hit_rates_df: pd.DataFrame,
        ortholog_df: pd.DataFrame,
        go_df: pd.DataFrame,
        output_dir: str,
    ) -> None:
        """Generate publication-ready comparative figures."""
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)

        StatisticalValidation.publication_style()
        sns.set_context("paper", font_scale=1.0)

        self._plot_hit_rate_bar(hit_rates_df, output / "hit_rate_barplot.png")
        self._plot_ip6_scatter(hit_rates_df, output / "ip6_vs_hit_rate_scatter.png")
        self._plot_ortholog_venn(ortholog_df, output / "ortholog_hits_venn.png")
        self._plot_functional_heatmap(go_df, output / "functional_category_heatmap.png")

    @staticmethod
    def _plot_hit_rate_bar(hit_rates_df: pd.DataFrame, outpath: Path) -> None:
        fig, ax = plt.subplots(figsize=(4.4, 3.0))
        x = np.arange(len(hit_rates_df))
        y = hit_rates_df["hit_rate_percent"].values
        yerr_lower = y - hit_rates_df["ci_low_percent"].values
        yerr_upper = hit_rates_df["ci_high_percent"].values - y

        ax.bar(x, y, color=["#4C78A8", "#F58518", "#54A24B"], width=0.7)
        ax.errorbar(x, y, yerr=[yerr_lower, yerr_upper], fmt="none", ecolor="black", capsize=3)
        ax.set_xticks(x)
        ax.set_xticklabels(hit_rates_df["organism"].values, rotation=20, ha="right")
        ax.set_ylabel("Hit rate (% proteins)")
        ax.set_title("Cryptic IP-binding prevalence by organism")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)

    @staticmethod
    def _plot_ip6_scatter(hit_rates_df: pd.DataFrame, outpath: Path) -> None:
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        sns.regplot(
            data=hit_rates_df,
            x="ip6_uM",
            y="hit_rate_percent",
            ax=ax,
            scatter_kws={"s": 45, "color": "#4C78A8"},
            line_kws={"color": "#333333", "lw": 1.2},
            ci=None,
        )
        for _, row in hit_rates_df.iterrows():
            ax.text(row["ip6_uM"], row["hit_rate_percent"], row["organism"], fontsize=6)
        ax.set_xlabel("Cellular IP6 concentration (µM)")
        ax.set_ylabel("Hit rate (%)")
        ax.set_title("IP6 concentration vs cryptic-site prevalence")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)

    @staticmethod
    def _plot_ortholog_venn(ortholog_df: pd.DataFrame, outpath: Path) -> None:
        try:
            from matplotlib_venn import venn3
        except ImportError as exc:
            raise ImportError(
                "matplotlib-venn is required for Venn diagrams. Install with `pip install matplotlib-venn`."
            ) from exc

        sets = {}
        for organism in sorted(ortholog_df["organism"].unique()):
            subset = ortholog_df[(ortholog_df["organism"] == organism) & (ortholog_df["is_hit"])]
            sets[organism] = set(subset["orthogroup"].astype(str))

        if len(sets) != 3:
            raise ValueError("Venn diagram requires exactly 3 organisms in ortholog data.")

        labels = list(sets.keys())
        fig, ax = plt.subplots(figsize=(4.2, 3.8))
        venn3([sets[labels[0]], sets[labels[1]], sets[labels[2]]], set_labels=labels, ax=ax)
        ax.set_title("Orthogroups with cryptic IP-binding hits")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)

    @staticmethod
    def _plot_functional_heatmap(go_df: pd.DataFrame, outpath: Path, top_n: int = 15) -> None:
        if go_df.empty:
            return

        sig = go_df.dropna(subset=["fdr_bh"]).copy()
        sig["neglog10_fdr"] = -np.log10(np.clip(sig["fdr_bh"], 1e-300, None))

        top_terms = (
            sig.sort_values("fdr_bh")
            .groupby("go_term", as_index=False)
            .first()
            .sort_values("fdr_bh")
            .head(top_n)["go_term"]
            .tolist()
        )
        plot_df = sig[sig["go_term"].isin(top_terms)]
        heat = plot_df.pivot_table(
            index="go_term",
            columns="organism",
            values="neglog10_fdr",
            aggfunc="max",
            fill_value=0,
        )

        fig, ax = plt.subplots(figsize=(5.6, max(3.0, 0.28 * len(heat))))
        sns.heatmap(heat, cmap="viridis", linewidths=0.3, linecolor="white", ax=ax)
        ax.set_xlabel("Organism")
        ax.set_ylabel("GO term")
        ax.set_title("Functional enrichment of cryptic IP-binding hits")
        fig.tight_layout()
        fig.savefig(outpath)
        plt.close(fig)

    def run_pipeline(
        self,
        hits_df: pd.DataFrame,
        ip6_map: Mapping[str, float],
        ortholog_df: pd.DataFrame,
        go_association_files: Mapping[str, str],
        go_obo_path: str,
        output_dir: str,
    ) -> ComparativeResult:
        """Run full comparative pipeline and write publication tables/figures."""
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)

        hit_rates = self.compute_hit_rates(hits_df, ip6_map)
        corr = self.spearman_ip6_correlation(hit_rates)
        pairwise = self.pairwise_organism_tests(hit_rates)

        ortholog_summary = self.ortholog_conservation(hits_df, ortholog_df)
        ortholog_hits = ortholog_df.merge(
            hits_df[["organism", "protein_id", "is_hit"]],
            on=["organism", "protein_id"],
            how="left",
        ).fillna({"is_hit": False})
        ortholog_hits["is_hit"] = self._to_bool(ortholog_hits["is_hit"])

        go_enrichment = self.run_go_enrichment(hits_df, go_association_files, go_obo_path)

        self.generate_figures(hit_rates, ortholog_hits, go_enrichment, str(output))

        stats_table = hit_rates.merge(corr, how="cross")
        stats_table.to_csv(output / "comparative_statistics_table.csv", index=False)
        hit_rates.to_csv(output / "hit_rates.csv", index=False)
        corr.to_csv(output / "ip6_hit_rate_spearman.csv", index=False)
        pairwise.to_csv(output / "pairwise_organism_tests.csv", index=False)
        ortholog_summary.to_csv(output / "ortholog_conservation.csv", index=False)
        go_enrichment.to_csv(output / "go_enrichment_results.csv", index=False)

        return ComparativeResult(
            hit_rate_table=hit_rates,
            correlation_table=corr,
            pairwise_tests_table=pairwise,
            ortholog_conservation_table=ortholog_summary,
        )
