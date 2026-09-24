"""Phase 3: call hits from a proteome screen and compare proteomes.

The screen records *every* pocket with its full descriptor set, not only the
pockets that pass. Hit calling therefore happens here, after the fact, which
has two consequences worth the storage:

* The strict filter, the score threshold and the denominators can be varied
  without re-running ~40,000 structures, so a hit rate can be reported across a
  range of thresholds rather than at one.
* A protein that fails is still ranked, so "no candidates" comes with the list
  of what came closest and which criterion stopped each one.

Denominators matter for the comparative question. A proteome rich in
intrinsically disordered proteins - Dictyostelium, with its long poly-N/Q
tracts - has many models with no confidently predicted pocket at all, which
dilutes a per-proteome hit rate without saying anything about IP binding. Hit
rates are therefore reported over all screened proteins *and* over proteins
with at least one confidently modelled pocket, and the organism comparison is
adjusted for length and model confidence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy import stats

from .scorer import ScoringParameters


@dataclass(frozen=True)
class HitCriteria:
    """The strict candidate definition, applied per pocket.

    Defaults follow the project plan and the pilot screen: composite score at
    least 0.75, lining-residue SASA at most 10 A^2, at least four basic
    residues, a cavity sized for an inositol phosphate, and a confidently
    predicted lining (mean pLDDT >= 70).

    The volume window is the scorer's calibrated **cavity** window
    (300-1600 A^3), not the 300-800 A^3 ligand volume the pilot's hard filter
    still used. fpocket measures the cavity, which is systematically larger
    than the ligand; the real ADAR2 InsP6 pocket measures ~1525 A^3, so the
    old gate rejected the paradigm positive outright.
    """

    min_score: float = 0.75
    max_sasa: float = 10.0
    min_basic: int = 4
    min_volume: float = ScoringParameters().volume_optimum_low
    max_volume: float = ScoringParameters().volume_optimum_high
    min_plddt: float = 70.0
    #: Minimum depth to the convex hull (A); no hull gate when ``None``.
    min_hull_depth: Optional[float] = None

    def gates(self, pockets: pd.DataFrame) -> pd.DataFrame:
        """Boolean pass/fail for each criterion, one column per gate.

        A gate whose threshold is ``None`` (or infinite, for SASA) always passes.
        """
        volume = pockets["volume"].astype(float)
        gates = {
            "score": pockets["composite_score"].astype(float) >= self.min_score,
            "sasa": (
                pockets["sasa"].astype(float) <= self.max_sasa
                if self.max_sasa is not None and np.isfinite(self.max_sasa)
                else pd.Series(True, index=pockets.index)
            ),
            "basic": pockets["basic_residues"].astype(float) >= self.min_basic,
            "volume": (volume >= self.min_volume) & (volume <= self.max_volume),
            "plddt": pockets["plddt_mean"].astype(float).fillna(0.0) >= self.min_plddt,
        }
        if self.min_hull_depth is not None:
            gates["hull"] = pockets["hull_depth"].astype(float).fillna(-np.inf) >= self.min_hull_depth
        else:
            gates["hull"] = pd.Series(True, index=pockets.index)
        return pd.DataFrame(gates, index=pockets.index)

    def replace(self, **changes) -> "HitCriteria":
        values = dict(self.__dict__)
        values.update(changes)
        return HitCriteria(**values)

    def passes(self, pockets: pd.DataFrame) -> pd.Series:
        return self.gates(pockets).all(axis=1)


#: The plan's strict candidate definition (sections 2 and 11): as the pilot
#: screen applied it, with the cavity volume window corrected.
PLAN_CRITERIA = HitCriteria()

#: Criteria recalibrated on the Phase 1 controls **as the screen sees them** -
#: AlphaFold models, no ligand. There, ADAR2's InsP6 site scores 0.575 against
#: 0.45-0.51 for the four PH-domain sites, sits 17.3 A inside the convex hull
#: against 5.7-8.2 A, and has a lining SASA of 31 A^2, so the plan's 0.75 score
#: and 10 A^2 SASA gates reject the paradigm positive itself. The thresholds
#: sit between ADAR2 and the highest negative; they rest on one positive and
#: four negatives, and are reported beside the plan's criteria, not instead.
#: A test fails if a scorer change stops them separating the controls.
CALIBRATED_CRITERIA = HitCriteria(min_score=0.54, max_sasa=float("inf"), min_hull_depth=10.0)


GATE_ORDER = ("plddt", "volume", "basic", "sasa", "hull", "score")


def clopper_pearson(k: int, n: int, confidence: float = 0.95) -> tuple:
    """Exact binomial interval for ``k`` successes in ``n`` trials."""
    if n <= 0:
        return (float("nan"), float("nan"))
    alpha = 1.0 - confidence
    low = 0.0 if k == 0 else float(stats.beta.ppf(alpha / 2, k, n - k + 1))
    high = 1.0 if k == n else float(stats.beta.ppf(1 - alpha / 2, k + 1, n - k))
    return (low, high)


def protein_table(
    pockets: pd.DataFrame,
    screened: pd.DataFrame,
    criteria: HitCriteria = HitCriteria(),
) -> pd.DataFrame:
    """One row per screened protein: best pocket, hit status, blocking gate.

    Args:
        pockets: All pockets, with ``uniprot_id`` and the descriptor columns.
        screened: One row per protein that was screened successfully, with
            ``uniprot_id`` and optionally ``length`` and ``mean_plddt``.
        criteria: Hit definition.
    """
    base = screened.drop_duplicates("uniprot_id").set_index("uniprot_id")
    table = pd.DataFrame(index=base.index)
    for column in ("organism_key", "length", "mean_plddt", "fraction_plddt_70"):
        if column in base.columns:
            table[column] = base[column]

    if pockets.empty:
        table["n_pockets"] = 0
        table["n_confident_pockets"] = 0
        table["best_score"] = np.nan
        table["best_confident_score"] = np.nan
        table["n_passing_pockets"] = 0
        table["is_hit"] = False
        table["eligible"] = False
        table["blocking_gate"] = "no pocket"
        return table.reset_index()

    gates = criteria.gates(pockets)
    passing = gates.all(axis=1)
    frame = pockets.assign(_pass=passing, _confident=gates["plddt"])
    grouped = frame.groupby("uniprot_id")

    table["n_pockets"] = grouped.size().reindex(table.index).fillna(0).astype(int)
    table["n_confident_pockets"] = grouped["_confident"].sum().reindex(table.index).fillna(0).astype(int)
    table["best_score"] = grouped["composite_score"].max().reindex(table.index)
    confident = frame[frame["_confident"]]
    table["best_confident_score"] = (
        confident.groupby("uniprot_id")["composite_score"].max().reindex(table.index)
    )
    table["n_passing_pockets"] = grouped["_pass"].sum().reindex(table.index).fillna(0).astype(int)
    table["is_hit"] = table["n_passing_pockets"] > 0
    table["eligible"] = table["n_confident_pockets"] > 0
    table["blocking_gate"] = _blocking_gate(frame, gates).reindex(table.index).fillna("no pocket")
    return table.reset_index()


def _blocking_gate(frame: pd.DataFrame, gates: pd.DataFrame) -> pd.Series:
    """For each protein, the first gate its closest pocket fails.

    "Closest" is the pocket failing the fewest gates, ties broken by score.
    Gates are checked in :data:`GATE_ORDER`, cheapest-to-explain first, so the
    answer reads as "stopped by low confidence" before "stopped by score".
    """
    n_failed = (~gates).sum(axis=1)
    ordered = frame.assign(_n_failed=n_failed).sort_values(
        ["uniprot_id", "_n_failed", "composite_score"], ascending=[True, True, False]
    )
    closest = ordered.groupby("uniprot_id").head(1)
    labels = {}
    for index, row in closest.iterrows():
        failed = [gate for gate in GATE_ORDER if not bool(gates.loc[index, gate])]
        labels[row["uniprot_id"]] = failed[0] if failed else "passes"
    return pd.Series(labels)


def hit_rates(proteins: pd.DataFrame, confidence: float = 0.95) -> pd.DataFrame:
    """Hit rate per organism, over all screened and over eligible proteins."""
    rows = []
    for organism, group in proteins.groupby("organism_key"):
        k = int(group["is_hit"].sum())
        n = int(len(group))
        eligible = group[group["eligible"]]
        n_e = int(len(eligible))
        k_e = int(eligible["is_hit"].sum())
        lo, hi = clopper_pearson(k, n, confidence)
        lo_e, hi_e = clopper_pearson(k_e, n_e, confidence)
        rows.append(
            {
                "organism_key": organism,
                "screened": n,
                "hits": k,
                "hit_rate": k / n if n else float("nan"),
                "ci_low": lo,
                "ci_high": hi,
                "eligible": n_e,
                "hit_rate_eligible": k_e / n_e if n_e else float("nan"),
                "ci_low_eligible": lo_e,
                "ci_high_eligible": hi_e,
            }
        )
    return pd.DataFrame(rows)


def threshold_sweep(
    pockets: pd.DataFrame,
    screened: pd.DataFrame,
    thresholds: Sequence[float],
    criteria: HitCriteria = HitCriteria(),
) -> pd.DataFrame:
    """Hit rate per organism across score thresholds, other gates fixed."""
    frames = []
    for threshold in thresholds:
        swept = criteria.replace(min_score=float(threshold))
        rates = hit_rates(protein_table(pockets, screened, swept))
        rates.insert(0, "min_score", float(threshold))
        frames.append(rates)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def pairwise_fisher(proteins: pd.DataFrame, reference: str = "dictyostelium") -> pd.DataFrame:
    """Unadjusted comparison of ``reference`` against each other organism."""
    rows = []
    ref = proteins[proteins["organism_key"] == reference]
    for organism, group in proteins.groupby("organism_key"):
        if organism == reference:
            continue
        for scope, a, b in (
            ("screened", ref, group),
            ("eligible", ref[ref["eligible"]], group[group["eligible"]]),
        ):
            table = [
                [int(a["is_hit"].sum()), int((~a["is_hit"]).sum())],
                [int(b["is_hit"].sum()), int((~b["is_hit"]).sum())],
            ]
            odds, p = stats.fisher_exact(table)
            rows.append(
                {
                    "comparison": f"{reference} vs {organism}",
                    "denominator": scope,
                    "hits_ref": table[0][0],
                    "n_ref": sum(table[0]),
                    "hits_other": table[1][0],
                    "n_other": sum(table[1]),
                    "odds_ratio": float(odds),
                    "p_value": float(p),
                }
            )
    return pd.DataFrame(rows)


def adjusted_comparison(
    proteins: pd.DataFrame,
    reference: str = "dictyostelium",
    *,
    score_column: str = "best_confident_score",
) -> Dict[str, object]:
    """Compare organisms after adjusting for protein length and model confidence.

    Two models, because hits may be too rare for the first:

    * **Hit status** (logistic): ``is_hit ~ organism + log(length) +
      fraction_plddt_70``. Only fitted when every organism has at least one hit
      and one non-hit; otherwise reported as not estimable, which is itself the
      honest answer when hits number in single digits.
    * **Best confident pocket score** (linear, eligible proteins only): the
      same covariates, on a continuous outcome that exists for every eligible
      protein. It asks whether one proteome's pockets look more like buried IP
      sites overall, which remains answerable when hits are rare.
    """
    try:
        import statsmodels.formula.api as smf
    except ImportError:  # pragma: no cover - optional at import time
        return {"available": False, "reason": "statsmodels not installed"}

    data = proteins.copy()
    data = data[data["length"].astype(float) > 0]
    data["log_length"] = np.log(data["length"].astype(float))
    data["fraction_plddt_70"] = data["fraction_plddt_70"].astype(float)
    data["is_hit_int"] = data["is_hit"].astype(int)
    formula_rhs = (
        f"C(organism_key, Treatment(reference='{reference}')) + log_length + fraction_plddt_70"
    )

    result: Dict[str, object] = {"available": True, "reference": reference}

    counts = data.groupby("organism_key")["is_hit_int"].agg(["sum", "count"])
    estimable = bool(((counts["sum"] > 0) & (counts["sum"] < counts["count"])).all())
    if estimable and data["organism_key"].nunique() > 1:
        try:
            fit = smf.logit(f"is_hit_int ~ {formula_rhs}", data=data).fit(disp=False)
            result["hit_logit"] = _coefficients(fit, exponentiate=True)
        except Exception as exc:  # separation, singular design
            result["hit_logit"] = {"estimable": False, "reason": f"{type(exc).__name__}: {exc}"}
    else:
        result["hit_logit"] = {
            "estimable": False,
            "reason": "an organism has no hits (or only hits); odds ratios are not identifiable",
            "hits_by_organism": {k: int(v) for k, v in counts["sum"].items()},
        }

    eligible = data[data["eligible"] & data[score_column].notna()]
    if eligible["organism_key"].nunique() > 1:
        fit = smf.ols(f"{score_column} ~ {formula_rhs}", data=eligible).fit(cov_type="HC3")
        result["best_score_ols"] = _coefficients(fit, exponentiate=False)
        result["best_score_ols"]["n"] = int(len(eligible))
    return result


def _coefficients(fit, *, exponentiate: bool) -> Dict[str, object]:
    params = fit.params
    conf = fit.conf_int()
    out = {"estimable": True, "terms": {}}
    for name in params.index:
        estimate, low, high = float(params[name]), float(conf.loc[name, 0]), float(conf.loc[name, 1])
        if exponentiate:
            estimate, low, high = math.exp(estimate), math.exp(low), math.exp(high)
        out["terms"][_clean_term(name)] = {
            "estimate": estimate,
            "ci_low": low,
            "ci_high": high,
            "p_value": float(fit.pvalues[name]),
        }
    out["scale"] = "odds ratio" if exponentiate else "difference in best score"
    return out


def _clean_term(name: str) -> str:
    if "[T." in name:
        return "organism=" + name.split("[T.", 1)[1].rstrip("]")
    return name


def keyword_enrichment(
    proteins: pd.DataFrame,
    annotations: pd.DataFrame,
    selected: pd.Series,
    keywords: Iterable[str],
) -> pd.DataFrame:
    """Fisher test of UniProt keywords in a selected set against the screened rest.

    Args:
        proteins: Protein table (must include ``uniprot_id``).
        annotations: UniProt table with ``uniprot_id`` and ``keywords``
            (semicolon-separated).
        selected: Boolean mask over ``proteins`` marking the set tested.
        keywords: Keywords to test, matched case-insensitively as whole entries.
    """
    kw = annotations.set_index("uniprot_id")["keywords"].reindex(proteins["uniprot_id"]).fillna("")
    kw_sets = kw.map(lambda s: {part.strip().lower() for part in str(s).split(";") if part.strip()})
    selected = selected.to_numpy(dtype=bool)
    rows = []
    for keyword in keywords:
        has = kw_sets.map(lambda parts: keyword.lower() in parts).to_numpy(dtype=bool)
        a = int((has & selected).sum())
        b = int((~has & selected).sum())
        c = int((has & ~selected).sum())
        d = int((~has & ~selected).sum())
        odds, p = stats.fisher_exact([[a, b], [c, d]]) if (a + b) and (c + d) else (float("nan"), float("nan"))
        rows.append(
            {
                "keyword": keyword,
                "in_selected": a,
                "selected": a + b,
                "in_background": c,
                "background": c + d,
                "odds_ratio": float(odds),
                "p_value": float(p),
            }
        )
    frame = pd.DataFrame(rows)
    if not frame.empty:
        frame["p_holm"] = _holm(frame["p_value"].to_numpy())
    return frame


def _holm(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(np.where(np.isnan(p), np.inf, p))
    adjusted = np.full_like(p, np.nan)
    running = 0.0
    m = int(np.sum(~np.isnan(p)))
    for rank, index in enumerate(order):
        if np.isnan(p[index]):
            continue
        running = max(running, min(1.0, (m - rank) * p[index]))
        adjusted[index] = running
    return adjusted


def rank_percentile(proteins: pd.DataFrame, uniprot_id: str, column: str = "best_confident_score") -> Optional[float]:
    """Fraction of the protein's own proteome it outranks (1.0 = top)."""
    row = proteins[proteins["uniprot_id"] == uniprot_id]
    if row.empty or pd.isna(row[column].iloc[0]):
        return None
    same = proteins[proteins["organism_key"] == row["organism_key"].iloc[0]][column].dropna()
    return float((same < row[column].iloc[0]).mean())


def known_ip_annotation(annotations: pd.DataFrame) -> pd.Series:
    """Proteins whose UniProt binding-site annotation names an inositol phosphate.

    The automated form of the plan's literature check: a candidate already
    annotated as binding an inositol phosphate validates the pipeline rather
    than being a new finding.
    """
    text = (
        annotations.get("binding_site", pd.Series("", index=annotations.index)).fillna("")
        + " "
        + annotations.get("function", pd.Series("", index=annotations.index)).fillna("")
    ).str.lower()
    flagged = text.str.contains(r"inositol (?:hexakis|pentakis|tetrakis|tris)phosphate|phytate|insp[3-8]\b", regex=True)
    return pd.Series(flagged.to_numpy(), index=annotations["uniprot_id"])
