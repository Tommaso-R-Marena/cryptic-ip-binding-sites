"""Cross-fitted re-ranking and its group-bootstrap evaluation (docs/RERANK_PLAN.md).

Each seed run is a pose list with Vina scores, RMSDs and electrostatic energies.
For a weight w the run's top pose is the argmin of ``vina + w * e_el`` (ties keep
Vina's order). A copy's success is the fraction of its runs whose top pose is
within 2 Å. The weight is chosen on four folds of strict homology groups and
applied to the fifth; the bootstrap repeats that choice inside every resample.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import pandas as pd

from ..benchmark import protocol
from ..docking.stats import MIN_GROUPS, per_copy_mean, per_group_mean

GRID = (0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
N_FOLDS = 5
SEED = 20260929
SUCCESS = 2.0
N_BOOTSTRAP = 2000
N_PERMUTATIONS = 100  # docs/RERANK_PLAN_AMENDMENT_1.md


@dataclass
class Run:
    copy_key: str
    group: str
    seed: int
    vina: np.ndarray
    rmsd: np.ndarray
    eel: np.ndarray


def top_index(run: Run, w: float) -> int:
    return int(np.argmin(run.vina + w * run.eel))  # argmin keeps the first of ties: Vina's order


def success_matrix(runs: Sequence[Run], copies: Sequence[str], grid: Sequence[float] = GRID,
                   cutoff: float = SUCCESS) -> np.ndarray:
    """(copies, grid): the fraction of each copy's runs whose top pose under w is within ``cutoff``."""
    index = {k: i for i, k in enumerate(copies)}
    hits = np.zeros((len(copies), len(grid)))
    counts = np.zeros(len(copies))
    for run in runs:
        i = index[run.copy_key]
        counts[i] += 1
        for j, w in enumerate(grid):
            hits[i, j] += float(run.rmsd[top_index(run, w)] <= cutoff)
    with np.errstate(invalid="ignore", divide="ignore"):
        return hits / counts[:, None]


def fold_of_groups(codes: np.ndarray, n_folds: int = N_FOLDS, seed: int = SEED) -> np.ndarray:
    """A fold per copy: a fixed random permutation of the groups, dealt round-robin."""
    n_groups = int(codes.max()) + 1 if len(codes) else 0
    perm = np.random.default_rng(seed).permutation(n_groups)
    fold = np.empty(n_groups, dtype=int)
    fold[perm] = np.arange(n_groups) % n_folds
    return fold[codes]


def crossfit(succ: np.ndarray, codes: np.ndarray, folds: np.ndarray, weights: Optional[np.ndarray] = None,
             n_folds: int = N_FOLDS):
    """Out-of-fold success per copy and the grid index chosen for each fold.

    Training maximises the group-equal mean success (ties: the smaller w). Rows of
    weight 0 do not train; the returned vector covers every row.
    """
    w = np.ones(len(succ)) if weights is None else np.asarray(weights, dtype=float)
    out = np.full(len(succ), np.nan)
    chosen: List[int] = []
    for k in range(n_folds):
        test = folds == k
        train = (~test) & (w > 0)
        if not train.any():
            j = 0
        else:
            objective = [per_group_mean(succ[train, g], codes[train], w[train]) for g in range(succ.shape[1])]
            objective = np.nan_to_num(np.asarray(objective), nan=-np.inf)
            j = int(np.argmax(objective))  # argmax keeps the first (smallest w) of ties
        chosen.append(j)
        out[test] = succ[test, j]
    return out, chosen


def f1_statistic(succ: np.ndarray, codes: np.ndarray, folds: np.ndarray, estimand: str = "group"):
    """A function of bootstrap weights: re-ranked minus Vina success, re-fitting w within the resample."""
    def stat(w: np.ndarray) -> float:
        oof, _ = crossfit(succ, codes, folds, w)
        diff = oof - succ[:, 0]
        return per_group_mean(diff, codes, w) if estimand == "group" else per_copy_mean(diff, w)
    return stat


def _mean_stat(y: np.ndarray, codes: np.ndarray, estimand: str):
    def stat(w: np.ndarray) -> float:
        return per_group_mean(y, codes, w) if estimand == "group" else per_copy_mean(y, w)
    return stat


def evaluate(runs: Sequence[Run], *, n_bootstrap: int = N_BOOTSTRAP, n_permutations: int = N_PERMUTATIONS,
             seed: int = SEED) -> Dict[str, object]:
    """F1 (with its decision, which uses F3), F2 and F3 on a set of seed runs."""
    runs = [r for r in runs if len(r.vina) and np.isfinite(r.rmsd).all()]
    copies = sorted({r.copy_key for r in runs})
    group_of = {r.copy_key: r.group for r in runs}
    groups = np.array([group_of[c] for c in copies])
    codes = pd.factorize(pd.Series(groups))[0]
    n_groups = int(len(set(groups)))
    out: Dict[str, object] = {"copies": len(copies), "groups": n_groups, "runs": len(runs), "grid": list(GRID)}
    if n_groups == 0:
        out["F1"] = {"decision": "not evaluable: no runs"}
        return out
    folds = fold_of_groups(codes)
    succ = success_matrix(runs, copies)
    oof, chosen = crossfit(succ, codes, folds)
    out["chosen_w_by_fold"] = [GRID[j] for j in chosen]
    out["success_by_w_in_sample"] = {str(w): per_group_mean(succ[:, j], codes, np.ones(len(copies)))
                                     for j, w in enumerate(GRID)}
    f1 = {}
    for estimand in ("group", "copy"):
        est = protocol.bootstrap_statistic(f1_statistic(succ, codes, folds, estimand), groups, null=0.0,
                                           n_bootstrap=n_bootstrap, seed=seed).as_dict()
        base = protocol.bootstrap_statistic(_mean_stat(succ[:, 0], codes, estimand), groups, null=0.5,
                                            n_bootstrap=n_bootstrap, seed=seed).as_dict()
        rer = protocol.bootstrap_statistic(_mean_stat(oof, codes, estimand), groups, null=0.5,
                                           n_bootstrap=n_bootstrap, seed=seed).as_dict()
        f1[estimand] = {"difference": est, "vina": base, "reranked": rer}
    out["F3"] = permutation_control(runs, copies, codes, folds, f1["group"]["difference"]["point"],
                                    n_permutations)
    out["F1"] = {**f1, "decision": decide(f1["group"]["difference"], out["F3"]["p_value"], n_groups)}
    out["F2"] = decompose(runs, copies, codes, chosen, folds, n_bootstrap, seed)
    return out


def decide(diff: Dict[str, float], p_permutation: float, n_groups: int) -> str:
    """F1 as amended (docs/RERANK_PLAN_AMENDMENT_1.md): the interval and the permutation control."""
    if n_groups < MIN_GROUPS:
        return f"not evaluable: {n_groups} groups (fewer than {MIN_GROUPS})"
    if diff["low"] > 0:
        return "improves" if p_permutation < 0.05 else "gain not specific to electrostatics"
    if diff["high"] < 0:
        return "worsens"
    return "no detectable difference"


def decompose(runs, copies, codes, chosen, folds, n_bootstrap, seed) -> Dict[str, object]:
    """Scoring versus sampling failures, per seed run, for Vina and for the out-of-fold weight."""
    index = {k: i for i, k in enumerate(copies)}
    rows = []
    for run in runs:
        i = index[run.copy_key]
        w = GRID[chosen[folds[i]]]
        sampled = bool(np.min(run.rmsd) <= SUCCESS)
        for arm, top in (("vina", top_index(run, 0.0)), ("reranked", top_index(run, w))):
            ok = bool(run.rmsd[top] <= SUCCESS)
            rows.append({"arm": arm, "ok": ok, "scoring_failure": (not ok) and sampled,
                         "sampling_failure": not sampled})
    frame = pd.DataFrame(rows)
    shares = {arm: {c: float(g[c].mean()) for c in ("ok", "scoring_failure", "sampling_failure")}
              for arm, g in frame.groupby("arm")}
    ceiling = np.zeros(len(copies))
    n = np.zeros(len(copies))
    for run in runs:
        ceiling[index[run.copy_key]] += float(np.min(run.rmsd) <= SUCCESS)
        n[index[run.copy_key]] += 1
    ceiling = ceiling / n
    groups = np.asarray(codes)
    est = protocol.bootstrap_statistic(_mean_stat(ceiling, codes, "group"), groups, null=0.5,
                                       n_bootstrap=n_bootstrap, seed=seed).as_dict()
    return {"per_run_shares": shares, "sampling_ceiling_group": est}


def permutation_control(runs, copies, codes, folds, observed: float, n_permutations: int) -> Dict[str, object]:
    """Shuffle E_el among each run's poses and recompute the F1 point estimate."""
    values = []
    for p in range(1, n_permutations + 1):
        rng = np.random.default_rng(p)
        shuffled = [Run(r.copy_key, r.group, r.seed, r.vina, r.rmsd, rng.permutation(r.eel)) for r in runs]
        succ = success_matrix(shuffled, copies)
        oof, _ = crossfit(succ, codes, folds)
        values.append(per_group_mean(oof - succ[:, 0], codes, np.ones(len(copies))))
    arr = np.asarray(values)
    at_least = int(np.sum(arr >= observed))
    return {"n": n_permutations, "values": [float(v) for v in arr], "mean": float(arr.mean()),
            "max": float(arr.max()), "fraction_at_least_observed": float(at_least / n_permutations),
            "p_value": float((1 + at_least) / (1 + n_permutations))}


def stratum(runs: Sequence[Run], keys: set, chosen_w: Dict[str, float]) -> Dict[str, object]:
    """Descriptive: Vina and out-of-fold re-ranked success on a subset of copies."""
    sub = [r for r in runs if r.copy_key in keys]
    copies = sorted({r.copy_key for r in sub})
    if not copies:
        return {"copies": 0, "groups": 0, "evidence": False}
    group_of = {r.copy_key: r.group for r in sub}
    codes = pd.factorize(pd.Series([group_of[c] for c in copies]))[0]
    index = {k: i for i, k in enumerate(copies)}
    vina, rer, n = np.zeros(len(copies)), np.zeros(len(copies)), np.zeros(len(copies))
    for r in sub:
        i = index[r.copy_key]
        vina[i] += float(r.rmsd[top_index(r, 0.0)] <= SUCCESS)
        rer[i] += float(r.rmsd[top_index(r, chosen_w[r.copy_key])] <= SUCCESS)
        n[i] += 1
    ones = np.ones(len(copies))
    n_groups = int(codes.max()) + 1
    return {"copies": len(copies), "groups": n_groups, "evidence": n_groups >= MIN_GROUPS,
            "vina": per_group_mean(vina / n, codes, ones), "reranked": per_group_mean(rer / n, codes, ones)}
