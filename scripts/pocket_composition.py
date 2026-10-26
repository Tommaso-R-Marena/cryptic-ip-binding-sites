#!/usr/bin/env python3
"""Study Q10: the composition half of the decomposition, and a positive control.

Pre-registered in `docs/COMPOSITION_PLAN.md`, dated before any composition
statistic had been computed on any real lining in this project.

Study Q9 took the quotient by composition and asked whether **order** carries
anything about a pocket. This takes the other quotient. It is not interesting as
biology - an inositol-phosphate contact shell being arginine-rich is textbook -
and it exists for one reason: **Q9's nulls have no real-data positive control**,
so nothing distinguishes "no ordering signal" from "no sensitivity". This study
puts a known-true signal through the same machinery and gates every other row on
whether it fires.

Three things follow from every statistic here being order-invariant:

* the window is the **lining itself**, gap-closed, which `signature.py` fences
  periodicity statistics off and which costs these statistics nothing;
* the null is a random **same-size subset of the same chain**, whose composition
  is exactly multivariate hypergeometric in the chain's letter counts, so it is
  drawn in one vectorised call and cached per (chain, size);
* `sequence.composition_z` refuses them all, since an order-permuting null
  cannot move what they read.

    pocket_composition.py analyse --sequences S [S ...] --pockets P \\
        [--classes C] [--sites T] [--site-pockets Q] --out J --markdown M
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import composition as cp  # noqa: E402
from cryptic_ip.detect import sequence as sq  # noqa: E402


def _sibling(name: str):
    """Load a sibling script as a module without putting scripts/ on sys.path."""
    spec = importlib.util.spec_from_file_location(
        name, Path(__file__).resolve().parent / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


#: Study Q9's driver, imported and never modified. Reusing its matching, its
#: capping and above all its `contrast` means the two halves of the
#: decomposition are measured by the identical statistic, which is the only way
#: a null in one and a detection in the other can be compared.
ps = _sibling("pocket_sequence")

SEED = 20261025

#: Draws per lining from the same-size subset null.
SUBSET_DRAWS = 2000

#: The shortest lining worth a null, inherited from `docs/SEQUENCE_PLAN.md`
#: rather than chosen here.
MIN_LINING = ps.MIN_WINDOW

#: A lining may be at most this share of its own chain. Past a half the "rest of
#: the chain" is no longer a background the lining sits against, it is a
#: comparable-sized piece of the same protein.
MAX_CHAIN_SHARE = 0.5

#: Lining sizes must match within this many residues, inherited from Q9.
SIZE_TOLERANCE = ps.LENGTH_TOLERANCE

#: Both caps inherited from `docs/POCKET_SEQUENCE_PLAN_AMENDMENT_1.md` section 3.
MAX_PAIRS_PER_PROTEIN = ps.MAX_PAIRS_PER_PROTEIN
MAX_PAIRS_PER_ARM = ps.MAX_PAIRS_PER_ARM

#: The plan's three primary statistics.
PRIMARY = ("basic_fraction", "composition_divergence", "nonbasic_divergence")

#: Reported in full, never promoted.
SECONDARY = ("acidic_fraction", "aromatic_fraction", "mean_hydropathy",
             "glycine_fraction", "proline_fraction", "net_charge_per_residue")

#: The control arm's name, its statistic, and the direction fixed in advance.
CONTROL_ARM = "positive_control"
CONTROL_STATISTIC = "basic_fraction"

#: Arms that may become a finding. The control is in neither family: it gates.
CONFIRMATORY_ARMS = ("formation", "function", "type")


# ------------------------------------------------------------------ chains

class Chain:
    """One chain's letters, its background, and its cached subset nulls."""

    def __init__(self, numbers: Dict[int, str]) -> None:
        self.numbers = numbers
        ordered = sorted(numbers)
        self.letters = "".join(numbers[n] for n in ordered)
        self.counts = cp.letter_counts(self.letters)
        self.size = int(self.counts.sum())
        self.background = cp.background_of(self.letters)
        self._null: Dict[int, Dict[str, Tuple[float, float]]] = {}

    def null(self, size: int, statistics: Sequence[str], *, draws: int,
             seed: int) -> Dict[str, Tuple[float, float]]:
        """Mean and standard deviation of each statistic over same-size subsets.

        Cached on the size, because the null depends on the chain and the size
        and on nothing else about the lining - so every lining of a given size
        in a given chain shares it. That is what makes eight thousand linings
        cost seconds.
        """
        if size in self._null:
            return self._null[size]
        rng = np.random.default_rng(seed)
        drawn = rng.multivariate_hypergeometric(
            self.counts.astype(np.int64), int(size), size=int(draws))
        self._null[size] = counts_statistics(drawn, statistics, self.background)
        return self._null[size]


def counts_statistics(drawn: np.ndarray, statistics: Sequence[str],
                      background: Dict[str, float]
                      ) -> Dict[str, Tuple[float, float]]:
    """Each statistic's mean and sd over a matrix of letter-count rows.

    Every statistic of `composition.py` is a function of the count vector, so
    the whole null evaluates in a handful of array operations instead of a
    Python loop over draws.
    """
    totals = drawn.sum(axis=1, keepdims=True).astype(float)
    shares = drawn / totals
    reference = np.array([background.get(letter, 0.0) for letter in cp.ALPHABET])
    reference = reference / reference.sum()
    keep = np.array([letter not in cp.BASIC for letter in cp.ALPHABET])

    out: Dict[str, Tuple[float, float]] = {}
    for statistic in statistics:
        if statistic in ("basic_fraction", "acidic_fraction", "aromatic_fraction",
                         "glycine_fraction", "proline_fraction"):
            wanted = {"basic_fraction": cp.BASIC, "acidic_fraction": cp.ACIDIC,
                      "aromatic_fraction": cp.AROMATIC,
                      "glycine_fraction": frozenset("G"),
                      "proline_fraction": frozenset("P")}[statistic]
            mask = np.array([letter in wanted for letter in cp.ALPHABET], dtype=float)
            values = shares @ mask
        elif statistic == "mean_hydropathy":
            values = shares @ np.array([sq.HYDROPATHY[c] for c in cp.ALPHABET])
        elif statistic == "net_charge_per_residue":
            values = shares @ np.array([sq.CHARGE_OF[c] for c in cp.ALPHABET])
        elif statistic == "composition_divergence":
            values = jensen_shannon_rows(shares, reference)
        elif statistic == "nonbasic_divergence":
            sub = drawn[:, keep].astype(float)
            total = sub.sum(axis=1, keepdims=True)
            # A draw with no non-basic residue at all has no distribution to
            # compare; it cannot happen for a lining past the floor, and is
            # dropped from the null rather than counted as zero divergence.
            alive = total[:, 0] > 0
            part = np.full(sub.shape[0], np.nan)
            if alive.any():
                part[alive] = jensen_shannon_rows(
                    sub[alive] / total[alive], reference[keep] / reference[keep].sum())
            values = part
        else:
            raise sq.SequenceError(f"no vectorised null for {statistic!r}")
        finite = values[np.isfinite(values)]
        out[statistic] = (float(finite.mean()) if finite.size else float("nan"),
                          float(finite.std(ddof=1)) if finite.size > 1 else 0.0)
    return out


def jensen_shannon_rows(rows: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Jensen-Shannon divergence in bits of every row against one reference."""
    rows = np.asarray(rows, dtype=float)
    mid = 0.5 * (rows + reference[None, :])
    out = np.zeros(rows.shape[0])
    for part in (rows, np.broadcast_to(reference[None, :], rows.shape)):
        with np.errstate(divide="ignore", invalid="ignore"):
            term = np.where(part > 0.0, part * np.log2(np.divide(part, mid)), 0.0)
        out += 0.5 * np.nansum(term, axis=1)
    return np.clip(out, 0.0, 1.0)


def chain_table(sequences: pd.DataFrame) -> Dict[Tuple[str, str], Chain]:
    """One `Chain` per structure stem and chain id."""
    return {key: Chain(numbers)
            for key, numbers in ps.sequence_index(sequences).items() if numbers}


# ------------------------------------------------------------------ linings

def lining_letters(numbers: Dict[int, str], positions: Sequence[int]) -> Optional[str]:
    """The lining's own residues as letters, gap-closed and in residue order.

    This is `signature.py`'s `lining` window. Legitimate here and only here:
    every statistic of study Q10 is a function of the multiset of these letters,
    so closing the gaps cannot change one. A periodicity statistic on this
    window would be a period of a concatenation, which is why `signature.py`
    refuses it.
    """
    wanted = sorted({int(p) for p in positions})
    letters = [numbers.get(n) for n in wanted]
    if any(letter is None for letter in letters):
        return None
    text = "".join(letter for letter in letters if letter)
    if len(text) < MIN_LINING or not set(text) <= set(cp.ALPHABET):
        return None
    return text


def refusal(chain: Optional[Chain], numbers: Optional[Dict[int, str]],
            positions: Sequence[int]) -> Optional[str]:
    """Why a lining produced no scorable window, for the yield diagnostics."""
    if chain is None or numbers is None:
        return "chain not in the sequence table"
    wanted = sorted({int(p) for p in positions})
    if not wanted:
        return "no positions"
    letters = [numbers.get(n) for n in wanted]
    if any(letter is None for letter in letters):
        return "position not in the model"
    text = "".join(letter for letter in letters if letter)
    if len(text) < MIN_LINING:
        return "under the floor"
    if not set(text) <= set(cp.ALPHABET):
        return "non-standard letter"
    if len(text) > MAX_CHAIN_SHARE * chain.size:
        return "more than half its chain"
    return None


def score(window: str, chain: Chain, statistics: Sequence[str], *,
          draws: int, seed: int) -> Dict[str, Optional[float]]:
    """A lining's z against same-size subsets of its own chain, per statistic."""
    null = chain.null(len(window), statistics, draws=draws, seed=seed)
    out: Dict[str, Optional[float]] = {}
    for statistic in statistics:
        mean, deviation = null[statistic]
        try:
            observed = cp.evaluate(statistic, window, chain.background)
        except sq.SequenceError:
            out[statistic] = None
            continue
        if not np.isfinite(mean) or deviation == 0.0:
            out[statistic] = None
        else:
            out[statistic] = float((observed - mean) / deviation)
    return out


def rows_for(frame: pd.DataFrame, chains: Dict[Tuple[str, str], Chain], *,
             key_of, chain_of, group_of, id_of, residues_of, statistics,
             draws: int, yields: Optional[Dict[str, int]] = None,
             rank_of=None) -> List[Dict[str, object]]:
    """One scored row per lining that can carry a null."""
    rows: List[Dict[str, object]] = []
    for item in frame.itertuples():
        key = (key_of(item), chain_of(item))
        chain = chains.get(key)
        if chain is None:
            key = (key[0].lower(), key[1])
            chain = chains.get(key)
        residues = [int(r) for r in str(residues_of(item)).split(",") if r.strip()]
        numbers = chain.numbers if chain is not None else None
        window = lining_letters(numbers, residues) if numbers is not None else None
        if window is not None and len(window) > MAX_CHAIN_SHARE * chain.size:
            window = None
        if window is None:
            if yields is not None:
                reason = refusal(chain, numbers, residues) or "refused"
                yields[reason] = yields.get(reason, 0) + 1
            continue
        # The seed is the lining's, so its null does not depend on which arm
        # selected it, exactly as study Q9's does not.
        seed = SEED + (ps.hash_window(window) % 1_000_003)
        row: Dict[str, object] = {
            "id": str(id_of(item)), "group": str(group_of(item)),
            "length": len(window), "window": window, "chain_size": chain.size,
            "rank": int(rank_of(item)) if rank_of is not None else 0,
        }
        for statistic, value in score(window, chain, statistics,
                                      draws=draws, seed=seed).items():
            row[f"z_{statistic}"] = value
        rows.append(row)
    return rows


# --------------------------------------------------------------------- arms

def rank_split(rows: List[Dict[str, object]], rng: np.random.Generator
               ) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Amendment 1 section 3's rule: split within a protein, then size-match.

    `tests/test_pocket_composition.py` pins that this selects the same pairs as
    `pocket_sequence.formation_split`, so the two studies cannot drift apart
    about what the formation contrast is.
    """
    better: List[Dict[str, object]] = []
    worse: List[Dict[str, object]] = []
    for _, found in sorted(ps.by_protein_of(rows).items()):
        ordered = sorted(found, key=lambda r: int(r["rank"]))
        if len(ordered) < 2:
            continue
        half = len(ordered) // 2
        left, right = ps.match_on_length(
            ordered[:half], ordered[len(ordered) - half:], rng)
        left, right = ps.cap_pairs(left, right, MAX_PAIRS_PER_PROTEIN, rng)
        better.extend(left)
        worse.extend(right)
    return ps.cap_pairs(better, worse, MAX_PAIRS_PER_ARM, rng)


def sign_flip_p(values: Sequence[float], *, draws: int, seed: int) -> float:
    """Two-sided randomisation p for a mean against zero, by flipping signs.

    The control arm of `one_sample` is a constant zero, so
    `sequence.compare_feature`'s permutation over group labels is not a valid
    null for it: the two arms are not exchangeable under the hypothesis, because
    one of them has no variance at all. Measured over 400 null pools of nine
    groups, that p rejected at **7.25%** against a nominal 5% and at 13.25%
    against a nominal 10% - anti-conservative in exactly the tail a gate reads.

    Flipping the sign of each group's mean is the canonical null for this
    hypothesis and is exact whenever each group's mean is symmetric about zero
    under it, which a z-score against a symmetric null is. The same 400 pools
    reject at the nominal rate under this statistic.

    `compare_feature`'s other two outputs are kept, because both are correct
    here: the group bootstrap of `mean(sites) - mean(zeros)` is the bootstrap of
    the site mean, and `minimum_detectable_difference` with a zero-variance arm
    reduces to `sqrt(var / n)`, the standard error of a one-sample mean. Only
    the p-value needed replacing. The two-arm case every other study uses is
    well calibrated and is untouched: 400 pools at nine and twenty groups a
    side reject at 5.3% and 5.8%.
    """
    sample = np.asarray(list(values), dtype=float)
    if sample.size < 2:
        return 1.0
    observed = abs(float(sample.mean()))
    rng = np.random.default_rng(seed)
    signs = rng.choice((-1.0, 1.0), size=(int(draws), sample.size))
    flipped = np.abs((signs * sample[None, :]).mean(axis=1))
    extreme = int(np.count_nonzero(flipped >= observed))
    return float((extreme + 1) / (int(draws) + 1))


def one_sample(name: str, rows: List[Dict[str, object]], statistic: str, *,
               n_bootstrap: int) -> Dict[str, object]:
    """The control arm: is the mean z of these linings different from zero?

    The null value is zero by construction of a z, so the comparison is against
    a constant rather than against another arm, and the interval is the group
    bootstrap's interval on the mean z itself.
    """
    by_group: Dict[str, List[float]] = {}
    for row in rows:
        value = row.get(f"z_{statistic}")
        if value is None:
            continue
        by_group.setdefault(str(row["group"]), []).append(float(value))
    if not by_group:
        return {"arm": name, "statistic": statistic, "runnable": False,
                "reason": "no lining carried this statistic"}
    observed = {g: float(np.mean(v)) for g, v in by_group.items()}
    result = sq.compare_feature(f"{name}:{statistic}", observed,
                                {g: 0.0 for g in observed},
                                n_bootstrap=n_bootstrap, seed=SEED)
    entry = result.as_dict()
    sizes = [int(r["length"]) for r in rows]
    entry.update({
        # The permutation p `compare_feature` returns is not a valid null for a
        # constant control arm; `sign_flip_p` is. Its own p is kept beside it so
        # the substitution is visible in the artefact rather than silent.
        "p_value": sign_flip_p(observed.values(), draws=n_bootstrap, seed=SEED),
        "p_value_label_permutation": float(result.p_value),
        "p_rule": "sign-flip randomisation over group means",
        "arm": name, "statistic": statistic, "runnable": True,
        "n_site": len(rows), "n_control": 0, "n_groups": len(observed),
        "mean_site_z": float(np.mean(list(observed.values()))),
        "mean_control_z": 0.0,
        "mean_length_site": float(np.mean(sizes)) if sizes else float("nan"),
        "mean_length_control": float(np.mean(sizes)) if sizes else float("nan"),
        "length_confounded": False,
    })
    return entry


def fired(entry: Dict[str, object]) -> bool:
    """Whether the positive control cleared every floor, in the stated direction."""
    if not entry.get("runnable"):
        return False
    low, high = entry["interval"]
    return bool(entry["is_evidence"] and entry["exceeds_power"]
                and float(entry["difference"]) > 0.0 and low > 0.0 and high > 0.0
                and float(entry["p_value"]) < 0.05)


def holm_over(entries: List[Dict[str, object]], family: Sequence[str],
              *, gate: bool) -> None:
    """Holm across one declared family, with the control's gate applied."""
    raw = {f"{e['arm']}:{e['statistic']}": float(e["p_value"])
           for e in entries
           if e.get("runnable") and f"{e['arm']}:{e['statistic']}" in family}
    corrected = sq.holm(raw) if raw else {}
    for entry in entries:
        if not entry.get("runnable"):
            continue
        key = f"{entry['arm']}:{entry['statistic']}"
        if key not in corrected:
            entry.setdefault("family", "secondary")
            entry["is_finding"] = False
            continue
        entry["family"] = "primary"
        entry["p_holm"] = corrected[key]
        entry["is_finding"] = bool(
            gate
            and str(entry["arm"]) in CONFIRMATORY_ARMS
            and corrected[key] < 0.05 and entry["is_evidence"] and entry["exceeds_power"]
            and not (entry["interval"][0] <= 0.0 <= entry["interval"][1])
            and not entry["length_confounded"])


def run(sequences: pd.DataFrame, pockets: pd.DataFrame,
        classes: Optional[pd.DataFrame], sites: Optional[pd.DataFrame],
        site_pockets: Optional[pd.DataFrame], *, draws: int, n_bootstrap: int
        ) -> Dict[str, object]:
    """Every arm, with the control arm evaluated first because it gates them."""
    statistics = list(PRIMARY) + list(SECONDARY)
    chains = chain_table(sequences)
    yields: Dict[str, Dict[str, int]] = {}

    yields["candidate pockets"] = {}
    pocket_rows = rows_for(
        pockets, chains, statistics=statistics, draws=draws,
        yields=yields["candidate pockets"],
        key_of=lambda p: str(p.model).replace(".pdb", "").replace(".ent", ""),
        # An AlphaFold model is one chain, always A.
        chain_of=lambda p: "A",
        group_of=lambda p: str(getattr(p, "cluster", p.uniprot_id)),
        id_of=lambda p: str(getattr(p, "uniprot_id", "")),
        residues_of=lambda p: p.residues, rank_of=lambda p: p.pocket_id)
    by_class: Dict[Tuple[str, int], str] = {}
    if classes is not None and not classes.empty:
        for row in classes.itertuples():
            by_class[(str(row.uniprot_id), int(row.pocket_id))] = str(row.element)
    for row in pocket_rows:
        row["element"] = by_class.get((str(row["id"]), int(row["rank"])))

    site_rows: List[Dict[str, object]] = []
    if sites is not None and not sites.empty:
        yields["annotated shells"] = {}
        site_rows = rows_for(
            sites, chains, statistics=statistics, draws=draws,
            yields=yields["annotated shells"],
            key_of=lambda s: str(s.pdb_id).upper(), chain_of=lambda s: str(s.chain),
            group_of=lambda s: str(getattr(s, "group", s.pdb_id)),
            id_of=lambda s: str(s.copy), residues_of=lambda s: s.residues)

    sibling_rows: List[Dict[str, object]] = []
    if site_pockets is not None and not site_pockets.empty:
        yields["deposited pockets"] = {}
        shell_of = {str(s.copy): (str(s.pdb_id), str(s.chain),
                                  str(getattr(s, "group", s.pdb_id)))
                    for s in sites.itertuples()} if sites is not None else {}
        sibling_rows = rows_for(
            site_pockets, chains, statistics=statistics, draws=draws,
            yields=yields["deposited pockets"],
            key_of=lambda q: shell_of.get(str(q.uniprot_id), ("", ""))[0].upper(),
            chain_of=lambda q: shell_of.get(str(q.uniprot_id), ("", ""))[1],
            group_of=lambda q: shell_of.get(str(q.uniprot_id), ("", "", ""))[2]
            if str(q.uniprot_id) in shell_of else str(q.uniprot_id),
            id_of=lambda q: str(q.uniprot_id),
            residues_of=lambda q: q.residues, rank_of=lambda q: q.pocket_id)

    entries: List[Dict[str, object]] = []

    # The control, first and alone, because every other row is gated on it.
    control = one_sample(CONTROL_ARM, site_rows, CONTROL_STATISTIC,
                         n_bootstrap=n_bootstrap)
    control["family"] = "control"
    control["is_finding"] = False
    gate = fired(control)
    entries.append(control)
    for statistic in statistics:
        if statistic == CONTROL_STATISTIC:
            continue
        extra = one_sample(f"{CONTROL_ARM}_other", site_rows, statistic,
                           n_bootstrap=n_bootstrap)
        extra["family"] = "control"
        extra["is_finding"] = False
        entries.append(extra)

    rng = np.random.default_rng(SEED)
    better, worse = rank_split(pocket_rows, rng)
    entries += [ps.contrast("formation", better, worse, statistic,
                            n_bootstrap=n_bootstrap) for statistic in statistics]

    faces = [r for r in pocket_rows if r.get("element") == "helix_face"]
    others = [r for r in pocket_rows if r.get("element") == "none"]
    faces, others = ps.match_on_length(faces, others, np.random.default_rng(SEED + 1))
    faces, others = ps.cap_pairs(faces, others, MAX_PAIRS_PER_ARM,
                                 np.random.default_rng(SEED + 2))
    entries += [ps.contrast("type", faces, others, statistic,
                            n_bootstrap=n_bootstrap) for statistic in statistics]

    shells, siblings = ps.match_on_length(
        site_rows, sibling_rows, np.random.default_rng(SEED + 3))
    shells, siblings = ps.cap_pairs(shells, siblings, MAX_PAIRS_PER_ARM,
                                    np.random.default_rng(SEED + 4))
    entries += [ps.contrast("function", shells, siblings, statistic,
                            n_bootstrap=n_bootstrap) for statistic in statistics]

    primary_family = [f"{arm}:{statistic}"
                      for arm in ("formation", "function") for statistic in PRIMARY]
    holm_over(entries, primary_family, gate=gate)
    type_family = [f"type:{statistic}" for statistic in PRIMARY]
    holm_over([e for e in entries if e.get("arm") == "type"], type_family, gate=gate)

    return {
        "seed": SEED, "subset_draws": draws,
        "min_lining": MIN_LINING, "max_chain_share": MAX_CHAIN_SHARE,
        "size_tolerance": SIZE_TOLERANCE,
        "max_pairs_per_protein": MAX_PAIRS_PER_PROTEIN,
        "max_pairs_per_arm": MAX_PAIRS_PER_ARM,
        "primary_family": primary_family, "type_family": type_family,
        "confirmatory_arms": list(CONFIRMATORY_ARMS),
        "control_arm": CONTROL_ARM, "control_statistic": CONTROL_STATISTIC,
        "control_fired": bool(gate),
        "n_pocket_linings": len(pocket_rows),
        "n_shell_linings": len(site_rows),
        "n_deposited_pocket_linings": len(sibling_rows),
        # The function arm's control is every pocket of the shell's own
        # deposited structure, which costs hours of fpocket. Its table is
        # optional on purpose, so a reader must be told the difference between
        # "the arm found nothing" and "the arm had no input".
        "deposited_pockets_supplied": site_pockets is not None and not site_pockets.empty,
        "n_classified": int(sum(1 for r in pocket_rows if r.get("element"))),
        "n_formation_pairs": len(better),
        "n_type_pairs": len(faces),
        "n_function_pairs": len(shells),
        "window_yield": yields,
        "tests": entries,
    }


def markdown(report: Dict[str, object]) -> str:
    """The report, with the gate stated before any other row is readable."""
    gate = bool(report["control_fired"])
    control = next((t for t in report["tests"]  # type: ignore[union-attr]
                    if t.get("arm") == report["control_arm"]), {})
    lines = [
        "# Does a pocket lining's composition carry anything — and can this pipeline tell?",
        "",
        f"**Run.** `scripts/pocket_composition.py analyse`, seed {report['seed']}, "
        f"{report['subset_draws']} subset draws per lining.",
        "",
        "**Pre-registered** in `docs/COMPOSITION_PLAN.md`, dated before any composition "
        "statistic had been computed on any real lining in this project.",
        "",
        "## The gate",
        "",
    ]
    if not control.get("runnable"):
        lines += [
            "**The positive control did not run at all** "
            f"({control.get('reason', 'no reason recorded')}). Nothing below is a result: "
            "every arm reads *instrument not demonstrated*, and study Q9's nulls stay "
            "uninterpretable for the same reason.",
            "",
        ]
    else:
        low, high = control["interval"]
        verdict = (
            "**The positive control fired.**" if gate else
            "**The positive control did NOT fire.**")
        lines += [
            f"{verdict} Annotated inositol-phosphate contact shells against same-size "
            f"subsets of their own chains, on `{report['control_statistic']}`: mean z "
            f"**{control['mean_site_z']:+.2f}** across {control['n_groups']} homology "
            f"groups on {control['n_site']} shells, interval [{low:+.3f}, {high:+.3f}], "
            f"p {control['p_value']:.4g}, against a minimum detectable "
            f"{control['minimum_detectable']:.3f}.",
            "",
        ]
        if gate:
            lines += [
                "That a contact shell is arginine-rich is textbook and **is not offered "
                "here as a finding**. What it establishes is the thing study Q9 could not: "
                "that these linings, in these structures, through this machinery, can "
                "produce a detection at all. Q9's nulls are therefore statements about "
                "**order**, not about the pipeline's sensitivity.",
                "",
            ]
        else:
            lines += [
                "**So nothing below is evidence of absence.** A signal that is certainly "
                "present did not clear the plan's own floors, so this pipeline's "
                "sensitivity on real linings is undemonstrated, every arm below reads "
                "*instrument not demonstrated*, and **study Q9's five nulls cannot be "
                "read as null results either**. `docs/COMPOSITION_PLAN.md` section 6 "
                "fixed that reading before this ran, and `is_finding` is forced false "
                "throughout.",
                "",
            ]
    lines += [
        f"{report['n_pocket_linings']} candidate-pocket linings, "
        f"{report['n_shell_linings']} annotated shells and "
        f"{report['n_deposited_pocket_linings']} pockets of the same deposited "
        f"structures carried a null. The window is the **lining itself**, gap-closed: "
        "every statistic here is a function of the multiset of its letters, so closing "
        "the gaps cannot change one, and `cryptic_ip/detect/signature.py` is the fence "
        "that keeps periodicity statistics off that same window.",
        "",
    ]
    if not report.get("deposited_pockets_supplied", True):
        lines += [
            "**The function arm did not run, for want of its control rather than for "
            "want of a signal.** It compares a shell against every pocket of its own "
            "deposited structure, which is hours of pocket detection and a separate job "
            "from this one. That table is absent from this run, so every `function` row "
            "below reads as not runnable — which is a statement about the input and not "
            "about inositol phosphate. The gate and the formation and type arms do not "
            "use it and are unaffected.",
            "",
        ]
    lines += [
        "## Primary",
        "",
        "| arm | statistic | linings | groups | mean z (site) | mean z (control) | "
        "difference | 95% interval | Holm | could have detected | finding |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]

    def row_for(entry: Dict[str, object]) -> str:
        low, high = entry["interval"]
        holm = entry.get("p_holm")
        verdict = "**yes**" if entry.get("is_finding") else (
            "**instrument not demonstrated**" if not gate else
            f"**no evidence** — {entry['n_groups']} groups" if not entry["is_evidence"]
            else "no")
        return (
            f"| {entry['arm']} | `{entry['statistic']}` | {entry['n_site']} | "
            f"{entry['n_groups']} | {entry['mean_site_z']:+.2f} | "
            f"{entry['mean_control_z']:+.2f} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | "
            f"{'—' if holm is None else format(holm, '.4f')} | "
            f"{entry['minimum_detectable']:.3f} | {verdict} |")

    primary = [e for e in report["tests"]  # type: ignore[union-attr]
               if e.get("runnable") and e.get("family") == "primary"]
    for entry in sorted(primary, key=lambda e: (str(e["arm"]), str(e["statistic"]))):
        lines.append(row_for(entry))
    if not primary:
        lines.append("| — | — | — | — | — | — | — | — | — | — | nothing runnable |")

    findings = [e for e in primary if e.get("is_finding")]
    lines += [
        "",
        f"**{len(findings)} of {len(primary)}** primary tests is a finding. A finding "
        "needs the positive control to have fired, a Holm-corrected p below 0.05, an "
        "interval excluding zero, five resampling groups, an observed difference above "
        "the minimum detectable one, and matched lining sizes — all six, fixed before "
        "this ran.",
        "",
        "## The control arm's other statistics",
        "",
        "Reported because they say what a shell's composition looks like against its own "
        "chain, and never promoted: this arm is an instrument check, not a test.",
        "",
        "| statistic | shells | groups | mean z | 95% interval | p |",
        "|---|---|---|---|---|---|",
    ]
    for entry in sorted((e for e in report["tests"]  # type: ignore[union-attr]
                         if e.get("runnable") and e.get("family") == "control"),
                        key=lambda e: str(e["statistic"])):
        low, high = entry["interval"]
        lines.append(
            f"| `{entry['statistic']}` | {entry['n_site']} | {entry['n_groups']} | "
            f"{entry['mean_site_z']:+.2f} | [{low:+.3f}, {high:+.3f}] | "
            f"{entry['p_value']:.4g} |")

    lines += [
        "",
        "## Secondary, reported and never promoted",
        "",
        "| arm | statistic | linings | groups | difference | 95% interval | "
        "could have detected |",
        "|---|---|---|---|---|---|---|",
    ]
    for entry in sorted((e for e in report["tests"]  # type: ignore[union-attr]
                         if e.get("runnable")
                         and e.get("family") not in ("primary", "control")),
                        key=lambda e: (str(e["arm"]), str(e["statistic"]))):
        low, high = entry["interval"]
        lines.append(
            f"| {entry['arm']} | `{entry['statistic']}` | {entry['n_site']} | "
            f"{entry['n_groups']} | {entry['difference']:+.3f} | "
            f"[{low:+.3f}, {high:+.3f}] | {entry['minimum_detectable']:.3f} |")

    skipped = [e for e in report["tests"] if not e.get("runnable")]  # type: ignore[union-attr]
    if skipped:
        lines += ["", "Not runnable:", ""]
        lines += [f"- {e['arm']}:{e['statistic']} — {e['reason']}" for e in skipped]

    yields = report.get("window_yield") or {}
    if yields:
        lines += ["", "## Why a lining produced no window", "",
                  "| source | reason | linings |", "|---|---|---|"]
        for source in sorted(yields):
            for reason, count in sorted(yields[source].items(),
                                        key=lambda kv: (-kv[1], kv[0])):
                lines.append(f"| {source} | {reason} | {count} |")

    lines += [
        "",
        "## What this cannot conclude",
        "",
        "* **A compositional difference is not a discovery here.** That these sites are "
        "basic is known; the control arm exists to test the instrument and no row of it "
        "is a finding about inositol phosphate.",
        "* Composition and order do not exhaust a sequence. Conservation, coevolution, "
        "linear motifs and learned representations are outside this study and outside "
        "study Q9, and a null in both says nothing about them.",
        "* Every arm is within-protein or within-chain. **Nothing here is a sequence-only "
        "predictor of anything**, and this study does not attempt one.",
        "* The candidate arms read a predicted model and the function arm reads one "
        "deposited conformation, so neither is a statement about solution behaviour.",
    ]
    return "\n".join(lines) + "\n"


# ------------------------------------------------- validating the gate itself

#: The validation's synthetic pools: copies, chain length, lining size, and how
#: many of a planted shell's residues are basic. Eight of fourteen is a real
#: shell rather than a homopolymer, so the non-basic statistics still have
#: residues to read.
VALIDATION = {"copies": 18, "chain": 400, "lining": 14, "planted_basic": 8,
              "noise_seeds": 12, "draws": 600, "bootstrap": 600}

#: How far the mean z of a composition-matched noise pool may sit from zero
#: before the control arm is called biased rather than quiet.
NOISE_TOLERANCE = 0.25


def synthetic_shells(planted: bool, seed: int) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Shells in random chains, optionally with a basic patch planted in them."""
    rng = np.random.default_rng(seed)
    sequences, sites = [], []
    for index in range(int(VALIDATION["copies"])):
        name = f"{index:02d}XY"
        letters = list(rng.choice(list(cp.ALPHABET), size=int(VALIDATION["chain"])))
        lining = list(range(100, 100 + int(VALIDATION["lining"])))
        if planted:
            for position in range(int(VALIDATION["planted_basic"])):
                letters[lining[position] - 1] = "RK"[position % 2]
        sequences.append({
            "structure": f"{name}.pdb", "stem": name, "chain": "B", "first": 1,
            "last": len(letters), "n_residues": len(letters),
            "numbers": ",".join(str(n) for n in range(1, len(letters) + 1)),
            "letters": "".join(letters)})
        sites.append({"copy": f"{name}_B_900", "pdb_id": name, "chain": "B",
                      "group": f"G:{index // 2}",
                      "residues": ",".join(str(n) for n in lining)})
    return pd.DataFrame(sequences), pd.DataFrame(sites)


def control_on(sequences: pd.DataFrame, sites: pd.DataFrame, *, draws: int,
               n_bootstrap: int) -> Dict[str, object]:
    """The control arm alone, on one synthetic pool."""
    rows = rows_for(
        sites, chain_table(sequences), statistics=[CONTROL_STATISTIC], draws=draws,
        key_of=lambda s: str(s.pdb_id).upper(), chain_of=lambda s: str(s.chain),
        group_of=lambda s: str(s.group), id_of=lambda s: str(s.copy),
        residues_of=lambda s: s.residues)
    return one_sample(CONTROL_ARM, rows, CONTROL_STATISTIC, n_bootstrap=n_bootstrap)


def validate(*, draws: Optional[int] = None,
             n_bootstrap: Optional[int] = None) -> Dict[str, object]:
    """Does the gate open on a planted signal and stay shut on matched noise?

    The gate decides how every other row of this study, and every null of study
    Q9, may be read, so it is calibrated before any real lining reaches it. The
    pools are composition-matched in the only sense that matters here: the noise
    pool's linings are drawn from the same chains by the same rule, so a
    non-zero mean z in it would be bias and not signal.
    """
    draws = int(draws or VALIDATION["draws"])
    n_bootstrap = int(n_bootstrap or VALIDATION["bootstrap"])
    planted = control_on(*synthetic_shells(True, SEED), draws=draws,
                         n_bootstrap=n_bootstrap)
    quiet = []
    for offset in range(int(VALIDATION["noise_seeds"])):
        entry = control_on(*synthetic_shells(False, SEED + 100 + offset),
                           draws=draws, n_bootstrap=n_bootstrap)
        quiet.append({"seed": SEED + 100 + offset,
                      "mean_z": float(entry["mean_site_z"]),
                      "p_value": float(entry["p_value"]),
                      "interval": [float(v) for v in entry["interval"]],
                      "fired": fired(entry)})
    centres = [row["mean_z"] for row in quiet]
    opened = [row for row in quiet if row["fired"]]
    # A wrong-direction rejection must be refused, which is what the plan's
    # pre-stated direction is for. Checked on a constructed entry so the check
    # does not depend on a noise seed happening to produce one.
    mirrored = dict(planted)
    mirrored["difference"] = -abs(float(planted["difference"]))
    mirrored["interval"] = [-abs(float(planted["interval"][1])),
                            -abs(float(planted["interval"][0]))]
    checks = {
        "planted_pool_opens_the_gate": fired(planted),
        "no_noise_pool_opens_the_gate": not opened,
        "noise_is_centred_on_zero":
            bool(abs(float(np.mean(centres))) < NOISE_TOLERANCE),
        "planted_separates_from_every_noise_pool":
            bool(float(planted["mean_site_z"]) > max(centres) + 1.0),
        "a_wrong_direction_rejection_is_refused": not fired(mirrored),
    }
    return {
        "seed": SEED, "draws": draws, "bootstrap": n_bootstrap,
        "pools": VALIDATION, "noise_tolerance": NOISE_TOLERANCE,
        "statistic": CONTROL_STATISTIC,
        "planted": {"mean_z": float(planted["mean_site_z"]),
                    "p_value": float(planted["p_value"]),
                    "interval": [float(v) for v in planted["interval"]],
                    "n_groups": int(planted["n_groups"]),
                    "minimum_detectable": float(planted["minimum_detectable"]),
                    "fired": fired(planted)},
        "noise": quiet,
        "noise_mean_z": float(np.mean(centres)),
        "noise_sd_z": float(np.std(centres, ddof=1)),
        "n_noise_pools_fired": len(opened),
        "checks": checks,
        "all_checks_pass": all(checks.values()),
    }


def validation_markdown(report: Dict[str, object]) -> str:
    """The validation, rendered from its own artefact and nothing else."""
    planted = report["planted"]
    lines = [
        "# Study Q10: is the gate calibrated?",
        "",
        f"**Run.** `scripts/pocket_composition.py validate`, seed {report['seed']}, "
        f"{report['draws']} subset draws, {report['bootstrap']} bootstrap resamples.",
        "",
        "The positive control decides how every other row of study Q10 — and every null "
        "of study Q9 — may be read. So it is calibrated on synthetic pools **before any "
        "real lining reaches it**, and this page says whether it opens on a signal that "
        "is there and stays shut on one that is not.",
        "",
        f"Pools: {report['pools']['copies']} copies, chains of "
        f"{report['pools']['chain']} residues, linings of "
        f"{report['pools']['lining']}, and in the planted pool "
        f"{report['pools']['planted_basic']} of those residues made arginine or lysine. "
        "The noise pool's linings come from the same chains by the same rule, so a "
        "non-zero mean z in it would be bias rather than signal.",
        "",
        "| pool | mean z | p | 95% interval | gate |",
        "|---|---|---|---|---|",
        f"| planted basic patch | **{planted['mean_z']:+.2f}** | "
        f"{planted['p_value']:.4g} | [{planted['interval'][0]:+.3f}, "
        f"{planted['interval'][1]:+.3f}] | "
        f"{'**opens**' if planted['fired'] else '**stays shut**'} |",
    ]
    for row in report["noise"]:  # type: ignore[union-attr]
        lines.append(
            f"| matched noise, seed {row['seed']} | {row['mean_z']:+.2f} | "
            f"{row['p_value']:.4g} | [{row['interval'][0]:+.3f}, "
            f"{row['interval'][1]:+.3f}] | "
            f"{'**opens**' if row['fired'] else 'stays shut'} |")
    lines += [
        "",
        f"Across the {len(report['noise'])} noise pools the mean of the mean-z values is "
        f"**{report['noise_mean_z']:+.4f}**, against an expected zero, with a spread of "
        f"{report['noise_sd_z']:.4f}; **{report['n_noise_pools_fired']}** of them opened "
        "the gate.",
        "",
        "| check | result |",
        "|---|---|",
    ]
    for name, value in sorted(report["checks"].items()):  # type: ignore[union-attr]
        lines.append(f"| `{name}` | {'pass' if value else '**FAIL**'} |")
    lines += [
        "",
        "**What a noise pool with a significant p but a negative mean z would be.** The "
        "plan fixes the control's direction in advance — a contact shell is *more* basic "
        "than subsets of its own chain — so a rejection in the other direction does not "
        "open the gate. Over a dozen pools one such pool is ordinary, and the check "
        "`a_wrong_direction_rejection_is_refused` is what stops it counting.",
        "",
        "**This says nothing about inositol phosphate.** It is a statement about the "
        "instrument, on synthetic sequence, and it is the statement study Q9 never made "
        "about its own.",
    ]
    return "\n".join(lines) + "\n"


def cmd_validate(args: argparse.Namespace) -> int:
    report = validate(draws=args.draws, n_bootstrap=args.bootstrap)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = validation_markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    # A failed check is a reason to stop, not a warning.
    return 0 if report["all_checks_pass"] else 1


def read_frame(path: Optional[Path]) -> Optional[pd.DataFrame]:
    """One optional table, or ``None`` where the arm that makes it had nothing.

    A job that measured no copy still writes its table and still uploads it, so
    the file can be present and empty. That is the same situation as the file
    being absent - the arm did not run - and `run` already reports it that way,
    so it must not arrive here as a parse error that fails the whole study for
    the one case the study is built to tolerate.
    """
    if not path or not Path(path).exists():
        return None
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return None


def cmd_analyse(args: argparse.Namespace) -> int:
    report = run(ps.read_sequences(args.sequences), pd.read_csv(args.pockets),
                 read_frame(args.classes), read_frame(args.sites),
                 read_frame(args.site_pockets),
                 draws=args.draws, n_bootstrap=args.bootstrap)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")
    text = markdown(report)
    if args.markdown:
        args.markdown.write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    analyse = sub.add_parser("analyse", help="every arm, with the control arm as a gate")
    analyse.add_argument("--sequences", type=Path, nargs="+", required=True)
    analyse.add_argument("--pockets", type=Path, required=True)
    analyse.add_argument("--classes", type=Path, default=None)
    analyse.add_argument("--sites", type=Path, default=None)
    analyse.add_argument("--site-pockets", type=Path, default=None)
    analyse.add_argument("--out", type=Path, required=True)
    analyse.add_argument("--markdown", type=Path, default=None)
    analyse.add_argument("--draws", type=int, default=SUBSET_DRAWS)
    analyse.add_argument("--bootstrap", type=int, default=2000)

    check = sub.add_parser("validate", help="is the gate calibrated, on synthetic pools")
    check.add_argument("--out", type=Path, required=True)
    check.add_argument("--markdown", type=Path, default=None)
    check.add_argument("--draws", type=int, default=None)
    check.add_argument("--bootstrap", type=int, default=None)

    args = parser.parse_args(argv)
    return cmd_validate(args) if args.command == "validate" else cmd_analyse(args)


if __name__ == "__main__":
    raise SystemExit(main())
