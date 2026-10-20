#!/usr/bin/env python3
"""Study Q8: is the candidate pocket different from its own protein's pockets?

Pre-registered in `docs/WITHIN_PROTEIN_PLAN.md`, dated before any structure was
fetched for this study.

Two amendments to the architecture plan name the same missing control - other
pockets from the same proteins - and record that no committed table carries
residue lists for them. That was verified, and it holds. So this study makes
the control instead of wishing for it: the same detector that produced the
candidates re-detects **every** pocket in the same proteins and keeps each
pocket's lining residues.

    within_protein.py pockets  --out pockets.csv.gz     # fetches models, needs fpocket
    within_protein.py analyse  --pockets pockets.csv.gz --out within_protein.json

The two commands are separate because only the first needs the network and
fpocket; the second is a pure function of two committed tables and is tested
offline.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.detect import architecture as ar  # noqa: E402
from cryptic_ip.detect import sequence as sq  # noqa: E402

SEED = 20261020

#: The candidate pocket is re-identified by overlap, because pocket indices do
#: not survive a re-run: a newer model release or a different fpocket build
#: renumbers them. Below this Jaccard the pocket being compared would not be
#: the pocket that was selected, so the protein is excluded and counted.
MIN_JACCARD = 0.5

#: Siblings must match the candidate's lining size to within this many
#: residues. Both halves of the architecture decomposition depend on n, so an
#: unmatched comparison would rediscover "the selected pocket is bigger".
SIZE_TOLERANCE = 2

#: Fewest matched siblings worth a quantile, and the same floor for the
#: *resolved* siblings of a half - a half resolved for only a handful of
#: siblings is unresolved for that protein, never a p of 1.
MIN_SIBLINGS = 8

#: Most siblings scored per protein, sampled with the seed above when more
#: qualify. A bound on compute, fixed in advance so it cannot be chosen to suit
#: a result.
MAX_SIBLINGS = 40

#: The halves, reported separately and never pooled.
HALVES = ("topology", "register")


# ------------------------------------------------------------ pocket table

def parse_residues(cell: object) -> List[int]:
    """Residue numbers from a delimited cell, sorted and de-duplicated."""
    if cell is None or (isinstance(cell, float) and np.isnan(cell)):
        return []
    text = str(cell).strip().strip("[]")
    if not text:
        return []
    out: Set[int] = set()
    for token in text.replace(";", ",").replace(" ", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            out.add(int(float(token)))
        except ValueError:
            continue
    return sorted(out)


def jaccard(left: Sequence[int], right: Sequence[int]) -> float:
    """Overlap of two residue sets, 0 when either is empty."""
    a, b = set(int(x) for x in left), set(int(x) for x in right)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def _fetch_model(uniprot_id: str, cache_dir: Path) -> Path:
    """The real AlphaFold fetch, injectable so the rest can be tested offline."""
    from cryptic_ip.database.alphafold_client import AlphaFoldClient

    return Path(AlphaFoldClient(cache_dir=cache_dir).fetch_structure(uniprot_id))


def _analyzer(path: Path, work: Path, timeout_s: float):
    """The real detector, injectable for the same reason."""
    from cryptic_ip.analysis import ProteinAnalyzer

    analyzer = ProteinAnalyzer(str(path), work_dir=str(work), skip_electrostatics=True)
    analyzer.fpocket_timeout_s = float(timeout_s)
    return analyzer


def detect_one(uniprot_id: str, cache_dir: Path, timeout_s: float = 900.0,
               fetch=_fetch_model, analyzer_factory=_analyzer) -> Dict[str, object]:
    """Fetch one AlphaFold model and return every pocket's lining residues.

    Errors are recorded rather than raised: one unavailable model must not end
    a run over 75 proteins, and a protein that failed has to be visible in the
    output instead of silently absent. The fetch and the detector are injected
    so that everything except the network and fpocket is tested offline.
    """
    record: Dict[str, object] = {"uniprot_id": uniprot_id, "model": "", "error": "",
                                 "n_pockets": 0, "pockets": []}
    work = Path(tempfile.mkdtemp(prefix=f"q8_{uniprot_id}_"))
    try:
        path = Path(fetch(uniprot_id, cache_dir))
        record["model"] = path.name
        # fpocket writes its output directory next to the structure it reads, so
        # the model is copied into this protein's own temporary directory first:
        # otherwise every run would leave a pocket directory beside the cached
        # model and the model cache would carry them from run to run.
        if path.exists():
            local = work / path.name
            shutil.copy(path, local)
            path = local
        analyzer = analyzer_factory(path, work, timeout_s)
        pockets = analyzer.detect_pockets()
        rows = []
        for pocket_id in pockets["pocket_id"]:
            try:
                residues = analyzer.get_pocket_residues(int(pocket_id))
            except Exception as error:  # noqa: BLE001 - one bad pocket, not one bad protein
                rows.append({"pocket_id": int(pocket_id), "residues": "",
                             "error": f"{type(error).__name__}: {error}"[:200]})
                continue
            rows.append({"pocket_id": int(pocket_id),
                         "residues": ",".join(str(int(r)) for r in sorted(set(residues))),
                         "error": ""})
        record["pockets"] = rows
        record["n_pockets"] = len(rows)
    except Exception as error:  # noqa: BLE001
        record["error"] = f"{type(error).__name__}: {error}"[:300]
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return record


def cmd_pockets(args: argparse.Namespace) -> int:
    frame = pd.read_csv(args.table)
    if "uniprot_id" not in frame.columns:
        raise SystemExit(f"{args.table} has no uniprot_id column")
    accessions = [str(a) for a in frame["uniprot_id"].tolist()]
    if args.limit:
        accessions = accessions[: int(args.limit)]
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    rows, status = [], []
    for index, accession in enumerate(accessions, start=1):
        record = detect_one(accession, args.cache_dir, timeout_s=args.fpocket_timeout)
        status.append({"uniprot_id": accession, "model": record["model"],
                       "n_pockets": record["n_pockets"], "error": record["error"]})
        for pocket in record["pockets"]:  # type: ignore[union-attr]
            rows.append({"uniprot_id": accession, "model": record["model"], **pocket})
        print(f"[{index}/{len(accessions)}] {accession}: {record['n_pockets']} pockets"
              f"{' ERROR ' + str(record['error']) if record['error'] else ''}", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    report = {"n_proteins": len(accessions),
              "n_with_pockets": int(sum(1 for s in status if s["n_pockets"])),
              "n_failed": int(sum(1 for s in status if s["error"])),
              "n_pockets": len(rows), "status": status}
    if args.report:
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "status"}, indent=2))
    return 0


# ------------------------------------------------- the ground-truth arm

#: A protein residue is in contact with the ligand when any heavy atom comes
#: this close to any of the ligand's. The census already uses it for protein
#: contacts; it is not tuned here.
CONTACT_ANGSTROMS = 4.5

#: The site pocket must cover at least this much of the annotated contact
#: shell. Recall rather than Jaccard, because a pocket may legitimately be
#: larger than a ligand's contact shell.
MIN_RECALL = 0.3


def contact_residues(arrays, chain: str, resseq: int, icode: str,
                     comp_id: str = "IHP",
                     cutoff: float = CONTACT_ANGSTROMS) -> Dict[str, object]:
    """The annotated site: polymer residues near one ligand copy, by chain.

    Returns the contacting residue numbers per chain, so the caller can apply
    the one-chain rule rather than pooling positions from several chains into a
    lining that is not a lining.
    """
    from scipy.spatial import cKDTree

    elements = np.asarray(arrays.elements)
    heavy = elements != "H"
    icodes = np.asarray([str(c).strip() for c in arrays.icodes])
    ligand = (np.asarray(arrays.resnames) == comp_id.upper()) & heavy \
        & (np.asarray(arrays.chain_ids) == chain) \
        & (np.asarray(arrays.resseqs).astype(int) == int(resseq)) \
        & (icodes == str(icode or "").strip())
    polymer = np.asarray(arrays.is_polymer) & heavy
    if not ligand.any() or not polymer.any():
        return {"found": False, "by_chain": {}}
    tree = cKDTree(np.asarray(arrays.coords)[polymer])
    hits: Set[int] = set()
    for index in tree.query_ball_point(np.asarray(arrays.coords)[ligand], r=float(cutoff)):
        hits.update(int(i) for i in index)
    if not hits:
        return {"found": True, "by_chain": {}}
    chains = np.asarray(arrays.chain_ids)[polymer]
    numbers = np.asarray(arrays.resseqs)[polymer].astype(int)
    by_chain: Dict[str, Set[int]] = {}
    for index in sorted(hits):
        by_chain.setdefault(str(chains[index]), set()).add(int(numbers[index]))
    return {"found": True, "by_chain": {k: sorted(v) for k, v in by_chain.items()}}


def recall(pocket: Sequence[int], site: Sequence[int]) -> float:
    """How much of the annotated site a pocket's lining covers."""
    a, b = set(int(x) for x in pocket), set(int(x) for x in site)
    if not b:
        return 0.0
    return len(a & b) / len(b)


def chain_aware_linings(pockets_dir: Path, pocket_id: int) -> Dict[str, List[int]]:
    """One pocket's lining residues **per chain**, from fpocket's own atom file.

    The detector's shared parser reads the residue number from columns 23-26 of
    exactly these records and drops the chain in column 22, which is why a
    homodimer's two residue 50s are indistinguishable downstream. Reading the
    one extra column gives a chain-aware lining with no heuristic and no second
    definition of what a lining is: it is still fpocket's own set of lining
    atoms. See `docs/WITHIN_PROTEIN_PLAN_AMENDMENT_2.md`.
    """
    path = pockets_dir / f"pocket{int(pocket_id)}_atm.pdb"
    out: Dict[str, Set[int]] = {}
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            chain = line[21:22].strip() or " "
            try:
                number = int(line[22:26].strip())
            except ValueError:
                continue
            out.setdefault(chain, set()).add(number)
    return {chain: sorted(numbers) for chain, numbers in out.items()}


def site_one(structure_path: Path, chain: str, resseq: int, icode: str, work: Path,
             comp_id: str = "IHP", timeout_s: float = 900.0,
             analyzer_factory=None) -> Dict[str, object]:
    """One deposited copy: its annotated site, and every pocket in the structure."""
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    factory = analyzer_factory or _analyzer
    record: Dict[str, object] = {"site": [], "pockets": [], "error": "",
                                 "chain": chain, "n_chains_in_contact": 0}
    try:
        arrays = load_structure_arrays(structure_path)
        contacts = contact_residues(arrays, chain, resseq, icode, comp_id=comp_id)
        if not contacts["found"]:
            record["error"] = "ligand copy not present in the structure"
            return record
        by_chain = contacts["by_chain"]  # type: ignore[assignment]
        record["n_chains_in_contact"] = len(by_chain)
        if chain not in by_chain:
            record["error"] = "no contact with the ligand's own chain"
            return record
        record["site"] = list(by_chain[chain])
        analyzer = factory(structure_path, work, timeout_s)
        pockets = analyzer.detect_pockets()
        # fpocket writes its output beside the structure it read, which is this
        # copy's own temporary directory.
        pockets_dir = Path(structure_path).parent / f"{Path(structure_path).stem}_out" / "pockets"
        rows, skipped = [], 0
        for pocket_id in pockets["pocket_id"]:
            by_pocket_chain = chain_aware_linings(pockets_dir, int(pocket_id))
            residues = by_pocket_chain.get(chain, [])
            if not residues:
                # A pocket with no residue in the ligand's chain is not a pocket
                # of that chain, and pooling chains would not give a lining.
                skipped += 1
                continue
            rows.append({"pocket_id": int(pocket_id),
                         "residues": ",".join(str(int(r)) for r in residues),
                         "error": ""})
        record["pockets"] = rows
        record["skipped_other_chain"] = skipped
    except Exception as error:  # noqa: BLE001
        record["error"] = f"{type(error).__name__}: {error}"[:300]
    return record


# --------------------------------------------------------------- analysis

def architecture_p(positions: Sequence[int], *, n_draws: int, seed: int) -> Dict[str, object]:
    """Both halves' Holm-corrected p for one lining, with None for unresolved.

    An unresolved half is not a p of 1. The architecture module separates
    "cannot be tested" from "no structure", and collapsing them produced a
    spurious result in the paired within-lining study; the distinction is kept
    here and pinned by this study's tests.
    """
    try:
        lining = ar.describe_architecture(positions, n_draws=n_draws, seed=seed)
    except ar.ArchitectureError:
        return {"topology": None, "register": None, "element": None, "size": len(set(positions))}
    return {"topology": lining.topology_p, "register": lining.register_p,
            "element": lining.element, "size": int(lining.n_positions),
            "architecture_class": lining.architecture_class}


def mid_p_quantile(observed: float, others: Sequence[float]) -> Tuple[float, float]:
    """Where `observed` sits among `others`, splitting ties, and the tied share.

    A **smaller** p means more structure, so the quantile counts how many
    siblings are *more* structured than the candidate: below 0.5 means the
    candidate is more structured than most of its protein's same-size pockets.
    Ties are split because permutation p-values are discrete and counting them
    as "at or below" pushes every quantile upward.
    """
    values = np.asarray(list(others), dtype=float)
    if values.size == 0:
        return float("nan"), float("nan")
    below = float(np.sum(values < observed - 1e-12))
    tied = float(np.sum(np.abs(values - observed) <= 1e-12))
    return float((below + 0.5 * tied) / values.size), float(tied / values.size)


def match_siblings(candidate_size: int, siblings: Sequence[Dict[str, object]],
                   rng: np.random.Generator) -> List[Dict[str, object]]:
    """Siblings within the size tolerance, capped by a seeded sample."""
    matched = [s for s in siblings
               if abs(int(s["size"]) - int(candidate_size)) <= SIZE_TOLERANCE]
    if len(matched) > MAX_SIBLINGS:
        index = rng.choice(len(matched), size=MAX_SIBLINGS, replace=False)
        matched = [matched[int(i)] for i in sorted(index)]
    return matched


#: How the committed lining is matched to a re-detected pocket, per arm. The
#: candidate arm matches a pocket to a pocket, so overlap is symmetric and
#: Jaccard is right; the ground-truth arm matches a pocket to a ligand's
#: contact shell, which a pocket may legitimately exceed, so it is recall.
MATCH_RULES = {"jaccard": (jaccard, MIN_JACCARD), "recall": (recall, MIN_RECALL)}


def measure(pockets: pd.DataFrame, table: pd.DataFrame, *, n_draws: int,
            n_bootstrap: int, rule: str = "jaccard") -> Dict[str, object]:
    """The selected lining against its own structure's same-size pockets.

    One function for both arms, because the arms' comparison means nothing
    unless they compute the same statistic. Only the rule that re-identifies
    the selected pocket differs, and it differs because the two arms match
    different kinds of object.
    """
    if rule not in MATCH_RULES:
        raise SystemExit(f"rule must be one of {sorted(MATCH_RULES)}, got {rule!r}")
    overlap_of, floor = MATCH_RULES[rule]
    rng = np.random.default_rng(SEED)
    per_protein: List[Dict[str, object]] = []
    excluded: Dict[str, int] = {}

    def exclude(reason: str) -> None:
        excluded[reason] = excluded.get(reason, 0) + 1

    groups = {str(k): v for k, v in pockets.groupby("uniprot_id")}
    for _index, row in table.iterrows():
        accession = str(row["uniprot_id"])
        committed = parse_residues(row.get("top_pocket_residues"))
        found = groups.get(accession)
        if found is None or not committed:
            exclude("no re-detected pockets for this protein")
            continue
        linings = [(int(r.pocket_id), parse_residues(r.residues))
                   for r in found.itertuples() if parse_residues(r.residues)]
        if not linings:
            exclude("every re-detected pocket has an empty lining")
            continue
        overlaps = [(overlap_of(lining, committed), pocket_id, lining)
                    for pocket_id, lining in linings]
        best_overlap, best_id, candidate_lining = max(overlaps, key=lambda t: (t[0], -t[1]))
        if best_overlap < floor:
            exclude(f"selected pocket not re-identified ({rule} below {floor})")
            continue

        candidate = architecture_p(candidate_lining, n_draws=n_draws, seed=SEED)
        siblings = []
        for pocket_id, lining in linings:
            if pocket_id == best_id:
                continue
            siblings.append({"pocket_id": pocket_id, "size": len(lining), "lining": lining})
        matched = match_siblings(int(candidate["size"]), siblings, rng)
        if len(matched) < MIN_SIBLINGS:
            exclude(f"fewer than {MIN_SIBLINGS} siblings within "
                    f"{SIZE_TOLERANCE} residues of the candidate")
            continue
        scored = [architecture_p(s["lining"], n_draws=n_draws, seed=SEED + 1 + i)
                  for i, s in enumerate(matched)]

        entry: Dict[str, object] = {
            "id": accession, "group": str(row.get("cluster", accession)),
            "overlap": float(best_overlap), "pocket_id": int(best_id),
            "candidate_size": int(candidate["size"]),
            "n_siblings": len(matched),
            "mean_sibling_size": float(np.mean([s["size"] for s in matched])),
            "size_difference": float(int(candidate["size"])
                                     - np.mean([s["size"] for s in matched])),
            "stratum": str(row.get("burial_class", "") or ""),
            "candidate_element": candidate.get("element"),
            "sibling_face_fraction": float(np.mean(
                [1.0 if s.get("element") in {"helix_face", "strand_face"} else 0.0
                 for s in scored])),
        }
        for half in HALVES:
            observed = candidate[half]
            resolved = [s[half] for s in scored if s[half] is not None]
            if observed is None or len(resolved) < MIN_SIBLINGS:
                continue
            quantile, tied = mid_p_quantile(float(observed), resolved)
            entry[f"{half}_quantile"] = quantile
            entry[f"{half}_tied_fraction"] = tied
            entry[f"{half}_resolved_siblings"] = len(resolved)
            entry[f"{half}_candidate_p"] = float(observed)
        per_protein.append(entry)

    tests, p_values = [], {}
    for half in HALVES:
        rows = [r for r in per_protein if f"{half}_quantile" in r]
        if not rows:
            tests.append({"half": half, "runnable": False,
                          "reason": "no protein has a resolved quantile for this half"})
            continue
        by_group: Dict[str, List[float]] = {}
        for row in rows:
            by_group.setdefault(str(row["group"]), []).append(float(row[f"{half}_quantile"]))
        observed = {g: float(np.mean(v)) for g, v in by_group.items()}
        # The null is a uniform quantile with mean 0.5, so the comparison arm is
        # a constant rather than a second sample.
        result = sq.compare_feature(half, observed, {g: 0.5 for g in observed},
                                    n_bootstrap=n_bootstrap, seed=SEED)
        entry = result.as_dict()
        entry.update({"half": half, "runnable": True,
                      "mean_quantile": float(np.mean(list(observed.values()))),
                      "n_proteins": len(rows), "n_groups": len(by_group),
                      "mean_tied_fraction": float(np.mean(
                          [r[f"{half}_tied_fraction"] for r in rows])),
                      "mean_size_difference": float(np.mean(
                          [r["size_difference"] for r in rows]))})
        tests.append(entry)
        p_values[half] = result.p_value

    corrected = sq.holm(p_values) if p_values else {}
    for entry in tests:
        if entry.get("runnable"):
            entry["p_holm"] = corrected.get(entry["half"])
            entry["is_finding"] = bool(
                entry["p_holm"] is not None and entry["p_holm"] < 0.05
                and entry["is_evidence"] and entry["exceeds_power"]
                and not (entry["interval"][0] <= 0.0 <= entry["interval"][1]))
    # Declared secondary: the same quantiles split by the census's burial class.
    # Reported, never promoted - a stratum with fewer than five groups gives no
    # evidence, and that floor is the project's standing rule.
    strata: Dict[str, Dict[str, object]] = {}
    for name in sorted({str(r.get("stratum", "")) for r in per_protein} - {""}):
        in_stratum = [r for r in per_protein if str(r.get("stratum", "")) == name]
        entry: Dict[str, object] = {
            "n_proteins": len(in_stratum),
            "n_groups": len({str(r["group"]) for r in in_stratum}),
        }
        for half in HALVES:
            values = [float(r[f"{half}_quantile"]) for r in in_stratum
                      if f"{half}_quantile" in r]
            entry[f"{half}_mean_quantile"] = float(np.mean(values)) if values else float("nan")
            entry[f"{half}_n"] = len(values)
        entry["is_evidence"] = bool(entry["n_groups"] >= 5)
        strata[name] = entry
    faces = [r for r in per_protein if r.get("candidate_element")]
    return {
        "seed": SEED,
        "rule": rule,
        "overlap_floor": float(floor),
        "min_jaccard": MIN_JACCARD,
        "size_tolerance": SIZE_TOLERANCE,
        "min_siblings": MIN_SIBLINGS,
        "max_siblings": MAX_SIBLINGS,
        "n_proteins_in_table": int(len(table)),
        "n_measured": len(per_protein),
        "excluded": dict(sorted(excluded.items())),
        "overlaps": {
            "median": float(np.median([r["overlap"] for r in per_protein]))
            if per_protein else float("nan"),
            "min": float(np.min([r["overlap"] for r in per_protein]))
            if per_protein else float("nan"),
        },
        "siblings": {
            "median": float(np.median([r["n_siblings"] for r in per_protein]))
            if per_protein else float("nan"),
        },
        "faces": {
            "candidate_face_fraction": float(np.mean(
                [1.0 if r["candidate_element"] in {"helix_face", "strand_face"} else 0.0
                 for r in faces])) if faces else float("nan"),
            "sibling_face_fraction": float(np.mean(
                [r["sibling_face_fraction"] for r in per_protein]))
            if per_protein else float("nan"),
        },
        "strata": strata,
        "per_protein": per_protein,
        "tests": tests,
        "holm_family_size": len(p_values),
    }


def _verdict(entry: Dict[str, object]) -> str:
    """The finding column, which must not read as a null when it is a power failure."""
    if entry.get("is_finding"):
        return "**yes**"
    if not entry.get("is_evidence"):
        return f"**no evidence** — {entry['n_groups']} groups"
    return "no"


def markdown(report: Dict[str, object]) -> str:
    """The report, with the design's one weak point stated where it is read."""
    excluded = report["excluded"]  # type: ignore[index]
    lines = [
        "# Is the candidate pocket different from its own protein's pockets?",
        "",
        "**Run.** `scripts/within_protein.py analyse`, seed "
        f"{report['seed']}, on pockets re-detected by `within_protein.py pockets`.",
        "",
        "**Pre-registered** in `docs/WITHIN_PROTEIN_PLAN.md`, dated before any structure "
        "was fetched for this study.",
        "",
        "## The control that did not exist",
        "",
        "`results/architecture/real/REAL.md` compares 75 candidate linings with random "
        "position sets matched on size and span, and says plainly that the comparison is "
        "nearly tautological: a random position set is not a protein, so any real lining "
        "beats that null. The control that would mean something is *other pockets from "
        "the same proteins*, and no committed table carries residue lists for them.",
        "",
        "So the control was built. The same detector re-detects every pocket in the same "
        "proteins and keeps each pocket's lining residues. Same protein, same fold, same "
        "model, same detector, same lining size - what varies is **which pocket**.",
        "",
        f"The candidate is re-identified by Jaccard overlap against the committed lining "
        f"(floor {report['min_jaccard']}), because pocket indices do not survive a "
        f"re-run. Median overlap {report['overlaps']['median']:.2f}, "  # type: ignore[index]
        f"lowest kept {report['overlaps']['min']:.2f}.",  # type: ignore[index]
        "",
        f"Siblings match the candidate's lining size to within "
        f"{report['size_tolerance']} residues, at least {report['min_siblings']} of them "
        f"and at most {report['max_siblings']}, sampled with the seed above. Median "
        f"{report['siblings']['median']:.0f} siblings per protein.",  # type: ignore[index]
        "",
        f"{report['n_measured']} of {report['n_proteins_in_table']} proteins measured.",
    ]
    measured = report.get("per_protein") or []
    if measured:
        sizes = sorted(int(r["candidate_size"]) for r in measured)
        lines += [
            "",
            f"**Which proteins survive, and why it matters.** The measured candidate "
            f"linings run {sizes[0]} to {sizes[-1]} residues (median "
            f"{sizes[len(sizes) // 2]}). The exclusions above are almost entirely a "
            "size effect: a large lining has few pockets of its own size to be compared "
            "with, so the proteins that clear the sibling floor are the smaller-lining "
            "end of the table. That is a property of the matched design rather than a "
            "result, it was not chosen after seeing anything, and loosening the "
            "tolerance afterwards to recover the rest would be exactly the tuning this "
            "plan fixed in advance to avoid.",
        ]
    if excluded:
        lines += ["", "Excluded, by cause:", ""]
        lines += [f"- {count} — {reason}" for reason, count in excluded.items()]  # type: ignore[union-attr]
    lines += [
        "",
        "## Result",
        "",
        "The per-protein statistic is the candidate's **mid-p quantile** among its "
        "matched siblings. A smaller architecture p means more structure, so the quantile "
        "counts how many siblings are *more* structured than the candidate: **below 0.5 "
        "means the candidate is more structured than its protein's same-size pockets**. "
        "An unresolved half is unresolved, never a p of 1, and ties are split because "
        "permutation p-values are discrete.",
        "",
        "| half | proteins | groups | ties | mean quantile | difference from 0.5 | "
        "95% interval | p | Holm | finding |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for entry in report["tests"]:  # type: ignore[union-attr]
        if not entry.get("runnable"):
            lines.append(f"| {entry['half']} | — | — | — | — | — | — | — | — | "
                         f"not runnable: {entry['reason']} |")
            continue
        low, high = entry["interval"]
        lines.append(
            f"| {entry['half']} | {entry['n_proteins']} | {entry['n_groups']} | "
            f"{entry['mean_tied_fraction']:.2f} | {entry['mean_quantile']:.3f} | "
            f"{entry['difference']:+.3f} | [{low:+.3f}, {high:+.3f}] | "
            f"{entry['p_value']:.4f} | {entry['p_holm']:.4f} | "
            f"{_verdict(entry)} |"
        )
    starved = [e for e in report["tests"]  # type: ignore[union-attr]
               if e.get("runnable") and not e.get("is_evidence")]
    if starved:
        lines += [
            "",
            "**This arm did not have the power to decide anything, and the table above "
            "is not a null result.** "
            + "; ".join(
                f"the {e['half']} half rests on {e['n_groups']} resampling groups and "
                f"could only have detected a difference of {e['minimum_detectable']:.2f}"
                for e in starved)
            + ". The project's standing rule is that a stratum with fewer than five "
            "groups yields no evidence, and that rule is applied here to the arm as a "
            "whole rather than quietly to a footnote. Whatever direction the means "
            "point in, this arm neither supports nor refutes the hypothesis.",
        ]
    faces = report["faces"]  # type: ignore[index]
    lines += [
        "",
        f"Holm across {report['holm_family_size']} primary tests, fixed by the plan "
        "before any of it was read.",
        "",
        "**Secondary.** Candidates called a helical or strand face: "
        f"{faces['candidate_face_fraction']:.2f}; their matched siblings: "  # type: ignore[index]
        f"{faces['sibling_face_fraction']:.2f}.",  # type: ignore[index]
        "",
        "**Residual size.** Matching is to within a tolerance, not exact, so the residual "
        "is reported: "
        + "; ".join(
            f"{e['half']} {e['mean_size_difference']:+.2f} residues"
            for e in report["tests"] if e.get("runnable"))  # type: ignore[union-attr]
        + ". A result accompanied by a residual difference in the same direction is "
        "confounded rather than a finding.",
    ]
    strata = report.get("strata") or {}
    if strata:
        lines += [
            "",
            "**By burial class**, declared secondary and never promoted: a stratum with "
            "fewer than five groups gives no evidence, which is the project's standing "
            "rule and not a judgement about this table.",
            "",
            "| class | proteins | groups | topology | register | evidence |",
            "|---|---|---|---|---|---|",
        ]
        for name, entry in strata.items():  # type: ignore[union-attr]
            lines.append(
                f"| {name} | {entry['n_proteins']} | {entry['n_groups']} | "
                f"{entry['topology_mean_quantile']:.3f} | "
                f"{entry['register_mean_quantile']:.3f} | "
                f"{'yes' if entry['is_evidence'] else 'no — too few groups'} |"
            )
    lines += [
        "",
        "## How to read a null here",
        "",
        "The transfer measurement in `results/architecture/ARCHITECTURE.md` stands over "
        "all of this: at the φ/ψ scatter real structures show, the method detects about "
        "**a third** of what it detects on idealised geometry. A null here is therefore "
        "weak evidence of absence, and the minimum detectable difference in the table is "
        "the honest statement of what this design could have seen.",
        "",
        "## What this cannot conclude",
        "",
        "* These are candidates, not annotated sites.",
        "* A sibling pocket is not a negative. The comparison is to the protein's own "
        "background, not to a labelled non-site.",
        "* AlphaFold models are predictions, and pocket detection inherits whatever the "
        "prediction got wrong — equally for candidate and sibling, which is the point of "
        "a within-protein design.",
        "* Nothing here involves an inositol phosphate.",
    ]
    return "\n".join(lines) + "\n"


def select_copies(census: pd.DataFrame, comp_id: str = "IHP") -> pd.DataFrame:
    """The copies the amendment names, by the census's own columns.

    The filter is the census's, not a new one: the right ligand, eligible,
    complete, not a symmetry contact, and the census's own selected
    representative per site. Kept in one function so the workflow's download
    list and the analysis cannot drift apart.
    """
    wanted = (census["comp_id"].astype(str).str.upper() == comp_id.upper())
    for column, value in (("status", "eligible"), ("complete", True),
                          ("selected", True)):
        if column in census.columns:
            wanted &= (census[column] == value)
    if "symmetry_contact" in census.columns:
        wanted &= (census["symmetry_contact"] != True)  # noqa: E712 - NaN must survive
    return census[wanted]


def cmd_sites(args: argparse.Namespace) -> int:
    """The ground-truth arm's data step: annotated sites and their structures' pockets.

    Writes the same two-column shape the candidate arm's `pockets` command
    writes, plus a table of annotated site linings, so that `analyse` runs the
    identical statistic on both arms.
    """
    copies = select_copies(pd.read_csv(args.census), args.comp_id)
    if args.limit:
        copies = copies.head(int(args.limit))
    rows, sites, status = [], [], []
    for index, row in enumerate(copies.itertuples(), start=1):
        pdb_id = str(row.pdb_id).lower()
        candidates = sorted(args.structures.glob(f"{pdb_id}.*")) + \
            sorted(args.structures.glob(f"{pdb_id.upper()}.*"))
        key = f"{row.pdb_id}_{row.chain}_{row.resseq}"
        if not candidates:
            status.append({"copy": key, "error": "structure not downloaded"})
            continue
        work = Path(tempfile.mkdtemp(prefix=f"q8site_{pdb_id}_"))
        try:
            local = work / candidates[0].name
            shutil.copy(candidates[0], local)
            icode = "" if pd.isna(getattr(row, "icode", "")) else str(getattr(row, "icode", ""))
            record = site_one(local, str(row.chain), int(row.resseq), icode, work,
                              comp_id=args.comp_id, timeout_s=args.fpocket_timeout)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        status.append({"copy": key, "error": record["error"],
                       "n_pockets": len(record["pockets"]),  # type: ignore[arg-type]
                       "n_site_residues": len(record["site"])})  # type: ignore[arg-type]
        if record["error"]:
            print(f"[{index}/{len(copies)}] {key}: {record['error']}", flush=True)
            continue
        sites.append({"copy": key, "pdb_id": row.pdb_id, "chain": row.chain,
                      "group": getattr(row, "homology_group_strict", row.pdb_id),
                      "burial_class": getattr(row, "burial_class", ""),
                      "residues": ",".join(str(r) for r in record["site"])})  # type: ignore[arg-type]
        for pocket in record["pockets"]:  # type: ignore[union-attr]
            rows.append({"uniprot_id": key, "model": candidates[0].name, **pocket})
        print(f"[{index}/{len(copies)}] {key}: {len(record['pockets'])} pockets, "
              f"{len(record['site'])} site residues", flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, index=False)
    pd.DataFrame(sites).to_csv(args.sites, index=False)
    report = {"n_copies": int(len(copies)), "n_measured": len(sites),
              "n_pockets": len(rows),
              "excluded": {},
              "status": status}
    for entry in status:
        if entry.get("error"):
            report["excluded"][entry["error"]] = report["excluded"].get(entry["error"], 0) + 1
    if args.report:
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    # The per-copy status goes in the artefact; the log gets the summary, so
    # that reading why copies were lost does not mean reading 137 entries.
    summary = {k: v for k, v in report.items() if k != "status"}
    if args.summary:
        args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


#: What each arm's committed table calls the three columns `measure` needs.
ARMS = {
    "candidate": {"rule": "jaccard",
                  "columns": {"uniprot_id": "uniprot_id",
                              "top_pocket_residues": "top_pocket_residues",
                              "cluster": "cluster"}},
    "site": {"rule": "recall",
             "columns": {"uniprot_id": "copy", "top_pocket_residues": "residues",
                         "cluster": "group"}},
}


def as_table(frame: pd.DataFrame, arm: str) -> pd.DataFrame:
    """One arm's table in the column names `measure` reads."""
    columns = ARMS[arm]["columns"]
    missing = [source for source in columns.values() if source not in frame.columns]
    if missing:
        raise SystemExit(f"the {arm} table is missing {missing}")
    out = frame.rename(columns={source: target for target, source in columns.items()})
    keep = [c for c in out.columns if c not in set(columns.values()) - set(columns)]
    return out[keep] if keep else out


def cmd_analyse(args: argparse.Namespace) -> int:
    pockets = pd.read_csv(args.pockets)
    if "uniprot_id" not in pockets.columns:
        raise SystemExit(f"{args.pockets} is missing uniprot_id")
    table = as_table(pd.read_csv(args.table), args.arm)
    report = measure(pockets, table, n_draws=args.draws, n_bootstrap=args.bootstrap,
                     rule=str(ARMS[args.arm]["rule"]))
    report["arm"] = args.arm
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

    detect = sub.add_parser("pockets", help="fetch models and re-detect every pocket")
    detect.add_argument("--table", type=Path,
                        default=Path("results/specificity/candidates.csv"))
    detect.add_argument("--out", type=Path, required=True)
    detect.add_argument("--report", type=Path, default=None)
    detect.add_argument("--cache-dir", type=Path, default=Path("structures/alphafold"))
    detect.add_argument("--fpocket-timeout", type=float, default=900.0)
    detect.add_argument("--limit", type=int, default=0)

    analyse = sub.add_parser("analyse", help="candidate against its own protein's pockets")
    analyse.add_argument("--pockets", type=Path, required=True)
    analyse.add_argument("--table", type=Path,
                         default=Path("results/specificity/candidates.csv"))
    analyse.add_argument("--out", type=Path, required=True)
    analyse.add_argument("--markdown", type=Path, default=None)
    analyse.add_argument("--arm", choices=sorted(ARMS), default="candidate")
    analyse.add_argument("--draws", type=int, default=800)
    analyse.add_argument("--bootstrap", type=int, default=2000)

    sites = sub.add_parser("sites", help="annotated sites and their structures' pockets")
    sites.add_argument("--census", type=Path, default=Path("results/redocking/census.csv"))
    sites.add_argument("--structures", type=Path, required=True)
    sites.add_argument("--comp-id", default="IHP")
    sites.add_argument("--out", type=Path, required=True)
    sites.add_argument("--sites", type=Path, required=True)
    sites.add_argument("--report", type=Path, default=None)
    sites.add_argument("--summary", type=Path, default=None)
    sites.add_argument("--fpocket-timeout", type=float, default=900.0)
    sites.add_argument("--limit", type=int, default=0)

    args = parser.parse_args(argv)
    commands = {"pockets": cmd_pockets, "sites": cmd_sites, "analyse": cmd_analyse}
    return commands[args.command](args)


if __name__ == "__main__":
    raise SystemExit(main())
