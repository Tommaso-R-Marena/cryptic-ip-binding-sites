#!/usr/bin/env python3
"""Study I (docs/ORTHOGONAL_PLAN.md): sequence-only evidence on the screen's candidates.

Two tools that measure something neither the screen nor the docking sees:

* **metapredict** v3 - per-residue disorder from sequence alone. A pocket predicted in a
  disordered region is an artefact, and AlphaFold pLDDT is only a proxy for that.
* **SHARK** (``bio-shark``) - alignment-free conservation. ``shark-capture`` gives a
  second opinion on study H's MAFFT criterion; ``shark-dive`` asks whether the candidates
  study H could not evaluate have remote relatives at all.

    python scripts/orthogonal.py disorder --candidates c.csv --proteins p.csv.gz --out disorder.json
    python scripts/orthogonal.py motifs   --candidates c.csv --proteins p.csv.gz --work w --out motifs.json
    python scripts/orthogonal.py targets  --proteins p.csv.gz --out targets.fasta
    python scripts/orthogonal.py dive     --conservation cons.json --targets targets.fasta --out dive.json
    python scripts/orthogonal.py report   --candidates c.csv --disorder d.json ... --out-dir results/orthogonal

The three arms - candidates, their pLDDT- and depth-matched controls, and the annotated
IP binders - are built by ``scripts/triage.py``'s own matching, imported rather than
reimplemented, so study I's arms are study H's arms.
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for _p in (ROOT, ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

LOGGER = logging.getLogger("orthogonal")

SEED = 20261002
N_BOOTSTRAP = 2000
DISORDER_CUT = 0.5          # metapredict's own disordered/ordered boundary
NO_DIFFERENCE = 0.05        # I1's "no difference" margin
K_MIN, K_MAX = 3, 8         # shark-capture motif lengths, fixed in the plan
MAX_CAPTURE = 25            # homologues given to shark-capture (runtime cap, fixed in the plan)
CAPTURE_TIMEOUT = 2400      # seconds; capture is superlinear, and the step is best-effort
TARGET_SAMPLE = 3000        # shark-dive target sequences, seeded sample of the screened proteomes
DIVE_TOP = 100              # remote homologues retained per query
BASIC = ("K", "R", "H")


# --------------------------------------------------------------------------- arms
def roles_frame(candidates: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    """Study H's three arms, through study H's matching."""
    from triage import matched_controls

    cands = candidates.copy()
    cands["role"] = "candidate"
    cands["matched_to"] = None
    positives = proteins[proteins.get("annotated", False).astype(bool)
                         & ~proteins.get("seen", False).astype(bool)].copy()
    positives["role"] = "positive"
    positives["matched_to"] = None
    controls = matched_controls(cands, proteins)
    if controls.empty:
        raise RuntimeError("no matched control could be built, so I1 would have no comparison arm; "
                           "refusing to run a study whose only hypothesis test is empty")
    # A protein can be both an annotated binder and a matched control. Study H keyed its
    # results by accession, so those proteins ended up counted as controls only; keeping the
    # arms commensurable means resolving the clash the same way, and saying how often it fell.
    clash = set(controls["uniprot_id"]) & set(positives["uniprot_id"])
    positives = positives[~positives["uniprot_id"].isin(clash)]
    if clash:
        LOGGER.info("%d annotated binders were also drawn as controls and count as controls only: %s",
                    len(clash), sorted(clash))
    cols = ["uniprot_id", "organism_key", "top_pocket_residues", "cluster", "role", "matched_to"]
    return pd.concat([cands[cols], positives[cols], controls[cols]], ignore_index=True)


def sequences_for(accessions: Sequence[str]) -> Dict[str, str]:
    from learned_screen import fetch_uniprot_fasta

    return fetch_uniprot_fasta(sorted({str(a) for a in accessions}))


def basic_pocket_positions(sequence: str, positions: Sequence[int]) -> List[int]:
    return [int(p) for p in positions if 0 < int(p) <= len(sequence) and sequence[int(p) - 1] in BASIC]


# ----------------------------------------------------------------------- disorder
def pocket_disorder(sequence: str, positions: Sequence[int]) -> Dict[str, object]:
    """Mean metapredict v3 disorder over the pocket's residues, and over the whole chain."""
    import metapredict

    scores = list(metapredict.predict_disorder(sequence))
    inside = [scores[int(p) - 1] for p in positions if 0 < int(p) <= len(scores)]
    if not inside:
        return {"error": "no pocket residue falls inside the sequence"}
    mean = float(sum(inside) / len(inside))
    return {"pocket_disorder": mean, "disordered": bool(mean > DISORDER_CUT),
            "residues_scored": len(inside), "length": len(scores),
            "protein_disorder": float(sum(scores) / len(scores))}


def cmd_disorder(args: argparse.Namespace) -> int:
    from triage import parse_positions

    frame = roles_frame(pd.read_csv(args.candidates), pd.read_csv(args.proteins))
    seqs = sequences_for(frame["uniprot_id"])
    LOGGER.info("%d proteins, %d sequences", len(frame), len(seqs))
    out: Dict[str, object] = {}
    for _, row in frame.iterrows():
        acc = str(row["uniprot_id"])
        entry: Dict[str, object] = {"role": row["role"], "cluster": row.get("cluster"),
                                    "matched_to": row.get("matched_to"),
                                    "organism_key": row.get("organism_key")}
        positions = parse_positions(row.get("top_pocket_residues"))
        if acc not in seqs:
            entry["error"] = "sequence unavailable"
        elif not positions:
            entry["error"] = "no pocket residues"
        else:
            entry.update(pocket_disorder(seqs[acc], positions))
        out[acc] = entry
        LOGGER.info("%s (%s): %s", acc, row["role"], entry.get("pocket_disorder", entry.get("error")))
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({"proteins": len(out)}))
    return 0


# -------------------------------------------------------------------------- SHARK
def write_fasta(records: Dict[str, str], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f">{n}\n{s}\n" for n, s in records.items()))
    return path


def run_capture(fasta: Path, out_dir: Path, timeout: int = CAPTURE_TIMEOUT) -> Path:
    """shark-capture on a homologue set; returns the occurrences table."""
    subprocess.run(["shark-capture", str(fasta), str(out_dir), "--k_min", str(K_MIN),
                    "--k_max", str(K_MAX), "--no_per_sequence_kmer_plots"],
                   check=True, capture_output=True, text=True, timeout=timeout)
    return out_dir / "outputs" / "occurrences" / "all_occurrences.tsv"


def motif_support(occurrences: Path, acc: str, positions: Sequence[int]) -> Dict[str, object]:
    """Fraction of the pocket's basic positions inside a captured motif's match in the target."""
    table = pd.read_csv(occurrences, sep="\t")
    mine = table[table["sequenceID"].astype(str) == acc]
    spans = [(int(a), int(b)) for a, b in zip(mine["start"], mine["end"])]
    covered = [p for p in positions if any(a <= p <= b for a, b in spans)]
    return {"motif_support": (len(covered) / len(positions)) if positions else None,
            "covered_positions": covered, "basic_positions": list(positions),
            "motifs": int(mine["reference_kmer"].nunique()) if len(mine) else 0}


def capture_one(acc: str, positions: Sequence[int], work: Path) -> Dict[str, object]:
    from arrestin import uniref50_members
    from learned_screen import fetch_uniprot_fasta

    cluster, members = uniref50_members(acc)
    if acc not in members:
        members.update(fetch_uniprot_fasta([acc]))
    if acc not in members:
        return {"error": "target sequence unavailable", "uniref50": cluster}
    basic = basic_pocket_positions(members[acc], positions)
    if not basic:
        return {"error": "no basic pocket position", "uniref50": cluster}
    names = [acc] + [n for n in sorted(members) if n != acc][:MAX_CAPTURE]
    if len(names) < 2:
        return {"error": "no homologue to compare against", "uniref50": cluster}
    fasta = write_fasta({n: members[n] for n in names}, work / acc / "homologues.fasta")
    occurrences = run_capture(fasta, work / acc / "capture")
    out = motif_support(occurrences, acc, basic)
    out.update({"uniref50": cluster, "homologues": len(names) - 1})
    return out


def cmd_motifs(args: argparse.Namespace) -> int:
    from triage import parse_positions

    frame = roles_frame(pd.read_csv(args.candidates), pd.read_csv(args.proteins))
    frame = frame[frame["role"].isin(("candidate", "control"))].sort_values("uniprot_id")
    if args.shards > 1:
        frame = frame.iloc[args.shard::args.shards]
        LOGGER.info("shard %d of %d: %d proteins", args.shard, args.shards, len(frame))
    out: Dict[str, object] = {}
    for _, row in frame.iterrows():
        acc = str(row["uniprot_id"])
        entry: Dict[str, object] = {"role": row["role"], "cluster": row.get("cluster"),
                                    "matched_to": row.get("matched_to")}
        try:
            entry.update(capture_one(acc, parse_positions(row.get("top_pocket_residues")), args.work))
        except Exception as exc:  # noqa: BLE001 - SHARK is best-effort, per the plan
            entry["error"] = f"{type(exc).__name__}: {exc}"[:200]
        out[acc] = entry
        LOGGER.info("%s (%s): %s", acc, row["role"], entry.get("motif_support", entry.get("error")))
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({"proteins": len(out)}))
    return 0


def cmd_targets(args: argparse.Namespace) -> int:
    """A seeded sample of the screened proteomes: shark-dive's target database."""
    proteins = pd.read_csv(args.proteins)
    accessions = sorted({str(a) for a in proteins["uniprot_id"].dropna()})
    rng = random.Random(SEED)
    if len(accessions) > args.n:
        accessions = sorted(rng.sample(accessions, args.n))
    records = sequences_for(accessions)
    write_fasta(records, args.out)
    print(json.dumps({"sampled": len(accessions), "sequences": len(records)}))
    return 0


def unevaluable_candidates(conservation: Dict[str, dict]) -> List[str]:
    return sorted(acc for acc, v in conservation.items()
                  if v.get("role") == "candidate" and str(v.get("decision", "")).startswith("not evaluable"))


def run_dive(query: Path, targets: Path, out_dir: Path, timeout: int = 7200) -> Optional[Path]:
    subprocess.run(["shark-dive", str(query), str(targets), "--output_dir", str(out_dir)],
                   check=True, capture_output=True, text=True, timeout=timeout)
    hits = sorted(out_dir.rglob("*.csv")) + sorted(out_dir.rglob("*.tsv"))
    return hits[0] if hits else None


def cmd_dive(args: argparse.Namespace) -> int:
    conservation = json.loads(args.conservation.read_text())
    accessions = unevaluable_candidates(conservation)
    seqs = sequences_for(accessions)
    out: Dict[str, object] = {}
    for acc in accessions:
        entry: Dict[str, object] = {"reason": conservation[acc].get("decision")}
        try:
            if acc not in seqs:
                raise RuntimeError("sequence unavailable")
            query = write_fasta({acc: seqs[acc]}, args.work / acc / "query.fasta")
            table = run_dive(query, args.targets, args.work / acc / "dive")
            if table is None:
                raise RuntimeError("shark-dive produced no table")
            frame = pd.read_csv(table, sep="\t" if table.suffix == ".tsv" else ",")
            entry["remote_homologues"] = int(min(len(frame), DIVE_TOP))
        except Exception as exc:  # noqa: BLE001 - best-effort, per the plan
            entry["error"] = f"{type(exc).__name__}: {exc}"[:200]
        out[acc] = entry
        LOGGER.info("%s: %s", acc, entry.get("remote_homologues", entry.get("error")))
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({"queries": len(out)}))
    return 0


def cmd_merge(args: argparse.Namespace) -> int:
    """One JSON object from the shards of a matrix job."""
    merged: Dict[str, object] = {}
    for path in sorted(args.inputs):
        for chunk in sorted(path.rglob("*.json")) if path.is_dir() else [path]:
            merged.update(json.loads(chunk.read_text()))
    args.out.write_text(json.dumps(merged, indent=2, default=str))
    print(json.dumps({"proteins": len(merged)}))
    return 0


# ------------------------------------------------------------------------- report
def _paired(records: Dict[str, dict], field: str, n_bootstrap: int) -> Tuple[List[dict], Optional[dict]]:
    """Candidate minus its matched control on ``field``, resampled as a pair."""
    from cryptic_ip.docking import stats

    pairs = []
    for acc, v in records.items():
        if v.get("role") != "control" or v.get("matched_to") is None:
            continue
        cand = records.get(str(v["matched_to"]))
        if cand is None or cand.get(field) is None or v.get(field) is None:
            continue
        pairs.append({"cluster": str(cand.get("cluster") or v["matched_to"]),
                      "cand": float(cand[field]), "ctrl": float(v[field])})
    if not pairs:
        return pairs, None
    frame = pd.DataFrame(pairs)
    est = stats.paired_difference(frame["cand"].to_numpy(), frame["ctrl"].to_numpy(),
                                  frame["cluster"].to_numpy(), n_bootstrap=n_bootstrap, seed=SEED)
    est["clusters"] = int(frame["cluster"].nunique())
    return pairs, est


def _rate(records: Dict[str, dict], role: str, field: str, n_bootstrap: int) -> Dict[str, object]:
    from cryptic_ip.docking import stats

    part = [v for v in records.values() if v.get("role") == role and v.get(field) is not None]
    out: Dict[str, object] = {"n": sum(1 for v in records.values() if v.get("role") == role),
                              "scored": len(part)}
    if part:
        out["estimate"] = stats.mean_estimates([float(v[field]) for v in part],
                                               [str(v.get("cluster") or "") for v in part],
                                               null=0.0, n_bootstrap=n_bootstrap, seed=SEED)
    return out


def decide_i1(est: Optional[dict], clusters: int) -> str:
    from cryptic_ip.docking import stats

    if est is None:
        return "not evaluable: no matched pair with both sides scored"
    if clusters < stats.MIN_GROUPS:
        return f"not evaluable: {clusters} clusters (fewer than {stats.MIN_GROUPS})"
    per_group = est["per_group"]
    if per_group["high"] < 0:
        return "candidates more ordered"
    if per_group["low"] >= -NO_DIFFERENCE and per_group["high"] <= NO_DIFFERENCE:
        return "no difference"
    return "inconclusive"


def build(candidates: pd.DataFrame, disorder: Dict[str, dict], motifs: Dict[str, dict],
          dive: Dict[str, dict], n_bootstrap: int = N_BOOTSTRAP) -> Tuple[Dict[str, object], pd.DataFrame]:
    from cryptic_ip.benchmark import protocol

    report: Dict[str, object] = {"plan": "docs/ORTHOGONAL_PLAN.md"}

    # I1: does sequence disorder separate candidates from pLDDT-matched controls?
    pairs, est = _paired(disorder, "pocket_disorder", n_bootstrap)
    clusters = int(est["clusters"]) if est else 0
    report["I1"] = {"question": "mean pocket disorder, candidates - matched controls",
                    "pairs": len(pairs), "clusters": clusters,
                    "estimate": (est or {}).get("per_group"), "per_copy": (est or {}).get("per_copy"),
                    "decision": decide_i1(est, clusters)}
    p_values = {"I1": est["per_group"]["p_value"]} if est else {}
    report["holm"] = protocol.holm(p_values) if p_values else {}

    # I2: the disorder QC, and whether it would also reject known binders.
    for entry in disorder.values():
        if entry.get("pocket_disorder") is not None:
            entry["disordered_flag"] = float(bool(entry["pocket_disorder"] > DISORDER_CUT))
    rates = {role: _rate(disorder, role, "disordered_flag", n_bootstrap)
             for role in ("candidate", "control", "positive")}
    pos = ((rates["positive"].get("estimate") or {}).get("per_group") or {}).get("point")
    informative = bool(pos is not None and pos <= 0.25)
    demoted = sorted(acc for acc, v in disorder.items()
                     if v.get("role") == "candidate" and v.get("disordered_flag") == 1.0)
    report["I2"] = {"rates": rates, "informative": informative,
                    "demoted": demoted if informative else [],
                    "note": ("annotated binders are rarely flagged, so a flagged candidate is demoted"
                             if informative else
                             "annotated binders are flagged at a comparable rate: the disorder filter is "
                             "not informative here, and no candidate is demoted on it")}

    # I3: alignment-free motif support, descriptive.
    mpairs, mest = _paired(motifs, "motif_support", n_bootstrap)
    report["I3"] = {"question": "shark-capture motif support, candidates - matched controls",
                    "pairs": len(mpairs), "clusters": int(mest["clusters"]) if mest else 0,
                    "estimate": (mest or {}).get("per_group"),
                    "failures": sum(1 for v in motifs.values() if v.get("error")),
                    "note": "descriptive: orthogonal to study H's MAFFT criterion, not combined with it"}

    # I4: remote homology for the candidates study H could not evaluate.
    found = {acc: v.get("remote_homologues") for acc, v in dive.items()}
    report["I4"] = {"queries": len(dive),
                    "with_remote_homologues": sum(1 for n in found.values() if n),
                    "counts": found,
                    "note": "a within-proteome remote hit is not an orthologue set: study H's "
                            "conservation criterion is not re-run on these"}

    merged = candidates.copy()
    merged["pocket_disorder"] = merged["uniprot_id"].map(
        {a: v.get("pocket_disorder") for a, v in disorder.items()})
    merged["disordered"] = merged["uniprot_id"].map(
        {a: v.get("disordered") for a, v in disorder.items()})
    merged["motif_support"] = merged["uniprot_id"].map(
        {a: v.get("motif_support") for a, v in motifs.items()})
    merged["remote_homologues"] = merged["uniprot_id"].map(found)
    report["notes"] = [
        "The controls are already pLDDT-matched, so I1 asks whether sequence disorder carries "
        "information beyond AlphaFold's confidence. A null is a real answer.",
        "No sequence tool speaks to ligand identity; study C measured that and found the "
        "descriptors too weak to change a ranking.",
    ]
    return report, merged


def markdown(r: Dict[str, object]) -> str:
    def fmt(e):
        return "-" if not e or "point" not in e else f"{e['point']:.3f} [{e['low']:.3f}, {e['high']:.3f}]"

    i1, i2, i3, i4 = (r.get(k, {}) for k in ("I1", "I2", "I3", "I4"))
    lines = ["## Orthogonal sequence evidence on the candidates (docs/ORTHOGONAL_PLAN.md)", "",
             f"**I1:** {i1.get('decision')} - {fmt(i1.get('estimate'))} over {i1.get('pairs')} matched "
             f"pairs in {i1.get('clusters')} clusters.", "",
             "| role | n | scored | disordered-pocket rate |", "|---|---|---|---|"]
    for role, v in (i2.get("rates") or {}).items():
        lines.append(f"| {role} | {v.get('n')} | {v.get('scored')} | "
                     f"{fmt((v.get('estimate') or {}).get('per_group'))} |")
    lines += ["", f"**I2 (disorder QC):** informative = {i2.get('informative')}. {i2.get('note')} "
              f"Demoted: {', '.join(i2.get('demoted') or []) or 'none'}.", "",
              f"**I3 (motif support):** {fmt(i3.get('estimate'))} over {i3.get('pairs')} pairs in "
              f"{i3.get('clusters')} clusters; {i3.get('failures')} proteins gave no result. "
              f"{i3.get('note')}", "",
              f"**I4 (remote homology):** {i4.get('with_remote_homologues')} of {i4.get('queries')} "
              f"unevaluable candidates have a remote homologue in the sampled proteomes. "
              f"{i4.get('note')}", "",
              *[f"- {n}" for n in r.get("notes", [])]]
    return "\n".join(lines) + "\n"


def cmd_report(args: argparse.Namespace) -> int:
    report, merged = build(pd.read_csv(args.candidates), json.loads(args.disorder.read_text()),
                           json.loads(args.motifs.read_text()) if args.motifs.exists() else {},
                           json.loads(args.dive.read_text()) if args.dive.exists() else {},
                           args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "orthogonal.json").write_text(json.dumps(report, indent=2, default=str))
    merged.to_csv(args.out_dir / "orthogonal_candidates.csv", index=False)
    text = markdown(report)
    (args.out_dir / "ORTHOGONAL.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    d = sub.add_parser("disorder")
    d.add_argument("--candidates", type=Path, required=True)
    d.add_argument("--proteins", type=Path, required=True)
    d.add_argument("--out", type=Path, required=True)

    m = sub.add_parser("motifs")
    m.add_argument("--candidates", type=Path, required=True)
    m.add_argument("--proteins", type=Path, required=True)
    m.add_argument("--work", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    m.add_argument("--shard", type=int, default=0)
    m.add_argument("--shards", type=int, default=1)

    t = sub.add_parser("targets")
    t.add_argument("--proteins", type=Path, required=True)
    t.add_argument("--n", type=int, default=TARGET_SAMPLE)
    t.add_argument("--out", type=Path, required=True)

    v = sub.add_parser("dive")
    v.add_argument("--conservation", type=Path, required=True)
    v.add_argument("--targets", type=Path, required=True)
    v.add_argument("--work", type=Path, required=True)
    v.add_argument("--out", type=Path, required=True)

    g = sub.add_parser("merge")
    g.add_argument("--inputs", type=Path, nargs="+", required=True)
    g.add_argument("--out", type=Path, required=True)

    r = sub.add_parser("report")
    r.add_argument("--candidates", type=Path, required=True)
    r.add_argument("--disorder", type=Path, required=True)
    r.add_argument("--motifs", type=Path, required=True)
    r.add_argument("--dive", type=Path, required=True)
    r.add_argument("--out-dir", type=Path, required=True)
    r.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)

    args = parser.parse_args(argv)
    return {"disorder": cmd_disorder, "motifs": cmd_motifs, "targets": cmd_targets,
            "dive": cmd_dive, "merge": cmd_merge, "report": cmd_report}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
