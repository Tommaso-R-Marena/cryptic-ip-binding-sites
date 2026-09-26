#!/usr/bin/env python3
"""Study H (docs/TRIAGE_PLAN.md): is a candidate's top pocket conserved-basic across its orthologues?

    python scripts/triage.py conserve --candidates candidates.csv --proteins proteins.csv.gz \
        --out-dir work --out conservation.json
    python scripts/triage.py report --candidates candidates.csv --conservation conservation.json \
        --out-dir results/triage

Conservation reuses ``scripts/arrestin.py``: the UniRef50 fetch, the MAFFT alignment and
the criterion (at least 3 basic pocket positions with a basic fraction >= 0.80) are
imported, not reimplemented, so study B's filter is applied unchanged. That script is
never edited here, because editing a finished study's script would re-run its workflow.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

LOGGER = logging.getLogger("triage")
SEED = 20261001
N_BOOTSTRAP = 2000
PLDDT_TOLERANCE = 5.0
DEPTH_TOLERANCE = 3.0
MAX_HOMOLOGUES = 500


def parse_positions(text: object) -> List[int]:
    if not isinstance(text, str) or not text.strip():
        return []
    return [int(t) for t in text.replace(";", ",").split(",") if t.strip().lstrip("-").isdigit()]


def matched_controls(candidates: pd.DataFrame, proteins: pd.DataFrame) -> pd.DataFrame:
    """One pLDDT- and depth-matched, low-ranking protein per candidate (the plan's rule)."""
    required = ("uniprot_id", "organism_key", "top_pocket_residues", "plddt_mean", "hull_depth", "combined")
    missing = [c for c in required if c not in proteins.columns]
    if missing:
        raise ValueError(
            f"the proteins table lacks {missing}; the matching needs study C's proteins_combined.csv.gz "
            f"(the specificity-rerank artifact), not the learned screen's proteins.csv.gz. "
            f"Columns present: {sorted(proteins.columns)[:20]}")
    pool = proteins[~proteins["uniprot_id"].isin(set(candidates["uniprot_id"]))].copy()
    for col in ("plddt_mean", "hull_depth", "combined"):
        pool[col] = pd.to_numeric(pool[col], errors="coerce")
    pool = pool.dropna(subset=["plddt_mean", "hull_depth", "combined"])
    rows, used = [], set()
    for _, cand in candidates.sort_values("uniprot_id").iterrows():
        same = pool[pool["organism_key"] == cand["organism_key"]]
        if same.empty:
            continue
        cutoff = same["combined"].quantile(0.5)
        near = same[(same["combined"] <= cutoff)
                    & (same["plddt_mean"] - cand["plddt_mean"]).abs().le(PLDDT_TOLERANCE)
                    & (same["hull_depth"] - cand["hull_depth"]).abs().le(DEPTH_TOLERANCE)
                    & ~same["uniprot_id"].isin(used)]
        if near.empty:
            continue
        pick = near.sort_values("uniprot_id").iloc[0]
        used.add(pick["uniprot_id"])
        rows.append({"uniprot_id": pick["uniprot_id"], "organism_key": pick["organism_key"],
                     "top_pocket_residues": pick.get("top_pocket_residues"), "cluster": pick.get("cluster"),
                     "role": "control", "matched_to": cand["uniprot_id"]})
    return pd.DataFrame(rows)


def conserve_one(acc: str, positions: Sequence[int], work: Path) -> Dict[str, object]:
    """UniRef50 homologues, MAFFT, then study B's criterion. Reused verbatim from arrestin.py."""
    from arrestin import conservation, conservation_decision, read_alignment, uniref50_members
    from learned_screen import fetch_uniprot_fasta

    if not positions:
        return {"decision": "not evaluable: no pocket residues"}
    cluster, seqs = uniref50_members(acc)
    unique: Dict[str, str] = {}
    for name in sorted(seqs):
        if seqs[name] not in unique.values() or name == acc:
            unique[name] = seqs[name]
    if acc not in unique:
        unique.update(fetch_uniprot_fasta([acc]))
    if acc not in unique:
        return {"decision": "not evaluable: target sequence unavailable", "uniref50": cluster}
    names = [acc] + [n for n in sorted(unique) if n != acc][:MAX_HOMOLOGUES]
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / f"{acc}.fasta"
    fasta.write_text("".join(f">{n}\n{unique[n]}\n" for n in names))
    if len(names) > 1:
        proc = subprocess.run(["mafft", "--auto", "--anysymbol", "--quiet", str(fasta)],
                              capture_output=True, text=True, timeout=3600)
        aligned = work / f"{acc}_aligned.fasta"
        aligned.write_text(proc.stdout)
        alignment = read_alignment(aligned)
    else:
        alignment = {acc: unique[acc]}
    scores = conservation(alignment, acc, positions)
    out = dict(conservation_decision(unique[acc], positions, scores, len(names) - 1))
    out["uniref50"] = cluster
    out["pocket_positions"] = len(positions)
    return out


def cmd_conserve(args: argparse.Namespace) -> int:
    cands = pd.read_csv(args.candidates)
    cands["role"] = "candidate"
    proteins = pd.read_csv(args.proteins)
    positives = proteins[proteins.get("annotated", False).astype(bool)
                         & ~proteins.get("seen", False).astype(bool)].copy()
    positives["role"] = "positive"
    controls = matched_controls(cands, proteins)
    if controls.empty:
        raise RuntimeError("no matched control could be built, so H1 would have no comparison arm; "
                           "refusing to run a triage whose only hypothesis test is empty")
    frame = pd.concat([cands[["uniprot_id", "organism_key", "top_pocket_residues", "cluster", "role"]],
                       positives[["uniprot_id", "organism_key", "top_pocket_residues", "cluster", "role"]],
                       controls], ignore_index=True)
    LOGGER.info("%d proteins: %s", len(frame), frame["role"].value_counts().to_dict())
    out: Dict[str, object] = {}
    for _, row in frame.iterrows():
        acc = str(row["uniprot_id"])
        entry: Dict[str, object] = {"role": row["role"], "organism_key": row.get("organism_key"),
                                    "cluster": row.get("cluster"), "matched_to": row.get("matched_to")}
        try:
            entry.update(conserve_one(acc, parse_positions(row.get("top_pocket_residues")), args.out_dir))
        except Exception as exc:  # noqa: BLE001 - recorded as this protein's reason
            entry["decision"] = f"not evaluable: {type(exc).__name__}: {exc}"[:200]
        out[acc] = entry
        LOGGER.info("%s (%s): %s", acc, row["role"], entry.get("decision"))
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps({"proteins": len(out)}))
    return 0


def rate(frame: pd.DataFrame) -> Dict[str, object]:
    from cryptic_ip.docking import stats

    evaluable = frame[~frame["decision"].astype(str).str.startswith("not evaluable")]
    y = (evaluable["decision"] == "conserved").to_numpy(dtype=float)
    groups = evaluable["cluster"].fillna(evaluable["uniprot_id"]).astype(str).to_numpy()
    out: Dict[str, object] = {"n": int(len(frame)), "evaluable": int(len(evaluable)),
                              "not_evaluable": int(len(frame) - len(evaluable))}
    if len(evaluable):
        out["estimate"] = stats.mean_estimates(y, groups, null=0.0, n_bootstrap=N_BOOTSTRAP, seed=SEED)
    return out


def build(cands: pd.DataFrame, conservation: Dict[str, dict], n_bootstrap: int) -> Dict[str, object]:
    from cryptic_ip.benchmark import protocol
    from cryptic_ip.docking import stats

    frame = pd.DataFrame([{"uniprot_id": k, **v} for k, v in conservation.items()])
    frame["decision"] = frame["decision"].astype(str)
    report: Dict[str, object] = {"plan": "docs/TRIAGE_PLAN.md",
                                 "roles": frame["role"].value_counts().to_dict()}
    report["rates"] = {role: rate(part) for role, part in frame.groupby("role")}

    # H1: paired candidate - matched control, resampling the pair together.
    pairs = []
    by_acc = frame.set_index("uniprot_id")
    for _, ctrl in frame[frame["role"] == "control"].iterrows():
        cand = by_acc.loc[ctrl["matched_to"]] if ctrl["matched_to"] in by_acc.index else None
        if cand is None or str(cand["decision"]).startswith("not evaluable") \
                or str(ctrl["decision"]).startswith("not evaluable"):
            continue
        pairs.append({"cluster": str(cand.get("cluster") or ctrl["matched_to"]),
                      "cand": float(cand["decision"] == "conserved"),
                      "ctrl": float(ctrl["decision"] == "conserved")})
    p_values: Dict[str, float] = {}
    if pairs:
        pf = pd.DataFrame(pairs)
        est = stats.paired_difference(pf["cand"].to_numpy(), pf["ctrl"].to_numpy(),
                                      pf["cluster"].to_numpy(), n_bootstrap=n_bootstrap, seed=SEED)
        n_groups = int(pf["cluster"].nunique())
        decision = (f"not evaluable: {n_groups} clusters (fewer than {stats.MIN_GROUPS})"
                    if n_groups < stats.MIN_GROUPS else
                    "enriched" if est["per_group"]["low"] > 0 else
                    "not enriched" if est["per_group"]["high"] < 0.05 else "inconclusive")
        report["H1"] = {"question": "conserved-basic rate, candidates - matched controls",
                        "pairs": len(pairs), "clusters": n_groups,
                        "estimate": est.get("per_group"), "per_copy": est.get("per_copy"),
                        "decision": decision}
        p_values["H1"] = est["per_group"]["p_value"]
    else:
        report["H1"] = {"decision": "not evaluable: no matched pair with both sides evaluable"}
    report["holm"] = protocol.holm(p_values) if p_values else {}

    # H2: does the filter keep known binders?
    pos = report["rates"].get("positive", {})
    pos_point = (pos.get("estimate") or {}).get("per_group", {}).get("point")
    report["H2"] = {"positive_controls": pos.get("n"), "evaluable": pos.get("evaluable"),
                    "conserved_rate": (pos.get("estimate") or {}).get("per_group"),
                    "informative": bool(pos_point is not None and pos_point >= 0.5),
                    "note": ("the filter keeps most annotated binders, so it can triage"
                             if (pos_point or 0) >= 0.5 else
                             "the filter rejects most annotated binders: it is not informative for triage, "
                             "and no candidate is promoted or demoted on it")}

    # H3: the triaged list.
    merged = cands.merge(frame[["uniprot_id", "decision", "homologues", "basic_positions",
                                "conserved_positions", "basic_fraction"]], on="uniprot_id", how="left")

    def label(row):
        if bool(row.get("explained")):
            return "explained"
        d = str(row.get("decision"))
        return "conserved basic pocket" if d == "conserved" else d
    merged["triage"] = merged.apply(label, axis=1)
    report["H3_counts"] = merged["triage"].value_counts().to_dict()
    keep = merged[merged["triage"] == "conserved basic pocket"].sort_values(["organism_key", "rank"])
    wanted = ["uniprot_id", "gene", "protein_name", "organism_key", "rank", "combined", "p_ip_given_site",
              "hull_depth", "plddt_mean", "homologues", "basic_positions", "conserved_positions"]
    report["H3_shortlist"] = [
        {k: (None if not isinstance(v, (list, dict)) and pd.isna(v) else v) for k, v in row.items()}
        for row in keep[[c for c in wanted if c in keep]].to_dict(orient="records")]
    report["notes"] = [
        "A conserved basic pocket is not evidence of binding: study C showed these descriptors do not "
        "separate IP from other polyanions well enough to change a ranking.",
        "The candidate list was read before this plan was written, so this is a filter, not a test of "
        "the screen.",
    ]
    return report, merged


def markdown(r: Dict[str, object]) -> str:
    def fmt(e):
        return "-" if not e or "point" not in e else f"{e['point']:.3f} [{e['low']:.3f}, {e['high']:.3f}]"
    h1 = r.get("H1", {})
    lines = ["## Triaging the candidates by pocket conservation (docs/TRIAGE_PLAN.md)", "",
             f"Roles: {json.dumps(r.get('roles'))}", "",
             f"**H1:** {h1.get('decision')} - {fmt(h1.get('estimate'))} over {h1.get('pairs')} matched pairs "
             f"in {h1.get('clusters')} clusters.", "",
             "| role | n | evaluable | conserved-basic rate |", "|---|---|---|---|"]
    for role, v in (r.get("rates") or {}).items():
        lines.append(f"| {role} | {v.get('n')} | {v.get('evaluable')} | "
                     f"{fmt((v.get('estimate') or {}).get('per_group'))} |")
    h2 = r.get("H2", {})
    lines += ["", f"**H2 (filter calibration):** annotated binders conserved-basic "
              f"{fmt(h2.get('conserved_rate'))}; informative = {h2.get('informative')}. {h2.get('note')}", "",
              f"**H3:** {json.dumps(r.get('H3_counts'))}", ""]
    short = r.get("H3_shortlist") or []
    if short:
        lines += ["| organism | rank | protein | homologues | basic positions | conserved |",
                  "|---|---|---|---|---|---|"]
        for c in short:
            name = c.get("gene") or c.get("uniprot_id")
            lines.append(f"| {c.get('organism_key')} | {c.get('rank')} | {name} ({c.get('uniprot_id')}) | "
                         f"{c.get('homologues')} | {len(c.get('basic_positions') or [])} | "
                         f"{len(c.get('conserved_positions') or [])} |")
    lines += ["", *[f"- {n}" for n in r.get("notes", [])]]
    return "\n".join(lines) + "\n"


def cmd_report(args: argparse.Namespace) -> int:
    cands = pd.read_csv(args.candidates)
    conservation = json.loads(args.conservation.read_text())
    report, merged = build(cands, conservation, args.n_bootstrap)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "triage.json").write_text(json.dumps(report, indent=2, default=str))
    merged.to_csv(args.out_dir / "triaged_candidates.csv", index=False)
    text = markdown(report)
    (args.out_dir / "TRIAGE.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    c = sub.add_parser("conserve")
    c.add_argument("--candidates", type=Path, required=True)
    c.add_argument("--proteins", type=Path, required=True)
    c.add_argument("--out-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    r = sub.add_parser("report")
    r.add_argument("--candidates", type=Path, required=True)
    r.add_argument("--conservation", type=Path, required=True)
    r.add_argument("--out-dir", type=Path, required=True)
    r.add_argument("--n-bootstrap", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args(argv)
    return cmd_conserve(args) if args.command == "conserve" else cmd_report(args)


if __name__ == "__main__":
    sys.exit(main())
