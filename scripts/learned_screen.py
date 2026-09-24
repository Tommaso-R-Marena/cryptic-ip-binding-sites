#!/usr/bin/env python3
"""Rank screened proteomes with the benchmark's learned model (docs/LEARNED_SCREEN_PLAN.md).

``train``
    Lock the ``ip_site`` model on the whole benchmark table.
``fasta``
    Fetch UniProt sequences for a list of accessions (accession-only headers).
``evaluate``
    Score every screened pocket; protein score is the best confident pocket;
    test L1/L2 against UniProt inositol phosphate annotations on proteins with
    no homologue among the benchmark's proteins; list candidates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.analysis.proteome_stats import known_ip_annotation  # noqa: E402
from cryptic_ip.benchmark import protocol  # noqa: E402
from cryptic_ip.benchmark.homology import MAX_EVALUE, MIN_IDENTITY, MIN_SHORTER_COVERAGE  # noqa: E402

MIN_PLDDT = 70.0
N_CANDIDATES = 25
MIN_BINDERS_FOR_DECISION = 5
SEED = 20260924
SEARCH_COLUMNS = ("query", "target", "fident", "qstart", "qend", "tstart", "tend", "qlen", "tlen", "evalue")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ------------------------------------------------------------------- model
def train_locked(table: pd.DataFrame, *, n_draws: int = 10, seed: int = 0) -> Dict[str, object]:
    rows = table[table["label_ip_site"] >= 0].reset_index(drop=True)
    features = list(protocol.ARMS["full"])
    X = rows[features].astype(float)
    y = rows["label_ip_site"].to_numpy(dtype=int)
    groups = rows["group_strict"].astype(str).to_numpy()
    specs = protocol.benchmark_specs()
    candidates = protocol.draw_candidates(n_draws, 7, specs)
    selection = protocol.select_candidate(X, y, groups, candidates, n_inner=3, seed=seed * 1000 + 999, specs=specs)
    model = protocol.fit_candidate(selection.candidate, X, y, seed, specs)
    return {"model": model, "selection": selection, "features": features,
            "candidate": selection.candidate.as_dict(), "rows": int(len(rows)), "positives": int(y.sum())}


def score_pockets(bundle: Mapping[str, object], pockets: pd.DataFrame) -> np.ndarray:
    X = pockets[list(bundle["features"])].astype(float)
    raw = protocol._scores(bundle["model"], X)
    return protocol._calibrate(bundle["selection"], raw)


# -------------------------------------------------------------- sequences
def parse_fasta(text: str) -> Dict[str, str]:
    """Accession -> sequence, from UniProt FASTA (``>sp|P12345|NAME ...``) or bare headers."""
    records: Dict[str, List[str]] = {}
    current = None
    for line in text.splitlines():
        if line.startswith(">"):
            head = line[1:].split()[0]
            current = head.split("|")[1] if head.count("|") >= 2 else head
            records[current] = []
        elif current is not None:
            records[current].append(line.strip())
    return {k: "".join(v) for k, v in records.items() if v}


def fetch_uniprot_fasta(accessions: Sequence[str], batch: int = 150) -> Dict[str, str]:
    out: Dict[str, str] = {}
    unique = sorted({a.strip().upper() for a in accessions if a and a.strip()})
    for start in range(0, len(unique), batch):
        chunk = unique[start:start + batch]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode({"format": "fasta", "query": query})
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=120) as response:
                    out.update(parse_fasta(response.read().decode()))
                break
            except Exception:  # transient network error: back off and retry, then fail loudly
                if attempt == 4:
                    raise
                time.sleep(2 ** attempt)
    return out


def write_fasta(records: Mapping[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for acc, seq in sorted(records.items()):
            fh.write(f">{acc}\n{seq}\n")


def seen_proteins(hits_path: Path) -> Set[str]:
    """Queries with a hit meeting the benchmark's homology criterion."""
    if not hits_path.exists() or hits_path.stat().st_size == 0:
        return set()
    hits = pd.read_csv(hits_path, sep="\t", header=None, names=list(SEARCH_COLUMNS))
    ident = hits["fident"].astype(float)
    ident = np.where(ident > 1.0, ident / 100.0, ident)
    shorter = np.minimum(hits["qlen"], hits["tlen"]).astype(float)
    coverage = np.minimum(hits["qend"] - hits["qstart"] + 1, hits["tend"] - hits["tstart"] + 1) / shorter
    keep = (ident >= MIN_IDENTITY) & (coverage >= MIN_SHORTER_COVERAGE) & (hits["evalue"] <= MAX_EVALUE)
    return set(hits.loc[keep, "query"].astype(str))


def read_clusters(path: Optional[Path]) -> Dict[str, str]:
    if path is None or not path.exists():
        return {}
    table = pd.read_csv(path, sep="\t", header=None, names=["rep", "member"], dtype=str)
    return dict(zip(table["member"], table["rep"]))


# -------------------------------------------------------------- evaluation
def protein_table(pockets: pd.DataFrame, min_plddt: float = MIN_PLDDT) -> pd.DataFrame:
    confident = pockets[pockets["plddt_mean"].astype(float).fillna(0) >= min_plddt]
    best = confident.sort_values("learned_score", ascending=False).drop_duplicates("uniprot_id")
    rule = confident.groupby("uniprot_id")["composite_score"].max().rename("rule_score")
    keep = [c for c in ("uniprot_id", "learned_score", "pocket_id", "pocket_residues", "hull_depth",
                        "plddt_mean", "volume") if c in best.columns]
    return best[keep].rename(columns={"pocket_id": "top_pocket", "hull_depth": "top_pocket_hull_depth",
                                      "plddt_mean": "top_pocket_plddt", "volume": "top_pocket_volume",
                                      "pocket_residues": "top_pocket_residues"}).merge(
        rule, on="uniprot_id", how="left")


def evaluate_ranking(proteins: pd.DataFrame, n_bootstrap: int) -> Dict[str, object]:
    y = proteins["annotated"].to_numpy(dtype=int)
    if y.sum() == 0 or y.sum() == len(y):
        return {"proteins": int(len(y)), "binders": int(y.sum()), "evaluable": False}
    groups = proteins["cluster"].astype(str).to_numpy()
    learned = protocol.Ranked(y, proteins["learned_score"].to_numpy(dtype=float))
    rule = protocol.Ranked(y, proteins["rule_score"].fillna(proteins["rule_score"].min()).to_numpy(dtype=float))
    l1 = protocol.bootstrap_statistic(learned.roc_auc, groups, null=0.5, n_bootstrap=n_bootstrap, seed=SEED)
    l2 = protocol.bootstrap_statistic(lambda w: learned.roc_auc(w) - rule.roc_auc(w), groups,
                                      n_bootstrap=n_bootstrap, seed=SEED)
    rule_auc = protocol.bootstrap_statistic(rule.roc_auc, groups, null=0.5, n_bootstrap=n_bootstrap, seed=SEED)
    return {"proteins": int(len(y)), "binders": int(y.sum()), "evaluable": True,
            "L1_learned_roc_auc": l1.as_dict(), "rule_roc_auc": rule_auc.as_dict(),
            "L2_learned_minus_rule": l2.as_dict()}


def candidates(proteins: pd.DataFrame, n: int = N_CANDIDATES) -> pd.DataFrame:
    ranked = proteins.sort_values("learned_score", ascending=False).reset_index(drop=True)
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    ranked["precision_at_rank"] = ranked["annotated"].cumsum() / ranked["rank"]
    return ranked[~ranked["annotated"]].head(n)


def cmd_train(args: argparse.Namespace) -> int:
    table = pd.read_csv(args.table, low_memory=False)
    bundle = train_locked(table, n_draws=args.n_draws)
    bundle["table_sha256"] = _sha256(args.table)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(bundle, args.output)
    print(json.dumps({k: bundle[k] for k in ("candidate", "rows", "positives", "table_sha256")}, default=str))
    return 0


def cmd_fasta(args: argparse.Namespace) -> int:
    accessions: List[str] = []
    for path in args.accessions_csv:
        frame = pd.read_csv(path, dtype=str)
        column = next(c for c in args.columns if c in frame.columns)
        for value in frame[column].dropna():
            accessions += [a for a in str(value).replace(",", ";").split(";") if a.strip()]
    records = fetch_uniprot_fasta(accessions)
    write_fasta(records, args.output)
    print(f"{len(set(accessions))} accessions -> {len(records)} sequences")
    return 0


def _read_many(paths: Iterable[Path], **kwargs) -> pd.DataFrame:
    frames = [pd.read_csv(p, **kwargs) for p in paths]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def cmd_evaluate(args: argparse.Namespace) -> int:
    bundle = joblib.load(args.model)
    pockets = _read_many(sorted(args.shards_dir.rglob("*_pockets_part*.csv.gz")), low_memory=False)
    catalog = _read_many(sorted(args.catalog_dir.rglob("*_catalog.csv")))
    annotations = _read_many(sorted(args.catalog_dir.rglob("*_uniprot.tsv")), sep="\t", dtype=str)
    if pockets.empty:
        raise SystemExit("no screened pockets")
    from cryptic_ip.analysis.scorer import PocketScorer

    pockets["composite_score"] = PocketScorer().score_frame(pockets)
    pockets["learned_score"] = score_pockets(bundle, pockets)
    proteins = protein_table(pockets)
    organism = catalog.drop_duplicates("uniprot_id").set_index("uniprot_id")["organism_key"]
    proteins["organism_key"] = proteins["uniprot_id"].map(organism)
    known = known_ip_annotation(annotations) if not annotations.empty else pd.Series(dtype=bool)
    proteins["annotated"] = proteins["uniprot_id"].map(known).fillna(False).astype(bool)
    seen = seen_proteins(args.seen_hits)
    proteins["seen"] = proteins["uniprot_id"].isin(seen)
    clusters = read_clusters(args.clusters)
    proteins["cluster"] = proteins["uniprot_id"].map(clusters).fillna(proteins["uniprot_id"])
    unseen = proteins[~proteins["seen"]]

    report: Dict[str, object] = {"plan": "docs/LEARNED_SCREEN_PLAN.md", "model": bundle["candidate"],
                                 "model_table_sha256": bundle.get("table_sha256"),
                                 "proteins_scored": int(len(proteins)), "proteins_seen": int(proteins["seen"].sum()),
                                 "organisms": {}}
    for org, part in unseen.groupby("organism_key"):
        result = evaluate_ranking(part, args.n_bootstrap)
        result["decision_eligible"] = result.get("binders", 0) >= MIN_BINDERS_FOR_DECISION
        report["organisms"][org] = result
    pooled = evaluate_ranking(unseen, args.n_bootstrap)
    report["pooled"] = pooled
    if pooled.get("evaluable"):
        adjusted = protocol.holm({"L1": pooled["L1_learned_roc_auc"]["p_value"],
                                  "L2": pooled["L2_learned_minus_rule"]["p_value"]})
        l1_ok = pooled["L1_learned_roc_auc"]["low"] > 0.5 and adjusted["L1"] < 0.05
        l2_ok = pooled["L2_learned_minus_rule"]["low"] > 0 and adjusted["L2"] < 0.05
        report["decisions"] = {
            "L1": "supported" if l1_ok else "not supported",
            "L2": "supported" if l2_ok else "not supported",
            "holm": adjusted,
        }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for org, part in unseen.groupby("organism_key"):
        frames.append(candidates(part))
    cand = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if not annotations.empty and not cand.empty:
        cand = cand.merge(annotations[["uniprot_id", "gene", "protein_name"]].drop_duplicates("uniprot_id"),
                          on="uniprot_id", how="left")
    cand.to_csv(args.output_dir / "candidates.csv", index=False)
    proteins.to_csv(args.output_dir / "proteins.csv.gz", index=False)
    (args.output_dir / "learned_screen.json").write_text(json.dumps(report, indent=2, default=float))
    lines = ["## Learned proteome ranking (docs/LEARNED_SCREEN_PLAN.md)", "",
             f"{report['proteins_scored']} proteins scored; {report['proteins_seen']} have a homologue among "
             "the benchmark's proteins and are excluded from evaluation.", "",
             "| organism | unseen proteins | annotated binders | learned ROC-AUC | rule ROC-AUC | learned − rule |",
             "|---|---|---|---|---|---|"]

    def ci(d):
        return f"{d['point']:.3f} [{d['low']:.3f}, {d['high']:.3f}]"

    for org, r in list(report["organisms"].items()) + [("pooled", pooled)]:
        if r.get("evaluable"):
            lines.append(f"| {org} | {r['proteins']} | {r['binders']} | {ci(r['L1_learned_roc_auc'])} | "
                         f"{ci(r['rule_roc_auc'])} | {ci(r['L2_learned_minus_rule'])} |")
        else:
            lines.append(f"| {org} | {r['proteins']} | {r['binders']} | – | – | – |")
    if "decisions" in report:
        lines += ["", f"Decisions (Holm): L1 {report['decisions']['L1']}, L2 {report['decisions']['L2']}."]
    if not cand.empty:
        cols = [c for c in ("organism_key", "rank", "uniprot_id", "gene", "protein_name", "learned_score",
                            "top_pocket_hull_depth", "top_pocket_plddt", "precision_at_rank") if c in cand.columns]
        lines += ["", "### Candidates (hypotheses, not findings)", "",
                  cand[cols].to_markdown(index=False, floatfmt=".3f")]
    text = "\n".join(lines) + "\n"
    (args.output_dir / "LEARNED_SCREEN.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    t = sub.add_parser("train")
    t.add_argument("--table", type=Path, required=True)
    t.add_argument("--output", type=Path, required=True)
    t.add_argument("--n-draws", type=int, default=10)
    f = sub.add_parser("fasta")
    f.add_argument("--accessions-csv", type=Path, nargs="+", required=True)
    f.add_argument("--columns", nargs="+", default=["uniprot_ids", "uniprot_id"])
    f.add_argument("--output", type=Path, required=True)
    e = sub.add_parser("evaluate")
    e.add_argument("--model", type=Path, required=True)
    e.add_argument("--shards-dir", type=Path, required=True)
    e.add_argument("--catalog-dir", type=Path, required=True)
    e.add_argument("--seen-hits", type=Path, required=True)
    e.add_argument("--clusters", type=Path, default=None)
    e.add_argument("--output-dir", type=Path, required=True)
    e.add_argument("--n-bootstrap", type=int, default=2000)
    args = parser.parse_args(argv)
    return {"train": cmd_train, "fasta": cmd_fasta, "evaluate": cmd_evaluate}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
