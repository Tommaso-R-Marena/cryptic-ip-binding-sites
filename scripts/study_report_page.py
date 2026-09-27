#!/usr/bin/env python3
"""Self-contained HTML pages for the redocking, α-arrestin, specificity and hull-gate studies.

Built from the benchmark page's parts (``Forest`` and ``PAGE_CSS`` in
``scripts/benchmark_report_page.py``): inline SVG forest plots, no external
assets. A page is a list of sections; each section holds decision rows,
forest plots and tables, all taken from a study's result JSON.
"""

from __future__ import annotations

import html
import sys
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

from benchmark_report_page import PAGE_CSS, Forest  # noqa: E402

GOOD = ("supported", "reliable", "trustworthy", "discriminates", "learnable", "remove", "relax",
        "buried easier", "chance", "pass", "improves")
BAD = ("not supported", "unreliable", "not trustworthy", "does not", "not learnable", "buried harder", "leak", "fail",
       "worsens")


def badge(decision: str) -> str:
    text = str(decision)
    klass = "meh"
    if any(text.startswith(b) for b in BAD):
        klass = "no"
    elif any(text.startswith(g) for g in GOOD):
        klass = "ok"
    return f'<span class="badge {klass}">{html.escape(text)}</span>'


def forest(rows: Sequence[Tuple[str, Optional[Mapping[str, float]], str]], *, null: float,
           caption: str = "") -> str:
    clean = [(label, dict(est or {}), note) for label, est, note in rows]
    svg = Forest(clean, null=null).svg()
    cap = f'<p class="axis">{html.escape(caption)}</p>' if caption else ""
    return f'<article class="card">{cap}{svg}</article>'


def table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    head = "".join(f"<th>{html.escape(str(h))}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c if str(c).startswith('<span') else html.escape(str(c))}</td>"
                                    for c in row) + "</tr>" for row in rows)
    return f'<article class="card"><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></article>'


def decisions(rows: Sequence[Tuple[str, str, str]]) -> str:
    """(name, decision, detail) rows."""
    return table(["question", "decision", "detail"], [(n, badge(d), detail) for n, d, detail in rows])


def page(title: str, plan: str, sections: Sequence[Tuple[str, List[str]]], footer: str = "") -> str:
    body = "".join(f"<h2>{html.escape(name)}</h2>{''.join(parts)}" for name, parts in sections)
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
{PAGE_CSS}</style></head>
<body><main>
<h1>{html.escape(title)}</h1>
<p class="sub">Pre-registered in <code>{html.escape(plan)}</code>. Intervals are 95 % percentile intervals from
resampling whole homology groups (or sequence clusters), never structures, copies or pockets.</p>
{body}
<footer>Green intervals exclude the reference value; grey intervals include it. {html.escape(footer)}</footer>
</main></body></html>
"""


def fmt(est: Optional[Mapping[str, float]], digits: int = 3) -> str:
    if not est or "point" not in est:
        return "–"
    return f"{est['point']:.{digits}f} [{est['low']:.{digits}f}, {est['high']:.{digits}f}]"


def evidence_note(entry: Mapping[str, object]) -> str:
    if entry.get("evidence") is False:
        return "fewer than 5 groups: not evidence"
    return ""


def as_rows(entries: Mapping[str, Mapping[str, object]], key: str = "per_group") -> List[Tuple[str, Dict, str]]:
    return [(name, dict(e.get(key) or {}), evidence_note(e)) for name, e in entries.items()]


def arrestin_page(report: Mapping[str, object]) -> str:
    """The α-arrestin study (results/arrestin/arrestin.json, built by ``scripts/arrestin.py decide``)."""
    b1 = report["B1"]
    validity = report["protocol_validity"]
    rows = [("B1 family ranks high among unseen proteins", b1["decision"],
             f"ROC-AUC {fmt(b1['roc_auc'])}, 5th percentile {b1['roc_auc']['p5']:.3f}; "
             f"{b1['members_unseen']} unseen α-arrestins in {b1['member_clusters']} clusters"),
            ("protocol validity (crystal redocks of arrestin–IP sites)",
             "valid" if validity["valid"] else "not valid",
             f"{validity['crystal_sites']} sites, mean top-pose success {validity['mean_top_pose_success']:.2f} "
             "(needs ≥ 0.5)")]
    for acc, p in report["proteins"].items():
        failing = ", ".join(p["failing"]) or "none"
        rows.append((f"B6 {p['gene']} ({acc})", p["verdict"], f"failing: {failing}"))
    organisms = [(org, e["roc_auc"], e["decision"] if e["decision"].startswith("not evaluable") else "")
                 for org, e in report["B1_per_organism_descriptive"].items()]
    members = table(["protein", "organism", "rank among unseen", "percentile", "cluster"],
                    [(m["uniprot_id"], m["organism_key"], m["rank"], f"{m['rank_percentile']:.2f}", m["cluster"])
                     for m in b1["members"]])
    criteria = []
    for acc, p in report["proteins"].items():
        for name, c in p["criteria"].items():
            detail = c.get("decision") or ""
            if "jaccard" in c:
                detail = f"Jaccard {c['jaccard']:.2f} vs {c['reference']} (TM-score {c['tm_score']:.3f})"
            criteria.append((f"{p['gene']} ({acc})", name, badge("pass" if c["pass"] else "fail"), detail))
    scores = []
    for d in report["docking"]:
        if d.get("kind") in ("lead_site", "negative") and d.get("ligand") == "IHP":
            scores.append((d["gene"], d["kind"], d["site"], f"{d['best_score']:.3f}"))
    scores.append(("controls", "weakest positive", "", f"{report['weakest_positive']:.3f}"))
    return page("The α-arrestin lead", str(report["plan"]), [
        ("Decisions", [decisions(rows)]),
        ("B1: the family among unseen proteins", [
            forest([("all organisms", b1["roc_auc"], "")] + organisms, null=0.5,
                   caption="ROC-AUC of α-arrestins against all other unseen proteins (per organism: descriptive)"),
            members]),
        ("B2–B6: the two leads", [table(["protein", "criterion", "result", "detail"], criteria)]),
        ("Docking scores (descriptive; the protocol failed its validity gate)",
         [table(["protein", "kind", "site", "best IP6 score (kcal/mol)"], scores)]),
    ], footer="Docking scores are not evidence here: no crystal arrestin–IP pose was reproduced.")


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(description="Rebuild a study page from its committed result JSON.")
    parser.add_argument("study", choices=["arrestin"])
    parser.add_argument("json", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(argv)
    args.output.write_text(arrestin_page(json.loads(args.json.read_text())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
