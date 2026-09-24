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
        "buried easier", "chance")
BAD = ("not supported", "unreliable", "not trustworthy", "does not", "not learnable", "buried harder", "leak")


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
