#!/usr/bin/env python3
"""Render the benchmark report as one self-contained HTML page.

Every estimate in the benchmark is an interval, not a number, and the plan's
decisions turn on where those intervals sit. A table of numbers hides that; a
forest plot shows it at a glance - which intervals clear zero, which straddle
it, and which are so wide that the evaluation could not have decided anything.

The page is a single file with inline SVG and no external assets, so it renders
from a CI artifact, offline, years from now.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

ROW_HEIGHT = 30
PLOT_WIDTH = 470
LABEL_WIDTH = 290
VALUE_WIDTH = 400
MARGIN = 14

DECISION_CLASS = {"supported": "ok", "refuted": "no", "inconclusive": "meh"}


def _finite(*values: float) -> List[float]:
    return [v for v in values if v is not None and isinstance(v, (int, float)) and math.isfinite(v)]


class Forest:
    """A forest plot: one interval per row, against a reference value."""

    def __init__(
        self,
        rows: Sequence[Tuple[str, Dict[str, float], str]],
        *,
        null: float = 0.0,
        margin: Optional[float] = None,
        domain: Optional[Tuple[float, float]] = None,
        excluded_class: str = "hit",
    ) -> None:
        self.rows = list(rows)
        self.excluded_class = excluded_class
        self.null = null
        self.margin = margin
        values: List[float] = []
        for _, estimate, _ in self.rows:
            values += _finite(estimate.get("low"), estimate.get("high"), estimate.get("point"))
        values.append(null)
        if margin:
            values += [null - margin, null + margin]
        low, high = (min(values), max(values)) if values else (0.0, 1.0)
        if domain:
            low, high = min(low, domain[0]), max(high, domain[1])
        pad = max((high - low) * 0.08, 1e-3)
        self.lo, self.hi = low - pad, high + pad

    def _x(self, value: float) -> float:
        return LABEL_WIDTH + (value - self.lo) / (self.hi - self.lo) * PLOT_WIDTH

    def _ticks(self) -> List[float]:
        span = self.hi - self.lo
        step = 10 ** math.floor(math.log10(span / 3.5))
        for multiple in (1, 2, 5, 10):
            if span / (step * multiple) <= 5:
                step *= multiple
                break
        first = math.ceil(self.lo / step) * step
        return [first + i * step for i in range(int(span / step) + 1) if first + i * step <= self.hi]

    def svg(self) -> str:
        height = len(self.rows) * ROW_HEIGHT + 46
        parts = [
            f'<svg class="forest" viewBox="0 0 {LABEL_WIDTH + PLOT_WIDTH + VALUE_WIDTH} {height}" '
            f'role="img" aria-label="interval plot">'
        ]
        if self.margin:
            x0, x1 = self._x(self.null - self.margin), self._x(self.null + self.margin)
            parts.append(
                f'<rect x="{x0:.1f}" y="6" width="{x1 - x0:.1f}" height="{height - 40}" class="band"/>'
                f'<text x="{(x0 + x1) / 2:.1f}" y="{height - 20}" class="tick mid">no effect of any useful size</text>'
            )
        for tick in self._ticks():
            x = self._x(tick)
            parts.append(
                f'<line x1="{x:.1f}" y1="6" x2="{x:.1f}" y2="{height - 40}" class="grid"/>'
                f'<text x="{x:.1f}" y="{height - 26}" class="tick mid">{tick:g}</text>'
            )
        null_x = self._x(self.null)
        parts.append(f'<line x1="{null_x:.1f}" y1="6" x2="{null_x:.1f}" y2="{height - 40}" class="nullline"/>')

        for index, (label, estimate, note) in enumerate(self.rows):
            y = 6 + index * ROW_HEIGHT + ROW_HEIGHT / 2
            point, low, high = estimate.get("point"), estimate.get("low"), estimate.get("high")
            parts.append(
                f'<text x="{LABEL_WIDTH - 12}" y="{y + 4:.1f}" class="rowlabel">{html.escape(label)}</text>'
            )
            if not _finite(point, low, high) or len(_finite(point, low, high)) < 3:
                parts.append(f'<text x="{LABEL_WIDTH + 8}" y="{y + 4:.1f}" class="missing">not measured</text>')
                continue
            excludes = low > self.null or high < self.null
            klass = self.excluded_class if excludes else "null"
            parts.append(
                f'<line x1="{self._x(low):.1f}" y1="{y:.1f}" x2="{self._x(high):.1f}" y2="{y:.1f}" '
                f'class="ci {klass}"/>'
                f'<circle cx="{self._x(point):.1f}" cy="{y:.1f}" r="4.5" class="pt {klass}"/>'
                f'<text x="{LABEL_WIDTH + PLOT_WIDTH + 10}" y="{y + 4:.1f}" class="value">'
                f'{point:.3f} <tspan class="dim">[{low:.3f}, {high:.3f}]</tspan>'
                f'{" " + html.escape(note) if note else ""}</text>'
            )
        parts.append("</svg>")
        return "".join(parts)


def _decision_badge(decision: str) -> str:
    klass = next((v for k, v in DECISION_CLASS.items() if decision.startswith(k)), "meh")
    return f'<span class="badge {klass}">{html.escape(decision)}</span>'


def hypothesis_section(report: Dict[str, object]) -> str:
    titles = {
        "H1": "H1 — hull depth improves finding <em>buried</em> inositol phosphate sites",
        "H2": "H2 — hull depth improves telling buried sites from surface sites",
    }
    blocks: List[str] = []
    for name, result in (report.get("hypotheses") or {}).items():
        rows: List[Tuple[str, Dict[str, float], str]] = []
        for grouping in ("sequence", "strict"):
            estimate = (result.get("development") or {}).get(grouping, {}).get("roc_auc")
            rows.append((f"development, {grouping} grouping", estimate or {}, ""))
        holdout = result.get("holdout") or {}
        if holdout:
            groups = holdout.get("positive_groups", 0)
            rows.append(("holdout (locked model)", holdout.get("roc_auc", {}), f"· {groups} positive groups"))
        holm = result.get("holm_p")
        holm_text = "n/a" if holm is None or not math.isfinite(holm) else f"{holm:.3g}"
        extra = ""
        if result.get("not_evaluable_runs"):
            extra = (
                '<p class="note">Runs that could not be evaluated: '
                + ", ".join(html.escape(s) for s in result["not_evaluable_runs"])
                + "</p>"
            )
        if result.get("unsupported_groupings"):
            extra += (
                '<p class="note">Grouping(s) where one homology group holds too many positives for a '
                "five-fold evaluation: " + ", ".join(html.escape(g) for g in result["unsupported_groupings"]) + "</p>"
            )
        blocks.append(
            f'<article class="card"><header><h3>{titles.get(name, html.escape(name))}</h3>'
            f'<div class="meta">task <code>{html.escape(str(result.get("task", "")))}</code>'
            f' · Holm-adjusted p {holm_text} · {_decision_badge(str(result.get("decision", "?")))}</div></header>'
            f'<p class="axis">ROC-AUC with hull depth minus without, paired on identical folds; '
            f"95 % intervals from resampling homology groups.</p>"
            f'{Forest(rows, null=0.0, margin=0.01).svg()}{extra}</article>'
        )
    return "\n".join(blocks)


def evaluations_section(report: Dict[str, object]) -> str:
    evaluations = report.get("evaluations") or {}
    by_split: Dict[str, List[Tuple[str, Dict[str, float], str]]] = {}
    ml_vs_rule: List[Tuple[str, Dict[str, float], str]] = []
    for key, evaluation in sorted(evaluations.items()):
        if "/permuted" in key:
            continue
        split = "holdout" if "/holdout" in key else "cross-validation"
        note = f"· {evaluation.get('positives', 0)} positives in {evaluation.get('positive_groups', 0)} groups"
        by_split.setdefault(split, []).append((key, evaluation.get("roc_auc", {}), note))
        by_split.setdefault(split, [])
        ml_vs_rule.append((key, evaluation.get("ml_minus_rule_roc_auc", {}), ""))
    blocks = []
    for split, rows in by_split.items():
        blocks.append(
            f'<article class="card"><h3>Discrimination ({html.escape(split)})</h3>'
            f'<p class="axis">ROC-AUC of the learned model. 0.5 is chance.</p>'
            f'{Forest(rows, null=0.5, domain=(0.5, 1.0)).svg()}</article>'
        )
    if ml_vs_rule:
        blocks.append(
            '<article class="card"><h3>Learned model minus rule-based score</h3>'
            '<p class="axis">Positive means the learned model ranks sites better than the '
            "transparent rule-based score on the same pockets.</p>"
            f"{Forest(ml_vs_rule, null=0.0).svg()}</article>"
        )
    return "\n".join(blocks)


def permutation_section(report: Dict[str, object]) -> str:
    permutation = report.get("permutation") or {}
    if not permutation:
        return ""
    rows = [
        (f"{task} (labels shuffled)", value.get("roc_auc", {}), "" if value.get("passes") else "· FAIL")
        for task, value in sorted(permutation.items())
    ]
    failed = [task for task, value in permutation.items() if not value.get("passes")]
    verdict = (
        '<p class="note bad">A control failed: with its labels shuffled, the pipeline still scored '
        "above chance for that task. That is either a chance excursion of the control or a leak; "
        "until a null distribution over many shuffles settles which, results for that task are not "
        "evidence.</p>"
        if failed
        else '<p class="note">Every control includes 0.5, as it must: with the labels shuffled the '
        "pipeline learns nothing, so no label information leaks through it.</p>"
    )
    return (
        '<article class="card"><h3>Leak control</h3>'
        '<p class="axis">The whole pipeline re-run with shuffled labels. Each interval must contain 0.5.</p>'
        f"{Forest(rows, null=0.5, excluded_class='bad').svg()}{verdict}</article>"
    )


def data_section(report: Dict[str, object]) -> str:
    data = report.get("data") or {}
    if not data:
        return ""
    cells = [
        ("entries", data.get("n_entries")),
        ("pockets", data.get("n_pockets")),
        ("holdout entries", data.get("holdout_entries")),
    ]
    for name, count in (data.get("groups") or {}).items():
        cells.append((f"{name} groups", count))
    tiles = "".join(
        f'<div class="tile"><div class="n">{html.escape(str(value))}</div>'
        f'<div class="k">{html.escape(str(label))}</div></div>'
        for label, value in cells
        if value is not None
    )
    rows = []
    for task, info in (data.get("tasks") or {}).items():
        development = info.get("development", {})
        holdout = info.get("holdout", {})
        share = info.get("largest_group_positive_share", {})
        rows.append(
            f"<tr><td><code>{html.escape(task)}</code></td>"
            f"<td>{development.get('pockets', 0):,}</td><td>{development.get('positives', 0):,}</td>"
            f"<td>{(development.get('positive_groups') or {}).get('sequence', 0)}"
            f" / {(development.get('positive_groups') or {}).get('strict', 0)}</td>"
            f"<td>{holdout.get('positives', 0):,}</td>"
            f"<td>{share.get('sequence', float('nan')):.0%} / {share.get('strict', float('nan')):.0%}</td></tr>"
        )
    table = (
        "<table><thead><tr><th>task</th><th>pockets</th><th>positives</th>"
        "<th>positive groups (seq / strict)</th><th>holdout positives</th>"
        "<th>largest group's share of positives</th></tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table>"
        if rows
        else ""
    )
    return f'<section class="card"><h3>What was measured</h3><div class="tiles">{tiles}</div>{table}</section>'


def render(report: Dict[str, object], title: str) -> str:
    smoke = bool(report.get("smoke"))
    banner = (
        '<p class="banner">Smoke run: every label was shuffled. Nothing on this page is a result.</p>'
        if smoke
        else ""
    )
    not_evaluable = report.get("not_evaluable") or {}
    not_evaluable_block = (
        '<article class="card"><h3>Not evaluable</h3><ul>'
        + "".join(
            f"<li><code>{html.escape(stem)}</code> — {html.escape(str(reason))}</li>"
            for stem, reason in not_evaluable.items()
        )
        + "</ul><p class=\"note\">Reported rather than dropped: an evaluation the data cannot support is "
        "a finding about the data.</p></article>"
        if not_evaluable
        else ""
    )
    return f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>
:root {{
  --bg: #fbfbfa; --card: #fff; --ink: #1c1b19; --dim: #6b6862; --line: #e4e1dc;
  --hit: #1a7f5a; --null: #8a8681; --bad: #b3261e; --band: #f0efeb; --accent: #2f5d8a;
}}
@media (prefers-color-scheme: dark) {{
  :root:not([data-theme="light"]) {{
    --bg: #17181a; --card: #1e2022; --ink: #e9e7e4; --dim: #9a978f; --line: #303336;
    --hit: #4cc38a; --null: #7d7a75; --bad: #f2645a; --band: #26292c; --accent: #7fb2e5;
  }}
}}
:root[data-theme="dark"] {{
  --bg: #17181a; --card: #1e2022; --ink: #e9e7e4; --dim: #9a978f; --line: #303336;
  --hit: #4cc38a; --null: #7d7a75; --bad: #f2645a; --band: #26292c; --accent: #7fb2e5;
}}
* {{ box-sizing: border-box; }}
body {{ margin: 0; background: var(--bg); color: var(--ink);
  font: 15px/1.5 ui-sans-serif, system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }}
main {{ max-width: 1120px; margin: 0 auto; padding: 32px 16px 64px; }}
h1 {{ font-size: 1.6rem; margin: 0 0 4px; letter-spacing: -0.01em; }}
h2 {{ font-size: 1.1rem; margin: 34px 0 12px; color: var(--dim); text-transform: uppercase;
  letter-spacing: 0.08em; font-weight: 600; }}
h3 {{ font-size: 1.02rem; margin: 0 0 6px; font-weight: 650; }}
.sub {{ color: var(--dim); margin: 0 0 8px; }}
.card {{ background: var(--card); border: 1px solid var(--line); border-radius: 12px;
  padding: 18px 20px; margin: 14px 0; overflow-x: auto; }}
.meta {{ color: var(--dim); font-size: 0.9rem; margin-bottom: 6px; }}
.axis {{ color: var(--dim); font-size: 0.9rem; margin: 2px 0 10px; }}
.note {{ color: var(--dim); font-size: 0.88rem; margin: 10px 0 0; }}
.note.bad {{ color: var(--bad); }}
.banner {{ background: var(--bad); color: #fff; padding: 10px 14px; border-radius: 10px; font-weight: 600; }}
code {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: 0.92em; }}
.badge {{ display: inline-block; padding: 2px 10px; border-radius: 999px; font-size: 0.82rem;
  font-weight: 650; border: 1px solid currentColor; }}
.badge.ok {{ color: var(--hit); }} .badge.no {{ color: var(--bad); }} .badge.meh {{ color: var(--dim); }}
.forest {{ width: 100%; min-width: 760px; height: auto; display: block; }}
.forest .rowlabel {{ fill: var(--ink); font-size: 12.5px; text-anchor: end; }}
.forest .value {{ fill: var(--ink); font-size: 12.5px; font-family: ui-monospace, Menlo, monospace; }}
.forest .dim {{ fill: var(--dim); }}
.forest .missing {{ fill: var(--dim); font-size: 12.5px; font-style: italic; }}
.forest .tick {{ fill: var(--dim); font-size: 11px; }}
.forest .mid {{ text-anchor: middle; }}
.forest .grid {{ stroke: var(--line); stroke-width: 1; }}
.forest .nullline {{ stroke: var(--dim); stroke-width: 1.5; stroke-dasharray: 4 3; }}
.forest .band {{ fill: var(--band); }}
.forest .ci {{ stroke-width: 2.5; stroke-linecap: round; }}
.forest .ci.hit, .forest .pt.hit {{ stroke: var(--hit); fill: var(--hit); }}
.forest .ci.bad, .forest .pt.bad {{ stroke: var(--bad); fill: var(--bad); }}
.forest .ci.null, .forest .pt.null {{ stroke: var(--null); fill: var(--null); }}
.tiles {{ display: flex; flex-wrap: wrap; gap: 10px; margin: 10px 0 16px; }}
.tile {{ border: 1px solid var(--line); border-radius: 10px; padding: 10px 14px; min-width: 104px; }}
.tile .n {{ font-size: 1.35rem; font-weight: 650; }}
.tile .k {{ color: var(--dim); font-size: 0.82rem; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.9rem; }}
th, td {{ text-align: left; padding: 7px 10px; border-bottom: 1px solid var(--line); }}
th {{ color: var(--dim); font-weight: 600; }}
footer {{ color: var(--dim); font-size: 0.85rem; margin-top: 30px; }}
@media (max-width: 640px) {{ main {{ padding: 20px 16px 48px; }} }}
</style></head>
<body><main>
<h1>{html.escape(title)}</h1>
<p class="sub">Pre-registered in <code>docs/ANALYSIS_PLAN.md</code>. Every interval is a 95 % percentile
interval from resampling whole homology groups, so proteins - not pockets - are the unit of uncertainty.</p>
{banner}
{data_section(report)}
<h2>Primary hypotheses</h2>
{hypothesis_section(report)}
<h2>How well anything works at all</h2>
{evaluations_section(report)}
<h2>Controls</h2>
{permutation_section(report)}
{not_evaluable_block}
<footer>Green intervals exclude the reference value; grey intervals include it. Generated from
<code>benchmark_report.json</code>.</footer>
</main></body></html>
"""


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-json", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="Cryptic inositol phosphate sites — benchmark")
    args = parser.parse_args(argv)

    report = json.loads(args.report_json.read_text(encoding="utf-8"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(render(report, args.title), encoding="utf-8")
    print(f"wrote {args.output} ({args.output.stat().st_size:,} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
