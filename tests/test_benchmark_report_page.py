"""The HTML report page: built from a real report, self-contained, honest."""

from __future__ import annotations

import json
import re

import pytest

from scripts import benchmark_report_page as page


def _estimate(point, low, high, p=0.01):
    return {"point": point, "low": low, "high": high, "p_value": p}


@pytest.fixture
def report():
    return {
        "data": {
            "n_entries": 409, "n_pockets": 51234, "holdout_entries": 80,
            "groups": {"sequence": 85, "strict": 40},
            "tasks": {
                "burial": {
                    "development": {"pockets": 300, "positives": 120, "positive_groups": {"sequence": 30, "strict": 12}},
                    "holdout": {"pockets": 60, "positives": 20, "positive_groups": {"sequence": 6, "strict": 3}},
                    "largest_group_positive_share": {"sequence": 0.2, "strict": 0.55},
                }
            },
        },
        "evaluations": {
            "burial/full/sequence/cv": {
                "pockets": 300, "positives": 120, "groups": 40, "positive_groups": 30, "repeats": 3,
                "roc_auc": _estimate(0.81, 0.74, 0.88), "pr_auc": _estimate(0.6, 0.5, 0.7),
                "rule_roc_auc": _estimate(0.7, 0.62, 0.78), "rule_pr_auc": _estimate(0.4, 0.3, 0.5),
                "ml_minus_rule_roc_auc": _estimate(0.11, 0.04, 0.18),
                "ml_minus_rule_pr_auc": _estimate(0.2, 0.1, 0.3),
                "mcc": 0.5, "rule_mcc": 0.3, "families_chosen": {"extra_trees": 9},
            },
            "burial/full/sequence/holdout": {
                "pockets": 60, "positives": 20, "groups": 8, "positive_groups": 6, "repeats": 1,
                "roc_auc": _estimate(0.77, 0.6, 0.9), "pr_auc": _estimate(0.5, 0.3, 0.7),
                "rule_roc_auc": _estimate(0.66, 0.5, 0.8), "rule_pr_auc": _estimate(0.35, 0.2, 0.5),
                "ml_minus_rule_roc_auc": _estimate(0.11, -0.05, 0.25),
                "ml_minus_rule_pr_auc": _estimate(0.15, -0.1, 0.4),
                "mcc": 0.4, "rule_mcc": 0.25, "families_chosen": {"extra_trees": 1},
            },
            "burial/full/sequence/cv/permuted": {
                "pockets": 300, "positives": 120, "groups": 40, "positive_groups": 30, "repeats": 1,
                "roc_auc": _estimate(0.5, 0.43, 0.57), "pr_auc": _estimate(0.4, 0.3, 0.5),
                "rule_roc_auc": _estimate(0.5, 0.42, 0.58), "rule_pr_auc": _estimate(0.4, 0.3, 0.5),
                "ml_minus_rule_roc_auc": _estimate(0.0, -0.08, 0.08),
                "ml_minus_rule_pr_auc": _estimate(0.0, -0.1, 0.1),
                "mcc": 0.0, "rule_mcc": 0.0, "families_chosen": {"logistic_regression": 3},
            },
        },
        "hypotheses": {
            "H1": {
                "task": "cryptic_ip_site",
                "development": {
                    "sequence": {"roc_auc": _estimate(0.002, -0.004, 0.008, 0.6), "pr_auc": _estimate(0.0, -0.02, 0.02)},
                    "strict": {"roc_auc": _estimate(0.001, -0.006, 0.007, 0.8), "pr_auc": _estimate(0.0, -0.03, 0.03)},
                },
                "holm_p": 0.8, "decision": "refuted",
            },
            "H2": {
                "task": "burial",
                "development": {
                    "sequence": {"roc_auc": _estimate(0.05, 0.01, 0.09), "pr_auc": _estimate(0.08, 0.02, 0.15)},
                },
                "holdout": {"roc_auc": _estimate(0.03, -0.02, 0.08), "positive_groups": 6},
                "holm_p": 0.02,
                "decision": "not evaluable: one group holds > 40% of positives under ['strict']",
                "unsupported_groupings": ["strict"],
                "not_evaluable_runs": ["burial__full__strict__r0"],
            },
        },
        "permutation": {"burial": {"roc_auc": _estimate(0.5, 0.43, 0.57), "passes": True}},
        "not_evaluable": {"burial__full__strict__r0": "no grouped split gives every training fold both classes"},
    }


def test_page_is_self_contained_and_complete(report, tmp_path):
    out = tmp_path / "report.html"
    json_path = tmp_path / "benchmark_report.json"
    json_path.write_text(json.dumps(report))
    assert page.main(["--report-json", str(json_path), "--output", str(out)]) == 0
    text = out.read_text()

    # No external assets: it must render from a CI artifact, offline.
    assert not re.search(r'(src|href)\s*=\s*"https?://', text)
    assert "<svg" in text and text.count("<svg") >= 4
    # Every decision, and the honest ones especially, reach the page.
    assert "refuted" in text and "not evaluable" in text
    assert "burial__full__strict__r0" in text
    assert "no grouped split" in text
    # Permuted evaluations are a control, not a result row.
    assert "/permuted" not in text.split("Leak control")[0]
    assert "51,234" in text or "51234" in text


def test_intervals_are_coloured_by_whether_they_exclude_the_reference(report, tmp_path):
    out = tmp_path / "r.html"
    json_path = tmp_path / "r.json"
    json_path.write_text(json.dumps(report))
    page.main(["--report-json", str(json_path), "--output", str(out)])
    text = out.read_text()
    assert 'class="ci hit"' in text and 'class="ci null"' in text


def test_missing_estimates_say_so_rather_than_plotting_zero():
    forest = page.Forest([("absent", {}, "")], null=0.0)
    assert "not measured" in forest.svg()
    assert "circle" not in forest.svg()


def test_smoke_run_is_labelled_as_no_result(report, tmp_path):
    report["smoke"] = True
    json_path, out = tmp_path / "s.json", tmp_path / "s.html"
    json_path.write_text(json.dumps(report))
    page.main(["--report-json", str(json_path), "--output", str(out)])
    assert "Nothing on this page is a result" in out.read_text()


def test_ticks_and_scale_cover_the_data():
    forest = page.Forest([("a", _estimate(0.9, 0.85, 0.95), ""), ("b", _estimate(0.6, 0.5, 0.7), "")], null=0.5)
    assert forest.lo < 0.5 and forest.hi > 0.95
    assert 2 <= len(forest._ticks()) <= 6
