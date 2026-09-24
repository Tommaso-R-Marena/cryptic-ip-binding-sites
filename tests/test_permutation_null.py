"""Diagnostic D1: the permutation-control null summary."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from scripts import permutation_null


def _write(directory, task, repeat, shift, seed):
    rng = np.random.default_rng(seed)
    n = 300
    label = (rng.random(n) < 0.2).astype(int)
    score = rng.normal(size=n) + shift * label
    pd.DataFrame({"split": "cv", "task": task, "repeat": repeat, "fold": np.arange(n) % 3,
                  "label": label, "score": score}).to_csv(
        directory / f"{task}__full__sequence__r{repeat}__permuted.csv.gz", index=False)
    (directory / f"{task}__full__sequence__r{repeat}__permuted.json").write_text(json.dumps({"task": task}))


def test_null_centred_at_half_is_chance_and_shifted_is_leak(tmp_path):
    for r in range(12):
        _write(tmp_path, "burial", r, 0.0, r)
        _write(tmp_path, "cryptic_ip_site", r, 1.0, 100 + r)
    (tmp_path / "burial__full__sequence__r99__permuted.json").write_text(
        json.dumps({"task": "burial", "repeat": 99, "not_evaluable": "no split"}))
    out = tmp_path / "null.json"
    assert permutation_null.main(["--predictions-dir", str(tmp_path), "--output", str(out)]) == 0
    report = json.loads(out.read_text())
    assert report["tasks"]["burial"]["verdict"] == "chance"
    assert report["tasks"]["cryptic_ip_site"]["verdict"] == "leak"
    assert report["tasks"]["burial"]["pooled"]["n"] == 12
    assert report["not_evaluable"] == [{"task": "burial", "repeat": 99, "reason": "no split"}]
    assert out.with_suffix(".md").read_text().startswith("## Permutation-control null")
