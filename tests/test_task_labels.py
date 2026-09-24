"""Every benchmark task's labels derive from one per-pocket record (docs/ANALYSIS_PLAN.md, section 4)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.analysis.labeling import LigandSite, assign_pocket_labels, task_label, task_labels

CASES = [
    # site label, burial class of the matched copy, expected (ip_site, cryptic_ip_site, burial)
    (1, "cryptic", (1, 1, 1)),
    (1, "surface", (1, -1, 0)),
    (1, "semi_cryptic", (1, -1, -1)),
    (1, "crystal_artifact", (-1, -1, -1)),  # a lattice contact is not a site
    (1, "unknown", (1, -1, -1)),
    (0, "", (0, 0, -1)),
    (-1, "cryptic", (-1, -1, -1)),  # ambiguous overlap stays excluded
]


@pytest.mark.parametrize("site_label,burial_class,expected", CASES)
def test_rules(site_label, burial_class, expected):
    got = tuple(task_label(task, site_label, burial_class) for task in ("ip_site", "cryptic_ip_site", "burial"))
    assert got == expected


def test_unknown_task_is_refused():
    with pytest.raises(ValueError):
        task_label("everything", 1, "cryptic")


def test_vectorised_matches_scalar():
    frame = pd.DataFrame({"label": [c[0] for c in CASES], "matched_burial_class": [c[1] for c in CASES]})
    for index, task in enumerate(("ip_site", "cryptic_ip_site", "burial")):
        assert task_labels(task, frame).tolist() == [c[2][index] for c in CASES]


def test_assignment_records_the_matched_copy_burial_class():
    ligand = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    sites = [LigandSite("IHP_A_1", "IHP", ligand, burial_class="cryptic")]
    pockets = [(1, np.zeros(3), ligand + 0.5), (2, np.array([30.0, 0, 0]), np.array([[30.0, 0, 0]]))]
    rows = [a.to_row() for a in assign_pocket_labels(pockets, sites)]
    assert rows[0]["label"] == 1 and rows[0]["matched_burial_class"] == "cryptic"
    assert rows[1]["label"] == 0
