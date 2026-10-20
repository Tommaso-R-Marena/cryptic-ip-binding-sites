"""Tests for study Q8, the within-protein control.

The design's whole claim is that the comparison holds protein, fold, model,
detector and lining size fixed and varies only which pocket. Three things have
to be true for that claim to survive, and each has a test here:

* the candidate must be re-identified by overlap, because pocket indices do not
  survive a re-run;
* siblings must be matched on lining size, or the test rediscovers that the
  selected pocket is bigger;
* an unresolved half must stay unresolved — the error that manufactured a
  spurious result in the paired within-lining study.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "within_protein", ROOT / "scripts" / "within_protein.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules["within_protein"] = script
SPEC.loader.exec_module(script)


# --------------------------------------------------------------- parsing

def test_residues_parse_from_the_committed_cell_format():
    assert script.parse_residues("1387,1390,1391") == [1387, 1390, 1391]
    assert script.parse_residues("[3, 7, 11]") == [3, 7, 11]
    assert script.parse_residues("") == []
    assert script.parse_residues(float("nan")) == []


def test_residues_are_sorted_and_deduplicated():
    assert script.parse_residues("9,3,9,5") == [3, 5, 9]


def test_jaccard_is_zero_against_nothing_and_one_against_itself():
    assert script.jaccard([1, 2, 3], []) == 0.0
    assert script.jaccard([1, 2, 3], [1, 2, 3]) == 1.0
    assert script.jaccard([1, 2, 3, 4], [3, 4, 5, 6]) == pytest.approx(2 / 6)


# ------------------------------------------------------- the three claims

def test_the_candidate_is_found_by_overlap_not_by_index():
    """A re-run renumbers pockets; the committed index must not be trusted."""
    committed = list(range(40, 60))
    pockets = pd.DataFrame([
        {"uniprot_id": "P1", "pocket_id": 1, "residues": "1,2,3,4,5,6,7,8,9,10"},
        # the candidate, renumbered and one residue different
        {"uniprot_id": "P1", "pocket_id": 77,
         "residues": ",".join(str(r) for r in committed[:-1] + [999])},
    ])
    table = pd.DataFrame([{"uniprot_id": "P1", "cluster": "c1",
                           "top_pocket_residues": ",".join(str(r) for r in committed)}])
    report = script.measure(pockets, table, n_draws=60, n_bootstrap=50)
    # Too few siblings to score, but the exclusion must be about siblings - not
    # about failing to find the candidate.
    assert any("siblings" in reason for reason in report["excluded"])
    assert not any("re-identified" in reason for reason in report["excluded"])


def test_a_candidate_that_cannot_be_re_identified_is_excluded_and_counted():
    pockets = pd.DataFrame([
        {"uniprot_id": "P1", "pocket_id": 1, "residues": "1,2,3,4,5,6"},
    ])
    table = pd.DataFrame([{"uniprot_id": "P1", "cluster": "c1",
                           "top_pocket_residues": "500,501,502,503,504,505"}])
    report = script.measure(pockets, table, n_draws=60, n_bootstrap=50)
    assert report["n_measured"] == 0
    assert any("re-identified" in reason for reason in report["excluded"])


def test_siblings_outside_the_size_tolerance_are_not_matched():
    rng = np.random.default_rng(0)
    siblings = [{"pocket_id": i, "size": size, "lining": list(range(size))}
                for i, size in enumerate([10, 11, 12, 13, 14, 30, 40])]
    matched = script.match_siblings(12, siblings, rng)
    assert sorted(s["size"] for s in matched) == [10, 11, 12, 13, 14]


def test_the_sibling_sample_is_capped_and_seeded():
    siblings = [{"pocket_id": i, "size": 12, "lining": list(range(12))}
                for i in range(script.MAX_SIBLINGS * 3)]
    first = script.match_siblings(12, siblings, np.random.default_rng(5))
    again = script.match_siblings(12, siblings, np.random.default_rng(5))
    assert len(first) == script.MAX_SIBLINGS
    assert [s["pocket_id"] for s in first] == [s["pocket_id"] for s in again]


def test_an_unresolved_half_is_not_scored_as_a_p_of_one():
    """The error that manufactured a result in the paired within-lining study."""
    contiguous = list(range(20, 34))          # every gap 1: one ordering
    out = script.architecture_p(contiguous, n_draws=60, seed=0)
    assert out["register"] is None
    assert out["size"] == len(contiguous)


def test_a_half_with_too_few_resolved_siblings_is_left_unscored():
    committed = [4, 8, 11, 15, 18, 22, 25, 29, 32, 36, 39, 43]
    rows = [{"uniprot_id": "P1", "pocket_id": 1,
             "residues": ",".join(str(r) for r in committed)}]
    # Contiguous siblings: register is unresolved for every one of them.
    for index in range(script.MIN_SIBLINGS + 2):
        start = 200 + 40 * index
        rows.append({"uniprot_id": "P1", "pocket_id": 10 + index,
                     "residues": ",".join(str(r) for r in range(start, start + 12))})
    pockets = pd.DataFrame(rows)
    table = pd.DataFrame([{"uniprot_id": "P1", "cluster": "c1",
                           "top_pocket_residues": ",".join(str(r) for r in committed)}])
    report = script.measure(pockets, table, n_draws=60, n_bootstrap=50)
    assert report["n_measured"] == 1
    protein = report["per_protein"][0]
    assert "register_quantile" not in protein
    register = next(t for t in report["tests"] if t["half"] == "register")
    assert register["runnable"] is False


# ------------------------------------------------------------- quantiles

def test_the_quantile_counts_siblings_that_are_more_structured():
    """Smaller p is more structure, so below 0.5 means the candidate leads."""
    quantile, _ = script.mid_p_quantile(0.01, [0.5, 0.6, 0.7, 0.8])
    assert quantile == 0.0
    quantile, _ = script.mid_p_quantile(0.9, [0.1, 0.2, 0.3, 0.4])
    assert quantile == 1.0


def test_ties_are_split_rather_than_counted_as_below():
    quantile, tied = script.mid_p_quantile(0.5, [0.5, 0.5, 0.5, 0.5])
    assert quantile == 0.5
    assert tied == 1.0


def test_an_empty_sibling_set_gives_no_quantile():
    quantile, tied = script.mid_p_quantile(0.5, [])
    assert np.isnan(quantile) and np.isnan(tied)


# ------------------------------------------------------- the whole report

def _synthetic(n_proteins: int = 6, n_siblings: int = 12):
    """A table and a pocket set with no planted difference between the arms."""
    rng = np.random.default_rng(11)
    rows, table = [], []
    for index in range(n_proteins):
        accession = f"P{index}"
        start = 50
        lining = sorted(rng.choice(np.arange(start, start + 120), size=14, replace=False))
        rows.append({"uniprot_id": accession, "pocket_id": 1,
                     "residues": ",".join(str(int(r)) for r in lining)})
        for sibling in range(n_siblings):
            other = sorted(rng.choice(np.arange(500 + 200 * sibling, 620 + 200 * sibling),
                                      size=14, replace=False))
            rows.append({"uniprot_id": accession, "pocket_id": 2 + sibling,
                         "residues": ",".join(str(int(r)) for r in other)})
        table.append({"uniprot_id": accession, "cluster": f"c{index}",
                      "top_pocket_residues": ",".join(str(int(r)) for r in lining)})
    return pd.DataFrame(rows), pd.DataFrame(table)


def test_the_report_measures_every_protein_and_states_its_rules():
    pockets, table = _synthetic()
    report = script.measure(pockets, table, n_draws=120, n_bootstrap=200)
    assert report["n_measured"] == 6
    assert report["min_jaccard"] == script.MIN_JACCARD
    assert report["size_tolerance"] == script.SIZE_TOLERANCE
    for entry in report["tests"]:
        if entry.get("runnable"):
            assert 0.0 <= entry["mean_quantile"] <= 1.0
            assert entry["n_groups"] == 6


def test_a_null_design_does_not_manufacture_a_finding():
    """Candidate and siblings are drawn the same way, so nothing should fire."""
    pockets, table = _synthetic()
    report = script.measure(pockets, table, n_draws=120, n_bootstrap=200)
    for entry in report["tests"]:
        if entry.get("runnable"):
            assert entry["is_finding"] is False


def test_the_markdown_states_the_control_and_the_ceiling():
    pockets, table = _synthetic()
    text = script.markdown(script.measure(pockets, table, n_draws=120, n_bootstrap=200))
    assert "nearly tautological" in text
    assert "a third" in text
    assert "not a negative" in text
    assert "below 0.5" in text.lower()


def test_the_markdown_reports_the_residual_size_difference():
    pockets, table = _synthetic()
    text = script.markdown(script.measure(pockets, table, n_draws=120, n_bootstrap=200))
    assert "Residual size" in text
    assert "confounded" in text


# ------------------------------------------- the detection half, offline

class _FakeAnalyzer:
    """Stands in for ProteinAnalyzer: the detector is not what is under test."""

    def __init__(self, pockets, residues, fail: set = frozenset()):
        self._pockets = pockets
        self._residues = residues
        self._fail = set(fail)

    def detect_pockets(self):
        return pd.DataFrame({"pocket_id": self._pockets})

    def get_pocket_residues(self, pocket_id, distance_cutoff: float = 8.0):
        if int(pocket_id) in self._fail:
            raise RuntimeError("no alpha spheres")
        return self._residues[int(pocket_id)]


def test_detection_records_every_pocket_with_its_lining(tmp_path):
    record = script.detect_one(
        "P1", tmp_path,
        fetch=lambda accession, cache: tmp_path / "AF-P1-F1-model_v4.pdb",
        analyzer_factory=lambda path, work, timeout: _FakeAnalyzer(
            [1, 2], {1: [5, 9, 12], 2: [40, 44, 47]}),
    )
    assert record["error"] == ""
    assert record["n_pockets"] == 2
    assert record["model"] == "AF-P1-F1-model_v4.pdb"
    assert record["pockets"][0]["residues"] == "5,9,12"


def test_one_unreadable_pocket_does_not_lose_the_protein(tmp_path):
    record = script.detect_one(
        "P1", tmp_path,
        fetch=lambda accession, cache: tmp_path / "model.pdb",
        analyzer_factory=lambda path, work, timeout: _FakeAnalyzer(
            [1, 2], {2: [40, 44, 47]}, fail={1}),
    )
    assert record["n_pockets"] == 2
    failed = [p for p in record["pockets"] if p["error"]]
    assert len(failed) == 1
    assert failed[0]["residues"] == ""


def test_a_model_that_cannot_be_fetched_is_recorded_not_raised(tmp_path):
    def refuse(accession, cache):
        raise ConnectionError("AlphaFold unreachable")

    record = script.detect_one("P1", tmp_path, fetch=refuse)
    assert record["n_pockets"] == 0
    assert "ConnectionError" in record["error"]


def test_a_detector_that_fails_is_recorded_not_raised(tmp_path):
    def explode(path, work, timeout):
        raise RuntimeError("fpocket exceeded 900 s")

    record = script.detect_one(
        "P1", tmp_path,
        fetch=lambda accession, cache: tmp_path / "model.pdb",
        analyzer_factory=explode,
    )
    assert "fpocket" in record["error"]


# ------------------------------------------------- the ground-truth arm

class _Arrays:
    """A two-chain structure with a ligand copy, as parallel arrays."""

    def __init__(self, shared_numbering: bool = False):
        coords, elements, resnames, chains, resseqs, icodes, polymer = [], [], [], [], [], [], []
        for index in range(10):                       # chain A, residues 1-10
            coords.append([float(index), 0.0, 0.0])
            elements.append("C")
            resnames.append("ALA")
            chains.append("A")
            resseqs.append(index + 1)
            icodes.append("")
            polymer.append(True)
        for index in range(10):                       # chain B, far away
            coords.append([float(index), 100.0, 0.0])
            elements.append("C")
            resnames.append("ALA")
            chains.append("B")
            resseqs.append(index + 1 if shared_numbering else index + 101)
            icodes.append("")
            polymer.append(True)
        for offset in range(3):                       # the ligand, beside chain A
            coords.append([float(offset) + 1.0, 2.0, 0.0])
            elements.append("P")
            resnames.append("IHP")
            chains.append("A")
            resseqs.append(900)
            icodes.append("")
            polymer.append(False)
        self.coords = np.asarray(coords, dtype=float)
        self.elements = np.asarray(elements)
        self.resnames = np.asarray(resnames)
        self.chain_ids = np.asarray(chains)
        self.resseqs = np.asarray(resseqs)
        self.icodes = np.asarray(icodes)
        self.is_polymer = np.asarray(polymer)


def test_the_annotated_site_is_the_ligand_s_contact_shell_by_chain():
    out = script.contact_residues(_Arrays(), "A", 900, "")
    assert out["found"] is True
    assert set(out["by_chain"]) == {"A"}
    # The ligand sits 2 A off chain A at x = 1, 2, 3, so every residue within
    # sqrt(4.5^2 - 2^2) = 4.03 along x is in contact: residues 1 through 8.
    assert out["by_chain"]["A"] == [1, 2, 3, 4, 5, 6, 7, 8]


def test_a_ligand_copy_that_is_not_there_is_reported_not_guessed():
    out = script.contact_residues(_Arrays(), "A", 999, "")
    assert out["found"] is False


def test_shared_chain_numbering_is_detected():
    assert script.chain_numbering_is_ambiguous(_Arrays(shared_numbering=True), "A") is True
    assert script.chain_numbering_is_ambiguous(_Arrays(), "A") is False


def test_recall_is_coverage_of_the_site_not_a_symmetric_overlap():
    """A pocket may legitimately be larger than a ligand's contact shell."""
    assert script.recall(list(range(1, 101)), [1, 2, 3, 4]) == 1.0
    assert script.jaccard(list(range(1, 101)), [1, 2, 3, 4]) < 0.05
    assert script.recall([1, 2], [1, 2, 3, 4]) == 0.5


def test_each_arm_matches_with_its_own_rule_and_floor():
    assert script.MATCH_RULES["jaccard"][1] == script.MIN_JACCARD
    assert script.MATCH_RULES["recall"][1] == script.MIN_RECALL
    assert set(script.ARMS) == {"candidate", "site"}


def test_an_unknown_rule_is_refused():
    pockets, table = _synthetic(n_proteins=1, n_siblings=10)
    with pytest.raises(SystemExit):
        script.measure(pockets, table, n_draws=60, n_bootstrap=50, rule="whatever")


def test_the_site_arm_reads_its_own_column_names():
    frame = pd.DataFrame([{"copy": "1ZY7_A_901", "residues": "3,7,11",
                           "group": "G:1ZY7", "burial_class": "cryptic"}])
    out = script.as_table(frame, "site")
    assert list(out.columns)[:3] == ["uniprot_id", "top_pocket_residues", "cluster"]
    assert out.iloc[0]["cluster"] == "G:1ZY7"


def test_the_recall_rule_keeps_a_pocket_the_jaccard_rule_would_reject():
    """The arms differ only here, and this is why."""
    site = [10, 11, 12, 13]
    lining = list(range(10, 40, 2))            # covers the site, much larger
    pockets = pd.DataFrame([{"uniprot_id": "c1", "pocket_id": 1,
                             "residues": ",".join(str(r) for r in lining)}])
    table = pd.DataFrame([{"uniprot_id": "c1", "cluster": "g1",
                           "top_pocket_residues": ",".join(str(r) for r in site)}])
    by_jaccard = script.measure(pockets, table, n_draws=60, n_bootstrap=50, rule="jaccard")
    by_recall = script.measure(pockets, table, n_draws=60, n_bootstrap=50, rule="recall")
    assert any("not re-identified" in r for r in by_jaccard["excluded"])
    assert not any("not re-identified" in r for r in by_recall["excluded"])


def test_the_site_arm_runs_the_same_statistic_as_the_candidate_arm():
    pockets, table = _synthetic(n_proteins=5, n_siblings=12)
    as_site = table.rename(columns={"uniprot_id": "copy",
                                    "top_pocket_residues": "residues",
                                    "cluster": "group"})
    candidate = script.measure(pockets, table, n_draws=120, n_bootstrap=200)
    site = script.measure(pockets, script.as_table(as_site, "site"),
                          n_draws=120, n_bootstrap=200, rule="recall")
    assert candidate["n_measured"] == site["n_measured"]
    for left, right in zip(candidate["tests"], site["tests"]):
        if left.get("runnable") and right.get("runnable"):
            assert left["mean_quantile"] == pytest.approx(right["mean_quantile"])


def test_the_burial_split_is_reported_and_never_promoted():
    """Declared secondary: a stratum under five groups is not evidence."""
    pockets, table = _synthetic(n_proteins=6, n_siblings=12)
    table = table.copy()
    table["burial_class"] = ["cryptic", "cryptic", "surface", "surface",
                             "surface", "surface"]
    report = script.measure(pockets, table, n_draws=120, n_bootstrap=200)
    strata = report["strata"]
    assert set(strata) == {"cryptic", "surface"}
    # Each protein is its own cluster here, so neither stratum reaches five
    # groups except the larger one.
    assert strata["cryptic"]["is_evidence"] is False
    assert strata["cryptic"]["n_groups"] == 2
    text = script.markdown(report)
    assert "By burial class" in text
    assert "too few groups" in text


def test_a_table_without_a_burial_class_reports_no_strata():
    pockets, table = _synthetic(n_proteins=3, n_siblings=10)
    report = script.measure(pockets, table, n_draws=120, n_bootstrap=200)
    assert report["strata"] == {}
    assert "By burial class" not in script.markdown(report)


def test_the_census_filter_is_the_one_the_amendment_names():
    """Pinned against the committed census: 137 copies, 135 entries, 17 groups."""
    census = pd.read_csv(ROOT / "results" / "redocking" / "census.csv")
    copies = script.select_copies(census)
    assert len(copies) == 137
    assert copies["pdb_id"].nunique() == 135
    assert copies["homology_group_strict"].nunique() == 17
    assert set(copies["comp_id"].astype(str).str.upper()) == {"IHP"}


def test_the_filter_drops_symmetry_contacts_and_keeps_unknown_ones():
    census = pd.DataFrame([
        {"comp_id": "IHP", "status": "eligible", "complete": True, "selected": True,
         "symmetry_contact": True, "pdb_id": "1AAA"},
        {"comp_id": "IHP", "status": "eligible", "complete": True, "selected": True,
         "symmetry_contact": float("nan"), "pdb_id": "2BBB"},
        {"comp_id": "IHP", "status": "excluded: crystal artefact", "complete": True,
         "selected": True, "symmetry_contact": False, "pdb_id": "3CCC"},
    ])
    kept = script.select_copies(census)
    assert list(kept["pdb_id"]) == ["2BBB"]


def test_a_power_failure_is_not_reported_as_a_null():
    """Five groups is the floor; below it the arm decides nothing, loudly."""
    starved = {"half": "topology", "runnable": True, "n_groups": 4, "n_proteins": 5,
               "is_evidence": False, "is_finding": False, "minimum_detectable": 0.37}
    assert "no evidence" in script._verdict(starved)
    assert script._verdict({"half": "topology", "runnable": True, "n_groups": 30,
                            "is_evidence": True, "is_finding": False}) == "no"
    assert script._verdict({"is_finding": True, "is_evidence": True}) == "**yes**"


def test_the_committed_ground_truth_arm_says_it_decided_nothing():
    """Pinned: the arm as run is underpowered and the report must say so."""
    path = ROOT / "results" / "within_protein" / "within_protein_sites.json"
    if not path.exists():
        pytest.skip("the ground-truth arm has not been run")
    report = json.loads(path.read_text())
    runnable = [t for t in report["tests"] if t.get("runnable")]
    assert runnable, "the arm produced no runnable test"
    assert all(not t["is_evidence"] for t in runnable)
    text = script.markdown(report)
    assert "did not have the power to decide anything" in text
    assert "neither supports nor refutes" in text
