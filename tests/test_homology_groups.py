"""Homology grouping: what links two entries, and how groups form."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from cryptic_ip.benchmark import homology
from scripts import homology_groups


def _mmseqs(path, rows):
    # query target fident qstart qend tstart tend qlen tlen evalue
    path.write_text("".join("\t".join(str(v) for v in row) + "\n" for row in rows))
    return path


def test_entry_names_from_both_tools():
    assert homology.entry_of("1zy7|A") == "1ZY7"
    assert homology.entry_of("1ZY7.cif_B") == "1ZY7"
    assert homology.entry_of("1zy7.pdb.gz_A") == "1ZY7"


class TestSequenceLinks:
    def test_thresholds(self, tmp_path):
        hits = _mmseqs(tmp_path / "h.tsv", [
            ("AAAA|A", "BBBB|A", 0.35, 1, 100, 1, 100, 120, 400, 1e-20),   # 83 % of the shorter: link
            ("AAAA|A", "CCCC|A", 0.29, 1, 120, 1, 120, 120, 120, 1e-20),   # identity too low
            ("AAAA|A", "DDDD|A", 0.90, 1, 50, 1, 50, 120, 120, 1e-20),     # 42 % coverage
            ("AAAA|A", "EEEE|A", 0.90, 1, 120, 1, 120, 120, 120, 1e-2),    # E-value too high
            ("AAAA|A", "AAAA|B", 1.00, 1, 120, 1, 120, 120, 120, 0.0),     # same entry
        ])
        assert homology.sequence_links(hits) == {("AAAA", "BBBB")}

    def test_domain_shared_with_a_long_protein_links(self, tmp_path):
        """A PH domain inside a 1,000-residue protein: coverage is of the shorter chain."""
        hits = _mmseqs(tmp_path / "h.tsv", [("LONG|A", "PHPH|A", 0.4, 400, 510, 1, 110, 1000, 115, 1e-10)])
        assert homology.sequence_links(hits) == {("LONG", "PHPH")}

    def test_percent_identity_is_accepted(self, tmp_path):
        hits = _mmseqs(tmp_path / "h.tsv", [("AAAA|A", "BBBB|A", 35.0, 1, 100, 1, 100, 100, 100, 1e-9)])
        assert homology.sequence_links(hits) == {("AAAA", "BBBB")}

    def test_malformed_table_raises(self, tmp_path):
        bad = tmp_path / "h.tsv"
        bad.write_text("AAAA|A\tBBBB|A\t0.5\n")
        with pytest.raises(ValueError):
            homology.sequence_links(bad)


def test_structure_links(tmp_path):
    hits = tmp_path / "f.tsv"
    hits.write_text(
        "AAAA.cif_A\tBBBB.cif_A\t0.62\t1\t100\t1\t100\t110\t300\t1e-8\n"
        "AAAA.cif_A\tCCCC.cif_A\t0.41\t1\t100\t1\t100\t110\t110\t1e-8\n"
    )
    assert homology.structure_links(hits) == {("AAAA", "BBBB")}


class TestGroups:
    def test_links_are_transitive(self):
        groups = homology.connected_groups(["C", "A", "B", "D"], {("A", "B"), ("B", "C")})
        assert groups["A"] == groups["B"] == groups["C"] == "G:A"
        assert groups["D"] == "G:D"

    def test_names_do_not_depend_on_link_order(self):
        one = homology.connected_groups(list("ABCD"), [("C", "D"), ("B", "C"), ("A", "B")])
        two = homology.connected_groups(list("DCBA"), [("A", "B"), ("B", "C"), ("C", "D")])
        assert one == two and set(one.values()) == {"G:A"}

    def test_unknown_entries_in_links_are_ignored(self):
        assert homology.connected_groups(["A"], [("A", "Z")]) == {"A": "G:A"}

    def test_report(self):
        report = homology.GroupReport.of({"A": "G:A", "B": "G:A", "C": "G:C"})
        assert (report.n_groups, report.largest_group, report.n_singletons) == (2, 2, 1)


def test_chain_sequences_from_a_structure(tmp_path):
    from cryptic_ip.testing.synthetic import default_benchmark_specs, write_synthetic_structure

    spec = default_benchmark_specs(n_buried=1, n_surface=0, n_decoy=0, seed=0)[0]
    sequences = homology.chain_sequences(write_synthetic_structure(spec, tmp_path))
    assert sequences
    assert all(len(seq) >= homology.MIN_CHAIN_LENGTH for seq in sequences.values())
    assert all(set(seq) <= set("ACDEFGHIKLMNPQRSTVWYXUO") for seq in sequences.values())


def test_script_writes_both_groupings(tmp_path):
    entries = tmp_path / "entries.csv"
    pd.DataFrame({"pdb_id": ["1AAA", "2BBB", "3CCC", "4DDD"]}).to_csv(entries, index=False)
    structures = tmp_path / "structures"
    structures.mkdir()
    mm = _mmseqs(tmp_path / "mm.tsv", [("1AAA|A", "2BBB|A", 0.5, 1, 100, 1, 100, 100, 100, 1e-30)])
    fs = tmp_path / "fs.tsv"
    fs.write_text("2BBB.cif_A\t3CCC.cif_A\t0.7\t1\t90\t1\t90\t100\t100\t1e-9\n")
    out, report = tmp_path / "out.csv", tmp_path / "report.json"
    assert homology_groups.main([
        "--entry-csv", str(entries), "--structures-dir", str(structures),
        "--output-csv", str(out), "--report-json", str(report),
        "--mmseqs-hits", str(mm), "--foldseek-hits", str(fs), "--work-dir", str(tmp_path / "w"),
    ]) == 0
    table = pd.read_csv(out).set_index("pdb_id")
    assert table.loc["1AAA", "homology_group"] == table.loc["2BBB", "homology_group"] == "G:1AAA"
    assert table.loc["3CCC", "homology_group"] == "G:3CCC"  # structure link only
    assert set(table.loc[["1AAA", "2BBB", "3CCC"], "homology_group_strict"]) == {"G:1AAA"}
    assert table.loc["4DDD", "homology_group_strict"] == "G:4DDD"
    summary = json.loads(report.read_text())
    assert summary["sequence"]["n_groups"] == 3 and summary["strict"]["n_groups"] == 2
