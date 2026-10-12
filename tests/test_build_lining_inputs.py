"""Tests for the inputs builder, which exists mostly to refuse things.

The two refusals are the point: a chain-blind residue number must not become a
window on a multi-chain source, and a sequence-inputs file must never be
written with letters nobody supplied. Both are checked here, as is the schema
of the committed table the builder claims to know.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "build_lining_inputs", ROOT / "scripts" / "build_lining_inputs.py")
script = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
# Registered before exec: a module-level @dataclass resolves its own module
# through sys.modules, and without this the class body raises on import.
sys.modules["build_lining_inputs"] = script
SPEC.loader.exec_module(script)


def table(**overrides):
    base = {
        "uniprot_id": ["P1", "P2", "P3"],
        "top_pocket_residues": ["10,14,17,21,24", "5,8,12,15", "3,7"],
        "cluster": ["C1", "C2", "C1"],
        "explained": [True, False, True],
        "plddt_mean": [88.0, 91.0, 72.0],
    }
    base.update(overrides)
    return pd.DataFrame(base)


COLUMNS = script.Columns(identifier="uniprot_id", residues="top_pocket_residues",
                         group="cluster", label="explained",
                         covariates=("plddt_mean",))


# ------------------------------------------------------------- parsing

def test_residues_parse_from_the_format_the_repo_writes():
    assert script.parse_residues("10,14,17") == [10, 14, 17]
    assert script.parse_residues("10 14 17") == [10, 14, 17]
    assert script.parse_residues("10;14;17") == [10, 14, 17]


def test_residues_are_sorted_and_deduplicated():
    assert script.parse_residues("17,10,14,10") == [10, 14, 17]


def test_an_empty_cell_is_no_residues_not_an_error():
    assert script.parse_residues("") == []
    assert script.parse_residues(None) == []


def test_an_unreadable_residue_number_stops_the_run():
    """Silently dropping it would shrink the study without saying so."""
    with pytest.raises(script.InputError) as caught:
        script.parse_residues("10,14A,17")
    assert "14A" in str(caught.value)


# ------------------------------------------------------------- refusals

def test_a_multi_chain_source_without_a_chain_column_is_refused():
    """The repo's residue numbers are chain-blind by their own docstring."""
    with pytest.raises(script.InputError) as caught:
        script.build_records(table(), COLUMNS, source="t", single_chain=False)
    assert "chain-blind" in str(caught.value)


def test_a_multi_chain_source_with_a_chain_column_is_allowed():
    frame = table(chain=["A", "B", "A"])
    columns = script.Columns(identifier="uniprot_id", residues="top_pocket_residues",
                             group="cluster", chain="chain")
    out = script.build_records(frame, columns, source="t", single_chain=False)
    assert out["n_records"] == 3
    assert {r["chain"] for r in out["records"]} == {"A", "B"}
    # Each record carries the chain its positions were read on, so a window
    # built downstream can never mix two of them.
    assert all(len({c for c, _p in r["residues"]}) == 1 for r in out["records"])
    assert all(r["chain_is_implicit"] is False for r in out["records"])


def test_a_single_chain_source_records_that_its_chain_is_implicit():
    out = script.build_records(table(), COLUMNS, source="t", single_chain=True)
    assert all(r["chain"] == script.IMPLICIT_CHAIN for r in out["records"])
    assert all(r["chain_is_implicit"] is True for r in out["records"])


def test_a_missing_column_stops_the_run():
    with pytest.raises(script.InputError) as caught:
        script.build_records(table().drop(columns=["cluster"]), COLUMNS,
                             source="t", single_chain=True)
    assert "cluster" in str(caught.value)


def test_rows_are_dropped_with_a_counted_reason():
    out = script.build_records(table(), COLUMNS, source="t", single_chain=True,
                               min_positions=5)
    assert out["n_records"] == 1
    assert sum(out["dropped"].values()) == 2


def test_a_row_without_a_group_is_dropped():
    frame = table(cluster=["C1", None, "C2"])
    out = script.build_records(frame, COLUMNS, source="t", single_chain=True)
    assert out["n_records"] == 2
    assert "no grouping key" in out["dropped"]


# ------------------------------------------------------------- sequences

def test_without_sequences_no_sequence_file_is_written(tmp_path):
    source = tmp_path / "t.csv"
    table().to_csv(source, index=False)
    out = tmp_path / "architecture.jsonl"
    assert script.main([
        "--table", str(source), "--identifier", "uniprot_id",
        "--residues", "top_pocket_residues", "--group", "cluster",
        "--single-chain", "--architecture-out", str(out),
    ]) == 0
    assert out.exists()
    assert not (tmp_path / "sequence_inputs.jsonl").exists()


def test_sequences_are_attached_when_supplied():
    payload = script.build_records(table(), COLUMNS, source="t", single_chain=True)
    out = script.attach_sequences(payload, {"P1": "A" * 40, "P2": "A" * 40, "P3": "A" * 40})
    assert out["n_records"] == 3
    assert all("sequence" in r for r in out["records"])


def test_a_record_whose_positions_fall_outside_its_sequence_is_dropped():
    """The signature of a numbering offset, and guessing it shifts a study."""
    payload = script.build_records(table(), COLUMNS, source="t", single_chain=True)
    out = script.attach_sequences(payload, {"P1": "A" * 5, "P2": "A" * 40, "P3": "A" * 40})
    assert out["n_records"] == 2
    assert "lining positions fall outside the sequence" in out["dropped"]


def test_a_record_without_a_sequence_is_dropped_not_padded():
    payload = script.build_records(table(), COLUMNS, source="t", single_chain=True)
    out = script.attach_sequences(payload, {"P1": "A" * 40})
    assert out["n_records"] == 1
    assert "no sequence for this identifier" in out["dropped"]


def test_fasta_reads_by_first_token(tmp_path):
    path = tmp_path / "s.fasta"
    path.write_text(">P1 some description\nACDE\nFGHI\n>P2\nKLMN\n")
    assert script.read_fasta(path) == {"P1": "ACDEFGHI", "P2": "KLMN"}


# ------------------------------------------------- the committed table

def test_the_known_table_still_has_the_columns_the_builder_names():
    """If the table's schema moves, this fails instead of the study."""
    spec = script.KNOWN["specificity-candidates"]
    frame = pd.read_csv(ROOT / spec["path"])
    script.verify_columns(frame, spec["columns"])


def test_the_known_table_is_declared_single_chain_for_a_reason():
    spec = script.KNOWN["specificity-candidates"]
    assert spec["single_chain"] is True
    assert "one chain" in spec["note"]


def test_the_committed_inputs_match_the_committed_table():
    inputs = ROOT / "results" / "architecture" / "real" / "architecture_inputs.jsonl"
    if not inputs.exists():
        pytest.skip("inputs not built")
    records = [json.loads(line) for line in inputs.read_text().splitlines() if line.strip()]
    frame = pd.read_csv(ROOT / script.KNOWN["specificity-candidates"]["path"])
    assert len(records) == len(frame)
    assert {r["id"] for r in records} == set(frame["uniprot_id"].astype(str))
    for record in records:
        assert record["chain_is_implicit"] is True
        assert record["n_positions"] >= 2
        assert len(record["residues"]) == record["n_positions"]
