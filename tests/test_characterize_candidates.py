"""Tests for the candidate-dossier builder (pure helpers, no network/fpocket)."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "characterize_candidates", ROOT / "scripts" / "characterize_candidates.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


MOD = _load_module()

HIT = {
    "uniprot_id": "P07264",
    "composite_score": 0.9546,
    "classification": "High confidence cryptic IP site",
    "pocket_id": 2,
    "volume": 481.4,
    "burial_depth": 2.84,
    "sasa": 7.30,
    "basic_residues": 7,
    "plddt_confidence": 93.47,
    "structure_path": "data/structures/yeast_pilot/AF-P07264-F1-model_v6.pdb",
}

UNIPROT = {
    "gene_name": "LEU1",
    "protein_name": "3-isopropylmalate dehydratase",
    "organism": "Saccharomyces cerevisiae",
    "sequence_length": 779,
    "function": "Catalyzes the isomerization of 2-isopropylmalate to 3-isopropylmalate.",
    "subcellular_location": ["Cytoplasm"],
    "go_terms": {"molecular_function": ["F:aconitate hydratase activity"]},
}


def test_build_candidate_record_core_fields():
    rec = MOD.build_candidate_record(HIT, rank=1)
    assert rec["rank"] == 1
    assert rec["uniprot_id"] == "P07264"
    assert rec["composite_score"] == 0.9546
    assert rec["pocket_id"] == 2
    assert rec["burial_depth"] == 2.84
    assert rec["basic_residues"] == 7
    # No annotation/structure lookups were provided.
    assert "gene_name" not in rec
    assert "coordinating_residues" not in rec


def test_build_candidate_record_with_annotation_and_residues():
    rec = MOD.build_candidate_record(
        HIT, rank=1, uniprot_info=UNIPROT, coordinating_residues=["R123", "K130", "H140"]
    )
    assert rec["gene_name"] == "LEU1"
    assert rec["organism"] == "Saccharomyces cerevisiae"
    assert rec["coordinating_residues"] == ["R123", "K130", "H140"]
    assert rec["go_molecular_function"] == ["F:aconitate hydratase activity"]


def test_build_candidate_record_handles_missing_values():
    sparse = {"uniprot_id": "X", "composite_score": "nan", "pocket_id": None}
    rec = MOD.build_candidate_record(sparse, rank=3)
    assert rec["composite_score"] is None
    assert rec["pocket_id"] is None
    assert rec["burial_depth"] is None


def test_format_candidate_markdown_contains_key_facts():
    rec = MOD.build_candidate_record(
        HIT, rank=1, uniprot_info=UNIPROT, coordinating_residues=["R123", "K130"]
    )
    md = MOD.format_candidate_markdown(rec)
    assert "P07264" in md
    assert "3-isopropylmalate dehydratase" in md
    assert "LEU1" in md
    assert "0.955" in md  # composite score, 3 dp
    assert "R123, K130" in md
    assert MOD.ANGSTROM_CU in md  # volume unit rendered


def test_format_candidate_markdown_truncates_long_function():
    rec = MOD.build_candidate_record(HIT, rank=1, uniprot_info={"function": "x" * 600})
    md = MOD.format_candidate_markdown(rec)
    assert "…" in md
