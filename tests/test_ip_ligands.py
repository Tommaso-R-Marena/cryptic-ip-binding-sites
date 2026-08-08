"""Tests for inositol phosphate ligand discovery and validation.

Discovery must be self-correcting: a mistyped or non-existent component
identifier has to be *detected and dropped*, never silently contribute nothing.
"""

import pytest

from cryptic_ip.database.ip_ligands import (
    IPLigand,
    build_ligand,
    comp_ids,
    discover_ip_ligands,
    ip_series_label,
    looks_like_inositol_phosphate,
    parse_formula,
    phosphorylated_comp_ids,
)


class FakeClient:
    """Minimal stand-in for the RCSB client, driven by a fixed component table."""

    def __init__(self, components, search_hits=None, fail_search=False):
        self.components = components
        self.search_hits = search_hits or {}
        self.fail_search = fail_search
        self.lookups = []

    def search_chemcomp_full_text(self, term):
        if self.fail_search:
            raise RuntimeError("search service unavailable")
        return list(self.search_hits.get(term, []))

    def fetch_chemcomp(self, comp_id):
        self.lookups.append(comp_id)
        payload = self.components.get(comp_id)
        return {"chem_comp": payload} if payload else None


def _component(name, formula, weight=660.0):
    return {"id": "X", "name": name, "formula": formula, "formula_weight": weight}


IHP = _component("INOSITOL HEXAKISPHOSPHATE", "C6 H18 O24 P6", 660.04)
I3P = _component("D-MYO-INOSITOL 1,4,5-TRISPHOSPHATE", "C6 H15 O15 P3", 420.10)
INS = _component("MYO-INOSITOL", "C6 H12 O6", 180.16)
PIP2 = _component("PHOSPHATIDYLINOSITOL 4,5-BISPHOSPHATE", "C41 H81 O19 P3", 1042.0)
GLUCOSE = _component("ALPHA-D-GLUCOSE", "C6 H12 O6", 180.16)


def test_parse_formula_handles_counts_and_charges():
    assert parse_formula("C6 H18 O24 P6") == {"C": 6, "H": 18, "O": 24, "P": 6}
    assert parse_formula("C6 H6 O24 P6 12-") == {"C": 6, "H": 6, "O": 24, "P": 6}
    assert parse_formula("N") == {"N": 1}
    assert parse_formula("") == {}
    assert parse_formula(None) == {}


def test_series_label_is_derived_from_phosphorus_count():
    assert ip_series_label(6) == "InsP6"
    assert ip_series_label(0) == "InsP0"
    # Inositol pyrophosphates carry more than six phosphorus atoms.
    assert ip_series_label(8) == "InsP8"
    with pytest.raises(ValueError):
        ip_series_label(-1)


def test_inositol_phosphate_recognition_requires_name_and_formula():
    assert looks_like_inositol_phosphate(IHP["name"], IHP["formula"])
    # Right name, no phosphorus.
    assert not looks_like_inositol_phosphate(INS["name"], INS["formula"])
    # Right formula shape, wrong molecule.
    assert not looks_like_inositol_phosphate(GLUCOSE["name"], GLUCOSE["formula"])
    # Lipid-linked species are excluded unless explicitly allowed.
    assert not looks_like_inositol_phosphate(PIP2["name"], PIP2["formula"])
    assert looks_like_inositol_phosphate(PIP2["name"], PIP2["formula"], allow_lipid=True)


def test_build_ligand_extracts_series_and_counts():
    ligand = build_ligand("IHP", {"chem_comp": IHP}, sources=["seed"])
    assert isinstance(ligand, IPLigand)
    assert ligand.n_phosphorus == 6
    assert ligand.series == "InsP6"
    assert ligand.n_heavy_atoms == 36
    assert ligand.is_phosphorylated
    assert not ligand.is_lipid_linked
    assert ligand.verified


def test_build_ligand_rejects_non_inositol_components():
    assert build_ligand("GLC", {"chem_comp": GLUCOSE}) is None


def test_discovery_accepts_valid_and_records_rejected_identifiers():
    client = FakeClient(
        components={"IHP": IHP, "I3P": I3P, "INS": INS, "PIP2": PIP2},
        search_hits={"inositol phosphate": ["IHP", "I3P", "PIP2"]},
    )
    ligands, provenance = discover_ip_ligands(
        client, search_terms=["inositol phosphate"], seed_comp_ids=["INS", "NOPE"]
    )

    assert comp_ids(ligands) == ("I3P", "IHP")
    assert phosphorylated_comp_ids(ligands) == ("I3P", "IHP")
    # Unknown, lipid-linked and unphosphorylated candidates are dropped with a reason.
    assert provenance["rejected"]["NOPE"] == "not_found"
    assert provenance["rejected"]["PIP2"] == "lipid_linked"
    assert provenance["rejected"]["INS"] == "unphosphorylated"
    assert provenance["n_accepted"] == 2


def test_discovery_can_include_free_inositol_as_a_reference():
    client = FakeClient(components={"IHP": IHP, "INS": INS})
    ligands, _ = discover_ip_ligands(
        client,
        search_terms=[],
        seed_comp_ids=["IHP", "INS"],
        include_unphosphorylated=True,
    )
    assert comp_ids(ligands) == ("IHP", "INS")
    assert phosphorylated_comp_ids(ligands) == ("IHP",)


def test_discovery_survives_a_failing_search_service():
    """A search outage must degrade to the seed list, not abort collection."""
    client = FakeClient(components={"IHP": IHP}, fail_search=True)
    ligands, provenance = discover_ip_ligands(
        client, search_terms=["inositol"], seed_comp_ids=["IHP"]
    )
    assert comp_ids(ligands) == ("IHP",)
    assert provenance["n_candidates"] == 1


def test_ligands_are_sorted_by_phosphorylation_state():
    client = FakeClient(components={"IHP": IHP, "I3P": I3P})
    ligands, _ = discover_ip_ligands(client, search_terms=[], seed_comp_ids=["I3P", "IHP"])
    assert [lig.comp_id for lig in ligands] == ["IHP", "I3P"]
