"""Learned proteome ranking (docs/LEARNED_SCREEN_PLAN.md)."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from cryptic_ip.benchmark import protocol
from scripts import benchmark, learned_screen as ls
from tests.test_benchmark_script import _inputs


def test_parse_fasta_uses_accessions():
    text = ">sp|P12345|ABC_HUMAN Something\nMKT\nLLV\n>Q99999\nAAA\n"
    assert ls.parse_fasta(text) == {"P12345": "MKTLLV", "Q99999": "AAA"}


def test_seen_uses_the_benchmark_homology_criterion(tmp_path):
    rows = [
        # query target fident qstart qend tstart tend qlen tlen evalue
        ("A", "X", 0.45, 1, 90, 1, 90, 100, 300, 1e-20),   # 90 % of the shorter: seen
        ("B", "X", 0.25, 1, 90, 1, 90, 100, 300, 1e-20),   # identity too low
        ("C", "X", 0.90, 1, 30, 1, 30, 100, 300, 1e-20),   # coverage too low
        ("D", "X", 0.90, 1, 90, 1, 90, 100, 300, 1e-2),    # E-value too high
        ("E", "X", 45.0, 1, 90, 1, 90, 100, 300, 1e-20),   # identity in percent
    ]
    path = tmp_path / "hits.tsv"
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)
    assert ls.seen_proteins(path) == {"A", "E"}


@pytest.fixture(scope="module")
def bundle(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("learned")
    pockets, entries = _inputs(tmp, n_entries=40)
    table = tmp / "t.csv.gz"
    assert benchmark.main(["prepare", "--pockets-csv", str(pockets), "--entry-csv", str(entries),
                           "--output", str(table), "--summary-json", str(tmp / "s.json")]) == 0
    model = tmp / "model.joblib"
    assert ls.main(["train", "--table", str(table), "--output", str(model), "--n-draws", "1"]) == 0
    return model


def _screen(tmp_path, n_proteins=60, seed=0):
    """Screen shards where annotated binders carry the synthetic site signal (enclosure)."""
    rng = np.random.default_rng(seed)
    shards, catalogs = tmp_path / "shards", tmp_path / "catalog"
    shards.mkdir()
    catalogs.mkdir()
    rows, ann = [], []
    for i in range(n_proteins):
        acc = f"P{i:05d}"
        binder = i % 6 == 0
        for p in range(8):
            row = {name: rng.normal() for name in protocol.BENCHMARK_FEATURES}
            row["enclosure"] += 3.0 if (binder and p == 0) else 0.0
            row.update({"uniprot_id": acc, "pocket_id": p + 1, "plddt_mean": 90.0, "volume": 500.0,
                        "pocket_residues": "1,2,3"})
            rows.append(row)
        ann.append({"uniprot_id": acc, "gene": f"G{i}", "protein_name": f"protein {i}",
                    "binding_site": "inositol hexakisphosphate" if binder else "", "function": ""})
    pd.DataFrame(rows).to_csv(shards / "yeast_0_pockets_part000.csv.gz", index=False)
    pd.DataFrame({"uniprot_id": [f"P{i:05d}" for i in range(n_proteins)], "organism_key": "yeast",
                  "fragment": 1}).to_csv(catalogs / "yeast_catalog.csv", index=False)
    pd.DataFrame(ann).to_csv(catalogs / "yeast_uniprot.tsv", sep="\t", index=False)
    hits = tmp_path / "hits.tsv"
    hit = [("P00006", "X", 0.9, 1, 90, 1, 90, 100, 100, 1e-30)]
    pd.DataFrame(hit).to_csv(hits, sep="\t", header=False, index=False)
    return shards, catalogs, hits


def test_evaluate_excludes_seen_proteins_and_lists_candidates(bundle, tmp_path):
    shards, catalogs, hits = _screen(tmp_path)
    out = tmp_path / "out"
    assert ls.main(["evaluate", "--model", str(bundle), "--shards-dir", str(shards), "--catalog-dir", str(catalogs),
                    "--seen-hits", str(hits), "--output-dir", str(out), "--n-bootstrap", "100"]) == 0
    report = json.loads((out / "learned_screen.json").read_text())
    assert report["proteins_seen"] == 1
    pooled = report["pooled"]
    assert pooled["proteins"] == 59 and pooled["binders"] == 9  # P00006 (a binder) is excluded
    assert pooled["L1_learned_roc_auc"]["point"] > 0.7
    cand = pd.read_csv(out / "candidates.csv")
    assert not cand["annotated"].any() and not cand["uniprot_id"].eq("P00006").any()
    assert (out / "LEARNED_SCREEN.md").read_text().startswith("## Learned proteome ranking")


def test_low_confidence_pockets_are_ignored():
    pockets = pd.DataFrame({"uniprot_id": ["A", "A", "B"], "learned_score": [0.9, 0.2, 0.5],
                            "composite_score": [0.9, 0.3, 0.4], "plddt_mean": [40.0, 80.0, 90.0],
                            "pocket_id": [1, 2, 1]})
    table = ls.protein_table(pockets).set_index("uniprot_id")
    assert table.loc["A", "learned_score"] == 0.2 and table.loc["A", "rule_score"] == 0.3


def test_fasta_validator_rejects_non_fasta():
    from cryptic_ip.database.async_fetch import ValidationError

    ls._validate_fasta(b">P1\nMKT\n")
    ls._validate_fasta(b"")  # every accession in the batch obsolete
    with pytest.raises(ValidationError):
        ls._validate_fasta(b"<html>error</html>")


def test_uniprot_fetch_isolates_a_bad_accession_and_drops_non_accessions():
    from types import SimpleNamespace
    from urllib.parse import parse_qs, urlsplit

    bad = "Q0BAD1"
    calls = []

    def fake_fetch(jobs):
        out = []
        for job in jobs:
            query = parse_qs(urlsplit(job.url).query)["query"][0]
            accs = [t.split(":")[1] for t in query.split(" OR ")]
            calls.append(len(accs))
            if bad in accs:
                out.append(SimpleNamespace(ok=False, status=400, error="HTTP 400", payload=None))
            else:
                body = "".join(f">sp|{a}|X\nMKT\n" for a in accs).encode()
                out.append(SimpleNamespace(ok=True, status=200, error="", payload=body))
        return out

    accessions = [f"P{i:05d}" for i in range(20)] + [bad, "nan", "not-an-id"]
    records = ls.fetch_uniprot_fasta(accessions, batch=8, fetch=fake_fetch)
    assert set(records) == {f"P{i:05d}" for i in range(20)}
    assert max(calls) <= 8


def test_uniprot_fetch_stops_on_a_non_400_failure():
    from types import SimpleNamespace

    def failing(jobs):
        return [SimpleNamespace(ok=False, status=503, error="HTTP 503", payload=None) for _ in jobs]

    with pytest.raises(RuntimeError):
        ls.fetch_uniprot_fasta(["P12345"], fetch=failing)
