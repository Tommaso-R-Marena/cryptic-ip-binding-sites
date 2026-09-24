#!/usr/bin/env python3
"""The α-arrestin lead (docs/ARRESTIN_PLAN.md).

``family``      B1: arrestin-fold family lists from UniProt Pfam cross-references,
                and the learned scores' family-level test
``references``  classic arrestin-IP complexes among the benchmark entries and their
                IP-contacting residues
``map``         B2/B3: US-align of AlphaFold targets onto the references, mapped sites,
                overlap with the top learned pocket
``conserve``    B4: UniRef50 homologues, MAFFT, basic-residue conservation
``tasks``       the docking tasks of B5 (sites, controls, decoy)
``dock``        one shard of docking tasks
``decide``      B6: the four criteria per protein, dossiers, report
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import shutil
import subprocess
import sys
import urllib.parse
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Set, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

LOGGER = logging.getLogger("arrestin")

ORGANISMS = {"human": "9606", "yeast": "559292", "dictyostelium": "44689"}
ARRESTIN_PFAM = ("PF00339", "PF02752")
LEADS = {"Q8TBH0": "ARRDC2", "P53244": "ART5"}
CLASSIC_HUMAN = {"P10523": "SAG", "P36575": "ARR3", "P49407": "ARRB1", "P32121": "ARRB2"}
AF_CONTROLS = ("P49407", "P32121", "P10523")  # plan B5 (b): ARRB1, ARRB2, SAG
SITE_DISTANCE = 4.0
MIN_TM = 0.5
OVERLAP_JACCARD = 0.20
NEGATIVE_JACCARD = 0.10
N_NEGATIVES = 5
MIN_HOMOLOGUES = 10
MAX_HOMOLOGUES = 500
BASIC = set("KRH")
CONSERVED_FRACTION = 0.80
MIN_CONSERVED = 3
CONVERGENCE = 0.50
TOP_PER_SEED = 5
CLUSTER_RMSD = 2.0
SEED_B1 = 20260927
SEED_NEG = 20260928
UNIPROT_FIELDS = "accession,gene_primary,protein_name,xref_pfam,organism_id,length"


# ------------------------------------------------------------ UniProt
def _validate_tsv(payload: bytes) -> None:
    from cryptic_ip.database.async_fetch import ValidationError

    text = payload.lstrip()
    if text and not text.startswith(b"Entry"):
        raise ValidationError(f"not a UniProt TSV: {payload[:80]!r}")


def _fetch_text(url: str, key: str, validator=None, fetch=None) -> str:
    from cryptic_ip.database.async_fetch import FetchJob, fetch_all

    fetch = fetch or (lambda jobs: fetch_all(jobs, concurrency=2, per_host=2))
    result = fetch([FetchJob(key=key, url=url, validator=validator)])[0]
    if not result.ok:
        raise RuntimeError(f"{key}: {result.error}")
    return (result.payload or b"").decode()


def uniprot_tsv(query: str, fields: str = UNIPROT_FIELDS, fetch=None) -> pd.DataFrame:
    import io

    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"format": "tsv", "fields": fields, "query": query})
    text = _fetch_text(url, query[:60], _validate_tsv, fetch)
    if not text.strip():
        return pd.DataFrame(columns=["Entry"])
    return pd.read_csv(io.StringIO(text), sep="\t", dtype=str)


def recommended_name(protein_names: str) -> str:
    """UniProt's 'Protein names' field starts with the recommended name; alternatives follow in parentheses."""
    return re.split(r" \(", str(protein_names or ""), maxsplit=1)[0].strip()


def is_arrestin_fold(pfam: str) -> bool:
    return any(pf in str(pfam or "") for pf in ARRESTIN_PFAM)


def is_classic(protein_names: str) -> bool:
    name = recommended_name(protein_names).lower()
    return "arrestin" in name and "domain-containing" not in name and "trafficking adapter" not in name


def classify(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.rename(columns={"Entry": "uniprot_id", "Gene Names (primary)": "gene", "Protein names": "name",
                                "Pfam": "pfam", "Organism (ID)": "organism_id", "Length": "length"}).copy()
    out["arrestin_fold"] = out["pfam"].map(is_arrestin_fold)
    out["classic"] = out["arrestin_fold"] & out["name"].map(is_classic)
    out["alpha"] = out["arrestin_fold"] & ~out["classic"]
    return out


# --------------------------------------------------------------- B1
def family_test(proteins: pd.DataFrame, family: Set[str], n_bootstrap: int = 2000) -> Dict[str, object]:
    """One-sided: does the family's learned score rank above other unseen proteins?"""
    from scipy.stats import mannwhitneyu

    from cryptic_ip.benchmark import protocol

    unseen = proteins[~proteins["seen"].astype(str).str.lower().isin(["true", "1"])].copy()
    unseen = unseen[pd.to_numeric(unseen["learned_score"], errors="coerce").notna()]
    y = unseen["uniprot_id"].isin(family).to_numpy(dtype=int)
    out: Dict[str, object] = {"unseen_proteins": int(len(unseen)), "members_unseen": int(y.sum())}
    if y.sum() == 0:
        out["decision"] = "not evaluable: no unseen family member"
        return out
    clusters = unseen["cluster"].fillna(unseen["uniprot_id"]).astype(str).to_numpy()
    member_clusters = int(pd.Series(clusters[y == 1]).nunique())
    out["member_clusters"] = member_clusters
    scores = unseen["learned_score"].astype(float).to_numpy()
    ranked = protocol.Ranked(y, scores)
    point = ranked.roc_auc(np.ones(len(y)))
    boot = np.array([ranked.roc_auc(w) for w in protocol.group_bootstrap_weights(clusters, n_bootstrap, SEED_B1)])
    boot = boot[np.isfinite(boot)]
    out["roc_auc"] = {"point": float(point), "p5": float(np.percentile(boot, 5)),
                      "low": float(np.percentile(boot, 2.5)), "high": float(np.percentile(boot, 97.5)),
                      "one_sided_p": float(np.mean(boot <= 0.5)), "n_bootstrap": int(boot.size)}
    mw = mannwhitneyu(scores[y == 1], scores[y == 0], alternative="greater")
    out["mann_whitney_one_sided_p_independence_assumed"] = float(mw.pvalue)
    percentile = unseen["learned_score"].rank(pct=True)
    unseen["rank"] = unseen["learned_score"].rank(ascending=False, method="min").astype(int)
    members = unseen[y == 1].assign(rank_percentile=percentile[y == 1] * 100)
    out["members"] = members[["uniprot_id", "organism_key", "learned_score", "rank", "rank_percentile", "cluster"]
                             ].sort_values("learned_score", ascending=False).to_dict(orient="records")
    if member_clusters < 5:
        out["decision"] = f"not evaluable: {member_clusters} clusters (fewer than 5)"
    elif out["roc_auc"]["p5"] > 0.5:
        out["decision"] = "supported"
    else:
        out["decision"] = "not supported"
    return out


def cmd_family(args: argparse.Namespace) -> int:
    proteins = pd.read_csv(args.proteins, low_memory=False)
    frames = []
    for org, taxon in ORGANISMS.items():
        query = f"(organism_id:{taxon}) AND ((xref:pfam-{ARRESTIN_PFAM[0]}) OR (xref:pfam-{ARRESTIN_PFAM[1]}))"
        f = uniprot_tsv(query)
        f["organism_key"] = org
        frames.append(f)
    fam = classify(pd.concat(frames, ignore_index=True))
    missing_leads = [a for a in LEADS if a not in set(fam["uniprot_id"])]
    fam["in_screen"] = fam["uniprot_id"].isin(set(proteins["uniprot_id"]))
    scored = proteins.set_index("uniprot_id")
    for col in ("learned_score", "seen", "cluster", "top_pocket_residues", "top_pocket"):
        if col in scored:
            fam[col] = fam["uniprot_id"].map(scored[col])
    args.out_dir.mkdir(parents=True, exist_ok=True)
    fam.to_csv(args.out_dir / "family.csv", index=False)
    alpha = set(fam.loc[fam["alpha"], "uniprot_id"])
    result = {"plan": "docs/ARRESTIN_PLAN.md", "family_counts": {
        org: {"arrestin_fold": int(g["arrestin_fold"].sum()), "classic": int(g["classic"].sum()),
              "alpha": int(g["alpha"].sum()), "alpha_scored": int((g["alpha"] & g["in_screen"]).sum())}
        for org, g in fam.groupby("organism_key")},
        "leads_missing_from_definition": missing_leads,
        "B1": family_test(proteins, alpha, args.n_bootstrap)}
    per_org = {}
    for org in ORGANISMS:
        sub = proteins[proteins["organism_key"] == org] if "organism_key" in proteins else proteins.iloc[0:0]
        if len(sub):
            r = family_test(sub, alpha, args.n_bootstrap)
            per_org[org] = {k: r.get(k) for k in ("members_unseen", "member_clusters", "roc_auc", "decision")}
    result["B1_per_organism_descriptive"] = per_org
    classic = fam[fam["classic"]]
    ranks = proteins["learned_score"].rank(ascending=False, method="min")
    pct = proteins["learned_score"].rank(pct=True) * 100
    result["positive_control_classic"] = [
        {"uniprot_id": r.uniprot_id, "gene": r.gene, "organism": r.organism_key,
         "learned_score": None if pd.isna(r.learned_score) else float(r.learned_score),
         "seen": None if pd.isna(r.seen) else bool(str(r.seen).lower() == "true"),
         "rank_among_all_scored": int(ranks[proteins["uniprot_id"] == r.uniprot_id].iloc[0])
         if (proteins["uniprot_id"] == r.uniprot_id).any() else None,
         "percentile_among_all_scored": float(pct[proteins["uniprot_id"] == r.uniprot_id].iloc[0])
         if (proteins["uniprot_id"] == r.uniprot_id).any() else None}
        for r in classic.itertuples()]
    (args.out_dir / "family.json").write_text(json.dumps(result, indent=2, default=str))
    print(json.dumps({k: v for k, v in result.items() if k != "B1"}, indent=1, default=str)[:4000])
    print(json.dumps(result["B1"], indent=1, default=str)[:4000])
    return 0


# -------------------------------------------------------- references
def ca_residues(path: Path, chain: Optional[str] = None) -> Tuple[str, List[int], np.ndarray, str]:
    """(sequence, residue numbers, CA coords, chain) of one chain (first polymer chain by default)."""
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from redocking import chain_sequence

    arrays = load_structure_arrays(path)
    if chain is None:
        chain = str(arrays.chain_ids[np.flatnonzero(arrays.is_polymer)[0]])
    seq, numbers = chain_sequence(arrays, chain)
    mask = (arrays.chain_ids == chain) & (arrays.atom_names == "CA") & arrays.is_polymer
    return seq, numbers, arrays.coords[mask], chain


def write_chain_pdb(src: Path, chain: str, out: Path) -> None:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking.receptor import Atom, _pdb_atom_line

    arrays = load_structure_arrays(src)
    lines = []
    for i in np.flatnonzero((arrays.chain_ids == chain) & arrays.is_polymer):
        a = Atom("ATOM", str(arrays.atom_names[i]), str(arrays.resnames[i]), "A", int(arrays.resseqs[i]),
                 str(arrays.icodes[i]), arrays.coords[i], str(arrays.elements[i]))
        lines.append(_pdb_atom_line(len(lines) + 1, a))
    out.write_text("\n".join(lines + ["END"]) + "\n")


def cmd_references(args: argparse.Namespace) -> int:
    """Classic-arrestin IP complexes among benchmark entries, with each copy's contacting residues."""
    from scipy.spatial import cKDTree

    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.validation.burial_metrics import compute_ligand_burial, find_ligand_instances
    from learned_screen import fetch_uniprot_fasta
    from redocking import align_sequences, chain_sequence, structure_path

    entries = pd.read_csv(args.entries, dtype=str)
    accs = sorted({a.strip() for v in entries["uniprot_ids"].dropna() for a in str(v).split(";") if a.strip()})
    frames = []
    for i in range(0, len(accs), 50):
        chunk = accs[i:i + 50]
        frames.append(uniprot_tsv(" OR ".join(f"accession:{a}" for a in chunk)))
    info = classify(pd.concat(frames, ignore_index=True))
    classic = set(info.loc[info["classic"], "uniprot_id"])
    refs = entries[entries["uniprot_ids"].fillna("").map(lambda v: bool(set(v.split(";")) & classic))]
    seqs = fetch_uniprot_fasta(sorted(classic))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "chains").mkdir(exist_ok=True)
    sites = []
    for row in refs.itertuples():
        path = structure_path(args.structures_dir, row.pdb_id.upper())
        if path is None:
            sites.append({"pdb_id": row.pdb_id, "error": "structure not fetched"})
            continue
        arrays = load_structure_arrays(path)
        burial = {b.instance_id: b for b in compute_ligand_burial(path, n_points=256)}
        acc_here = [a for a in row.uniprot_ids.split(";") if a in classic]
        chains = sorted(set(arrays.chain_ids[arrays.is_polymer].tolist()))
        chain_acc = {}
        for ch in chains:
            seq, _ = chain_sequence(arrays, ch)
            for acc in acc_here:
                if acc in seqs and len(seq) >= 50:
                    _, ident, cov = align_sequences(seq, seqs[acc])
                    if ident >= 0.90 and cov >= 0.5:
                        chain_acc[ch] = acc
        for key, comp_id, atoms in find_ligand_instances(arrays):
            inst = f"{comp_id}_{key[1]}_{f'{key[2]}{key[3]}'.strip()}"
            b = burial.get(inst)
            if b is None or b.burial_class == "crystal_artifact":
                continue
            heavy = atoms[arrays.elements[atoms] != "H"]
            tree = cKDTree(arrays.coords[heavy])
            by_chain: Dict[str, Set[int]] = {}
            for i in np.flatnonzero(arrays.is_polymer):
                ch = str(arrays.chain_ids[i])
                if ch in chain_acc and tree.query_ball_point(arrays.coords[i], SITE_DISTANCE):
                    by_chain.setdefault(ch, set()).add(int(arrays.resseqs[i]))
            for ch, residues in by_chain.items():
                chain_file = args.out_dir / "chains" / f"{row.pdb_id.upper()}_{ch}.pdb"
                if not chain_file.exists():
                    write_chain_pdb(path, ch, chain_file)
                sites.append({"pdb_id": row.pdb_id.upper(), "chain": ch, "accession": chain_acc[ch],
                              "copy": inst, "comp_id": comp_id, "burial_class": b.burial_class,
                              "copy_chain": key[1], "copy_resseq": int(key[2]), "copy_icode": key[3] or "",
                              "residues": sorted(residues), "chain_file": chain_file.name,
                              "centroid": arrays.coords[heavy].mean(axis=0).tolist(),
                              "homology_group_strict": getattr(row, "homology_group_strict", None)})
    payload = {"classic_accessions_in_benchmark": sorted(classic),
               "reference_entries": sorted(refs["pdb_id"].str.upper()), "sites": sites}
    (args.out_dir / "references.json").write_text(json.dumps(payload, indent=2, default=str))
    print(json.dumps({"classic": sorted(classic), "entries": len(refs), "sites": len(sites)}))
    return 0


# ------------------------------------------------------------- mapping
def parse_usalign(text: str) -> Dict[str, object]:
    tm = re.findall(r"TM-score= ([0-9.]+) \(normalized by length of Structure_(\d)", text)
    out: Dict[str, object] = {f"tm_norm_structure{k}": float(v) for v, k in tm}
    m = re.search(r"Aligned length=\s*(\d+), RMSD=\s*([0-9.]+)", text)
    if m:
        out["aligned_length"], out["rmsd"] = int(m.group(1)), float(m.group(2))
    lines = text.splitlines()
    idx = next((i for i, ln in enumerate(lines) if ln.startswith('(":" denotes')), None)
    if idx is None or idx + 3 >= len(lines):
        raise ValueError("no alignment block in US-align output")
    out["seq1"], out["markers"], out["seq2"] = lines[idx + 1], lines[idx + 2], lines[idx + 3]
    return out


def mapping_from_alignment(seq1: str, markers: str, seq2: str, numbers1: Sequence[int],
                           numbers2: Sequence[int]) -> Dict[int, int]:
    """Structure-2 residue number -> structure-1 residue number, for pairs marked ':' (< 5 Å)."""
    i = j = 0
    out: Dict[int, int] = {}
    for a, mk, b in zip(seq1, markers, seq2):
        if a != "-" and b != "-" and mk == ":":
            out[int(numbers2[j])] = int(numbers1[i])
        if a != "-":
            i += 1
        if b != "-":
            j += 1
    return out


def usalign(target: Path, reference: Path) -> Dict[str, object]:
    exe = shutil.which("USalign") or shutil.which("usalign")
    if exe is None:
        raise RuntimeError("US-align is not installed")
    proc = subprocess.run([exe, str(target), str(reference)], capture_output=True, text=True, timeout=600)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-300:])
    return parse_usalign(proc.stdout)


def parse_residue_list(text: object) -> Set[int]:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return set()
    out = set()
    for token in re.split(r"[;,\s]+", str(text)):
        m = re.search(r"(-?\d+)[A-Za-z]?$", token.strip())
        if m:
            out.add(int(m.group(1)))
    return out


def jaccard(a: Set[int], b: Set[int]) -> float:
    return len(a & b) / len(a | b) if (a or b) else 0.0


def cmd_map(args: argparse.Namespace) -> int:
    refs = json.loads((args.refs_dir / "references.json").read_text())
    family = pd.read_csv(args.family)
    targets = family[(family["alpha"]) & (family["organism_key"].isin(["human", "yeast"]))]["uniprot_id"].tolist()
    targets = list(dict.fromkeys(targets + list(LEADS) + list(CLASSIC_HUMAN)))
    top = family.set_index("uniprot_id")
    by_chain: Dict[str, List[dict]] = {}
    for s in refs["sites"]:
        if "chain_file" in s:
            by_chain.setdefault(s["chain_file"], []).append(s)
    results = []
    for acc in targets:
        model = next(iter(sorted(args.af_dir.glob(f"AF-{acc}-F1-*.pdb"))), None)
        entry: Dict[str, object] = {"uniprot_id": acc, "gene": top["gene"].get(acc) if acc in top.index else None}
        if model is None:
            entry["error"] = "no AlphaFold model"
            results.append(entry)
            continue
        _, tnum, _, _ = ca_residues(model)
        alignments = []
        for chain_file in by_chain:
            ref_path = args.refs_dir / "chains" / chain_file
            try:
                al = usalign(model, ref_path)
            except Exception as exc:  # noqa: BLE001
                alignments.append({"reference": chain_file, "error": str(exc)[:200]})
                continue
            _, rnum, _, _ = ca_residues(ref_path)
            al["reference"] = chain_file
            al["map"] = mapping_from_alignment(al["seq1"], al["markers"], al["seq2"], tnum, rnum)
            alignments.append(al)
        ok = [a for a in alignments if "tm_norm_structure2" in a]
        keys = ("reference", "tm_norm_structure1", "tm_norm_structure2", "rmsd", "aligned_length", "error")
        entry["alignments"] = [{k: a.get(k) for k in keys} for a in alignments]
        if not ok:
            entry["error"] = "no alignment"
            results.append(entry)
            continue
        best = max(ok, key=lambda a: a["tm_norm_structure2"])
        entry["reference"] = best["reference"]
        entry["tm_reference_norm"] = best["tm_norm_structure2"]
        entry["tm_target_norm"] = best["tm_norm_structure1"]
        pocket = parse_residue_list(top["top_pocket_residues"].get(acc)) if acc in top.index else set()
        entry["top_pocket_residues"] = sorted(pocket)
        sites = []
        if best["tm_norm_structure2"] >= MIN_TM:
            for s in by_chain[best["reference"]]:
                mapping = {int(r): best["map"][int(r)] for r in s["residues"] if int(r) in best["map"]}
                mapped = set(mapping.values())
                sites.append({"copy": s["copy"], "pdb_id": s["pdb_id"], "reference_residues": s["residues"],
                              "mapping": mapping, "mapped_residues": sorted(mapped),
                              "fraction_mapped": len(mapping) / max(1, len(s["residues"])),
                              "jaccard_top_pocket": jaccard(mapped, pocket) if mapped else 0.0})
        entry["sites"] = sites
        entry["has_mapped_site"] = bool([s for s in sites if s["mapped_residues"]])
        if sites and pocket:
            lead = max(sites, key=lambda s: s["jaccard_top_pocket"])
            entry["lead_site"] = lead["copy"]
            entry["lead_jaccard"] = lead["jaccard_top_pocket"]
            entry["overlap"] = lead["jaccard_top_pocket"] >= OVERLAP_JACCARD
        elif sites:
            lead = max(sites, key=lambda s: len(s["mapped_residues"]))
            entry["lead_site"] = lead["copy"]
            entry["overlap"] = False
        else:
            entry["overlap"] = False
        results.append(entry)
    args.out.write_text(json.dumps(results, indent=2, default=str))
    for r in results:
        print(r["uniprot_id"], r.get("gene"), r.get("tm_reference_norm"), r.get("lead_jaccard"), r.get("overlap"),
              r.get("error", ""))
    return 0


# --------------------------------------------------------- conservation
def uniref50_members(acc: str, fetch=None) -> Tuple[Optional[str], Dict[str, str]]:
    from learned_screen import parse_fasta

    url = "https://rest.uniprot.org/uniref/search?" + urllib.parse.urlencode(
        {"query": f"uniprot_id:{acc} AND identity:0.5", "fields": "id", "format": "tsv"})
    text = _fetch_text(url, f"uniref:{acc}", None, fetch)
    ids = [ln.strip() for ln in text.splitlines()[1:] if ln.strip()]
    if not ids:
        return None, {}
    cluster = ids[0]
    url = "https://rest.uniprot.org/uniprotkb/stream?" + urllib.parse.urlencode(
        {"query": f"uniref_cluster_50:{cluster}", "format": "fasta"})
    return cluster, parse_fasta(_fetch_text(url, f"members:{cluster}", None, fetch))


def conservation(alignment: Mapping[str, str], target: str, positions: Sequence[int]) -> Dict[int, float]:
    """Target residue number -> fraction of other sequences with K/R/H in that column (gaps count as not basic)."""
    tseq = alignment[target]
    col_of: Dict[int, int] = {}
    k = 0
    for col, ch in enumerate(tseq):
        if ch != "-":
            k += 1
            col_of[k] = col
    others = [s for name, s in alignment.items() if name != target]
    out = {}
    for pos in positions:
        col = col_of.get(int(pos))
        if col is None or not others:
            continue
        out[int(pos)] = float(np.mean([s[col].upper() in BASIC for s in others]))
    return out


def conservation_decision(target_seq: str, lead_positions: Sequence[int], scores: Mapping[int, float],
                          n_homologues: int) -> Dict[str, object]:
    basic_positions = [p for p in lead_positions if 0 < p <= len(target_seq) and target_seq[p - 1] in BASIC]
    conserved = [p for p in basic_positions if scores.get(p, 0.0) >= CONSERVED_FRACTION]
    out = {"homologues": n_homologues, "basic_positions": basic_positions,
           "basic_fraction": {p: scores.get(p) for p in basic_positions}, "conserved_positions": conserved}
    if n_homologues < MIN_HOMOLOGUES:
        out["decision"] = f"not evaluable: {n_homologues} homologues (fewer than {MIN_HOMOLOGUES})"
    elif len(basic_positions) >= MIN_CONSERVED and len(conserved) >= MIN_CONSERVED:
        out["decision"] = "conserved"
    else:
        out["decision"] = "not conserved"
    return out


def read_alignment(path: Path) -> Dict[str, str]:
    from learned_screen import parse_fasta

    return parse_fasta(path.read_text())


def cmd_conserve(args: argparse.Namespace) -> int:
    mapping = {m["uniprot_id"]: m for m in json.loads(args.mapping.read_text())}
    out = {}
    for acc in LEADS:
        m = mapping.get(acc, {})
        lead = next((s for s in m.get("sites", []) if s["copy"] == m.get("lead_site")), None)
        entry: Dict[str, object] = {"gene": LEADS[acc]}
        if lead is None:
            entry["decision"] = "not evaluable: no mapped site"
            out[acc] = entry
            continue
        try:
            cluster, seqs = uniref50_members(acc)
        except Exception as exc:  # noqa: BLE001
            entry["decision"] = f"not evaluable: homologue fetch failed ({exc})"[:200]
            out[acc] = entry
            continue
        entry["uniref50"] = cluster
        unique: Dict[str, str] = {}
        for name in sorted(seqs):
            if seqs[name] not in unique.values() or name == acc:
                unique[name] = seqs[name]
        if acc not in unique:
            from learned_screen import fetch_uniprot_fasta

            unique.update(fetch_uniprot_fasta([acc]))
        names = [acc] + [n for n in sorted(unique) if n != acc][:MAX_HOMOLOGUES]
        fasta = args.out_dir / f"{acc}_homologues.fasta"
        args.out_dir.mkdir(parents=True, exist_ok=True)
        fasta.write_text("".join(f">{n}\n{unique[n]}\n" for n in names))
        aligned = args.out_dir / f"{acc}_aligned.fasta"
        if len(names) > 1:
            proc = subprocess.run(["mafft", "--auto", "--anysymbol", "--quiet", str(fasta)], capture_output=True,
                                  text=True, timeout=3600)
            aligned.write_text(proc.stdout)
            alignment = read_alignment(aligned)
        else:
            alignment = {acc: unique[acc]}
        scores = conservation(alignment, acc, lead["mapped_residues"])
        entry.update(conservation_decision(unique[acc], lead["mapped_residues"], scores, len(names) - 1))
        out[acc] = entry
    args.out.write_text(json.dumps(out, indent=2, default=str))
    print(json.dumps(out, indent=1, default=str)[:3000])
    return 0


# --------------------------------------------------------------- docking
def residue_ca_centroid(model: Path, residues: Sequence[int]) -> Optional[List[float]]:
    from cryptic_ip.analysis.structure_arrays import load_structure_arrays

    arrays = load_structure_arrays(model)
    mask = (arrays.atom_names == "CA") & np.isin(arrays.resseqs, list(residues))
    if not mask.any():
        return None
    return arrays.coords[mask].mean(axis=0).tolist()


def cmd_tasks(args: argparse.Namespace) -> int:
    """Every docking task of B5, fixed before any docking."""
    mapping = {m["uniprot_id"]: m for m in json.loads(args.mapping.read_text())}
    refs = json.loads((args.refs_dir / "references.json").read_text())
    pockets = pd.read_csv(args.lead_pockets)
    tasks: List[dict] = []
    rng = np.random.default_rng(SEED_NEG)
    for acc, gene in LEADS.items():
        m = mapping.get(acc, {})
        model = next(iter(sorted(args.af_dir.glob(f"AF-{acc}-F1-*.pdb"))), None)
        if model is None:
            continue
        sites = {s["copy"]: s for s in m.get("sites", [])}
        site_list = [sites[m["lead_site"]]] if m.get("lead_site") in sites else []
        if not site_list:
            site_list = [s for s in sites.values() if s["mapped_residues"]]
        exclude: Set[int] = set(m.get("top_pocket_residues", []))
        for s in site_list:
            exclude |= set(s["mapped_residues"])
            c = residue_ca_centroid(model, s["mapped_residues"])
            for lig in ("IHP", "ATP"):
                if c is not None:
                    tasks.append({"protein": acc, "gene": gene, "kind": "lead_site" if s["copy"] == m.get("lead_site")
                                  else "mapped_site", "site": s["copy"], "ligand": lig, "receptor": model.name,
                                  "centre": c, "residues": s["mapped_residues"], "seeds": [1, 2, 3]})
        top = m.get("top_pocket_residues") or []
        c = residue_ca_centroid(model, top) if top else None
        if c is not None:
            for lig in ("IHP", "ATP"):
                tasks.append({"protein": acc, "gene": gene, "kind": "top_pocket", "site": "top", "ligand": lig,
                              "receptor": model.name, "centre": c, "residues": top, "seeds": [1, 2, 3]})
        mine = pockets[(pockets["uniprot_id"] == acc) & (pd.to_numeric(pockets["plddt_mean"], errors="coerce") >= 70)]
        lead_res = set(site_list[0]["mapped_residues"]) if site_list else set()
        eligible = []
        for p in mine.sort_values("pocket_id").itertuples():
            res = parse_residue_list(p.pocket_residues)
            if res and jaccard(res, lead_res) < NEGATIVE_JACCARD and jaccard(res, set(top)) < NEGATIVE_JACCARD:
                eligible.append((int(p.pocket_id), sorted(res)))
        chosen = [eligible[i] for i in sorted(rng.choice(len(eligible), size=min(N_NEGATIVES, len(eligible)),
                                                         replace=False))] if eligible else []
        for pid, res in chosen:
            c = residue_ca_centroid(model, res)
            if c is not None:
                tasks.append({"protein": acc, "gene": gene, "kind": "negative", "site": f"pocket{pid}",
                              "ligand": "IHP", "receptor": model.name, "centre": c, "residues": res,
                              "seeds": [1, 2, 3]})
    for s in refs["sites"]:
        if "chain_file" not in s:
            continue
        tasks.append({"protein": s["accession"], "kind": "positive_crystal", "site": f"{s['pdb_id']}:{s['copy']}",
                      "ligand": s["comp_id"], "pdb_id": s["pdb_id"], "copy_chain": s["copy_chain"],
                      "copy_resseq": s["copy_resseq"], "copy_icode": s["copy_icode"],
                      "homology_group_strict": s.get("homology_group_strict"), "seeds": [1, 2, 3]})
    for acc in AF_CONTROLS:
        m = mapping.get(acc, {})
        model = next(iter(sorted(args.af_dir.glob(f"AF-{acc}-F1-*.pdb"))), None)
        for s in m.get("sites", []):
            c = residue_ca_centroid(model, s["mapped_residues"]) if model and s["mapped_residues"] else None
            if c is not None:
                tasks.append({"protein": acc, "gene": CLASSIC_HUMAN[acc], "kind": "positive_af", "site": s["copy"],
                              "ligand": "IHP", "receptor": model.name, "centre": c,
                              "residues": s["mapped_residues"], "seeds": [1, 2, 3]})
    for i, t in enumerate(tasks):
        t["task_id"] = i
    args.out.write_text(json.dumps(tasks, indent=2, default=str))
    print(pd.DataFrame(tasks)[["protein", "kind", "site", "ligand"]].to_string())
    return 0


def dock_af_task(task: dict, af_dir: Path, ccd_dir: Path, work: Path, exhaustiveness: int) -> dict:
    from rdkit import Chem

    from cryptic_ip.analysis.structure_arrays import load_structure_arrays
    from cryptic_ip.docking import engine
    from cryptic_ip.docking.ligand import poses_from_pdbqt, protonate, start_pose, to_pdbqt
    from cryptic_ip.docking.receptor import prepare_receptor
    from redocking import box_side, ccd_template

    model = af_dir / task["receptor"]
    arrays = load_structure_arrays(model)
    wd = work / Path(task["receptor"]).stem
    cached = (wd / "receptor_h.pdb", wd / "receptor.pqr")
    prepared = prepare_receptor(arrays, wd, protonated=cached if all(p.exists() for p in cached) else None)
    receptor = engine.Receptor(prepared.pdbqt)
    parent, _, _ = ccd_template(ccd_dir / f"{task['ligand']}.cif")
    ligand = protonate(parent, "primary")
    side = box_side(ligand)
    centre = np.asarray(task["centre"], dtype=float)
    runs, poses = [], []
    for seed in task["seeds"]:
        pose = start_pose(ligand, centre, seed, crystal=None)
        res = engine.dock(receptor, to_pdbqt(pose.mol), centre, [side] * 3, seed=seed, exhaustiveness=exhaustiveness)
        mols = poses_from_pdbqt(res.poses_pdbqt)
        scores = [float(s) for s in res.scores[: len(mols)]]
        runs.append({"seed": seed, "scores": scores, "top_score": scores[0] if scores else None})
        poses += [(seed, rank, scores[rank], mols[rank]) for rank in range(min(TOP_PER_SEED, len(mols)))]
    out = {"runs": runs, "best_score": min(r["top_score"] for r in runs if r["top_score"] is not None),
           "box_side": side}
    poses.sort(key=lambda p: p[2])
    from cryptic_ip.docking.rmsd import rmsd_matrix

    mat = rmsd_matrix([p[3] for p in poses])
    clusters = engine.greedy_clusters(mat, CLUSTER_RMSD)
    largest = max(clusters, key=len) if clusters else []
    out["convergence"] = {"poses": len(poses), "largest_cluster": len(largest),
                          "fraction": len(largest) / len(poses) if poses else 0.0,
                          "converged": bool(poses) and len(largest) / len(poses) >= CONVERGENCE}
    best = poses[0][3] if poses else None
    if best is not None:
        from scipy.spatial import cKDTree

        tree = cKDTree(Chem.RemoveHs(best).GetConformer().GetPositions())
        contacts: Dict[int, str] = {}
        for i in np.flatnonzero(arrays.is_polymer & (arrays.elements != "H")):
            if tree.query_ball_point(arrays.coords[i], SITE_DISTANCE):
                contacts[int(arrays.resseqs[i])] = str(arrays.resnames[i])
        out["best_pose_contacts"] = [f"{contacts[k]}{k}" for k in sorted(contacts)]
        out["best_pose_basic_contacts"] = [f"{contacts[k]}{k}" for k in sorted(contacts)
                                           if contacts[k] in ("LYS", "ARG", "HIS")]
    return out


def dock_crystal_task(task: dict, structures_dir: Path, ccd_dir: Path, work: Path, exhaustiveness: int) -> dict:
    import redocking

    redocking.EXHAUSTIVENESS[0] = exhaustiveness
    row = {"pdb_id": task["pdb_id"], "chain": task["copy_chain"], "resseq": task["copy_resseq"],
           "icode": task.get("copy_icode") or ""}
    ctx = redocking.Context(row, structures_dir, ccd_dir, work)
    receptor, _ = ctx.receptor()
    runs = [redocking.run_arm(ctx, receptor, seed=s, arm=f"vina_s{s}") for s in task["seeds"]]
    success = [float(r["top_rmsd"] <= 2.0) for r in runs if r.get("top_rmsd") is not None]
    return {"runs": [{k: r.get(k) for k in ("seed", "top_score", "top_rmsd", "best_rmsd")} for r in runs],
            "best_score": min(r["top_score"] for r in runs), "success": float(np.mean(success)) if success else None}


def cmd_dock(args: argparse.Namespace) -> int:
    tasks = json.loads(args.tasks.read_text())
    mine = [t for t in tasks if t["task_id"] % args.shard_count == args.shard_index]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / f"dock_{args.shard_index}.jsonl"
    with open(out, "w") as fh:
        for t in mine:
            rec = {"task_id": t["task_id"]}
            try:
                if t["kind"] == "positive_crystal":
                    rec.update(dock_crystal_task(t, args.structures_dir, args.ccd_dir, args.work_dir,
                                                 args.exhaustiveness))
                else:
                    rec.update(dock_af_task(t, args.af_dir, args.ccd_dir, args.work_dir, args.exhaustiveness))
            except Exception as exc:  # noqa: BLE001
                import traceback

                rec["error"] = f"{type(exc).__name__}: {exc}"[:400]
                rec["traceback"] = traceback.format_exc()[-1200:]
            fh.write(json.dumps(rec, default=float) + "\n")
            fh.flush()
            LOGGER.info("task %s %s %s %s: %s", t["task_id"], t.get("gene", t["protein"]), t["kind"], t["ligand"],
                        rec.get("error", rec.get("best_score")))
    return 0


# ---------------------------------------------------------------- decide
def decide(tasks: Sequence[dict], results: Mapping[int, dict], mapping: Mapping[str, dict],
           conserve: Mapping[str, dict]) -> Dict[str, object]:
    by = pd.DataFrame([{**t, **{k: v for k, v in results.get(t["task_id"], {}).items() if k != "runs"}}
                       for t in tasks])
    crystal = by[by["kind"] == "positive_crystal"]
    validity_values = pd.to_numeric(crystal.get("success"), errors="coerce").dropna() if len(crystal) else []
    validity = float(np.mean(validity_values)) if len(validity_values) else float("nan")
    protocol_valid = bool(np.isfinite(validity) and validity >= 0.5)
    positives = by[by["kind"].isin(["positive_crystal", "positive_af"]) & (by["ligand"] == "IHP")]
    pos_scores = pd.to_numeric(positives.get("best_score"), errors="coerce").dropna()
    weakest_positive = float(pos_scores.max()) if len(pos_scores) else float("nan")
    out: Dict[str, object] = {"protocol_validity": {"crystal_sites": int(len(crystal)),
                                                    "mean_top_pose_success": validity, "valid": protocol_valid},
                              "positive_control_best_scores": pos_scores.tolist(),
                              "weakest_positive": weakest_positive, "proteins": {}}
    for acc, gene in LEADS.items():
        m = mapping.get(acc, {})
        mine = by[(by["protein"] == acc)]
        lead = mine[(mine["kind"] == "lead_site") & (mine["ligand"] == "IHP")]
        neg = pd.to_numeric(mine[mine["kind"] == "negative"].get("best_score"), errors="coerce").dropna()
        crit: Dict[str, Dict[str, object]] = {}
        crit["1_overlap"] = {"pass": bool(m.get("overlap")), "jaccard": m.get("lead_jaccard"),
                             "tm_score": m.get("tm_reference_norm"), "reference": m.get("reference"),
                             "detail": "no mapped site" if not m.get("has_mapped_site") else ""}
        c = conserve.get(acc, {})
        crit["2_conservation"] = {"pass": c.get("decision") == "conserved", "decision": c.get("decision"),
                                  "conserved_positions": c.get("conserved_positions"),
                                  "basic_positions": c.get("basic_positions")}
        if not protocol_valid:
            crit["3_convergence"] = {"pass": False, "decision": "not evaluable: protocol validity failed"}
            crit["4_scores"] = {"pass": False, "decision": "not evaluable: protocol validity failed"}
        elif lead.empty or lead.iloc[0].get("error") or pd.isna(lead.iloc[0].get("best_score")):
            crit["3_convergence"] = {"pass": False, "decision": "not evaluable: no lead-site docking"}
            crit["4_scores"] = {"pass": False, "decision": "not evaluable: no lead-site docking"}
        else:
            conv = lead.iloc[0].get("convergence") or {}
            crit["3_convergence"] = {"pass": bool(conv.get("converged")), **conv}
            best = float(lead.iloc[0]["best_score"])
            ok_pos = np.isfinite(weakest_positive) and best <= weakest_positive
            ok_neg = len(neg) > 0 and best < float(neg.min())
            crit["4_scores"] = {"pass": bool(ok_pos and ok_neg), "lead_best_score": best,
                                "weakest_positive": weakest_positive,
                                "best_negative": float(neg.min()) if len(neg) else None,
                                "negatives": neg.tolist()}
        failing = [k for k, v in crit.items() if not v.get("pass")]
        verdict = "supported for experiment" if not failing else "not supported"
        extra = {k: mine[(mine["kind"] == k2) & (mine["ligand"] == lig)]["best_score"].tolist()
                 for k, k2, lig in (("top_pocket_IP6", "top_pocket", "IHP"), ("top_pocket_ATP", "top_pocket", "ATP"),
                                    ("lead_site_ATP", "lead_site", "ATP"))}
        contacts = lead.iloc[0].get("best_pose_contacts") if len(lead) else None
        out["proteins"][acc] = {"gene": gene, "verdict": verdict, "failing": failing, "criteria": crit,
                                "descriptive_scores": extra, "best_pose_contacts": contacts}
    return out


def dossier(acc: str, d: Mapping[str, object], mapping: Mapping[str, dict], conserve: Mapping[str, dict]) -> str:
    m = mapping.get(acc, {})
    lead = next((s for s in m.get("sites", []) if s["copy"] == m.get("lead_site")), None)
    lines = [f"### {d['gene']} ({acc}): **{d['verdict']}**", ""]
    if d["failing"]:
        lines += [f"Failing criteria: {', '.join(d['failing'])}.", ""]
    lines += ["| criterion | pass | detail |", "|---|---|---|"]
    for k, v in d["criteria"].items():
        detail = {kk: vv for kk, vv in v.items() if kk != "pass"}
        lines.append(f"| {k} | {'yes' if v.get('pass') else 'no'} | {json.dumps(detail, default=str)[:300]} |")
    if lead:
        pairs = ", ".join(f"{r}→{t}" for r, t in sorted(lead["mapping"].items(), key=lambda x: int(x[0])))
        lines += ["", f"Mapped site (reference {m.get('reference')}, copy {lead['copy']}; TM-score "
                  f"{m.get('tm_reference_norm')}): reference residue → {d['gene']} residue: {pairs}."]
    c = conserve.get(acc, {})
    if c.get("basic_fraction"):
        fractions = ", ".join(f"{p}: {f:.2f}" for p, f in c["basic_fraction"].items() if f is not None)
        lines += ["", "Basic residues at the mapped site and their conservation across UniRef50 homologues "
                  f"({c.get('homologues')} sequences): {fractions}."]
    if d.get("best_pose_contacts"):
        lines += ["", "Residues within 4 Å of the best docked IP6 pose: " + ", ".join(d["best_pose_contacts"]) + "."]
    lines += ["", "What would test it: isothermal titration calorimetry or a fluorescence binding assay of the "
              "purified protein with IP6 (and ATP as a polyanion control); charge-reversal mutants (K/R to E) of "
              "the mapped basic residues, which should abolish binding if the site is real. Docking scores for a "
              "-9 polyanion from a scoring function without electrostatics are weak evidence whatever they say.", ""]
    return "\n".join(lines)


def cmd_decide(args: argparse.Namespace) -> int:
    tasks = json.loads(args.tasks.read_text())
    results: Dict[int, dict] = {}
    for path in sorted(args.dock_dir.rglob("dock_*.jsonl")):
        for line in path.read_text().splitlines():
            if line.strip():
                r = json.loads(line)
                results[int(r["task_id"])] = r
    mapping = {m["uniprot_id"]: m for m in json.loads(args.mapping.read_text())}
    conserve = json.loads(args.conserve.read_text())
    family = json.loads(args.family.read_text())
    decision = decide(tasks, results, mapping, conserve)
    report = {"plan": "docs/ARRESTIN_PLAN.md", "B1": family.get("B1"),
              "B1_per_organism_descriptive": family.get("B1_per_organism_descriptive"),
              "family_counts": family.get("family_counts"),
              "positive_control_classic": family.get("positive_control_classic"),
              "mapping": [{k: v for k, v in m.items() if k not in ("sites",)} for m in mapping.values()],
              "conservation": conserve, "docking_tasks": len(tasks),
              "docking_failures": {t: r.get("error") for t, r in results.items() if r.get("error")},
              "docking": [{**t, **{k: v for k, v in results.get(t["task_id"], {}).items() if k != "traceback"}}
                          for t in tasks], **decision}
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "arrestin.json").write_text(json.dumps(report, indent=2, default=str))
    b1 = report["B1"] or {}
    auc = b1.get("roc_auc") or {}
    lines = ["## The α-arrestin lead (docs/ARRESTIN_PLAN.md)", "",
             f"**B1 (family ranks high among unseen proteins):** {b1.get('decision')} - ROC-AUC "
             f"{auc.get('point', float('nan')):.3f}, 95 % [{auc.get('low', float('nan')):.3f}, "
             f"{auc.get('high', float('nan')):.3f}], 5th percentile {auc.get('p5', float('nan')):.3f}, "
             f"{b1.get('members_unseen')} unseen α-arrestins in {b1.get('member_clusters')} clusters.", "",
             f"Protocol validity (crystal redocks of arrestin-IP sites): {decision['protocol_validity']}.", ""]
    lines += ["| protein | overlap | conservation | convergence | scores | verdict |", "|---|---|---|---|---|---|"]
    for acc, d in decision["proteins"].items():
        c = d["criteria"]
        marks = " | ".join("yes" if c[k].get("pass") else "no" for k in sorted(c))
        lines.append(f"| {d['gene']} ({acc}) | {marks} | **{d['verdict']}** |")
    lines += ["", "## Dossiers", ""]
    for acc, d in decision["proteins"].items():
        lines.append(dossier(acc, d, mapping, conserve))
    text = "\n".join(lines) + "\n"
    (args.out_dir / "ARRESTIN.md").write_text(text)
    print(text)
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("family")
    f.add_argument("--proteins", type=Path, required=True)
    f.add_argument("--out-dir", type=Path, required=True)
    f.add_argument("--n-bootstrap", type=int, default=2000)
    r = sub.add_parser("references")
    r.add_argument("--entries", type=Path, required=True)
    r.add_argument("--structures-dir", type=Path, required=True)
    r.add_argument("--out-dir", type=Path, required=True)
    m = sub.add_parser("map")
    m.add_argument("--refs-dir", type=Path, required=True)
    m.add_argument("--family", type=Path, required=True)
    m.add_argument("--af-dir", type=Path, required=True)
    m.add_argument("--out", type=Path, required=True)
    c = sub.add_parser("conserve")
    c.add_argument("--mapping", type=Path, required=True)
    c.add_argument("--out-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    t = sub.add_parser("tasks")
    t.add_argument("--mapping", type=Path, required=True)
    t.add_argument("--refs-dir", type=Path, required=True)
    t.add_argument("--af-dir", type=Path, required=True)
    t.add_argument("--lead-pockets", type=Path, required=True)
    t.add_argument("--out", type=Path, required=True)
    d = sub.add_parser("dock")
    d.add_argument("--tasks", type=Path, required=True)
    d.add_argument("--af-dir", type=Path, required=True)
    d.add_argument("--structures-dir", type=Path, required=True)
    d.add_argument("--ccd-dir", type=Path, required=True)
    d.add_argument("--work-dir", type=Path, default=Path("/tmp/arrestin_work"))
    d.add_argument("--out-dir", type=Path, required=True)
    d.add_argument("--shard-index", type=int, default=0)
    d.add_argument("--shard-count", type=int, default=1)
    d.add_argument("--exhaustiveness", type=int, default=32, help="the plan fixes 32; tests only")
    e = sub.add_parser("decide")
    e.add_argument("--tasks", type=Path, required=True)
    e.add_argument("--dock-dir", type=Path, required=True)
    e.add_argument("--mapping", type=Path, required=True)
    e.add_argument("--conserve", type=Path, required=True)
    e.add_argument("--family", type=Path, required=True)
    e.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    return {"family": cmd_family, "references": cmd_references, "map": cmd_map, "conserve": cmd_conserve,
            "tasks": cmd_tasks, "dock": cmd_dock, "decide": cmd_decide}[args.command](args)


if __name__ == "__main__":
    sys.exit(main())
