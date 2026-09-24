"""Phase 1 success criteria, measured exactly as the project plan states them.

The plan's decision point for proteome screening is a table of criteria on the
controls (sections 5-6):

=============================  =====================================  =========
criterion                      target (ADAR2 / Pds5B / HDAC1-3)       measured by
=============================  =====================================  =========
IP site pocket rank            top 3 / top 5 / top 5                  fpocket + composite
SASA at the site               < 5 / < 10 / < 15 A^2                  per-residue SASA
electrostatic potential        > +5 / > +4 / > +3 kT/e                APBS at pocket centre
basic residues                 >= 6 / >= 4 / >= 4 within 5 A          distance
AlphaFold vs crystal           RMSD < 2 A over the binding region     superposition
negative controls              IP site pocket in the bottom half      rank
separation                     no overlap, buried vs surface scores   composite
=============================  =====================================  =========

Every control is evaluated the way a proteome target is: the ligand, waters
and every non-polymer atom are removed first (the plan's "clean structure"
step), then pockets are detected and scored. Scoring the deposited holo
structure instead lets the bound ligand occlude its own pocket - measured on
this panel, it makes lining residues read ~0 A^2 that read hundreds of A^2 once
the ligand is gone - which is a signal no AlphaFold model can show. The holo
ligand is used only to say *which* pocket is the site; apo and holo share a
coordinate frame, so that identification is exact.

The same measurements are then repeated on each control's **AlphaFold model**,
which is the actual screening input. The crystal ligand is carried into the
model's frame by superposing the binding-region C-alpha atoms, and that
superposition's RMSD is the plan's AlphaFold-versus-crystal criterion.
"""

from __future__ import annotations

import json
import logging
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.PDB.Polypeptide import three_to_index, index_to_one
from Bio.PDB.SASA import ShrakeRupley

from ..analysis.structure_arrays import is_polymer_residue, load_structure_arrays
from .burial_metrics import find_ligand_instances

LOGGER = logging.getLogger(__name__)

BASIC = {"ARG", "LYS", "HIS"}
#: Heavy-atom distance defining a residue as coordinating the ligand.
CONTACT_DISTANCE = 4.0
#: The plan's "basic residues within 5 A".
BASIC_DISTANCE = 5.0
#: C-alpha atoms within this distance of the ligand define the binding region
#: for the AlphaFold-versus-crystal RMSD.
BINDING_REGION_RADIUS = 10.0


@dataclass(frozen=True)
class Control:
    """One control structure and the plan's expectation for it."""

    name: str
    pdb_id: str
    role: str  # "positive" or "negative"
    max_rank: Optional[int] = None  # positives: site must rank at or above this
    max_site_sasa: Optional[float] = None
    min_potential: Optional[float] = None
    min_basic: Optional[int] = None
    note: str = ""


#: The plan's panel (sections 3.3, 4 and 6). 1BTK is the identifier the plan
#: lists for Btk; if that entry carries no inositol phosphate, 1BWN - the Btk
#: PH domain with Ins(1,3,4,5)P4 - measures the same site.
CONTROLS: Tuple[Control, ...] = (
    Control("ADAR2", "1ZY7", "positive", 3, 5.0, 5.0, 6, "buried InsP6, folding cofactor"),
    Control("Pds5B", "5HDT", "positive", 5, 10.0, 4.0, 4, "InsP6 from insect-cell expression"),
    Control("HDAC1", "5ICN", "positive", 5, 15.0, 3.0, 4, "InsP4 at the MTA1 interface"),
    Control("HDAC3", "4A69", "positive", 5, 15.0, 3.0, 4, "InsP4 at the SMRT interface"),
    Control("PLCd1_PH", "1MAI", "negative", note="InsP3, surface signalling site"),
    Control("Btk_PH", "1BTK", "negative", note="plan's identifier for Btk"),
    Control("Btk_PH_IP4", "1BWN", "negative", note="Btk PH domain with InsP4"),
    Control("DAPP1_PH", "1FAO", "negative", note="InsP4, surface signalling site"),
    Control("Grp1_PH", "1FGY", "negative", note="InsP4, surface signalling site"),
)


# ------------------------------------------------------------------ structure
class _ChainPolymer(Select):
    """Keep one model, one chain, polymer residues only."""

    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return is_polymer_residue(residue.id[0], residue.get_resname())


def write_apo_chain(path: Path, chain_id: str, out_path: Path) -> Path:
    """Write one chain of ``path`` with every non-polymer atom removed.

    A monomer screen sees one chain. Restricting a crystal control to the chain
    that holds the ligand keeps the pocket count - and therefore the rank
    criterion - comparable with an AlphaFold monomer, and removes interface
    partners (MTA1 at HDAC1's site, for example) that a monomer model lacks.
    """
    structure = PDBParser(QUIET=True).get_structure("control", str(path))
    io = PDBIO()
    io.set_structure(structure)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_path), _ChainPolymer(chain_id))
    return out_path


class _ChainWithLigand(_ChainPolymer):
    """One chain's polymer residues plus one ligand copy; no waters."""

    def __init__(self, chain_id: str, ligand_chain: str, ligand_resseq: int):
        super().__init__(chain_id)
        self.ligand = (ligand_chain, ligand_resseq)

    def accept_chain(self, chain):
        return chain.id in (self.chain_id, self.ligand[0])

    def accept_residue(self, residue):
        chain = residue.get_parent().id
        if chain == self.chain_id and super().accept_residue(residue):
            return True
        return residue.id[0].startswith("H_") and (chain, int(residue.id[1])) == self.ligand


def write_holo_site(path: Path, chain_id: str, ligand_key: Sequence[str], out_path: Path) -> Path:
    """Write one chain with its ligand copy and without waters or other hetero groups.

    Used for the crystallographic view of site accessibility: the ligand
    occludes its coordinating residues, as in the complex, but crystal waters -
    which would occlude them too - are not counted as protein burial.
    """
    structure = PDBParser(QUIET=True).get_structure("control", str(path))
    io = PDBIO()
    io.set_structure(structure)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out_path), _ChainWithLigand(chain_id, str(ligand_key[1]), int(ligand_key[2])))
    return out_path


def select_ligand_copy(path: Path) -> Optional[Dict[str, Any]]:
    """The inositol phosphate copy with the most protein contacts, and its chain."""
    arrays = load_structure_arrays(path)
    instances = find_ligand_instances(arrays)
    if not instances:
        return None
    protein = np.where(arrays.is_polymer)[0]
    best = None
    for key, comp_id, atom_indices in instances:
        ligand = arrays.coords[atom_indices]
        distances = np.linalg.norm(
            arrays.coords[protein][:, None, :] - ligand[None, :, :], axis=2
        ).min(axis=1)
        near = protein[distances <= CONTACT_DISTANCE]
        if len(near) == 0:
            continue
        chains, counts = np.unique(arrays.chain_ids[near], return_counts=True)
        chain = str(chains[np.argmax(counts)])
        n_contacts = int(len(near))
        if best is None or n_contacts > best["n_contacts"]:
            best = {
                "residue_key": [str(k) for k in key],
                "comp_id": comp_id,
                "coords": ligand,
                "chain": chain,
                "n_contacts": n_contacts,
            }
    return best


def _residues(structure, chain_id: Optional[str] = None):
    model = next(iter(structure))
    for chain in model:
        if chain_id is not None and chain.id != chain_id:
            continue
        for residue in chain:
            if is_polymer_residue(residue.id[0], residue.get_resname()):
                yield chain.id, residue


def contact_residues(
    structure_path: Path, ligand: np.ndarray, chain_id: str, cutoff: float
) -> List[Tuple[str, int, str]]:
    """Polymer residues of ``chain_id`` with a heavy atom within ``cutoff`` of the ligand."""
    structure = PDBParser(QUIET=True).get_structure("s", str(structure_path))
    found = []
    for chain, residue in _residues(structure, chain_id):
        coords = np.array([a.coord for a in residue if a.element != "H"])
        if len(coords) and np.linalg.norm(coords[:, None] - ligand[None], axis=2).min() <= cutoff:
            found.append((chain, int(residue.id[1]), residue.get_resname()))
    return found


def residue_sasa(structure_path: Path, keys: Sequence[Tuple[str, int]]) -> Dict[Tuple[str, int], float]:
    """Per-residue SASA (Shrake-Rupley) for selected residues of a structure."""
    structure = PDBParser(QUIET=True).get_structure("s", str(structure_path))
    ShrakeRupley().compute(structure[0], level="R")
    wanted = set(keys)
    return {
        (chain, int(residue.id[1])): float(residue.sasa)
        for chain, residue in _residues(structure)
        if (chain, int(residue.id[1])) in wanted
    }


# ------------------------------------------------------------ superposition
def _ca_table(structure_path: Path, chain_id: Optional[str]) -> Tuple[str, List[Tuple[str, int]], np.ndarray]:
    structure = PDBParser(QUIET=True).get_structure("s", str(structure_path))
    sequence, keys, coords = [], [], []
    for chain, residue in _residues(structure, chain_id):
        if "CA" not in residue:
            continue
        try:
            letter = index_to_one(three_to_index(residue.get_resname()))
        except (KeyError, ValueError):
            letter = "X"
        sequence.append(letter)
        keys.append((chain, int(residue.id[1])))
        coords.append(residue["CA"].coord)
    return "".join(sequence), keys, np.asarray(coords, dtype=float)


def kabsch(mobile: np.ndarray, target: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Optimal rotation taking ``mobile`` onto ``target``.

    Returns ``(rotation, mobile_centroid, target_centroid, rmsd)``; a point
    ``x`` maps as ``(x - mobile_centroid) @ rotation + target_centroid``.
    """
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    h = (mobile - mc).T @ (target - tc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(u @ vt))
    rotation = u @ np.diag([1.0, 1.0, d]) @ vt
    moved = (mobile - mc) @ rotation + tc
    rmsd = float(np.sqrt(np.mean(np.sum((moved - target) ** 2, axis=1))))
    return rotation, mc, tc, rmsd


def binding_region_superposition(
    crystal_path: Path,
    chain_id: str,
    model_path: Path,
    ligand: np.ndarray,
    radius: float = BINDING_REGION_RADIUS,
) -> Dict[str, Any]:
    """Superpose the crystal binding region onto the AlphaFold model.

    Residues are paired by sequence alignment, not by number: crystal
    constructs are often renumbered, truncated or taken from another isoform or
    species, and pairing by number would silently compare different residues.
    """
    c_seq, c_keys, c_ca = _ca_table(crystal_path, chain_id)
    m_seq, m_keys, m_ca = _ca_table(model_path, None)
    if not c_seq or not m_seq:
        return {"ok": False, "reason": "no C-alpha atoms"}

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2.0, -1.0
    aligner.open_gap_score, aligner.extend_gap_score = -5.0, -0.5
    aligner.end_gap_score = 0.0
    alignment = aligner.align(c_seq, m_seq)[0]
    pairs = []
    for (c0, c1), (m0, m1) in zip(*alignment.aligned):
        pairs.extend(zip(range(c0, c1), range(m0, m1)))
    if not pairs:
        return {"ok": False, "reason": "sequences do not align"}
    pairs_arr = np.asarray(pairs)
    identity = float(np.mean([c_seq[i] == m_seq[j] for i, j in pairs]))

    near = np.linalg.norm(c_ca[:, None] - ligand[None], axis=2).min(axis=1) <= radius
    region = pairs_arr[near[pairs_arr[:, 0]]]
    if len(region) < 3:
        return {"ok": False, "reason": "binding region not covered by the model", "identity": identity}

    rotation, mc, tc, rmsd = kabsch(c_ca[region[:, 0]], m_ca[region[:, 1]])
    return {
        "ok": True,
        "rmsd": rmsd,
        "n_region_residues": int(len(region)),
        "aligned_residues": int(len(pairs)),
        "sequence_identity": identity,
        "crystal_residue_range": [c_keys[0][1], c_keys[-1][1]],
        "ligand_in_model_frame": (ligand - mc) @ rotation + tc,
    }


# ------------------------------------------------------------------ pockets
def hull_depth(protein_coords: np.ndarray, point: Sequence[float]) -> float:
    """Distance from ``point`` inward to the protein's convex hull (A).

    See :func:`cryptic_ip.analysis.geometry.hull_depths` for why this, rather
    than distance to the nearest exposed atom, measures burial on an apo site.
    """
    from ..analysis.geometry import hull_depths

    return float(hull_depths(protein_coords, np.asarray(point, dtype=float))[0])


def site_apbs_potential(structure_path: Path, center: Sequence[float], work: Path) -> Optional[float]:
    """APBS potential (kT/e) at ``center``, fine grid focused there.

    pdb2pqr at pH 7.0 with the AMBER force field, as the project plan
    specifies. ``None`` when APBS or pdb2pqr is unavailable or fails.
    """
    from ..analysis.electrostatics import ElectrostaticsCalculator

    calculator = ElectrostaticsCalculator()
    try:
        pqr = calculator.generate_pqr(structure_path, ph=7.0, output_dir=work, forcefield="AMBER")
        _, dx = calculator.run_apbs_with_map(pqr, work, focus_center=center)
        return float(calculator.sample_potential_at_point(dx, center))
    except Exception as exc:
        LOGGER.warning("APBS at %s failed: %s", structure_path.name, exc)
        return None


def evaluate_site(
    structure_path: Path,
    ligand: np.ndarray,
    *,
    use_apbs: bool = True,
) -> Dict[str, Any]:
    """Detect and score pockets on an apo structure; locate and rank the site."""
    from ..analysis import ProteinAnalyzer
    from .site_selection import select_ligand_pocket

    work = Path(tempfile.mkdtemp(prefix="phase1_"))
    try:
        local = work / structure_path.name
        shutil.copy(structure_path, local)
        # Pockets are scored exactly as the proteome screen scores them - the
        # screened-Coulomb surrogate, no APBS - so the rank measured here is the
        # rank the screen would assign. APBS is then run once and sampled at
        # the site, for the plan's potential criterion.
        analyzer = ProteinAnalyzer(str(local), work_dir=str(work / "work"), skip_electrostatics=True)
        scored = analyzer.run_pipeline(include_electrostatics=False)
        if scored.empty:
            return {"ok": False, "reason": "fpocket detected no pockets"}
        scored = scored.sort_values("composite_score", ascending=False).reset_index(drop=True)
        row, overlap = select_ligand_pocket(scored, analyzer, ligand)
        rank = int(scored.index[scored["pocket_id"] == row["pocket_id"]][0]) + 1
        n = int(len(scored))

        center = tuple(float(c) for c in row["center"])
        apbs = site_apbs_potential(local, center, work / "apbs") if use_apbs else None
        arrays = load_structure_arrays(local)
        depth_to_hull = hull_depth(arrays.coords[arrays.is_polymer], center)
        return {
            "ok": True,
            "n_pockets": n,
            "site_pocket_id": int(row["pocket_id"]),
            "ligand_overlap": float(overlap) if overlap == overlap else None,
            "site_rank": rank,
            "site_rank_fraction": rank / n,
            "site_score": float(row["composite_score"]),
            "top_score": float(scored["composite_score"].iloc[0]),
            "apbs_potential_kT": apbs,
            "coulomb_potential_kT": float(row.get("coulomb_potential_kt", np.nan)),
            "pocket_volume": float(row.get("volume", np.nan)),
            "burial_depth": float(row.get("burial_depth", np.nan)),
            "hull_depth": depth_to_hull,
            "enclosure": float(row.get("enclosure", np.nan)),
            "sasa_mean": float(row.get("sasa_mean", np.nan)),
            "n_basic_residues_pocket": int(row.get("n_basic_residues", 0)),
            "plddt_mean": float(row.get("plddt_mean", np.nan)),
        }
    finally:
        shutil.rmtree(work, ignore_errors=True)


def site_residue_measurements(holo: Path, apo: Path, ligand: np.ndarray, chain: str) -> Dict[str, Any]:
    """Coordinating-residue SASA (with and without the ligand) and basic counts."""
    contacts = contact_residues(holo, ligand, chain, CONTACT_DISTANCE)
    basic5 = [r for r in contact_residues(holo, ligand, chain, BASIC_DISTANCE) if r[2] in BASIC]
    keys = [(c, n) for c, n, _ in contacts]
    # "Holo" SASA keeps the ligand as an occluder: the crystallographic view.
    holo_sasa = residue_sasa(holo, keys)
    apo_sasa = residue_sasa(apo, keys)
    basic_keys = [(c, n) for c, n, name in contacts if name in BASIC]
    return {
        "coordinating_residues": [f"{name}{n}" for _, n, name in contacts],
        "basic_within_5A": [f"{name}{n}" for _, n, name in basic5],
        "n_basic_within_5A": len(basic5),
        "coordinating_sasa_holo_mean": _mean([holo_sasa.get(k) for k in keys]),
        "coordinating_sasa_apo_mean": _mean([apo_sasa.get(k) for k in keys]),
        "basic_coordinating_sasa_holo_mean": _mean([holo_sasa.get(k) for k in basic_keys]),
        "basic_coordinating_sasa_apo_mean": _mean([apo_sasa.get(k) for k in basic_keys]),
    }


def _mean(values) -> Optional[float]:
    vals = [v for v in values if v is not None]
    return float(np.mean(vals)) if vals else None


# ------------------------------------------------------------------ fetching
def _get_json(urls: Sequence[str]) -> List[Optional[Any]]:
    from ..database.async_fetch import FetchJob, fetch_all, validate_json

    results = fetch_all([FetchJob(u, u, validator=validate_json) for u in urls], concurrency=8)
    return [json.loads(r.payload) if r.ok and r.payload else None for r in results]


def uniprot_for_entry(pdb_id: str, chain: str) -> Optional[str]:
    """UniProt accession of the polymer entity carrying ``chain`` (RCSB API)."""
    (entry,) = _get_json([f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"])
    if not entry:
        LOGGER.warning("RCSB entry lookup failed for %s", pdb_id)
        return None
    entity_ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
    entities = _get_json(
        [f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{e}" for e in entity_ids]
    )
    for entity in entities:
        ids = (entity or {}).get("rcsb_polymer_entity_container_identifiers", {})
        if chain in (ids.get("auth_asym_ids") or []):
            accessions = ids.get("uniprot_ids") or []
            return accessions[0] if accessions else None
    return None


def fetch_alphafold_model(uniprot_id: str, out_dir: Path) -> Optional[Path]:
    """Download the current AlphaFold model for ``uniprot_id``.

    Resolved through the prediction API by accession and fragment (not the
    first entry listed, which can be an isoform), verified and cached.
    """
    from ..database.async_fetch import fetch_alphafold_models

    result = fetch_alphafold_models([uniprot_id], Path(out_dir), concurrency=2)[uniprot_id.upper()]
    if not result.ok:
        LOGGER.warning("AlphaFold model unavailable for %s: %s", uniprot_id, result.error)
        return None
    return Path(result.path)


# ---------------------------------------------------------------- criteria
def judge(control: Control, crystal: Dict[str, Any], model: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Apply the plan's criteria for ``control`` to its measurements.

    Each entry is ``{"target", "value", "passes"}``; ``passes`` is ``None``
    where the measurement could not be made, so an unmeasured criterion is never
    reported as met.
    """
    out: Dict[str, Any] = {}

    def item(target, value, ok):
        return {"target": target, "value": value, "passes": None if value is None else bool(ok)}

    site = crystal.get("site", {})
    residues = crystal.get("residues", {})
    if control.role == "positive":
        rank = site.get("site_rank")
        out["pocket_rank"] = item(f"top {control.max_rank}", rank, rank is not None and rank <= control.max_rank)
        sasa = residues.get("basic_coordinating_sasa_holo_mean")
        out["site_sasa_holo"] = item(f"< {control.max_site_sasa} A^2", sasa, sasa is not None and sasa < control.max_site_sasa)
        potential = site.get("apbs_potential_kT")
        out["apbs_potential"] = item(
            f"> +{control.min_potential} kT/e", potential, potential is not None and potential > control.min_potential
        )
        basic = residues.get("n_basic_within_5A")
        out["basic_residues"] = item(f">= {control.min_basic}", basic, basic is not None and basic >= control.min_basic)
        score = site.get("site_score")
        out["composite_score"] = item("> 0.7", score, score is not None and score > 0.7)
    else:
        fraction = site.get("site_rank_fraction")
        out["pocket_rank_bottom_half"] = item("bottom half", fraction, fraction is not None and fraction > 0.5)
        score = site.get("site_score")
        out["composite_score"] = item("< 0.4", score, score is not None and score < 0.4)

    if model and model.get("superposition", {}).get("ok"):
        rmsd = model["superposition"]["rmsd"]
        out["alphafold_vs_crystal_rmsd"] = item("< 2 A (binding region)", rmsd, rmsd < 2.0)
        m_site = model.get("site", {})
        if control.role == "positive" and m_site.get("ok"):
            rank = m_site.get("site_rank")
            out["alphafold_pocket_rank"] = item(
                f"top {control.max_rank}", rank, rank is not None and rank <= control.max_rank
            )
        elif m_site.get("ok"):
            fraction = m_site.get("site_rank_fraction")
            out["alphafold_pocket_rank_bottom_half"] = item(
                "bottom half", fraction, fraction is not None and fraction > 0.5
            )
    return out


def critical_test(results: Sequence[Dict[str, Any]], source: str = "crystal") -> Dict[str, Any]:
    """The plan's three conditions for proceeding to proteome screening."""

    def site(name):
        for r in results:
            if r["control"]["name"] == name and r.get(source, {}).get("site", {}).get("ok"):
                return r[source]["site"]
        return None

    adar2 = site("ADAR2")
    negatives = [
        (r["control"]["name"], r[source]["site"])
        for r in results
        if r["control"]["role"] == "negative" and r.get(source, {}).get("site", {}).get("ok")
    ]
    positives = [
        (r["control"]["name"], r[source]["site"])
        for r in results
        if r["control"]["role"] == "positive" and r.get(source, {}).get("site", {}).get("ok")
    ]
    plc, btk = site("PLCd1_PH"), site("Btk_PH_IP4") or site("Btk_PH")
    lowest_positive = min((s["site_score"] for _, s in positives), default=None)
    highest_negative = max((s["site_score"] for _, s in negatives), default=None)
    adar_neg_max = max((s["site_score"] for _, s in negatives), default=None)
    return {
        "source": source,
        "adar2_site_in_top_3": None if adar2 is None else adar2["site_rank"] <= 3,
        "plc_and_btk_sites_in_bottom_half": (
            None if plc is None or btk is None
            else plc["site_rank_fraction"] > 0.5 and btk["site_rank_fraction"] > 0.5
        ),
        "no_overlap_all_positives_vs_negatives": (
            None if lowest_positive is None or highest_negative is None
            else lowest_positive > highest_negative
        ),
        "no_overlap_adar2_vs_negatives": (
            None if adar2 is None or adar_neg_max is None else adar2["site_score"] > adar_neg_max
        ),
        "lowest_positive_site_score": lowest_positive,
        "highest_negative_site_score": highest_negative,
        "separation_adar2_minus_best_negative": (
            None if adar2 is None or adar_neg_max is None else adar2["site_score"] - adar_neg_max
        ),
    }
