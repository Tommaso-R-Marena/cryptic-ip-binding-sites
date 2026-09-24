"""Phosphate-dense ligands for the transfer dataset (docs/TRANSFER_PLAN.md).

The class is defined by a rule on the component's formula, not by a list: at
least :data:`MIN_PHOSPHORUS` phosphorus atoms, at least
:data:`MIN_PHOSPHORUS_PER_HEAVY_ATOM` phosphorus per heavy atom, not an inositol
phosphate, not lipid-linked. The seed identifiers below only propose
candidates; each is resolved against the chemical component dictionary and
kept only if its formula passes the rule, so a mistyped or unexpected
identifier is rejected on the record rather than silently included.
"""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

from cryptic_ip.database.ip_ligands import (
    LIPID_NAME_PATTERN,
    looks_like_inositol_phosphate,
    parse_formula,
)

LOGGER = logging.getLogger(__name__)

MIN_PHOSPHORUS = 2
MIN_PHOSPHORUS_PER_HEAVY_ATOM = 0.07

#: Candidates only (nucleotide di/triphosphates and analogues, phosphoribosyl
#: and sugar bisphosphates, isoprenoid and inorganic pyrophosphates, thiamine
#: diphosphate, adenosine bisphosphates). Resolved and filtered at run time.
SEED_COMP_IDS: Tuple[str, ...] = (
    "ATP", "ADP", "ANP", "ACP", "AGS", "APC", "GTP", "GDP", "GNP", "GCP", "GSP",
    "UTP", "UDP", "CTP", "CDP", "TTP", "DTP", "DGT", "DCP", "DUT", "DGP", "TYD",
    "TPP", "PRP", "FBP", "BPG", "DG2", "13P", "POP", "PPV", "DPO", "PPK", "PI",
    "IPE", "DMA", "GPP", "FPP", "GRG", "A3P", "PAP", "PPS", "5GP", "ADX",
)


@dataclass(frozen=True)
class PolyanionLigand:
    """A component that passes the class rule."""

    comp_id: str
    name: str
    formula: str
    n_phosphorus: int
    n_heavy_atoms: int
    phosphorus_per_heavy_atom: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def classify(comp_id: str, payload: Mapping[str, Any]) -> Tuple[Optional[PolyanionLigand], str]:
    """Apply the class rule to a chemical component payload.

    Returns:
        ``(ligand, reason)``: the ligand when it passes, else ``None`` and why not.
    """
    chem_comp = payload.get("chem_comp", payload) or {}
    name = str(chem_comp.get("name") or "")
    formula = str(chem_comp.get("formula") or "")
    counts = parse_formula(formula)
    if not counts:
        return None, f"unparseable formula {formula!r}"
    heavy = sum(n for element, n in counts.items() if element not in {"H", "D"})
    n_p = counts.get("P", 0)
    ratio = n_p / heavy if heavy else 0.0
    if n_p < MIN_PHOSPHORUS:
        return None, f"{n_p} phosphorus atoms"
    if ratio < MIN_PHOSPHORUS_PER_HEAVY_ATOM:
        return None, f"phosphorus per heavy atom {ratio:.3f}"
    if looks_like_inositol_phosphate(name, formula, allow_lipid=True, require_phosphate=False):
        return None, "inositol phosphate"
    if LIPID_NAME_PATTERN.search(name):
        return None, "lipid-linked"
    return PolyanionLigand(comp_id.upper(), name, formula, n_p, heavy, round(ratio, 4)), "accepted"


def resolve_polyanion_ligands(
    client: Any, seed_comp_ids: Sequence[str] = SEED_COMP_IDS
) -> Tuple[List[PolyanionLigand], Dict[str, Any]]:
    """Resolve the seeds against the component dictionary and apply the rule.

    Args:
        client: Exposes ``fetch_chemcomp(comp_id)``.

    Returns:
        ``(ligands, provenance)``; provenance records every decision.
    """
    ligands: List[PolyanionLigand] = []
    decisions: Dict[str, str] = {}
    for comp_id in dict.fromkeys(c.upper() for c in seed_comp_ids):
        payload = client.fetch_chemcomp(comp_id)
        if not payload:
            decisions[comp_id] = "unknown identifier"
            continue
        ligand, reason = classify(comp_id, payload)
        decisions[comp_id] = reason
        if ligand is not None:
            ligands.append(ligand)
    LOGGER.info("Class ligands: %s", ", ".join(lig.comp_id for lig in ligands))
    rejected = {k: v for k, v in decisions.items() if v != "accepted"}
    if rejected:
        LOGGER.info("Rejected: %s", rejected)
    provenance = {
        "rule": {"min_phosphorus": MIN_PHOSPHORUS, "min_phosphorus_per_heavy_atom": MIN_PHOSPHORUS_PER_HEAVY_ATOM},
        "decisions": decisions,
    }
    return sorted(ligands, key=lambda lig: lig.comp_id), provenance
