"""Registry and runtime discovery of inositol phosphate chemical components.

Rationale
---------
Earlier versions of the pipeline hard-coded eight PDB chemical component
identifiers (``IP3``, ``IP4``, ``IP5``, ``IP6``, ``IHP``, ``I3P``, ``4IP``,
``6A0``). That approach has two defects that matter for a systematic screen:

1. **Incompleteness.** The PDB chemical component dictionary contains many more
   inositol phosphate species than any hand-written list, including phosphate
   regioisomers, inositol pyrophosphates (InsP7/InsP8) and non-*myo* inositol
   stereoisomers. Anything missing from the list is silently absent from the
   ground-truth dataset, which biases every downstream benchmark.
2. **Unverifiability.** A mistyped identifier produces no error; it simply
   contributes nothing, so the dataset silently shrinks.

This module therefore treats the ligand vocabulary as *data to be discovered and
validated at run time*, not as a constant:

* :func:`discover_ip_ligands` performs a full-text search of the RCSB chemical
  component dictionary and keeps every component whose name and formula are
  consistent with an inositol phosphate.
* Every candidate identifier - discovered or seeded - is resolved against the
  chemical component API, so unknown identifiers are dropped with a warning
  rather than silently ignored.
* The inositol phosphate series (InsP1 ... InsP8) is derived by counting
  phosphorus atoms in the *reported formula*, never assumed from the identifier.

The curated seed identifiers below only widen the search; they are not trusted.
Each is validated before use, and :attr:`IPLigand.verified` records whether the
entry was confirmed against the chemical component dictionary.
"""

from __future__ import annotations

import logging
import re
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

LOGGER = logging.getLogger(__name__)

#: Identifiers used to seed discovery. These are *candidates only*: each one is
#: resolved against the RCSB chemical component dictionary before use, and any
#: identifier that does not exist, or whose formula is not consistent with an
#: inositol phosphate, is discarded. Seeding matters because full-text search
#: depends on curator-entered names, which are not perfectly uniform.
SEED_COMP_IDS: Tuple[str, ...] = (
    "IHP",  # myo-inositol hexakisphosphate (InsP6); ligand of ADAR2 (1ZY7)
    "I3P",  # D-myo-inositol 1,4,5-trisphosphate
    "4IP",  # inositol tetrakisphosphate
    "2IP",
    "5IP",
    "6IP",
    "IP3",
    "IP4",
    "IP5",
    "IP6",
    "I6P",
    "INS",  # myo-inositol (unphosphorylated; retained as a negative reference)
    "6A0",
    "MYI",
    "IPD",
)

#: Matches any name containing an inositol core. Deliberately not anchored on
#: word boundaries: curators name components in many ways ("myo-inositol",
#: "scyllo-inositol", "2-O-methyl-inositol", "phosphatidylinositol"), and
#: requiring a boundary would silently miss whole families. Discrimination is
#: left to the formula constraints and :data:`LIPID_NAME_PATTERN`, which is where
#: it can be applied precisely.
INOSITOL_NAME_PATTERN = re.compile(r"inositol", re.IGNORECASE)

#: Lipid-linked inositols (phosphatidylinositol phosphates) bind at membranes
#: and lipid-binding modules rather than forming buried structural cofactor
#: sites. They are excluded by default but can be requested explicitly.
LIPID_NAME_PATTERN = re.compile(
    r"phosphatidyl|diacyl|dioctanoyl|dibutanoyl|glycerol|ceramide|acyl", re.IGNORECASE
)

#: Upper bound on inositol-core carbon count. The inositol ring is C6; small
#: substituted derivatives are allowed, lipids are not.
MAX_INOSITOL_CARBONS = 12

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)\s*(\d*)")


def parse_formula(formula: Optional[str]) -> Dict[str, int]:
    """Parse a PDB chemical component formula into element counts.

    PDB formulae are space-separated element/count pairs with an optional
    trailing charge, e.g. ``"C6 H18 O24 P6"`` or ``"C6 H6 O24 P6 12-"``.

    Args:
        formula: Formula string as reported by the chemical component dictionary.

    Returns:
        Mapping from element symbol to atom count. Unparseable input yields ``{}``.
    """
    if not formula:
        return {}
    counts: Dict[str, int] = {}
    for token in str(formula).split():
        # Skip charge tokens such as "12-", "2+", "-".
        if re.fullmatch(r"\d*[+-]", token):
            continue
        match = _FORMULA_TOKEN.fullmatch(token.strip())
        if not match:
            continue
        element, digits = match.group(1), match.group(2)
        counts[element] = counts.get(element, 0) + (int(digits) if digits else 1)
    return counts


def ip_series_label(n_phosphorus: int) -> str:
    """Return the inositol phosphate series label for a phosphorus count.

    Counting phosphorus atoms in the formula is stereochemistry-agnostic and
    therefore robust: ``InsP6`` has six phosphorus atoms whether the component is
    named ``IHP`` or something else. Components with more than six phosphorus
    atoms are inositol pyrophosphates (InsP7/InsP8 and beyond).

    Args:
        n_phosphorus: Number of phosphorus atoms in the component.

    Returns:
        A label such as ``"InsP0"``, ``"InsP6"`` or ``"InsP8"``.
    """
    if n_phosphorus < 0:
        raise ValueError("n_phosphorus must be non-negative")
    return f"InsP{int(n_phosphorus)}"


@dataclass(frozen=True)
class IPLigand:
    """A validated inositol phosphate chemical component.

    Attributes:
        comp_id: PDB chemical component identifier (e.g. ``"IHP"``).
        name: Curator-assigned chemical name.
        formula: Reported formula string.
        n_phosphorus: Phosphorus atom count parsed from ``formula``.
        n_carbon: Carbon atom count parsed from ``formula``.
        n_heavy_atoms: Non-hydrogen atom count parsed from ``formula``.
        series: Series label from :func:`ip_series_label`.
        is_phosphorylated: ``True`` when ``n_phosphorus >= 1``.
        is_lipid_linked: ``True`` for phosphatidylinositol-type components.
        formula_weight: Formula weight in Da when reported.
        verified: ``True`` when resolved against the chemical component API.
        sources: Provenance tags describing how the component was found.
    """

    comp_id: str
    name: str
    formula: str
    n_phosphorus: int
    n_carbon: int
    n_heavy_atoms: int
    series: str
    is_phosphorylated: bool
    is_lipid_linked: bool
    formula_weight: Optional[float] = None
    verified: bool = False
    sources: Tuple[str, ...] = field(default_factory=tuple)

    def to_dict(self) -> Dict[str, Any]:
        """Return a JSON-serialisable representation for provenance manifests."""
        payload = asdict(self)
        payload["sources"] = list(self.sources)
        return payload


def looks_like_inositol_phosphate(
    name: Optional[str],
    formula: Optional[str],
    *,
    allow_lipid: bool = False,
    require_phosphate: bool = True,
) -> bool:
    """Decide whether a chemical component is an inositol phosphate.

    The test combines a name check with formula constraints so that neither
    signal alone can admit a false positive:

    * the name must mention an inositol stereoisomer;
    * the formula must contain phosphorus (unless ``require_phosphate`` is off);
    * the carbon count must be consistent with an inositol core rather than a
      lipid-linked species.

    Args:
        name: Chemical component name.
        formula: Chemical component formula.
        allow_lipid: Accept phosphatidylinositol-type components.
        require_phosphate: Require at least one phosphorus atom.

    Returns:
        ``True`` when the component should be treated as an inositol phosphate.
    """
    if not name or not INOSITOL_NAME_PATTERN.search(name):
        return False
    if not allow_lipid and LIPID_NAME_PATTERN.search(name):
        return False

    counts = parse_formula(formula)
    if not counts:
        return False
    n_p = counts.get("P", 0)
    if require_phosphate and n_p < 1:
        return False

    n_c = counts.get("C", 0)
    if n_c < 6:
        return False
    if not allow_lipid and n_c > MAX_INOSITOL_CARBONS:
        return False
    return True


def build_ligand(
    comp_id: str,
    payload: Mapping[str, Any],
    *,
    sources: Sequence[str] = (),
    verified: bool = True,
) -> Optional[IPLigand]:
    """Construct an :class:`IPLigand` from an RCSB ``chemcomp`` payload.

    Args:
        comp_id: Chemical component identifier.
        payload: ``chem_comp`` sub-document from the RCSB chemical component API.
        sources: Provenance tags to record.
        verified: Whether the payload came from the live API.

    Returns:
        The ligand record, or ``None`` when the payload is not inositol-like.
    """
    chem_comp = payload.get("chem_comp", payload) or {}
    name = chem_comp.get("name") or ""
    formula = chem_comp.get("formula") or ""
    counts = parse_formula(formula)
    if not counts:
        LOGGER.debug("Component %s has unparseable formula %r", comp_id, formula)
        return None

    is_lipid = bool(LIPID_NAME_PATTERN.search(name))
    if not looks_like_inositol_phosphate(
        name, formula, allow_lipid=is_lipid, require_phosphate=False
    ):
        return None

    n_p = counts.get("P", 0)
    heavy = sum(count for element, count in counts.items() if element not in {"H", "D"})
    weight = chem_comp.get("formula_weight")
    return IPLigand(
        comp_id=comp_id.upper(),
        name=name,
        formula=formula,
        n_phosphorus=n_p,
        n_carbon=counts.get("C", 0),
        n_heavy_atoms=heavy,
        series=ip_series_label(n_p),
        is_phosphorylated=n_p >= 1,
        is_lipid_linked=is_lipid,
        formula_weight=float(weight) if weight is not None else None,
        verified=verified,
        sources=tuple(sources),
    )


def discover_ip_ligands(
    client: Any,
    *,
    search_terms: Sequence[str] = ("inositol phosphate", "inositol", "inositol hexakisphosphate"),
    seed_comp_ids: Sequence[str] = SEED_COMP_IDS,
    allow_lipid: bool = False,
    include_unphosphorylated: bool = False,
    max_candidates: int = 2000,
) -> Tuple[List[IPLigand], Dict[str, Any]]:
    """Discover and validate every inositol phosphate component in the PDB.

    Args:
        client: An object exposing ``search_chemcomp_full_text(term, rows)`` and
            ``fetch_chemcomp(comp_id)`` - see
            :class:`cryptic_ip.database.rcsb_client.RcsbClient`.
        search_terms: Full-text queries issued against the component dictionary.
        seed_comp_ids: Extra candidate identifiers to resolve.
        allow_lipid: Keep phosphatidylinositol-type components.
        include_unphosphorylated: Keep free inositols (InsP0) as references.
        max_candidates: Safety bound on the number of components resolved.

    Returns:
        ``(ligands, provenance)`` where ``ligands`` is sorted by descending
        phosphorylation state then identifier, and ``provenance`` records the
        queries issued and the identifiers accepted and rejected.
    """
    candidates: Dict[str, Set[str]] = {}

    def note(comp_id: str, source: str) -> None:
        candidates.setdefault(comp_id.upper(), set()).add(source)

    for term in search_terms:
        try:
            hits = client.search_chemcomp_full_text(term)
        except Exception as exc:  # noqa: BLE001 - network/API failures must not abort
            LOGGER.warning("Chemical component search failed for %r: %s", term, exc)
            continue
        for comp_id in hits:
            note(comp_id, f"full_text:{term}")

    for comp_id in seed_comp_ids:
        note(comp_id, "seed")

    accepted: List[IPLigand] = []
    rejected: Dict[str, str] = {}
    for comp_id in sorted(candidates)[:max_candidates]:
        try:
            payload = client.fetch_chemcomp(comp_id)
        except Exception as exc:  # noqa: BLE001
            rejected[comp_id] = f"lookup_failed: {exc}"
            continue
        if not payload:
            rejected[comp_id] = "not_found"
            continue

        ligand = build_ligand(comp_id, payload, sources=sorted(candidates[comp_id]))
        if ligand is None:
            rejected[comp_id] = "not_inositol_like"
            continue
        if ligand.is_lipid_linked and not allow_lipid:
            rejected[comp_id] = "lipid_linked"
            continue
        if not ligand.is_phosphorylated and not include_unphosphorylated:
            rejected[comp_id] = "unphosphorylated"
            continue
        accepted.append(ligand)

    accepted.sort(key=lambda lig: (-lig.n_phosphorus, lig.comp_id))
    provenance = {
        "search_terms": list(search_terms),
        "seed_comp_ids": list(seed_comp_ids),
        "n_candidates": len(candidates),
        "n_accepted": len(accepted),
        "accepted_comp_ids": [lig.comp_id for lig in accepted],
        "rejected": rejected,
        "allow_lipid": allow_lipid,
        "include_unphosphorylated": include_unphosphorylated,
    }
    LOGGER.info(
        "Resolved %d inositol phosphate components from %d candidates",
        len(accepted),
        len(candidates),
    )
    return accepted, provenance


def comp_ids(ligands: Iterable[IPLigand]) -> Tuple[str, ...]:
    """Return the sorted, upper-case identifiers of a ligand collection."""
    return tuple(sorted({lig.comp_id.upper() for lig in ligands}))


def phosphorylated_comp_ids(ligands: Iterable[IPLigand]) -> Tuple[str, ...]:
    """Return identifiers of phosphorylated components only."""
    return tuple(sorted({lig.comp_id.upper() for lig in ligands if lig.is_phosphorylated}))
