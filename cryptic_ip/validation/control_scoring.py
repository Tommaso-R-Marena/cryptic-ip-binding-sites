"""Burial-aware scoring for Phase 1 positive/negative control benchmarks."""

from __future__ import annotations

from typing import Dict, Optional

BURIAL_CLASSES = ("cryptic", "semi_cryptic", "surface", "crystal_artifact", "unknown")


#: Relative-SASA anchors for the burial component, taken from the measured
#: controls (``scripts/calibrate_controls.py``): ADAR2 1ZY7 measures 0.093 and
#: PLC-delta-1 PH 1MAI measures 0.373. Full credit is given at or below the
#: cryptic boundary and none at or above the surface anchor.
RELATIVE_BURIAL_FULL_CREDIT = 0.12
RELATIVE_BURIAL_ZERO_CREDIT = 0.40


def burial_component(ligand_sasa: float) -> float:
    """Return 1.0 for cryptic ligands and 0.0 for surface ligands.

    Absolute-SASA interface, retained for callers and stored results that only
    have Å². Prefer :func:`burial_component_relative`, which is invariant to
    ligand size and copy count.
    """
    if ligand_sasa <= 5.0:
        return 1.0
    if ligand_sasa >= 50.0:
        return 0.0
    return float((50.0 - ligand_sasa) / 45.0)


def burial_component_relative(relative_sasa: Optional[float]) -> float:
    """Burial credit from the size-normalised relative SASA.

    Absolute thresholds cannot serve here: the same 50 Å² means "mostly buried"
    for InsP6 and "fully exposed" for InsP3, and it scales with the number of
    ligand copies a crystal happened to contain. The relative measure removes
    both dependencies, and its anchors come from the measured controls rather
    than from a guess.

    Args:
        relative_sasa: Fraction of the ligand surface accessible in the complex.

    Returns:
        Credit in ``[0, 1]``; ``0.5`` (neutral) when the measure is unavailable.
    """
    if relative_sasa is None or relative_sasa != relative_sasa:
        return 0.5
    value = float(relative_sasa)
    if value <= RELATIVE_BURIAL_FULL_CREDIT:
        return 1.0
    if value >= RELATIVE_BURIAL_ZERO_CREDIT:
        return 0.0
    span = RELATIVE_BURIAL_ZERO_CREDIT - RELATIVE_BURIAL_FULL_CREDIT
    return float((RELATIVE_BURIAL_ZERO_CREDIT - value) / span)


def electrostatic_component(pocket_potential: Optional[float]) -> float:
    """Normalize pocket electrostatic potential to 0-1 (higher = more positive)."""
    if pocket_potential is None or pocket_potential != pocket_potential:
        return 0.5
    # Typical APBS sampled values for basic pockets are positive single digits.
    return float(max(0.0, min(1.0, (pocket_potential + 2.0) / 12.0)))


def cryptic_likeness(
    pocket_composite: float,
    ligand_sasa: Optional[float],
    *,
    pocket_potential: Optional[float] = None,
    use_electrostatics: bool = False,
    relative_sasa: Optional[float] = None,
) -> float:
    """Unified score: higher values indicate more cryptic/buried binding.

    When ``relative_sasa`` is supplied it takes precedence over ``ligand_sasa``,
    because it is invariant to ligand size and copy count.
    """
    if relative_sasa is not None and relative_sasa == relative_sasa:
        burial = burial_component_relative(relative_sasa)
        base = float(0.55 * burial + 0.35 * pocket_composite)
        if use_electrostatics:
            elec = electrostatic_component(pocket_potential)
            return float(0.75 * base + 0.25 * elec)
        return base
    if ligand_sasa is not None:
        burial = burial_component(ligand_sasa)
        base = float(0.55 * burial + 0.35 * pocket_composite)
        if use_electrostatics:
            elec = electrostatic_component(pocket_potential)
            return float(0.75 * base + 0.25 * elec)
        return base
    return float(pocket_composite)


def validation_score(
    control_type: str,
    pocket_composite: float,
    ligand_sasa: Optional[float],
    basic_residues: int,
    *,
    pocket_potential: Optional[float] = None,
    use_electrostatics: bool = False,
    relative_sasa: Optional[float] = None,
) -> float:
    """Return cryptic-likeness on a single scale for separation analysis."""
    _ = control_type, basic_residues
    return cryptic_likeness(
        pocket_composite,
        ligand_sasa,
        pocket_potential=pocket_potential,
        use_electrostatics=use_electrostatics,
        relative_sasa=relative_sasa,
    )


def positive_passed(
    validation_value: float,
    ligand_sasa: Optional[float],
    basic_residues: int,
    pocket_composite: float,
    *,
    burial_class: str = "cryptic",
    score_threshold: float = 0.55,
    relative_sasa: Optional[float] = None,
) -> bool:
    """Phase 1 positive-control pass criteria.

    When ``relative_sasa`` is available the burial gate uses it, since the
    absolute-SASA cutoffs below were calibrated against a copy-summed measure
    that no longer exists. On the measured controls a buried site clears the
    cryptic boundary (ADAR2 0.093) while a surface site does not (PLC-delta-1
    0.373).
    """
    if burial_class == "crystal_artifact":
        return False
    if relative_sasa is not None and relative_sasa == relative_sasa:
        if burial_class == "semi_cryptic":
            return (
                basic_residues >= 4
                and pocket_composite >= 0.45
                and validation_value >= 0.40
            )
        return (
            float(relative_sasa) <= RELATIVE_BURIAL_FULL_CREDIT
            and basic_residues >= 4
            and pocket_composite >= 0.50
            and validation_value >= score_threshold
        )
    if burial_class == "semi_cryptic":
        return (
            ligand_sasa is not None
            and 10.0 <= ligand_sasa < 50.0
            and basic_residues >= 4
            and pocket_composite >= 0.45
            and validation_value >= 0.40
        )
    if ligand_sasa is not None:
        return (
            ligand_sasa < 10.0
            and basic_residues >= 4
            and pocket_composite >= 0.50
            and validation_value >= score_threshold
        )
    return basic_residues >= 3 and validation_value >= score_threshold


def negative_passed(
    validation_value: float,
    ligand_sasa: Optional[float],
    pocket_composite: float,
    *,
    score_threshold: float = 0.45,
    decoy_mode: bool = False,
    relative_sasa: Optional[float] = None,
) -> bool:
    """Phase 1 negative-control pass criteria (low cryptic-likeness)."""
    threshold = 0.55 if decoy_mode else score_threshold
    if relative_sasa is not None and relative_sasa == relative_sasa:
        # A surface ligand must not reach the cryptic-likeness of a buried one.
        return validation_value <= threshold
    if ligand_sasa is not None:
        if ligand_sasa >= 50.0:
            return validation_value <= threshold
        if ligand_sasa >= 20.0:
            return validation_value <= 0.65
        return validation_value <= threshold and pocket_composite < 0.55
    return validation_value <= threshold


def separation_quality(
    positive_scores,
    negative_scores,
    *,
    exclude_artifact_positives: bool = True,
    positive_classes: Optional[list] = None,
) -> Dict[str, object]:
    """Summarize whether positive and negative controls separate cleanly."""
    pos_scores = [float(score) for score in positive_scores if score == score]
    neg_scores = [float(score) for score in negative_scores if score == score]

    if positive_classes and exclude_artifact_positives:
        filtered = [
            score
            for score, cls in zip(positive_scores, positive_classes)
            if score == score and cls not in {"crystal_artifact"}
        ]
        if filtered:
            pos_scores = [float(s) for s in filtered]

    if not pos_scores or not neg_scores:
        return {}

    pos_mean = float(sum(pos_scores) / len(pos_scores))
    neg_mean = float(sum(neg_scores) / len(neg_scores))
    separation = pos_mean - neg_mean

    tier1_pos = [score for score in pos_scores[:1]]
    tier1_neg = [score for score in neg_scores[:1]]
    tier1_separation = float(tier1_pos[0] - tier1_neg[0]) if tier1_pos and tier1_neg else separation

    return {
        "positive_mean": pos_mean,
        "negative_mean": neg_mean,
        "separation": separation,
        "tier1_separation": tier1_separation,
        "clear_separation": separation > 0.30,
        "tier1_gate_passed": tier1_separation > 0.50,
        "phase1_ready": tier1_separation > 0.50,
    }
