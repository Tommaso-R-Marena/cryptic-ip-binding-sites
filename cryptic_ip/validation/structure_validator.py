"""Structure-level validation for PDB/CIF inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import mean
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from Bio.PDB import MMCIFParser, PDBParser, is_aa


@dataclass
class ValidationIssue:
    """Single validation issue with guidance."""

    level: str
    check: str
    message: str
    suggested_fix: str


@dataclass
class ValidationReport:
    """Structured validation report."""

    valid: bool
    issues: List[ValidationIssue]
    metrics: Dict[str, float]


class StructureValidator:
    """Validate structural quality and coordinate sanity for protein structures."""

    def __init__(
        self,
        min_resolution: float = 4.0,
        min_completeness: float = 0.7,
        min_chain_length: int = 50,
        low_confidence_threshold: float = 70.0,
        origin_tolerance: float = 1e-6,
        max_reasonable_b_factor: float = 300.0,
        pocket_residue_ids: Optional[Sequence[int]] = None,
    ) -> None:
        self.min_resolution = min_resolution
        self.min_completeness = min_completeness
        self.min_chain_length = min_chain_length
        self.low_confidence_threshold = low_confidence_threshold
        self.origin_tolerance = origin_tolerance
        self.max_reasonable_b_factor = max_reasonable_b_factor
        self.pocket_residue_ids = set(pocket_residue_ids or [])

    def validate(self, structure_path: str) -> ValidationReport:
        """Run all structure checks for a PDB/mmCIF file."""
        structure_file = Path(structure_path)
        issues: List[ValidationIssue] = []
        metrics: Dict[str, float] = {}

        structure = self._parse_structure(structure_file, issues)
        if structure is None:
            return ValidationReport(valid=False, issues=issues, metrics=metrics)

        model = next(structure.get_models(), None)
        if model is None:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="model_presence",
                    message="Structure contains no models.",
                    suggested_fix="Use a full coordinate file with at least one model.",
                )
            )
            return ValidationReport(valid=False, issues=issues, metrics=metrics)

        residues = [res for res in model.get_residues() if is_aa(res, standard=True)]
        atoms = list(model.get_atoms())
        metrics["residue_count"] = float(len(residues))
        metrics["atom_count"] = float(len(atoms))

        self._check_resolution(structure_file, issues, metrics)
        self._check_completeness(residues, issues, metrics)
        self._check_chain_lengths(model, issues, metrics)
        self._check_missing_pocket_residues(residues, issues)
        self._check_alphafold_confidence(atoms, issues, metrics)
        self._check_coordinates_and_bfactors(atoms, issues, metrics)

        is_valid = not any(issue.level == "error" for issue in issues)
        return ValidationReport(valid=is_valid, issues=issues, metrics=metrics)

    def _parse_structure(self, structure_file: Path, issues: List[ValidationIssue]):
        suffix = structure_file.suffix.lower()
        try:
            if suffix == ".cif":
                parser = MMCIFParser(QUIET=True)
            else:
                parser = PDBParser(QUIET=True)
            return parser.get_structure(structure_file.stem, str(structure_file))
        except Exception as exc:  # BioPython parser exceptions vary
            issues.append(
                ValidationIssue(
                    level="error",
                    check="readable_by_biopython",
                    message=f"BioPython failed to parse {structure_file.name}: {exc}",
                    suggested_fix="Re-download the structure file or convert it to valid PDB/mmCIF format.",
                )
            )
            return None

    def _check_resolution(
        self, structure_file: Path, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        resolution = None
        if structure_file.suffix.lower() != ".cif":
            try:
                with structure_file.open("r", encoding="utf-8", errors="ignore") as handle:
                    for line in handle:
                        if line.startswith("REMARK   2 RESOLUTION."):
                            parts = line.replace("ANGSTROMS.", "").split()
                            for token in parts:
                                try:
                                    resolution = float(token)
                                    break
                                except ValueError:
                                    continue
                            if resolution is not None:
                                break
            except OSError:
                pass

        if resolution is not None:
            metrics["resolution"] = resolution
            if resolution > self.min_resolution:
                issues.append(
                    ValidationIssue(
                        level="warning",
                        check="resolution",
                        message=f"Low structural resolution: {resolution:.2f} Å (threshold {self.min_resolution:.2f} Å).",
                        suggested_fix="Prefer higher-resolution crystal structures when available.",
                    )
                )

    def _check_completeness(
        self, residues, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        seq_ids = sorted({res.id[1] for res in residues})
        if not seq_ids:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="completeness",
                    message="No amino-acid residues found in the structure.",
                    suggested_fix="Use a protein structure (not ligand-only or empty model).",
                )
            )
            return

        expected = seq_ids[-1] - seq_ids[0] + 1
        completeness = len(seq_ids) / expected if expected else 0.0
        metrics["completeness"] = completeness
        if completeness < self.min_completeness:
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="completeness",
                    message=f"Structure completeness is {completeness:.2%}, below threshold {self.min_completeness:.0%}.",
                    suggested_fix="Choose a structure with fewer unresolved regions or missing loops.",
                )
            )

    def _check_chain_lengths(
        self, model, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        chain_lengths = {}
        for chain in model:
            chain_res = [res for res in chain.get_residues() if is_aa(res, standard=True)]
            chain_lengths[chain.id] = len(chain_res)
            if 0 < len(chain_res) < self.min_chain_length:
                issues.append(
                    ValidationIssue(
                        level="warning",
                        check="chain_length",
                        message=f"Chain {chain.id} has only {len(chain_res)} residues (< {self.min_chain_length}).",
                        suggested_fix="Filter out short fragments/chains or use complete chain coordinates.",
                    )
                )
        if chain_lengths:
            metrics["min_chain_length"] = float(min(chain_lengths.values()))
            metrics["max_chain_length"] = float(max(chain_lengths.values()))

    def _check_missing_pocket_residues(self, residues, issues: List[ValidationIssue]) -> None:
        if not self.pocket_residue_ids:
            return
        present = {res.id[1] for res in residues}
        missing = sorted(self.pocket_residue_ids - present)
        if missing:
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="predicted_pocket_missing_residues",
                    message=f"Missing residues in predicted pocket region: {missing[:10]}",
                    suggested_fix="Use an alternative structure/model that covers the predicted pocket residues.",
                )
            )

    def _check_alphafold_confidence(
        self, atoms, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        if not atoms:
            return
        b_values = [atom.bfactor for atom in atoms if atom.bfactor is not None]
        if not b_values:
            return

        low_conf = sum(1 for value in b_values if value < self.low_confidence_threshold)
        low_fraction = low_conf / len(b_values)
        metrics["low_confidence_fraction"] = low_fraction
        metrics["mean_b_factor"] = mean(b_values)

        if low_fraction > 0.2:
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="alphafold_low_confidence",
                    message=(
                        f"{low_fraction:.1%} of atoms have pLDDT/B-factor below "
                        f"{self.low_confidence_threshold:.0f}."
                    ),
                    suggested_fix="Inspect low-confidence regions before interpreting pocket features.",
                )
            )

    def _check_coordinates_and_bfactors(
        self, atoms, issues: List[ValidationIssue], metrics: Dict[str, float]
    ) -> None:
        if not atoms:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="coordinates_present",
                    message="No atoms found in structure model.",
                    suggested_fix="Confirm coordinate section exists in PDB/mmCIF file.",
                )
            )
            return

        origin_atoms = 0
        unreasonable_b = 0
        for atom in atoms:
            x, y, z = atom.coord
            if (
                abs(x) < self.origin_tolerance
                and abs(y) < self.origin_tolerance
                and abs(z) < self.origin_tolerance
            ):
                origin_atoms += 1
            if atom.bfactor is not None and atom.bfactor > self.max_reasonable_b_factor:
                unreasonable_b += 1

        metrics["origin_atom_fraction"] = origin_atoms / len(atoms)
        metrics["unreasonable_b_factor_fraction"] = unreasonable_b / len(atoms)

        if origin_atoms > 0:
            issues.append(
                ValidationIssue(
                    level="error",
                    check="atom_origin",
                    message=f"Found {origin_atoms} atoms at (0,0,0), indicating likely coordinate corruption.",
                    suggested_fix="Re-download the file or regenerate structure with correct coordinate transforms.",
                )
            )

        if unreasonable_b > 0:
            issues.append(
                ValidationIssue(
                    level="warning",
                    check="b_factor_sanity",
                    message=f"Found {unreasonable_b} atoms with unusually high B-factors (> {self.max_reasonable_b_factor}).",
                    suggested_fix="Check for parsing errors or use cleaned structures before analysis.",
                )
            )
