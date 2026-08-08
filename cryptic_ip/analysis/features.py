"""Physically interpretable pocket descriptors for cryptic IP-site prediction.

Why the feature set was rebuilt
-------------------------------
The original model used six features, and two of them were broken in ways that
made learning impossible:

* ``pocket_depth`` was populated from fpocket's *mean local hydrophobic density*,
  a composition statistic that has nothing to do with depth. The pipeline
  separately computed a genuine geometric burial depth and then discarded it.
* ``electrostatic_potential`` required APBS. When APBS was unavailable - the
  common case, and the case in the shipped training run - the column was entirely
  missing, so a sixth of the feature space was constant ``NaN``.

Beyond those defects, six features cannot express the hypothesis. A cryptic
inositol phosphate site is defined by a *conjunction*: a cavity of the right size
and shape, enclosed rather than merely concave, lined by a cluster of basic side
chains whose nitrogens point inward, electrostatically positive, and confidently
modelled. This module computes 40 descriptors organised in seven blocks:

=====================  ===========================================================
Block                  What it captures
=====================  ===========================================================
Geometry               volume (fpocket and convex hull), radius of gyration,
                       asphericity, extent, alpha-sphere statistics
Burial                 geometric burial depth, enclosure, buried-residue fraction
Accessibility          absolute and relative residue SASA (mean/median/min/max)
Composition            basic, acidic, aromatic, hydroxyl, polar fractions;
                       Kyte-Doolittle hydropathy
Charge geometry        coordinating nitrogen count and their distance
                       distribution; charge density; net formal charge
Electrostatics         screened Coulomb potential at the pocket centre (always
                       available) plus the APBS value when computed
Confidence             pLDDT mean/min/fraction above cutoff for AlphaFold models
=====================  ===========================================================

Every descriptor is computed from chain-aware residue keys, so multi-chain
structures are handled correctly, and every one is documented with the physical
reason it belongs in the set. Features are returned in a stable order so that a
serialised model and a fresh extraction can never disagree about column meaning.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from .geometry import (
    DEFAULT_N_POINTS,
    burial_depth,
    convex_hull_volume,
    enclosure_fraction,
    point_cloud_shape,
    shrake_rupley_sasa,
)
from .structure_arrays import (
    ACIDIC_RESIDUES,
    AROMATIC_RESIDUES,
    BASIC_NITROGEN_ATOMS,
    BASIC_RESIDUES,
    FORMAL_CHARGES,
    HYDROXYL_RESIDUES,
    KYTE_DOOLITTLE,
    MAX_RESIDUE_SASA,
    STRONG_BASIC_RESIDUES,
    ResidueKey,
    StructureArrays,
)

LOGGER = logging.getLogger(__name__)

#: Radius around the pocket centre used to define the pocket-lining shell (Å).
#: Chosen to match the reach of an inositol phosphate plus its first coordination
#: shell: InsP6 spans roughly 11 Å, so residues within 8 Å of the centre are the
#: ones that could contact it.
POCKET_SHELL_RADIUS = 8.0

#: Tighter radius used for the "inner" coordination shell (Å).
POCKET_CORE_RADIUS = 5.0

#: Relative solvent accessibility below which a residue counts as buried.
BURIED_RSA_THRESHOLD = 0.20

#: Debye screening length at physiological ionic strength (Å). At 150 mM
#: monovalent salt the Debye length is about 8 Å; the screened Coulomb sum below
#: uses it so that distant charges contribute negligibly, as they do in solution.
DEBYE_LENGTH = 8.0

#: Relative permittivity used in the screened Coulomb surrogate. A value between
#: bulk water (78) and protein interior (4) reflects that a partially buried
#: pocket is neither; 20 is a standard compromise for continuum estimates of
#: protein electrostatics.
RELATIVE_PERMITTIVITY = 20.0

#: Conversion so the screened Coulomb sum is reported in units of kT/e at 298 K,
#: matching the APBS convention: e^2 / (4 pi eps_0 k_B T) = 560.7 Å for eps_r = 1.
COULOMB_KT_PREFACTOR_ANGSTROM = 560.7

#: pLDDT threshold for AlphaFold confidence gating.
PLDDT_CONFIDENCE_CUTOFF = 70.0

#: Canonical feature order. A stable order is essential: a serialised model
#: indexes columns positionally, so silently reordering them would produce
#: confident nonsense rather than an error.
FEATURE_NAMES: Tuple[str, ...] = (
    # --- geometry
    "pocket_volume",
    "hull_volume",
    "n_alpha_spheres",
    "alpha_sphere_density",
    "radius_of_gyration",
    "asphericity",
    "max_extent",
    # --- burial
    "burial_depth",
    "enclosure",
    "buried_residue_fraction",
    "mean_relative_sasa",
    # --- accessibility
    "sasa_mean",
    "sasa_median",
    "sasa_min",
    "sasa_max",
    "sasa_total",
    # --- composition (shell)
    "n_residues",
    "n_basic_residues",
    "n_strong_basic_residues",
    "n_acidic_residues",
    "n_aromatic_residues",
    "n_hydroxyl_residues",
    "basic_fraction",
    "acidic_fraction",
    "aromatic_fraction",
    "hydroxyl_fraction",
    "polar_fraction",
    "hydropathy_mean",
    # --- composition (core shell)
    "n_basic_residues_core",
    "n_acidic_residues_core",
    # --- charge geometry
    "n_basic_nitrogens",
    "basic_nitrogen_min_distance",
    "basic_nitrogen_mean_distance",
    "basic_nitrogen_dispersion",
    "net_formal_charge",
    "positive_charge_density",
    "charge_balance",
    # --- electrostatics
    "coulomb_potential_kt",
    "electrostatic_potential",
    # --- confidence
    "plddt_mean",
    "plddt_min",
    "plddt_fraction_above_cutoff",
)

#: Legacy six-feature schema kept so that models serialised by earlier versions
#: can still be loaded and scored.
LEGACY_FEATURE_NAMES: Tuple[str, ...] = (
    "pocket_depth",
    "sasa",
    "electrostatic_potential",
    "n_basic_residues",
    "pocket_volume",
    "plddt_confidence",
)


@dataclass
class PocketFeatures:
    """Descriptors for one pocket, with identifying metadata.

    Attributes:
        pocket_id: Pocket identifier from the detector.
        centre: Pocket centre coordinate.
        features: Descriptor name to value; ``nan`` marks a value that could not
            be computed rather than a real zero.
        residue_keys: Chain-aware keys of pocket-lining residues.
    """

    pocket_id: int
    centre: Tuple[float, float, float]
    features: Dict[str, float]
    residue_keys: List[ResidueKey]

    def vector(self, names: Sequence[str] = FEATURE_NAMES) -> np.ndarray:
        """Return the descriptor values in the requested order.

        Args:
            names: Feature names to extract.

        Returns:
            A float array; missing descriptors become ``nan``.
        """
        return np.asarray([self.features.get(name, np.nan) for name in names], dtype=float)

    def to_row(self) -> Dict[str, float]:
        """Return a flat mapping suitable for a DataFrame row."""
        row: Dict[str, float] = {"pocket_id": self.pocket_id}
        row.update({name: self.features.get(name, np.nan) for name in FEATURE_NAMES})
        row["center_x"], row["center_y"], row["center_z"] = self.centre
        return row


class PocketFeatureExtractor:
    """Compute the descriptor suite for pockets of one structure.

    SASA over the whole structure is computed once and reused for every pocket,
    which is what makes a proteome-scale screen tractable: the per-atom
    calculation dominates the cost and does not depend on the pocket.

    Args:
        arrays: Parsed structure arrays.
        n_points: SASA sample points per atom.
        shell_radius: Radius defining pocket-lining residues (Å).
        core_radius: Radius defining the inner coordination shell (Å).
    """

    def __init__(
        self,
        arrays: StructureArrays,
        *,
        n_points: int = DEFAULT_N_POINTS,
        shell_radius: float = POCKET_SHELL_RADIUS,
        core_radius: float = POCKET_CORE_RADIUS,
    ) -> None:
        self.arrays = arrays
        self.n_points = int(n_points)
        self.shell_radius = float(shell_radius)
        self.core_radius = float(core_radius)

        self._context_indices = np.flatnonzero(~arrays.is_solvent)
        self._context_coords = arrays.coords[self._context_indices]
        self._context_radii = arrays.radii[self._context_indices]
        self._atom_sasa = shrake_rupley_sasa(
            self._context_coords, self._context_radii, n_points=self.n_points
        )
        self._protein_positions = np.flatnonzero(arrays.is_polymer[self._context_indices])
        self._protein_coords = self._context_coords[self._protein_positions]
        self._protein_radii = self._context_radii[self._protein_positions]
        self._protein_sasa = self._atom_sasa[self._protein_positions]
        self._protein_atom_indices = self._context_indices[self._protein_positions]
        self._residue_sasa = self._aggregate_residue_sasa()

    def _aggregate_residue_sasa(self) -> Dict[ResidueKey, float]:
        """Sum per-atom SASA into chain-aware per-residue totals."""
        totals: Dict[ResidueKey, float] = {}
        arrays = self.arrays
        for position, atom_index in enumerate(self._context_indices):
            key = (
                int(arrays.model_ids[atom_index]),
                str(arrays.chain_ids[atom_index]),
                int(arrays.resseqs[atom_index]),
                str(arrays.icodes[atom_index]),
            )
            totals[key] = totals.get(key, 0.0) + float(self._atom_sasa[position])
        return totals

    def shell_residues(
        self, centre: Sequence[float], radius: float
    ) -> List[Tuple[ResidueKey, str]]:
        """Return polymer residues with any atom within ``radius`` of ``centre``.

        Args:
            centre: Pocket centre.
            radius: Shell radius in Å.

        Returns:
            ``(residue_key, residue_name)`` pairs, de-duplicated.
        """
        from scipy.spatial import cKDTree

        if self._protein_coords.size == 0:
            return []
        tree = cKDTree(self._protein_coords)
        hits = tree.query_ball_point(np.asarray(centre, dtype=float), float(radius))
        arrays = self.arrays
        found: Dict[ResidueKey, str] = {}
        for local_index in hits:
            atom_index = self._protein_atom_indices[local_index]
            key = (
                int(arrays.model_ids[atom_index]),
                str(arrays.chain_ids[atom_index]),
                int(arrays.resseqs[atom_index]),
                str(arrays.icodes[atom_index]),
            )
            found.setdefault(key, str(arrays.resnames[atom_index]))
        return sorted(found.items())

    def coulomb_potential(self, centre: Sequence[float]) -> float:
        """Screened Coulomb potential at a point from formal side-chain charges.

        This surrogate exists because APBS is frequently unavailable, and a
        feature that is ``NaN`` for every row carries no information. The estimate
        sums Debye-Huckel screened contributions from formally charged side-chain
        groups:

        ``phi = prefactor / eps_r * sum_i q_i * exp(-r_i / lambda) / r_i``

        It is a continuum approximation - no explicit solvent, no titration
        shifts, no dielectric boundary - so it is *not* a substitute for a
        Poisson-Boltzmann calculation. It is a monotone, always-computable proxy
        that ranks pockets by how cationic their environment is, in kT/e so the
        scale matches the APBS feature.

        Args:
            centre: Query point.

        Returns:
            Screened potential in kT/e.
        """
        arrays = self.arrays
        centre_arr = np.asarray(centre, dtype=float)

        charged_atoms: List[Tuple[float, np.ndarray]] = []
        for atom_index in self._protein_atom_indices:
            resname = str(arrays.resnames[atom_index])
            charge = FORMAL_CHARGES.get(resname)
            if charge is None:
                continue
            atom_name = str(arrays.atom_names[atom_index]).upper()
            # Place the charge on the atoms that actually carry it.
            if resname in {"ARG"} and atom_name not in {"NE", "NH1", "NH2"}:
                continue
            if resname == "LYS" and atom_name != "NZ":
                continue
            if resname == "HIS" and atom_name not in {"ND1", "NE2"}:
                continue
            if resname == "ASP" and atom_name not in {"OD1", "OD2"}:
                continue
            if resname == "GLU" and atom_name not in {"OE1", "OE2"}:
                continue
            # Split the formal charge across the atoms bearing it.
            divisor = {"ARG": 3.0, "LYS": 1.0, "HIS": 2.0, "ASP": 2.0, "GLU": 2.0}[resname]
            charged_atoms.append((charge / divisor, arrays.coords[atom_index]))

        if not charged_atoms:
            return float("nan")

        charges = np.asarray([q for q, _ in charged_atoms], dtype=float)
        coords = np.asarray([c for _, c in charged_atoms], dtype=float)
        distances = np.linalg.norm(coords - centre_arr, axis=1)
        # Floor the distance so an atom essentially at the query point cannot
        # produce a singular contribution.
        distances = np.maximum(distances, 1.5)
        screened = np.exp(-distances / DEBYE_LENGTH) / distances
        return float(
            COULOMB_KT_PREFACTOR_ANGSTROM / RELATIVE_PERMITTIVITY * np.sum(charges * screened)
        )

    def extract(
        self,
        pocket_id: int,
        centre: Sequence[float],
        *,
        alpha_sphere_coords: Optional[np.ndarray] = None,
        fpocket_volume: Optional[float] = None,
        apbs_potential: Optional[float] = None,
    ) -> PocketFeatures:
        """Compute all descriptors for one pocket.

        Args:
            pocket_id: Pocket identifier.
            centre: Pocket centre coordinate.
            alpha_sphere_coords: Alpha-sphere centres from the detector, used for
                hull volume and shape descriptors.
            fpocket_volume: Volume reported by fpocket (Å³).
            apbs_potential: APBS potential at the centre (kT/e), when available.

        Returns:
            The pocket descriptors.
        """
        centre_arr = np.asarray(centre, dtype=float)
        features: Dict[str, float] = {}

        # ---------------------------------------------------------- geometry
        features["pocket_volume"] = (
            float(fpocket_volume) if fpocket_volume is not None else float("nan")
        )
        if alpha_sphere_coords is not None and len(alpha_sphere_coords) >= 4:
            spheres = np.asarray(alpha_sphere_coords, dtype=float)
            shape = point_cloud_shape(spheres)
            hull = convex_hull_volume(spheres)
            features["hull_volume"] = hull
            features["n_alpha_spheres"] = float(spheres.shape[0])
            features["alpha_sphere_density"] = (
                float(spheres.shape[0] / hull) if hull and np.isfinite(hull) and hull > 0 else np.nan
            )
            features["radius_of_gyration"] = shape.radius_of_gyration
            features["asphericity"] = shape.asphericity
            features["max_extent"] = shape.max_extent
        else:
            for name in (
                "hull_volume",
                "n_alpha_spheres",
                "alpha_sphere_density",
                "radius_of_gyration",
                "asphericity",
                "max_extent",
            ):
                features[name] = float("nan")
            if alpha_sphere_coords is not None:
                features["n_alpha_spheres"] = float(len(alpha_sphere_coords))

        # ------------------------------------------------------------ burial
        features["burial_depth"] = burial_depth(
            centre_arr, self._protein_coords, self._protein_sasa
        )
        features["enclosure"] = enclosure_fraction(
            centre_arr, self._protein_coords, self._protein_radii
        )

        # ------------------------------------- accessibility and composition
        shell = self.shell_residues(centre_arr, self.shell_radius)
        core = self.shell_residues(centre_arr, self.core_radius)
        residue_keys = [key for key, _ in shell]

        if shell:
            sasa_values = np.asarray(
                [self._residue_sasa.get(key, 0.0) for key, _ in shell], dtype=float
            )
            relative = np.asarray(
                [
                    self._residue_sasa.get(key, 0.0) / MAX_RESIDUE_SASA[name]
                    for key, name in shell
                    if name in MAX_RESIDUE_SASA
                ],
                dtype=float,
            )
            features["sasa_mean"] = float(np.mean(sasa_values))
            features["sasa_median"] = float(np.median(sasa_values))
            features["sasa_min"] = float(np.min(sasa_values))
            features["sasa_max"] = float(np.max(sasa_values))
            features["sasa_total"] = float(np.sum(sasa_values))
            features["mean_relative_sasa"] = (
                float(np.mean(relative)) if relative.size else float("nan")
            )
            features["buried_residue_fraction"] = (
                float(np.mean(relative < BURIED_RSA_THRESHOLD)) if relative.size else float("nan")
            )

            names = [name for _, name in shell]
            n_res = float(len(names))
            n_basic = float(sum(1 for name in names if name in BASIC_RESIDUES))
            n_strong = float(sum(1 for name in names if name in STRONG_BASIC_RESIDUES))
            n_acidic = float(sum(1 for name in names if name in ACIDIC_RESIDUES))
            n_aromatic = float(sum(1 for name in names if name in AROMATIC_RESIDUES))
            n_hydroxyl = float(sum(1 for name in names if name in HYDROXYL_RESIDUES))
            polar = {"SER", "THR", "TYR", "ASN", "GLN", "HIS", "CYS", "TRP"}
            n_polar = float(sum(1 for name in names if name in polar))

            features["n_residues"] = n_res
            features["n_basic_residues"] = n_basic
            features["n_strong_basic_residues"] = n_strong
            features["n_acidic_residues"] = n_acidic
            features["n_aromatic_residues"] = n_aromatic
            features["n_hydroxyl_residues"] = n_hydroxyl
            features["basic_fraction"] = n_basic / n_res
            features["acidic_fraction"] = n_acidic / n_res
            features["aromatic_fraction"] = n_aromatic / n_res
            features["hydroxyl_fraction"] = n_hydroxyl / n_res
            features["polar_fraction"] = n_polar / n_res
            hydropathy = [KYTE_DOOLITTLE[name] for name in names if name in KYTE_DOOLITTLE]
            features["hydropathy_mean"] = (
                float(np.mean(hydropathy)) if hydropathy else float("nan")
            )

            charge = sum(FORMAL_CHARGES.get(name, 0.0) for name in names)
            features["net_formal_charge"] = float(charge)
            volume = features.get("hull_volume")
            if volume is None or not np.isfinite(volume) or volume <= 0:
                volume = features.get("pocket_volume")
            features["positive_charge_density"] = (
                float((n_basic) / volume * 1000.0)
                if volume is not None and np.isfinite(volume) and volume > 0
                else float("nan")
            )
            features["charge_balance"] = float((n_basic - n_acidic) / max(n_basic + n_acidic, 1.0))
        else:
            for name in (
                "sasa_mean", "sasa_median", "sasa_min", "sasa_max", "sasa_total",
                "mean_relative_sasa", "buried_residue_fraction", "n_residues",
                "n_basic_residues", "n_strong_basic_residues", "n_acidic_residues",
                "n_aromatic_residues", "n_hydroxyl_residues", "basic_fraction",
                "acidic_fraction", "aromatic_fraction", "hydroxyl_fraction",
                "polar_fraction", "hydropathy_mean", "net_formal_charge",
                "positive_charge_density", "charge_balance",
            ):
                features[name] = float("nan")

        core_names = [name for _, name in core]
        features["n_basic_residues_core"] = float(
            sum(1 for name in core_names if name in BASIC_RESIDUES)
        )
        features["n_acidic_residues_core"] = float(
            sum(1 for name in core_names if name in ACIDIC_RESIDUES)
        )

        # -------------------------------------------------- charge geometry
        nitrogen_distances = self._basic_nitrogen_distances(centre_arr)
        features["n_basic_nitrogens"] = float(nitrogen_distances.size)
        if nitrogen_distances.size:
            features["basic_nitrogen_min_distance"] = float(np.min(nitrogen_distances))
            features["basic_nitrogen_mean_distance"] = float(np.mean(nitrogen_distances))
            # Dispersion separates a genuine cluster converging on one site from
            # basic residues scattered around a large, shallow surface patch.
            features["basic_nitrogen_dispersion"] = (
                float(np.std(nitrogen_distances)) if nitrogen_distances.size > 1 else 0.0
            )
        else:
            features["basic_nitrogen_min_distance"] = float("nan")
            features["basic_nitrogen_mean_distance"] = float("nan")
            features["basic_nitrogen_dispersion"] = float("nan")

        # ---------------------------------------------------- electrostatics
        features["coulomb_potential_kt"] = self.coulomb_potential(centre_arr)
        features["electrostatic_potential"] = (
            float(apbs_potential) if apbs_potential is not None else float("nan")
        )

        # -------------------------------------------------------- confidence
        plddt = self._plddt_summary(residue_keys)
        features.update(plddt)

        return PocketFeatures(
            pocket_id=int(pocket_id),
            centre=(float(centre_arr[0]), float(centre_arr[1]), float(centre_arr[2])),
            features=features,
            residue_keys=residue_keys,
        )

    def _basic_nitrogen_distances(self, centre: np.ndarray) -> np.ndarray:
        """Distances from the pocket centre to basic side-chain nitrogens.

        Only nitrogens within :attr:`shell_radius` are considered, so the count
        reflects groups that could actually coordinate a ligand at this site.
        """
        arrays = self.arrays
        distances: List[float] = []
        for atom_index in self._protein_atom_indices:
            if str(arrays.resnames[atom_index]) not in BASIC_RESIDUES:
                continue
            if str(arrays.atom_names[atom_index]).upper() not in BASIC_NITROGEN_ATOMS:
                continue
            distance = float(np.linalg.norm(arrays.coords[atom_index] - centre))
            if distance <= self.shell_radius:
                distances.append(distance)
        return np.asarray(sorted(distances), dtype=float)

    def _plddt_summary(self, residue_keys: Sequence[ResidueKey]) -> Dict[str, float]:
        """Summarise pLDDT (CA B-factors) over pocket-lining residues.

        For AlphaFold models the B-factor column holds pLDDT. For experimental
        structures it holds a true B-factor, so these descriptors are only
        meaningful for predicted models; the training code therefore treats them
        as optional and the imputer handles their absence.
        """
        arrays = self.arrays
        wanted = set(residue_keys)
        values: List[float] = []
        for atom_index in self._protein_atom_indices:
            if str(arrays.atom_names[atom_index]).upper() != "CA":
                continue
            key = (
                int(arrays.model_ids[atom_index]),
                str(arrays.chain_ids[atom_index]),
                int(arrays.resseqs[atom_index]),
                str(arrays.icodes[atom_index]),
            )
            if key in wanted:
                values.append(float(arrays.bfactors[atom_index]))
        if not values:
            return {
                "plddt_mean": float("nan"),
                "plddt_min": float("nan"),
                "plddt_fraction_above_cutoff": float("nan"),
            }
        arr = np.asarray(values, dtype=float)
        return {
            "plddt_mean": float(np.mean(arr)),
            "plddt_min": float(np.min(arr)),
            "plddt_fraction_above_cutoff": float(np.mean(arr >= PLDDT_CONFIDENCE_CUTOFF)),
        }


def legacy_feature_row(features: Mapping[str, float]) -> Dict[str, float]:
    """Map the descriptor suite onto the legacy six-feature schema.

    Provided so models trained by earlier versions remain usable. Note the
    correction: legacy ``pocket_depth`` is now filled from the genuine geometric
    ``burial_depth`` rather than from fpocket's hydrophobic density, so legacy
    models are fed the quantity their column name always claimed.

    Args:
        features: A descriptor mapping.

    Returns:
        A mapping keyed by :data:`LEGACY_FEATURE_NAMES`.
    """
    electrostatic = features.get("electrostatic_potential", np.nan)
    if electrostatic is None or not np.isfinite(electrostatic):
        electrostatic = features.get("coulomb_potential_kt", np.nan)
    return {
        "pocket_depth": features.get("burial_depth", np.nan),
        "sasa": features.get("sasa_mean", np.nan),
        "electrostatic_potential": electrostatic,
        "n_basic_residues": features.get("n_basic_residues", np.nan),
        "pocket_volume": features.get("pocket_volume", np.nan),
        "plddt_confidence": features.get("plddt_mean", np.nan),
    }


def feature_documentation() -> Dict[str, str]:
    """Return a one-line rationale for each descriptor, for methods reporting."""
    return {
        "pocket_volume": "fpocket alpha-sphere volume; InsP3-InsP6 need roughly 300-800 A^3.",
        "hull_volume": "Convex hull volume of alpha-sphere centres; detector-independent size.",
        "n_alpha_spheres": "Alpha-sphere count; a raw measure of cavity size.",
        "alpha_sphere_density": "Spheres per 1000 A^3; high density marks a tight, enclosed cavity.",
        "radius_of_gyration": "Spread of alpha spheres; compact sites match the compact ligand.",
        "asphericity": "Anisotropy of the cavity; 0 is globular, 1 is a slot.",
        "max_extent": "Longest cavity dimension; must accommodate an ~11 A ligand.",
        "burial_depth": "Distance from centre to nearest solvent-exposed atom; the depth criterion.",
        "enclosure": "Fraction of directions blocked by protein; distinguishes buried from grooved.",
        "buried_residue_fraction": "Share of lining residues with relative accessibility < 0.20.",
        "mean_relative_sasa": "Mean residue accessibility normalised by residue type.",
        "sasa_mean": "Mean absolute residue SASA of the lining shell.",
        "sasa_median": "Median lining-residue SASA; robust to one exposed outlier.",
        "sasa_min": "Most buried lining residue.",
        "sasa_max": "Most exposed lining residue; large values reveal an opening.",
        "sasa_total": "Total lining-shell SASA.",
        "n_residues": "Lining-shell residue count; normalises the composition counts.",
        "n_basic_residues": "Arg/Lys/His in the shell; phosphate coordination requires several.",
        "n_strong_basic_residues": "Arg/Lys only, excluding pH-dependent His.",
        "n_acidic_residues": "Asp/Glu in the shell; electrostatically unfavourable for a polyanion.",
        "n_aromatic_residues": "Aromatics that stack against the inositol ring.",
        "n_hydroxyl_residues": "Ser/Thr/Tyr hydroxyls that hydrogen-bond to phosphates.",
        "basic_fraction": "Basic residues per lining residue; size-independent.",
        "acidic_fraction": "Acidic residues per lining residue.",
        "aromatic_fraction": "Aromatic residues per lining residue.",
        "hydroxyl_fraction": "Hydroxyl residues per lining residue.",
        "polar_fraction": "Polar residues per lining residue.",
        "hydropathy_mean": "Mean Kyte-Doolittle hydropathy; IP sites are polar, not greasy.",
        "n_basic_residues_core": "Basic residues within 5 A of the centre.",
        "n_acidic_residues_core": "Acidic residues within 5 A of the centre.",
        "n_basic_nitrogens": "Coordinating side-chain nitrogens; counts groups, not residues.",
        "basic_nitrogen_min_distance": "Closest basic nitrogen; direct coordination distance.",
        "basic_nitrogen_mean_distance": "Mean basic-nitrogen distance from the centre.",
        "basic_nitrogen_dispersion": "Spread of basic-nitrogen distances; low means a real cluster.",
        "net_formal_charge": "Sum of formal side-chain charges at pH 7.4.",
        "positive_charge_density": "Basic residues per 1000 A^3 of cavity.",
        "charge_balance": "(basic - acidic) / (basic + acidic); bounded charge asymmetry.",
        "coulomb_potential_kt": "Debye-screened Coulomb potential at the centre (kT/e); always available.",
        "electrostatic_potential": "APBS Poisson-Boltzmann potential at the centre (kT/e); optional.",
        "plddt_mean": "Mean pLDDT of lining residues (AlphaFold models only).",
        "plddt_min": "Minimum lining-residue pLDDT.",
        "plddt_fraction_above_cutoff": "Fraction of lining residues with pLDDT >= 70.",
    }
