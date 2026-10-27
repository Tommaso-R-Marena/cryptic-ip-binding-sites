"""Finding a pocket that is not there yet: cavity detection across conformers.

Why this exists
---------------
This is the part of the package that is actually new, and the part the rest is
built to support.

Every candidate in this project is a single apo AlphaFold conformer, and the
screen measures how buried a pocket is in it. A *cryptic* site is one that is
closed in the apo state. So the deployed method is looking for the thing it is
named after in the one structure where, by definition, it is hardest to see.
The project's own known-binder check shows the consequence: four annotated
binders are blocked by the gates.

The signal this module computes is not burial. It is **persistence and
crypticity across a deterministic conformer ensemble**:

* Displace the structure along the softest elastic-network modes, both
  directions, at a few amplitudes. This is a spanning set of the protein's
  cheap motions, not a Boltzmann sample, and it costs one eigendecomposition.
* Detect cavities in every conformer with the in-process detector.
* Match cavities across conformers by position, so one physical site becomes
  one track.
* For each track report **persistence** - the fraction of conformers in which
  it exists - and **crypticity** - how much larger it gets away from the apo
  structure.

A site with high persistence and low crypticity is an ordinary open pocket.
A site with low-to-middling persistence and high crypticity is the thing this
project set out to find and has never had a way to score. A site with very low
persistence is noise in the detector.

What would make this wrong
--------------------------
Elastic-network modes are global and side-chain-free, so they open cavities by
moving domains rather than by rotating a tryptophan out of the way, and
side-chain cryptic sites are a large share of the real ones. The displacement
amplitudes are in angstrom along an unweighted mode, which is not a physical
temperature. And the matching tolerance is a free parameter that trades
merging two neighbouring sites against splitting one moving site in two.

So crypticity here is a *hypothesis generator*, and the pre-registered
question in docs/DETECTOR_PLAN.md is narrow: does it rank the known cryptic
IP-binding sites above matched controls? If it does not, it is reported as a
negative result and the method is not used.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from cryptic_ip.detect.cavity import Cavity, CavityMap, detect_cavities
from cryptic_ip.detect.openability import NormalModes, perturbed_conformers
from cryptic_ip.detect.rotamer import (
    RotamerError, lining_residues, rotamer_conformers, rotatable_side_chains,
)

#: How far two cavity centres can be apart and still be the same site, in
#: angstrom. A site whose lining moves by 2 Å moves its centroid by less than
#: that, so 4.0 tracks a moving site without merging neighbours 8 Å apart.
DEFAULT_MATCH_TOLERANCE = 4.0

#: Displacement amplitudes along each mode, in angstrom of largest single-residue
#: motion. Both signs, because a mode is traversed in both directions, and a
#: pocket that opens one way usually closes the other.
DEFAULT_AMPLITUDES: Tuple[float, ...] = (-2.0, -1.0, 1.0, 2.0)

#: Which modes to displace along. The three softest are collective in almost
#: every protein; past about six they become local and the ensemble stops being
#: a spanning set of cheap motions.
DEFAULT_MODES: Tuple[int, ...] = (0, 1, 2)


class EnsembleError(ValueError):
    """An ensemble from which no persistence can be computed."""


@dataclass(frozen=True)
class CavityTrack:
    """One physical site, followed across the conformer ensemble."""

    centre: Tuple[float, float, float]
    persistence: float
    crypticity: float
    apo_volume: float
    max_volume: float
    mean_volume: float
    volume_spread: float
    n_conformers: int
    n_present: int
    mean_enclosure: float
    best_conformer: int

    @property
    def classification(self) -> str:
        """A label a reader can argue with, from two numbers.

        The boundaries are stated here rather than fitted, and they are the
        ones docs/DETECTOR_PLAN.md pre-registers.
        """
        if self.persistence >= 0.8 and self.crypticity < 0.25:
            return "open"
        if self.crypticity >= 0.5 and self.persistence >= 0.25:
            return "cryptic"
        if self.crypticity >= 0.25:
            return "semi_cryptic"
        if self.persistence < 0.25:
            return "transient"
        return "open"

    def as_dict(self) -> Dict[str, object]:
        return {
            "centre": list(self.centre),
            "persistence": float(self.persistence),
            "crypticity": float(self.crypticity),
            "apo_volume": float(self.apo_volume),
            "max_volume": float(self.max_volume),
            "mean_volume": float(self.mean_volume),
            "volume_spread": float(self.volume_spread),
            "n_conformers": int(self.n_conformers),
            "n_present": int(self.n_present),
            "mean_enclosure": float(self.mean_enclosure),
            "best_conformer": int(self.best_conformer),
            "classification": self.classification,
        }


@dataclass
class EnsembleReport:
    """Every tracked site, and what the ensemble was."""

    tracks: Tuple[CavityTrack, ...]
    n_conformers: int
    amplitudes: Tuple[float, ...]
    modes: Tuple[int, ...]
    match_tolerance: float

    def cryptic(self) -> List[CavityTrack]:
        return [t for t in self.tracks if t.classification in ("cryptic", "semi_cryptic")]

    def most_cryptic(self) -> Optional[CavityTrack]:
        return max(self.tracks, key=lambda t: t.crypticity) if self.tracks else None

    def as_dict(self) -> Dict[str, object]:
        return {
            "n_conformers": self.n_conformers,
            "amplitudes": list(self.amplitudes),
            "modes": list(self.modes),
            "match_tolerance": self.match_tolerance,
            "n_tracks": len(self.tracks),
            "n_cryptic": len(self.cryptic()),
            "tracks": [t.as_dict() for t in self.tracks],
        }


def match_cavities(
    maps: Sequence[CavityMap],
    *,
    tolerance: float = DEFAULT_MATCH_TOLERANCE,
) -> List[List[Optional[Cavity]]]:
    """Group cavities across conformers into tracks, one row per track.

    Greedy nearest-centre matching, seeded by the first conformer and extended
    by later ones. Greedy rather than optimal because the alternative is an
    assignment problem per conformer pair for a tolerance that is itself
    approximate, and because a site that moves more than the tolerance between
    adjacent amplitudes is one this method should not claim to track.

    ``maps[0]`` is taken as the apo conformer, and every track's row has one
    entry per conformer - ``None`` where the site was not detected.
    """
    if not maps:
        raise EnsembleError("no conformers to match")
    if tolerance <= 0:
        raise EnsembleError(f"tolerance must be positive, got {tolerance}")
    n = len(maps)
    tracks: List[List[Optional[Cavity]]] = []
    centres: List[np.ndarray] = []

    for conformer_index, cavity_map in enumerate(maps):
        taken = set()
        # Existing tracks claim their nearest unclaimed cavity.
        for track_index, centre in enumerate(centres):
            best, best_distance = None, tolerance
            for cavity_index, cavity in enumerate(cavity_map.cavities):
                if cavity_index in taken:
                    continue
                distance = float(np.linalg.norm(np.asarray(cavity.centre) - centre))
                if distance < best_distance:
                    best, best_distance = cavity_index, distance
            if best is not None:
                taken.add(best)
                tracks[track_index][conformer_index] = cavity_map.cavities[best]
        # Anything unclaimed starts a new track, absent from earlier conformers.
        for cavity_index, cavity in enumerate(cavity_map.cavities):
            if cavity_index in taken:
                continue
            row: List[Optional[Cavity]] = [None] * n
            row[conformer_index] = cavity
            tracks.append(row)
            centres.append(np.asarray(cavity.centre, dtype=float))
    return tracks


def _track_report(row: Sequence[Optional[Cavity]], n_conformers: int) -> CavityTrack:
    present = [(index, cavity) for index, cavity in enumerate(row) if cavity is not None]
    if not present:
        raise EnsembleError("a track with no detections cannot be summarised")
    volumes = np.array([cavity.volume for _, cavity in present], dtype=float)
    apo = row[0].volume if row[0] is not None else 0.0
    largest = float(volumes.max())
    # Crypticity: how much of the site's best volume is absent in the apo
    # structure. 0 when the apo is already its largest, 1 when the apo has
    # nothing at all.
    crypticity = float(1.0 - apo / largest) if largest > 0 else 0.0
    best_index = int(present[int(np.argmax(volumes))][0])
    # The centre is the mean over detections, which is the site's position
    # rather than any one conformer's.
    centre = np.mean([np.asarray(cavity.centre) for _, cavity in present], axis=0)
    return CavityTrack(
        centre=(float(centre[0]), float(centre[1]), float(centre[2])),
        persistence=float(len(present) / n_conformers),
        crypticity=crypticity,
        apo_volume=float(apo),
        max_volume=largest,
        mean_volume=float(volumes.mean()),
        volume_spread=float(volumes.std(ddof=1)) if volumes.size > 1 else 0.0,
        n_conformers=int(n_conformers),
        n_present=len(present),
        mean_enclosure=float(np.mean([cavity.mean_enclosure for _, cavity in present])),
        best_conformer=best_index,
    )


def analyse_ensemble(
    maps: Sequence[CavityMap],
    *,
    tolerance: float = DEFAULT_MATCH_TOLERANCE,
    amplitudes: Sequence[float] = DEFAULT_AMPLITUDES,
    modes: Sequence[int] = DEFAULT_MODES,
) -> EnsembleReport:
    """Persistence and crypticity for every site in a conformer ensemble."""
    if not maps:
        raise EnsembleError("no conformers")
    rows = match_cavities(maps, tolerance=tolerance)
    tracks = tuple(
        sorted(
            (_track_report(row, len(maps)) for row in rows),
            key=lambda t: -t.max_volume,
        )
    )
    return EnsembleReport(
        tracks=tracks,
        n_conformers=len(maps),
        amplitudes=tuple(float(a) for a in amplitudes),
        modes=tuple(int(m) for m in modes),
        match_tolerance=float(tolerance),
    )


def scan_structure(
    coords: np.ndarray,
    elements: Sequence[str],
    modes: NormalModes,
    node_coords: np.ndarray,
    node_of_atom: Sequence[int],
    *,
    amplitudes: Sequence[float] = DEFAULT_AMPLITUDES,
    which_modes: Sequence[int] = DEFAULT_MODES,
    tolerance: float = DEFAULT_MATCH_TOLERANCE,
    spacing: float = 0.8,
    min_enclosure: int = 5,
) -> EnsembleReport:
    """The ensemble scan for one structure. An alias of :func:`scan_with_nodes`.

    Kept as a separate name because "scan this structure" is what a caller
    wants to write, and because an earlier version of this function tried to
    infer the atom-to-residue map from nearest-neighbour distances. That
    inference is wrong in exactly the case that matters: inside a closed pocket
    the nearest alpha carbon to a side-chain atom is often a residue across the
    cavity, so displacing along a mode would pull the two walls together
    instead of apart. The map is required rather than guessed.
    """
    return scan_with_nodes(
        coords, elements, node_coords, node_of_atom, modes,
        amplitudes=amplitudes, which_modes=which_modes, tolerance=tolerance,
        spacing=spacing, min_enclosure=min_enclosure,
    )


def scan_two_tier(
    arrays,
    modes: NormalModes,
    node_coords: np.ndarray,
    node_of_atom: Sequence[int],
    centre: Optional[Sequence[float]] = None,
    *,
    amplitudes: Sequence[float] = DEFAULT_AMPLITUDES,
    which_modes: Sequence[int] = DEFAULT_MODES,
    tolerance: float = DEFAULT_MATCH_TOLERANCE,
    spacing: float = 0.8,
    min_enclosure: int = 5,
    lining_radius: float = 8.0,
    max_rotamers: Optional[int] = 60,
    rotamers: bool = True,
) -> Tuple[EnsembleReport, Dict[str, object]]:
    """Both tiers of the ensemble: backbone modes and side-chain rotamers.

    The two tiers open different things and neither subsumes the other. A
    backbone mode opens a site by moving a domain or a loop; a chi rotation
    opens one that a single side chain occludes, which no amount of
    low-frequency backbone motion will do. On a 6.5 angstrom helix-barrel
    channel plugged by tryptophans, the rotamer tier alone gives a crypticity
    of 0.45, and the figure falls monotonically as the barrel widens - the
    tighter the channel, the more one side chain decides.

    Returns the report and a dictionary describing what the rotamer tier did,
    including the side chains whose input geometry was already strained and
    the ones whose torsion axis was not a covalent bond at all.
    """
    coords = np.asarray(arrays.coords, dtype=float)
    elements = [str(e) for e in np.asarray(arrays.elements, dtype=object)]
    node_coords = np.asarray(node_coords, dtype=float)
    node_of_atom = np.asarray(node_of_atom, dtype=int)

    maps: List[CavityMap] = [
        detect_cavities(coords, elements, spacing=spacing, min_enclosure=min_enclosure)
    ]
    # Tier one: backbone modes.
    for displaced in perturbed_conformers(
        node_coords, modes, amplitudes=amplitudes, which_modes=which_modes,
    ):
        shift = displaced - node_coords
        maps.append(detect_cavities(
            coords + shift[node_of_atom], elements,
            spacing=spacing, min_enclosure=min_enclosure,
        ))
    n_backbone = len(maps) - 1

    # Tier two: side-chain rotamers around the pocket.
    rotamer_note: Dict[str, object] = {"enabled": bool(rotamers), "n_conformers": 0}
    if rotamers:
        if centre is None:
            largest = maps[0].largest()
            centre = largest.centre if largest is not None else coords.mean(axis=0)
        residues = lining_residues(arrays, centre, radius=lining_radius)
        try:
            chains = rotatable_side_chains(arrays, residues=residues)
        except RotamerError as exc:
            # Unphysical geometry is reported, never silently treated as a
            # tightly packed pocket - which is what it looks like otherwise.
            rotamer_note["unavailable"] = str(exc)
            chains = []
        if chains:
            ensemble = rotamer_conformers(arrays, chains, max_conformers=max_rotamers)
            for conformer in ensemble.conformers:
                maps.append(detect_cavities(
                    conformer, elements, spacing=spacing, min_enclosure=min_enclosure,
                ))
            rotamer_note.update(ensemble.as_dict())
            rotamer_note["n_conformers"] = ensemble.n_accepted
        rotamer_note["n_lining_residues"] = len(residues)
        rotamer_note["n_torsions"] = len(chains)

    report = analyse_ensemble(
        maps, tolerance=tolerance, amplitudes=amplitudes, modes=which_modes,
    )
    rotamer_note["n_backbone_conformers"] = n_backbone
    return report, rotamer_note


def scan_with_nodes(
    atom_coords: np.ndarray,
    elements: Sequence[str],
    node_coords: np.ndarray,
    node_of_atom: Sequence[int],
    modes: NormalModes,
    *,
    amplitudes: Sequence[float] = DEFAULT_AMPLITUDES,
    which_modes: Sequence[int] = DEFAULT_MODES,
    tolerance: float = DEFAULT_MATCH_TOLERANCE,
    spacing: float = 0.8,
    min_enclosure: int = 5,
) -> EnsembleReport:
    """The ensemble scan with an explicit atom-to-residue map.

    This is the path a real structure takes: ``node_of_atom[i]`` is the index
    of atom ``i``'s residue in the elastic network, so each atom moves with its
    own residue rather than with whichever Cα happens to be nearest.
    """
    atom_coords = np.asarray(atom_coords, dtype=float)
    node_coords = np.asarray(node_coords, dtype=float)
    node_of_atom = np.asarray(node_of_atom, dtype=int)
    if node_of_atom.shape != (len(atom_coords),):
        raise EnsembleError("one node index per atom")
    if node_coords.shape != (modes.n_nodes, 3):
        raise EnsembleError("node coordinates and modes describe different structures")
    if node_of_atom.min() < 0 or node_of_atom.max() >= modes.n_nodes:
        raise EnsembleError("a node index is outside the network")

    maps: List[CavityMap] = [
        detect_cavities(atom_coords, elements, spacing=spacing, min_enclosure=min_enclosure)
    ]
    for displaced in perturbed_conformers(
        node_coords, modes, amplitudes=amplitudes, which_modes=which_modes,
    ):
        shift = displaced - node_coords
        maps.append(detect_cavities(
            atom_coords + shift[node_of_atom], elements,
            spacing=spacing, min_enclosure=min_enclosure,
        ))
    return analyse_ensemble(maps, tolerance=tolerance, amplitudes=amplitudes, modes=which_modes)
