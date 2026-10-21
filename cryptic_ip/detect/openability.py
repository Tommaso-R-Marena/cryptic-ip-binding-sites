"""How hard a closed pocket is to open, from an elastic network model.

Why this exists
---------------
This project's biggest standing weakness, in its own words, is that every
candidate is a single apo AlphaFold conformer: no induced fit, no alternative
rotamers, no evidence the pocket exists as anything but a static prediction.
That weakness is worse than it sounds, because it is self-contradicting. A
*cryptic* site is by definition one that is not open in the apo structure. A
pipeline that measures how buried a pocket is in a single apo conformer is
measuring the wrong variable: it will find sites that are already open and
lined with basic residues, and it will miss exactly the class it is named
after.

The deployed screen's own known-binder check shows this happening. ADAR2 and
ADAR1 pass, and ADAT1, PDS5B and HDAC1 are blocked at the score gate with
HDAC3 blocked at the volume gate - all of them annotated binders.

So this module scores a different quantity: not *is the pocket open* but *what
would it cost to open it*. The answer has to be cheap enough for a proteome,
deterministic, and computable from one structure, which rules out molecular
dynamics and rules in a normal-mode argument.

The method, stated so it can be attacked
----------------------------------------
1. Build an anisotropic network model: every Cα pair within a cutoff is a
   spring, and the Hessian of that harmonic energy is a 3N x 3N matrix whose
   low eigenvectors are the protein's softest collective motions. This is a
   standard construction and nothing here is new.
2. What is new is the projection. For each pocket, take the lining residues and
   an outward direction for each - away from the pocket centre - and ask how
   much of each low mode moves those residues *outward together*. That is a
   single number per mode, the **opening overlap**, and it says whether the
   protein's cheap motions happen to be the motion that opens this pocket.
3. Weight each mode's overlap by the inverse of its eigenvalue, because the
   thermal amplitude of a harmonic mode goes as 1/lambda. The sum is an
   **openability**: large when a soft mode opens the pocket, small when
   opening it requires a stiff one.
4. Report the implied **opening cost** in kT for a stated opening distance,
   from the harmonic energy of displacing along the weighted mode combination.

What this is not
----------------
An elastic network model knows nothing about chemistry. It cannot tell a
hinge that opens a real cavity from one that opens a crevice, it has no side
chains, and its spring constant is a single global parameter with no units
anyone has calibrated against a measured opening free energy. The openability
number is therefore *ordinal at best*: it may rank pockets, and it certainly
does not predict a free energy. The pre-registered question in
docs/DETECTOR_PLAN.md asks only whether it ranks the known cryptic binders
above matched controls - and if it does not, that is the answer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

#: Spring cutoff in angstrom. 15 Å is the long-range convention for Cα
#: anisotropic network models; shorter cutoffs make the network too sparse to
#: carry collective motion and longer ones make every mode global.
DEFAULT_CUTOFF = 15.0

#: Spring constant, in arbitrary energy per angstrom squared. It scales every
#: eigenvalue linearly, so it cancels from any ranking and it is kept explicit
#: only so that the cost numbers have a stated origin rather than an implied
#: calibration.
DEFAULT_GAMMA = 1.0

#: How many non-trivial modes to use. The first six eigenvalues of a free
#: elastic network are zero by construction - three translations, three
#: rotations - so they are skipped. 20 is where the modes stop being
#: collective in most proteins.
DEFAULT_N_MODES = 20

#: An eigenvalue at or below this is treated as a rigid-body mode rather than a
#: real one. Scaled by the largest eigenvalue, so it is a conditioning test
#: rather than an absolute one.
RIGID_BODY_TOLERANCE = 1e-8

#: The opening distance the cost is quoted for, in angstrom. An inositol
#: phosphate is about 9-10 Å across, so a site that has to admit one needs an
#: opening of a few angstrom at least. 2.0 is a deliberately modest number and
#: the cost scales as its square.
DEFAULT_OPENING = 2.0

#: The cost this module reports is in units of the spring constant times
#: angstrom squared - *not* kT. Converting it to kT needs a gamma calibrated
#: against measured opening free energies, which nobody here has done, and the
#: field was called `opening_cost_kt` until that was noticed. A name that
#: implies a unit the number does not have is how this pipeline previously came
#: to read fpocket's mean local hydrophobic density as a burial depth in
#: angstrom, so the name now carries its real units.
COST_UNITS = "gamma*A^2"


class OpenabilityError(ValueError):
    """A structure or pocket from which no openability can be computed."""


def alpha_carbons(arrays) -> Tuple[np.ndarray, List[str]]:
    """Cα coordinates and residue keys, polymer only.

    An elastic network built on every atom is 20 times larger for no gain in
    the low modes, which are what this module uses.
    """
    polymer = np.asarray(arrays.is_polymer, dtype=bool)
    names = np.asarray(arrays.atom_names, dtype=object)
    is_ca = polymer & (np.char.upper(names.astype(str)) == "CA")
    index = np.nonzero(is_ca)[0]
    if len(index) < 4:
        raise OpenabilityError(
            f"an elastic network needs at least 4 residues, found {len(index)} alpha carbons"
        )
    chains = np.asarray(arrays.chain_ids, dtype=object)
    resseqs = np.asarray(arrays.resseqs, dtype=object)
    keys = [f"{chains[i]}:{resseqs[i]}" for i in index]
    return arrays.coords[index].astype(float), keys


def anm_hessian(
    coords: np.ndarray,
    *,
    cutoff: float = DEFAULT_CUTOFF,
    gamma: float = DEFAULT_GAMMA,
) -> coo_matrix:
    """The 3N x 3N anisotropic network model Hessian.

    For each pair within ``cutoff`` the off-diagonal 3x3 block is
    ``-gamma * outer(u, u)`` with ``u`` the unit vector between them, and the
    diagonal block is minus the sum of that residue's off-diagonal blocks. That
    construction is what makes the six rigid-body modes exactly zero, and
    `test_the_hessian_has_exactly_six_rigid_body_modes` is the check that the
    implementation has it right.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise OpenabilityError("coordinates must have shape (n, 3)")
    if cutoff <= 0:
        raise OpenabilityError(f"cutoff must be positive, got {cutoff}")
    n = len(coords)
    if n < 4:
        raise OpenabilityError(f"need at least 4 nodes, got {n}")

    rows: List[int] = []
    cols: List[int] = []
    data: List[float] = []
    diagonal = np.zeros((n, 3, 3))
    cutoff_sq = cutoff * cutoff

    for i in range(n):
        delta = coords[i + 1:] - coords[i]
        distance_sq = np.sum(delta * delta, axis=1)
        near = np.nonzero((distance_sq <= cutoff_sq) & (distance_sq > 0))[0]
        for offset in near:
            j = i + 1 + int(offset)
            unit = delta[offset] / np.sqrt(distance_sq[offset])
            block = -gamma * np.outer(unit, unit)
            for a in range(3):
                for b in range(3):
                    value = block[a, b]
                    if value == 0.0:
                        continue
                    rows.append(3 * i + a)
                    cols.append(3 * j + b)
                    data.append(value)
                    rows.append(3 * j + b)
                    cols.append(3 * i + a)
                    data.append(value)
            diagonal[i] -= block
            diagonal[j] -= block

    for i in range(n):
        for a in range(3):
            for b in range(3):
                value = diagonal[i][a, b]
                if value == 0.0:
                    continue
                rows.append(3 * i + a)
                cols.append(3 * i + b)
                data.append(value)
    return coo_matrix((data, (rows, cols)), shape=(3 * n, 3 * n))


@dataclass(frozen=True)
class NormalModes:
    """Low-frequency modes of an elastic network."""

    eigenvalues: np.ndarray        # (n_modes,) ascending, rigid body removed
    eigenvectors: np.ndarray       # (3N, n_modes)
    n_nodes: int
    n_rigid_found: int
    cutoff: float

    @property
    def n_modes(self) -> int:
        return int(len(self.eigenvalues))

    def displacement(self, mode: int) -> np.ndarray:
        """Mode ``mode`` reshaped to a per-node displacement, ``(N, 3)``."""
        if not 0 <= mode < self.n_modes:
            raise OpenabilityError(f"mode {mode} is outside 0..{self.n_modes - 1}")
        return self.eigenvectors[:, mode].reshape(self.n_nodes, 3)


def normal_modes(
    coords: np.ndarray,
    *,
    cutoff: float = DEFAULT_CUTOFF,
    gamma: float = DEFAULT_GAMMA,
    n_modes: int = DEFAULT_N_MODES,
) -> NormalModes:
    """The softest non-trivial modes of the network.

    Six eigenvalues are zero by construction; they are found rather than
    assumed, and a structure that produces a different number of them is a
    disconnected network - two chains with no contact between them, usually -
    which is reported rather than silently analysed.
    """
    hessian = anm_hessian(coords, cutoff=cutoff, gamma=gamma)
    n = len(coords)
    wanted = min(n_modes + 6, 3 * n - 1)
    if wanted < 7:
        raise OpenabilityError(f"a network of {n} nodes is too small for mode analysis")
    # `sigma=0` with 'LM' converges on the smallest eigenvalues, which is what
    # a shift-invert solve is for; the alternative ('SM' directly) is slow and
    # often fails to converge on a Hessian with a six-fold null space.
    matrix = hessian.tocsc()
    values, vectors = eigsh(matrix, k=wanted, sigma=0.0, which="LM")
    order = np.argsort(values)
    values = values[order]
    vectors = vectors[:, order]

    # The scale has to come from the whole Hessian, not from the eigenvalues
    # returned: those are deliberately the smallest ones, so scaling by their
    # maximum made every zero mode look large relative to its peers and the
    # disconnection guard below never fired. The mean diagonal entry is a
    # cheap, robust stand-in for the spectrum's scale.
    scale = float(matrix.diagonal().mean()) or 1.0
    rigid = int(np.sum(values <= RIGID_BODY_TOLERANCE * scale))
    if rigid > 6:
        raise OpenabilityError(
            f"the network has {rigid} zero modes, not 6: it is disconnected, so no "
            "collective motion is defined across it. Check for chains with no contacts."
        )
    keep = slice(6, 6 + n_modes)
    return NormalModes(
        eigenvalues=values[keep].copy(),
        eigenvectors=vectors[:, keep].copy(),
        n_nodes=n,
        n_rigid_found=rigid,
        cutoff=cutoff,
    )


def outward_directions(coords: np.ndarray, centre: Sequence[float]) -> np.ndarray:
    """Unit vectors from ``centre`` to each of ``coords``.

    These are the directions in which a residue has to move for the pocket to
    get bigger, which is what "opening" means here.
    """
    coords = np.asarray(coords, dtype=float)
    centre = np.asarray(centre, dtype=float).reshape(3)
    delta = coords - centre[None, :]
    norms = np.linalg.norm(delta, axis=1)
    if np.any(norms == 0):
        raise OpenabilityError("a lining residue sits exactly on the pocket centre")
    return delta / norms[:, None]


@dataclass(frozen=True)
class OpenabilityReport:
    """What the elastic network says about opening one pocket."""

    openability: float
    opening_cost_gamma_a2: float
    best_mode: int
    best_overlap: float
    per_mode_overlap: np.ndarray
    eigenvalues: np.ndarray
    n_lining: int
    opening_angstrom: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "openability": float(self.openability),
            "opening_cost_gamma_a2": float(self.opening_cost_gamma_a2),
            "best_mode": int(self.best_mode),
            "best_overlap": float(self.best_overlap),
            "per_mode_overlap": [float(v) for v in self.per_mode_overlap],
            "eigenvalues": [float(v) for v in self.eigenvalues],
            "n_lining": int(self.n_lining),
            "opening_angstrom": float(self.opening_angstrom),
        }


def opening_overlap(modes: NormalModes, lining: Sequence[int], outward: np.ndarray) -> np.ndarray:
    """How much each mode moves the lining residues outward together.

    The overlap of a mode is the mean projection of its displacement on the
    outward directions, normalised by the mode's displacement magnitude over
    those residues. It is signed, and the sign is taken positive because a mode
    and its negative are the same motion: a harmonic mode is traversed in both
    directions.
    """
    lining = np.asarray(list(lining), dtype=int)
    if lining.size == 0:
        raise OpenabilityError("a pocket with no lining residues has no opening motion")
    if outward.shape != (lining.size, 3):
        raise OpenabilityError("one outward direction per lining residue")
    out = np.zeros(modes.n_modes)
    # Cauchy-Schwarz bounds the summed projection by ||displacement|| times the
    # norm of the stacked outward directions, which is sqrt(n_lining) because
    # each is a unit vector. Dividing by ||displacement|| alone - as this did
    # at first - produced "overlaps" above 1.5, which is not a cosine of
    # anything and would have gone into a ranking as if it were.
    denominator_scale = float(np.sqrt(lining.size))
    for mode in range(modes.n_modes):
        displacement = modes.displacement(mode)[lining]
        magnitude = float(np.linalg.norm(displacement))
        if magnitude == 0.0:
            continue
        projection = float(np.sum(displacement * outward))
        out[mode] = abs(projection) / (magnitude * denominator_scale)
    return out


def describe_pocket(
    arrays,
    centre: Sequence[float],
    lining_keys: Optional[Sequence[str]] = None,
    *,
    lining_radius: float = 8.0,
    cutoff: float = DEFAULT_CUTOFF,
    n_modes: int = DEFAULT_N_MODES,
    opening: float = DEFAULT_OPENING,
) -> OpenabilityReport:
    """Openability of the pocket at ``centre``.

    The lining is given either explicitly as residue keys or, by default, as
    every residue whose Cα is within ``lining_radius`` of the centre. The
    default is deliberately generous: a closed cryptic pocket has no cavity to
    line, so a definition that requires contact with a cavity surface would
    return nothing for exactly the pockets this module exists to score.
    """
    coords, keys = alpha_carbons(arrays)
    centre = np.asarray(centre, dtype=float).reshape(3)
    if lining_keys is not None:
        wanted = {str(key) for key in lining_keys}
        lining = [i for i, key in enumerate(keys) if key in wanted]
        if not lining:
            raise OpenabilityError(
                f"none of the {len(wanted)} lining residue keys matched this structure; "
                "check the chain and residue numbering"
            )
    else:
        distance = np.linalg.norm(coords - centre[None, :], axis=1)
        lining = list(np.nonzero(distance <= lining_radius)[0])
        if not lining:
            raise OpenabilityError(
                f"no alpha carbon within {lining_radius} A of the pocket centre"
            )

    modes = normal_modes(coords, cutoff=cutoff, n_modes=n_modes)
    outward = outward_directions(coords[lining], centre)
    overlap = opening_overlap(modes, lining, outward)

    # Thermal amplitude of a harmonic mode goes as 1/lambda, so a soft mode
    # that opens the pocket counts far more than a stiff one that does.
    weights = 1.0 / modes.eigenvalues
    openability = float(np.sum(overlap * weights) / np.sum(weights))

    # The cost of an opening of `opening` angstrom along the best-overlapping
    # mode: harmonic energy is 0.5 * lambda * amplitude^2, and the amplitude
    # needed scales inversely with how much of the mode is the opening motion.
    best = int(np.argmax(overlap))
    best_overlap = float(overlap[best])
    if best_overlap <= 0.0:
        cost = float("inf")
    else:
        amplitude = opening / best_overlap
        cost = float(0.5 * modes.eigenvalues[best] * amplitude * amplitude)
    return OpenabilityReport(
        openability=openability,
        opening_cost_gamma_a2=cost,
        best_mode=best,
        best_overlap=best_overlap,
        per_mode_overlap=overlap,
        eigenvalues=modes.eigenvalues.copy(),
        n_lining=len(lining),
        opening_angstrom=float(opening),
    )


def perturbed_conformers(
    coords: np.ndarray,
    modes: NormalModes,
    *,
    amplitudes: Sequence[float] = (-2.0, -1.0, 1.0, 2.0),
    which_modes: Sequence[int] = (0, 1, 2),
) -> List[np.ndarray]:
    """Deterministic conformers displaced along the softest modes.

    This is the ensemble a cryptic-pocket search needs and a single apo model
    cannot give. It is not molecular dynamics and makes no claim to sample a
    Boltzmann distribution: it is a spanning set of the cheap motions, chosen
    so that pocket detection can be run on more than one structure without
    anybody waiting for a trajectory.

    The displacements are in angstrom along the mass-unweighted mode, scaled so
    that the largest single-residue displacement equals the requested
    amplitude - otherwise a mode spread over a big protein moves every residue
    by almost nothing and the conformer is indistinguishable from the input.
    """
    coords = np.asarray(coords, dtype=float)
    if coords.shape != (modes.n_nodes, 3):
        raise OpenabilityError("coordinates and modes describe different structures")
    out: List[np.ndarray] = []
    for mode in which_modes:
        displacement = modes.displacement(int(mode))
        largest = float(np.max(np.linalg.norm(displacement, axis=1)))
        if largest == 0.0:
            continue
        unit = displacement / largest
        for amplitude in amplitudes:
            out.append(coords + amplitude * unit)
    if not out:
        raise OpenabilityError("no usable modes to perturb along")
    return out
