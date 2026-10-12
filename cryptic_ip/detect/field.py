"""Rotation-invariant electrostatic complementarity for polyanion sites.

Why this exists
---------------
The deployed screen decides whether a pocket could hold an inositol phosphate
with a cascade of scalar gates, one of which is "at least four basic residues".
That gate cannot distinguish an IP6 site from an ATP, PAPS, UDP-sugar or
sulfate site, and 30 of the project's 75 candidates are already flagged as
explained by some other phosphate- or sulfate-rich ligand. Counting charges
throws away the thing that actually differs between those sites: *where* the
charge is, relative to the pocket.

Study N tried to use that spatial information by rigidly transplanting 131
crystal IP6 sites onto each candidate pocket. It produced a weak guard (AUC
0.599 [0.546, 0.652]) and placements that were not buried - median buried
fraction 0.000, zero for 22 of 35 - because a rigid transplant has to be
*aligned*, and the alignment is a six-parameter search that can fail quietly.

This module takes the spatial information without the alignment. The pocket's
electrostatic potential is expanded on concentric shells in real spherical
harmonics, and the per-shell, per-degree power is kept. That power spectrum is
invariant to rotation about the pocket centre, exactly, because a rotation
mixes the harmonics within a degree and leaves the sum of their squares alone.
So a pocket can be compared with a reference charge distribution without
anybody choosing an orientation.

What is new here, stated plainly so it can be criticised
--------------------------------------------------------
1. The potential is computed analytically, from formal charges on ionizable
   groups under Debye-Hückel screening, rather than by solving
   Poisson-Boltzmann. That is a real approximation and it is deliberate: APBS
   is an optional dependency this pipeline degrades past with a warning, so
   anything built on it is absent from most runs. An analytic field is
   deterministic, needs no external binary, and costs microseconds, which is
   what a proteome-scale screen can actually afford.
2. The descriptor is the *shell-resolved* power spectrum, not the usual single
   set of moments. A monopole tells you the pocket is positive; the radial
   profile of the dipole and quadrupole power tells you whether the positive
   charge forms a shell at the radius an IP6 would occupy, which is the shape
   of a sequestering site rather than a surface patch.
3. Complementarity is a correlation between spectra, so it is scale-free in
   the charge magnitude, and the absolute depth of the well is kept as a
   separate feature instead of being folded in.

None of this is validated. No number produced here has been compared with a
measured affinity, and the pre-registered questions are in
docs/DETECTOR_PLAN.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dataclass_field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.special import sph_harm_y

#: Formal charges at pH 7.4 on the atoms that carry them. Side-chain carboxyl
#: and amine groups are split across their equivalent oxygens or nitrogens so
#: the charge sits where the electron density is rather than on a single atom.
FORMAL_CHARGES: Dict[Tuple[str, str], float] = {
    ("ARG", "NH1"): 0.5, ("ARG", "NH2"): 0.5,
    ("LYS", "NZ"): 1.0,
    ("HIS", "ND1"): 0.05, ("HIS", "NE2"): 0.05,   # ~10% protonated at pH 7.4
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("TYR", "OH"): 0.0,
}

#: Backbone termini carry charge too, and a pocket formed by a loop end can
#: get a meaningful contribution from them.
TERMINAL_CHARGES: Dict[str, float] = {"OXT": -0.5, "N_TERM": 1.0}

#: Debye length in water at 150 mM monovalent salt, 310 K, in angstrom. This is
#: physiological ionic strength, not a fitted parameter.
DEBYE_LENGTH = 7.9

#: Relative permittivity used for the screened Coulomb field. 78.5 is bulk
#: water; a buried pocket is lower, which this deliberately does not model -
#: see the limitations note in docs/DETECTOR_PLAN.md.
PERMITTIVITY = 78.5

#: kT/e at 310 K in volts, so potentials come out in units of kT per unit
#: charge and are comparable with thermal energy without a second conversion.
KT_OVER_E = 0.0267

#: 332 kcal·Å/(mol·e²) is the Coulomb constant in the units structural biology
#: uses; dividing by kT at 310 K (0.616 kcal/mol) puts the field in kT/e.
COULOMB_KCAL = 332.0637
KT_KCAL = 0.6157

#: Shell radii in angstrom. An inositol hexakisphosphate is about 9-10 Å across,
#: so its phosphates sit near 4-5 Å from the ring centre; the shells bracket
#: that and extend far enough to see whether the charge keeps going (a surface
#: patch) or stops (a closed site).
DEFAULT_SHELLS: Tuple[float, ...] = (2.0, 3.5, 5.0, 6.5, 8.0, 10.0)

#: Harmonic degrees kept. l=0 is the monopole, l=1 the dipole, and so on. Past
#: l=4 the power is dominated by atomic granularity rather than by pocket shape.
DEFAULT_MAX_DEGREE = 4

#: Points per shell for the quadrature. The Fibonacci sphere is used instead of
#: a true Lebedev grid because it is exact to no order at all but is
#: near-uniform at every count, which matters more here: the descriptor is a
#: comparison between two fields sampled the same way.
#:
#: 768 is chosen from a measured convergence/cost curve, not by eye. The
#: discrepancy column is `spectrum_discrepancy` under four rotations - the
#: largest change in the power spectrum as a fraction of the spectrum's own
#: scale - and the cost column is single-core wall clock for the 78,920 pockets
#: of the deployed screen, on a 600-charge receptor with the default cutoff:
#:
#:   points   rotation   cost
#:     192   2.4e-02     0.6 h
#:     384   5.6e-03     1.6 h
#:     768   2.7e-03     3.0 h
#:    1536   1.6e-03     5.5 h
#:    3072   6.0e-04    10.6 h
#:
#: Sharded 272 ways, as study J was, 768 points is under a minute a shard, so
#: the accuracy is bought for nothing in wall clock.
DEFAULT_POINTS_PER_SHELL = 768

#: How many Debye lengths past the outermost shell a charge can sit and still
#: be included. Measured discrepancy against the uncut field, on the same
#: 600-charge receptor at 768 points:
#:
#:   cutoff   discrepancy   cost
#:      2.0      4.9e-02    33 ms
#:      3.0      2.5e-02    57 ms
#:      5.0      3.0e-04   125 ms
#:      8.0      2.5e-05   179 ms
#:     none      0           182 ms
#:
#: 5.0 is the first value whose error is an order of magnitude below the
#: quadrature error it sits beside, so the approximation is not the dominant
#: one. It saves 31% of the cost here and more on a larger receptor, because
#: what it removes grows with the protein while the pocket does not.
SCREENING_CUTOFF_DEBYE = 5.0


class FieldError(ValueError):
    """An input from which no field can be computed."""


def fibonacci_sphere(n_points: int) -> np.ndarray:
    """``(n_points, 3)`` near-uniform unit vectors.

    The spiral is deterministic, so two fields sampled at the same count are
    sampled at the same directions and their spectra are comparable.
    """
    if n_points < 4:
        raise FieldError(f"need at least 4 points on a sphere, got {n_points}")
    indices = np.arange(n_points, dtype=float) + 0.5
    # The golden-angle spiral: equal-area in z, golden-ratio advance in phi.
    cos_theta = 1.0 - 2.0 * indices / n_points
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    phi = np.pi * (1.0 + 5.0 ** 0.5) * indices
    return np.column_stack([
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta),
    ])


def real_spherical_harmonics(directions: np.ndarray, max_degree: int) -> List[np.ndarray]:
    """Orthonormal real spherical harmonics up to ``max_degree``.

    Returns one ``(n_directions, 2l+1)`` array per degree, ordered
    ``m = -l .. +l``.

    Built from `scipy.special.sph_harm_y` by the standard real combination
    rather than from hand-written Cartesian polynomials. The first version here
    *was* hand-written, and the l=4 forms were wrong in a way that is worth
    recording: every lower degree was exactly orthonormal, so the error was
    invisible in the obvious checks, and it showed up only as a rotation
    invariance that refused to converge - 130% error at l=4 whatever the
    sampling density, against 3e-4 at l=0. An invariance that does not improve
    with more points is not a quadrature problem, it is a basis problem.
    `test_the_basis_is_orthonormal_within_each_degree` now pins it.
    """
    if max_degree < 0:
        raise FieldError(f"max_degree cannot be negative, got {max_degree}")
    directions = np.asarray(directions, dtype=float)
    if directions.ndim != 2 or directions.shape[1] != 3:
        raise FieldError("directions must have shape (n, 3)")
    norms = np.linalg.norm(directions, axis=1)
    if np.any(norms == 0):
        raise FieldError("a direction vector cannot be zero")
    unit = directions / norms[:, None]
    # Polar angle from +z, azimuth in the xy plane.
    theta = np.arccos(np.clip(unit[:, 2], -1.0, 1.0))
    phi = np.arctan2(unit[:, 1], unit[:, 0])

    out: List[np.ndarray] = []
    for degree in range(max_degree + 1):
        columns = np.empty((len(unit), 2 * degree + 1))
        for order in range(-degree, degree + 1):
            m = abs(order)
            value = sph_harm_y(degree, m, theta, phi)
            if order == 0:
                real = value.real
            elif order > 0:
                real = np.sqrt(2.0) * ((-1.0) ** m) * value.real
            else:
                real = np.sqrt(2.0) * ((-1.0) ** m) * value.imag
            columns[:, order + degree] = real
        out.append(columns)
    return out


@dataclass(frozen=True)
class ChargeSet:
    """Point charges and where they sit."""

    coords: np.ndarray
    charges: np.ndarray
    labels: Tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.coords.ndim != 2 or self.coords.shape[1] != 3:
            raise FieldError("charge coordinates must have shape (n, 3)")
        if self.charges.shape != (self.coords.shape[0],):
            raise FieldError("one charge per coordinate")

    @property
    def net_charge(self) -> float:
        return float(np.sum(self.charges))


def protein_charges(arrays, *, include_termini: bool = True) -> ChargeSet:
    """Formal charges on a structure's ionizable groups.

    ``arrays`` is a :class:`cryptic_ip.analysis.structure_arrays.StructureArrays`.
    Only polymer atoms contribute: a bound ligand's own charge is not part of
    the pocket's field, and counting it would let a structure crystallised with
    a sulfate look like a site that attracts one.
    """
    coords: List[np.ndarray] = []
    charges: List[float] = []
    labels: List[str] = []
    polymer = np.asarray(arrays.is_polymer, dtype=bool)
    resnames = np.asarray(arrays.resnames, dtype=object)
    names = np.asarray(arrays.atom_names, dtype=object)
    for index in np.nonzero(polymer)[0]:
        resname = str(resnames[index]).upper()
        name = str(names[index]).upper()
        charge = FORMAL_CHARGES.get((resname, name))
        if charge is None and include_termini:
            charge = TERMINAL_CHARGES.get(name)
        if charge is None or charge == 0.0:
            continue
        coords.append(arrays.coords[index])
        charges.append(float(charge))
        labels.append(f"{resname}:{name}")
    if not coords:
        return ChargeSet(np.zeros((0, 3)), np.zeros(0), ())
    return ChargeSet(np.asarray(coords, dtype=float), np.asarray(charges, dtype=float), tuple(labels))


def screened_potential(
    charges: ChargeSet,
    points: np.ndarray,
    *,
    debye_length: float = DEBYE_LENGTH,
    permittivity: float = PERMITTIVITY,
    min_distance: float = 1.0,
) -> np.ndarray:
    """Debye-Hückel potential at ``points``, in units of kT per unit charge.

    The ``min_distance`` floor keeps a sample point that lands on top of a
    charge from producing an infinite value. A field with a singularity in it
    is not a descriptor; the floor is the atom's own radius, below which the
    continuum picture has no meaning anyway.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise FieldError("points must have shape (n, 3)")
    if debye_length <= 0:
        raise FieldError(f"debye_length must be positive, got {debye_length}")
    if len(charges.charges) == 0:
        return np.zeros(len(points))
    delta = points[:, None, :] - charges.coords[None, :, :]
    distance = np.sqrt(np.sum(delta * delta, axis=2))
    distance = np.maximum(distance, min_distance)
    screening = np.exp(-distance / debye_length)
    scale = COULOMB_KCAL / (permittivity * KT_KCAL)
    return scale * np.sum(charges.charges[None, :] * screening / distance, axis=1)


@dataclass(frozen=True)
class FieldSpectrum:
    """Shell-resolved power spectrum of a potential around one centre."""

    centre: Tuple[float, float, float]
    shells: Tuple[float, ...]
    power: np.ndarray                    # (n_shells, max_degree + 1)
    shell_mean: np.ndarray               # (n_shells,) the l=0 term, in kT/e
    max_degree: int

    @property
    def monopole(self) -> float:
        """Mean potential on the innermost shell: is the pocket positive at all."""
        return float(self.shell_mean[0])

    @property
    def well_depth(self) -> float:
        """Mean potential over all shells, the absolute attractiveness."""
        return float(np.mean(self.shell_mean))

    def feature_vector(self) -> np.ndarray:
        """The spectrum flattened, for use as a descriptor."""
        return self.power.reshape(-1)

    def as_dict(self) -> Dict[str, object]:
        return {
            "centre": list(self.centre),
            "shells": list(self.shells),
            "max_degree": self.max_degree,
            "monopole_kt": self.monopole,
            "well_depth_kt": self.well_depth,
            "shell_mean_kt": [float(v) for v in self.shell_mean],
            "power": [[float(v) for v in row] for row in self.power],
        }


def field_spectrum(
    charges: ChargeSet,
    centre: Sequence[float],
    *,
    shells: Sequence[float] = DEFAULT_SHELLS,
    max_degree: int = DEFAULT_MAX_DEGREE,
    points_per_shell: int = DEFAULT_POINTS_PER_SHELL,
    debye_length: float = DEBYE_LENGTH,
    cutoff_debye: Optional[float] = SCREENING_CUTOFF_DEBYE,
) -> FieldSpectrum:
    """Expand the potential around ``centre`` into rotation-invariant power.

    For each shell radius the potential is sampled on a near-uniform sphere,
    projected onto the real spherical harmonics of each degree, and the squared
    coefficients are summed within the degree. That within-degree sum is what a
    rotation cannot change: rotating the structure mixes the 2l+1 coefficients
    of a degree orthogonally, so their sum of squares is fixed.
    """
    centre = np.asarray(centre, dtype=float).reshape(3)
    shells = tuple(float(r) for r in shells)
    if not shells:
        raise FieldError("at least one shell radius is required")
    if any(r <= 0 for r in shells):
        raise FieldError("shell radii must be positive")
    directions = fibonacci_sphere(points_per_shell)
    harmonics = real_spherical_harmonics(directions, max_degree)
    weight = 4.0 * np.pi / points_per_shell   # quadrature weight of one point

    # Screening makes this a local calculation: a charge far enough outside the
    # outermost shell contributes less than the quadrature error. Dropping
    # those turns a cost that grows with the whole protein into one that grows
    # with the pocket's neighbourhood, which is what makes a 2,352-residue
    # candidate cost the same as a small one.
    if cutoff_debye is not None and len(charges.charges):
        reach = max(shells) + cutoff_debye * debye_length
        offsets = charges.coords - centre[None, :]
        near = np.sum(offsets * offsets, axis=1) <= reach * reach
        if not np.all(near):
            charges = ChargeSet(
                charges.coords[near],
                charges.charges[near],
                tuple(np.asarray(charges.labels, dtype=object)[near]) if charges.labels else (),
            )

    power = np.zeros((len(shells), max_degree + 1))
    shell_mean = np.zeros(len(shells))
    for shell_index, radius in enumerate(shells):
        points = centre[None, :] + radius * directions
        potential = screened_potential(charges, points, debye_length=debye_length)
        shell_mean[shell_index] = float(np.mean(potential))
        for degree, basis in enumerate(harmonics):
            coefficients = weight * (basis * potential[:, None]).sum(axis=0)
            power[shell_index, degree] = float(np.sum(coefficients * coefficients))
    return FieldSpectrum(
        centre=(float(centre[0]), float(centre[1]), float(centre[2])),
        shells=shells,
        power=power,
        shell_mean=shell_mean,
        max_degree=max_degree,
    )


@dataclass(frozen=True)
class LigandTemplate:
    """A reference ligand as a charge distribution, with no coordinates kept.

    Only the radial charge profile is retained, because that is all the
    spectrum comparison can use. Keeping coordinates would invite someone to
    align them, which is the trap study N fell into.
    """

    name: str
    net_charge: float
    shell_charge: Tuple[float, ...]       # charge within each shell radius
    radius_of_gyration: float
    note: str = ""

    def as_charge_set(self, *, points_per_shell: int = 64) -> ChargeSet:
        """A spherically-averaged realisation of the template.

        The charge in each shell is spread over a sphere of that radius. This
        is a *model* of the ligand, not the ligand: an inositol phosphate is
        not spherically symmetric. It is the right model for this descriptor
        because the descriptor is itself rotation-averaged within each degree.
        """
        coords: List[np.ndarray] = []
        charges: List[float] = []
        previous = 0.0
        for radius, cumulative in zip(DEFAULT_SHELLS, self.shell_charge):
            increment = cumulative - previous
            previous = cumulative
            if increment == 0.0:
                continue
            directions = fibonacci_sphere(points_per_shell)
            for direction in directions:
                coords.append(radius * direction)
                charges.append(increment / points_per_shell)
        if not coords:
            return ChargeSet(np.zeros((0, 3)), np.zeros(0), ())
        return ChargeSet(np.asarray(coords), np.asarray(charges), (self.name,))


#: The competitor panel. Every entry is a polyanion this project's own
#: specificity work found indistinguishable from IP6 by the deployed
#: descriptors, plus IP6 itself and two smaller anions. The numbers are formal
#: charge at pH 7.4 and the cumulative charge inside each DEFAULT_SHELLS
#: radius, computed from the ligand's own geometry - they are inputs to the
#: comparison, not results of it.
LIGAND_TEMPLATES: Dict[str, LigandTemplate] = {
    "IP6": LigandTemplate(
        name="IP6", net_charge=-9.0,
        # Inositol ring inside ~2.5 Å, six phosphates centred near 4.5 Å.
        shell_charge=(0.0, -1.0, -6.5, -9.0, -9.0, -9.0),
        radius_of_gyration=4.3,
        note="inositol hexakisphosphate at pH 7.4; the charge is a shell, not a point",
    ),
    "IP3": LigandTemplate(
        name="IP3", net_charge=-5.0,
        shell_charge=(0.0, -0.8, -3.6, -5.0, -5.0, -5.0),
        radius_of_gyration=3.9,
        note="inositol trisphosphate; same ring, half the shell occupied",
    ),
    "ATP": LigandTemplate(
        name="ATP", net_charge=-4.0,
        # Extended, not a shell: charge keeps increasing with radius.
        shell_charge=(0.0, -0.4, -1.6, -3.0, -4.0, -4.0),
        radius_of_gyration=5.6,
        note="adenosine triphosphate; the charge is at one end of a long molecule",
    ),
    "PAPS": LigandTemplate(
        name="PAPS", net_charge=-4.0,
        shell_charge=(0.0, -0.3, -1.4, -2.8, -4.0, -4.0),
        radius_of_gyration=6.1,
        note="the hardest competitor: a sulfate and two phosphates on a nucleotide",
    ),
    "UDP_GLUCOSE": LigandTemplate(
        name="UDP_GLUCOSE", net_charge=-2.0,
        shell_charge=(0.0, -0.2, -0.9, -1.6, -2.0, -2.0),
        radius_of_gyration=5.8,
        note="a sugar nucleotide; charge diffuse and low",
    ),
    "SULFATE": LigandTemplate(
        name="SULFATE", net_charge=-2.0,
        shell_charge=(-1.0, -2.0, -2.0, -2.0, -2.0, -2.0),
        radius_of_gyration=1.5,
        note="a crystallisation additive that occupies real anion sites",
    ),
    "PHOSPHATE": LigandTemplate(
        name="PHOSPHATE", net_charge=-2.0,
        shell_charge=(-1.2, -2.0, -2.0, -2.0, -2.0, -2.0),
        radius_of_gyration=1.6,
        note="inorganic phosphate; the smallest thing a basic pocket will bind",
    ),
    "HEPARIN_TETRA": LigandTemplate(
        name="HEPARIN_TETRA", net_charge=-8.0,
        # Comparable net charge to IP6 but extended, which is the point: net
        # charge alone cannot separate these two and the shape term must.
        shell_charge=(0.0, -0.5, -2.0, -4.0, -6.5, -8.0),
        radius_of_gyration=7.8,
        note="same order of charge as IP6, spread along a chain instead of a shell",
    ),
}


def complementarity(
    pocket: FieldSpectrum,
    template: LigandTemplate,
    *,
    points_per_shell: Optional[int] = None,
) -> Dict[str, float]:
    """How well a pocket's field shape matches what a ligand's charge wants.

    Two numbers, deliberately kept apart:

    * ``shape`` - the correlation between the pocket's shell-resolved power
      spectrum and the spectrum the template's own charge distribution would
      produce. Scale-free: a weakly but correctly shaped pocket scores high.
    * ``depth`` - the pocket's mean potential weighted by where the template
      puts its charge, in kT. This one has units and a sign, and it is what
      says whether binding is favourable at all.

    Folding them into one number is what the deployed composite score does, and
    it is why a shallow surface patch with many arginines can outscore a real
    buried site.
    """
    reference = field_spectrum(
        template.as_charge_set(points_per_shell=points_per_shell or 64),
        centre=(0.0, 0.0, 0.0),
        shells=pocket.shells,
        max_degree=pocket.max_degree,
        points_per_shell=len(fibonacci_sphere(max(4, pocket.power.shape[0] * 32))),
    )
    a = pocket.feature_vector()
    b = reference.feature_vector()
    shape = _correlation(a, b)

    # The template's charge increment per shell, as weights.
    increments = np.diff((0.0,) + tuple(template.shell_charge))
    weights = -increments            # positive where the ligand puts anion charge
    total = float(np.sum(np.abs(weights)))
    depth = float(np.sum(weights * pocket.shell_mean) / total) if total else 0.0
    return {
        "shape": float(shape),
        "depth_kt": depth,
        "net_charge": float(template.net_charge),
        "template_rg": float(template.radius_of_gyration),
    }


def spectrum_discrepancy(a: np.ndarray, b: np.ndarray) -> float:
    """How far apart two power spectra are, as a fraction of the larger one's scale.

    The obvious metric - the largest *relative* difference, component by
    component - is the wrong one here, and measuring with it was a mistake
    worth leaving on the record. A spectrum spans orders of magnitude between
    the monopole and l=4, and several of its components are near-cancelling, so
    a 10% relative error on a component whose value is 0.0098 against a
    spectrum whose scale is 35 reads as a failure while being nothing at all.

    So the difference is normalised by the scale of the spectrum rather than by
    each component: ``max|a - b| / max(a)``. That is a signal-to-noise measure
    and it is what the convergence figures in this module's constants are
    quoted in.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise FieldError("spectra must have the same shape to be compared")
    scale = float(np.max(np.abs(a)))
    if scale == 0.0:
        return float(np.max(np.abs(b)))
    return float(np.max(np.abs(a - b)) / scale)


def _correlation(a: np.ndarray, b: np.ndarray) -> float:
    """Pearson correlation on log-power, or 0.0 when either side is flat.

    Power spans orders of magnitude between the monopole and l=4, so a
    correlation on raw power is a correlation on the monopole alone. The log
    puts the degrees on comparable footing, which is the whole reason for
    keeping the higher degrees.
    """
    if a.shape != b.shape or a.size == 0:
        raise FieldError("spectra must have the same shape to be compared")
    floor = 1e-12
    la = np.log10(np.maximum(a, floor))
    lb = np.log10(np.maximum(b, floor))
    la = la - la.mean()
    lb = lb - lb.mean()
    denominator = float(np.linalg.norm(la) * np.linalg.norm(lb))
    if denominator == 0.0:
        return 0.0
    return float(np.dot(la, lb) / denominator)


@dataclass
class PocketFieldReport:
    """Everything the field module says about one pocket."""

    spectrum: FieldSpectrum
    per_ligand: Dict[str, Dict[str, float]] = dataclass_field(default_factory=dict)

    def best_ligand(self) -> Tuple[str, float]:
        """Which template the pocket's field shape fits best, and by how much."""
        if not self.per_ligand:
            raise FieldError("no ligands scored")
        ranked = sorted(self.per_ligand.items(), key=lambda item: -item[1]["shape"])
        return ranked[0][0], float(ranked[0][1]["shape"])

    def as_dict(self) -> Dict[str, object]:
        return {"spectrum": self.spectrum.as_dict(), "ligands": self.per_ligand}


def describe_pocket(
    arrays,
    centre: Sequence[float],
    *,
    templates: Optional[Mapping[str, LigandTemplate]] = None,
    shells: Sequence[float] = DEFAULT_SHELLS,
    max_degree: int = DEFAULT_MAX_DEGREE,
    points_per_shell: int = DEFAULT_POINTS_PER_SHELL,
) -> PocketFieldReport:
    """The field description of one pocket, against the whole competitor panel."""
    charges = protein_charges(arrays)
    spectrum = field_spectrum(
        charges, centre, shells=shells, max_degree=max_degree, points_per_shell=points_per_shell,
    )
    panel = dict(templates if templates is not None else LIGAND_TEMPLATES)
    report = PocketFieldReport(spectrum=spectrum)
    for name, template in sorted(panel.items()):
        report.per_ligand[name] = complementarity(spectrum, template)
    return report


def rotate(coords: np.ndarray, axis: Sequence[float], angle: float) -> np.ndarray:
    """Rotate ``coords`` about ``axis`` by ``angle`` radians.

    Here rather than in a test helper because the invariance this module claims
    is a property of the descriptor, and the property is worth being able to
    check from anywhere.
    """
    axis = np.asarray(axis, dtype=float)
    norm = np.linalg.norm(axis)
    if norm == 0:
        raise FieldError("rotation axis cannot be zero")
    axis = axis / norm
    cos, sin = np.cos(angle), np.sin(angle)
    cross = np.array([
        [0.0, -axis[2], axis[1]],
        [axis[2], 0.0, -axis[0]],
        [-axis[1], axis[0], 0.0],
    ])
    matrix = cos * np.eye(3) + sin * cross + (1.0 - cos) * np.outer(axis, axis)
    return np.asarray(coords, dtype=float) @ matrix.T


def iter_shell_labels(shells: Iterable[float], max_degree: int) -> List[str]:
    """Column names for a flattened spectrum, so a feature table is readable."""
    return [f"shell{radius:g}_l{degree}" for radius in shells for degree in range(max_degree + 1)]
