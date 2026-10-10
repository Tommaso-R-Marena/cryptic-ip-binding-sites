"""The next-generation pocket detector.

Six modules, each answering one of the findings that limit the deployed screen.
The pre-registered questions that would turn any of this into a result are in
`docs/DETECTOR_PLAN.md`; until those run, every output is a hypothesis ordering.

* `cavity`      - LIGSITE-style detection in process, so fpocket is optional
                  and hundreds of conformers are affordable.
* `ensemble`    - persistence and crypticity across elastic-network conformers.
                  The only part that can find a site closed in the apo model,
                  which is what "cryptic" means.
* `openability` - what it would cost to open a pocket, rather than whether it
                  is already open.
* `field`       - a rotation-invariant electrostatic power spectrum and its
                  complementarity with an IP6 charge shell, replacing a count
                  of basic residues that cannot tell IP6 from ATP.
* `specificity` - IP6 against a competitor panel, with `polyanion_generic` as a
                  reportable outcome rather than a forced candidate call.
* `monotone`    - a monotone additive score that abstains on missing input
                  instead of failing a protein for a measurement nobody took.
* `pu`          - positive-unlabelled learning, group-held-out calibration,
                  conformal sets, and the minimum detectable effect.
"""

from cryptic_ip.detect.cavity import Cavity, CavityMap, detect_cavities
from cryptic_ip.detect.ensemble import CavityTrack, EnsembleReport, analyse_ensemble, scan_with_nodes
from cryptic_ip.detect.field import (
    LIGAND_TEMPLATES, FieldSpectrum, LigandTemplate, complementarity, describe_pocket, field_spectrum,
)
from cryptic_ip.detect.monotone import DEFAULT_FEATURES, Feature, MonotoneScorer, ScoreInterval
from cryptic_ip.detect.openability import NormalModes, normal_modes, perturbed_conformers
from cryptic_ip.detect.pu import (
    ConformalSets, LabelFrequency, PUFit, conformal_sets, fit_pu, group_auc, minimum_detectable_auc,
)
from cryptic_ip.detect.specificity import SpecificityVerdict, summarise, verdict

__all__ = [
    "Cavity",
    "CavityMap",
    "CavityTrack",
    "ConformalSets",
    "DEFAULT_FEATURES",
    "EnsembleReport",
    "Feature",
    "FieldSpectrum",
    "LIGAND_TEMPLATES",
    "LabelFrequency",
    "LigandTemplate",
    "MonotoneScorer",
    "NormalModes",
    "PUFit",
    "ScoreInterval",
    "SpecificityVerdict",
    "analyse_ensemble",
    "complementarity",
    "conformal_sets",
    "describe_pocket",
    "detect_cavities",
    "field_spectrum",
    "fit_pu",
    "group_auc",
    "minimum_detectable_auc",
    "normal_modes",
    "perturbed_conformers",
    "scan_with_nodes",
    "summarise",
    "verdict",
]
