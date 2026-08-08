"""Testing utilities: deterministic synthetic structures and benchmark builders.

These helpers exist so that every numerical stage of the pipeline - SASA, burial,
enclosure, pocket detection, labelling, feature extraction, nested
cross-validation - can be exercised and verified without network access to the
PDB, and against ground truth that is known by construction.

They are a *machinery* test, not a biology test. Synthetic structures verify that
the code measures what it claims to measure; only real structures can validate
the biological conclusions.
"""

from .synthetic import (
    SyntheticStructureSpec,
    build_synthetic_benchmark,
    write_synthetic_structure,
)

__all__ = [
    "SyntheticStructureSpec",
    "build_synthetic_benchmark",
    "write_synthetic_structure",
]
