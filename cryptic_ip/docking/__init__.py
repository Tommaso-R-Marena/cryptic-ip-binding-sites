"""Docking of inositol phosphate ligands (docs/REDOCKING_PLAN.md, docs/ARRESTIN_PLAN.md).

``rmsd``      symmetry-corrected ligand RMSD without superposition
``ligand``    ligands built from the Chemical Component Dictionary, protonation
              states, configuration checks, start poses and PDBQT
``receptor``  receptor preparation: non-polymer atoms stripped, PDB2PQR/PROPKA
              protonation, AutoDock typing
``engine``    AutoDock Vina runs (Vina, Vinardo, AD4) and pose bookkeeping
"""
