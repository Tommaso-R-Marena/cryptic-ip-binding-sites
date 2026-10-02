# Study N, amendment 2: the template library was larger than the plan stated

Dated 2026-10-01. This corrects a factual error in `docs/TEMPLATE_PLAN.md` and in every
place the figure was repeated. It changes no result of study N.

## The error

The plan described the template library as the eligible, complete IHP copies with
`symmetry_contact == False`, and reported **131 copies across 77 PDB entries and 14 strict
homology groups**.

The code (`scripts/template_fit.py::cmd_templates`) and the workflow both filter on
`symmetry_contact` **not equal to true**, which keeps copies whose flag is empty. The library
that actually ran is therefore **222 copies across 137 PDB entries in 17 strict homology
groups**.

## Why the code's filter is the right one

The 91 additional copies are **all ELECTRON MICROSCOPY** entries. The symmetry-contact check
asks whether a ligand touches a polymer atom of a crystallographic symmetry mate; a cryo-EM
structure has no unit cell, so the test is *inapplicable* to it rather than failed. Excluding
those copies would discard sites for a reason that does not apply to them — and would discard
exactly the large assemblies in which this project expects IP6 to act as a structural
cofactor, which `docs/DATA_ACQUISITION.md` gives as the reason the pipeline accepts all
experimental methods rather than X-ray alone.

So the implementation was correct and the plan's wording was wrong and narrower than both the
code and the intent.

## What was corrected

- `results/template/PROVENANCE.md` and `formal/PROVENANCE.md`, where the 131/77/14 figures
  were repeated.
- `docs/TEMPLATE2_PLAN.md` (study O) states the filter explicitly, in words that match the
  code, and gives the corrected counts.

Study N's reported results stand: they were computed with the 222-copy library throughout.
No figure in `results/template/template.json` changes, because nothing about the run changes
— only the description of it.
