# Manuscript submission package

This folder holds the manuscript draft and links to auto-generated artifacts in `results/publication/`.

## Quick start

```bash
# Regenerate all publication artifacts
python scripts/run_publication_package.py --output-dir results/publication --skip-dataset-build

# Export supplementary tables only
python scripts/export_supplementary_tables.py --publication-dir results/publication
```

## Artifact map

| Manuscript section | Source in repo |
|--------------------|----------------|
| Results & Discussion (draft) | `manuscript/RESULTS_DRAFT.md` |
| Methods | `results/publication/METHODS_AUTO.md`, `docs/METHODS.md` |
| Results — controls | `results/publication/RESULTS_SUMMARY.md` |
| Candidate dossiers | `results/candidates/candidate_dossiers.md` |
| Figure 1 | `results/publication/figures/Figure1_Overview.pdf` |
| Figure 2 | `results/publication/figures/Figure2_Comparative_Proteomics.pdf` |
| Figure 3 | `results/publication/figures/Figure3_Top_Candidates.pdf` |
| Figure legends | `results/publication/figures/figure_legends.txt` |
| Table S1–S6 | `results/publication/supplementary/` |
| Provenance | `results/publication/provenance_manifest.jsonld` |

## Submission checklist

### Done (automation)

- [x] Tier-1 validation gate (ADAR2 vs PLCδ1)  
- [x] RCSB validation dataset (136 structures)  
- [x] Control benchmark CSVs  
- [x] ML vs threshold comparison  
- [x] Auto-methods and provenance manifest  
- [x] Figure generation pipeline (PDF/PNG)  
- [x] Supplementary table export  

### In progress

- [x] Real strict-filter yeast pilot screen (499 structures → 1 candidate, 0.2%)  
- [x] Candidate dossier generator (`scripts/characterize_candidates.py`) + pilot candidate P07264/LEU1  
- [x] Results & Discussion working draft (`manuscript/RESULTS_DRAFT.md`)  

### Remaining (science + writing)

- [ ] Full yeast / human / dicty proteome screens (pilot proven; compute-bound)  
- [ ] Phase-5 comparative statistics (Wilson CI, Spearman, GO enrichment)  
- [ ] Manual review + experimental follow-up of candidates  
- [ ] Introduction section + final Discussion polish  
- [ ] Zenodo archive + DOI in `CITATION.cff`  
- [ ] MD validation on top candidates (optional)  

## Suggested narrative structure

1. **Introduction** — cryptic IP sites as structural cofactors; ADAR2 precedent (Macbeth et al., 2005).  
2. **Results** — tier-1 gate → RCSB benchmark → proteome screen → comparative analysis.  
3. **Discussion** — threshold vs ML; crystal artifacts; evolutionary IP6 concentration hypothesis.  
4. **Methods** — paste from `METHODS_AUTO.md` + software versions from provenance.  

## Citation

Update `CITATION.cff` with Zenodo DOI before submission. Software citation:

```bibtex
@software{cryptic_ip_binding_sites,
  author = {Marena, Tommaso R.},
  title = {Cryptic IP Binding Sites: Systematic Discovery Pipeline},
  year = {2026},
  url = {https://github.com/Tommaso-R-Marena/cryptic-ip-binding-sites}
}
```
