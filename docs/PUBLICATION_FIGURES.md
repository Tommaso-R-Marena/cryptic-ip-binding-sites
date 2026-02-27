# Publication Figure Automation

Use the script below to generate publication-ready composite figures for manuscripts.

```bash
python scripts/generate_publication_figures.py --config scripts/figure_config_template.yaml --figures all
```

## Inputs
The generator expects CSV files with the following columns:

- `figure1.validation_roc_csv`:
  - `fpr`, `tpr`, `model`
- `figure2.hit_rates_csv`:
  - `organism`, `hit_rate`
  - Optional significance columns: `group1`, `group2`, `stars`
- `figure2.correlation_csv`:
  - `ip6_concentration`, `hit_rate`
- `figure2.go_heatmap_csv`:
  - first column = GO term labels, remaining columns = enrichment scores
- `figure2.phylo_conservation_csv`:
  - either (`species`, `conservation_score`, `candidate`) for point plots
  - or matrix format for heatmap fallback
- `figure3.top_candidates_csv`:
  - `candidate`, `structure_image`, `electrostatic_image`

## Outputs
For each figure, exports are written to `output_dir`:

- Vector PDF (`FigureX_*.pdf`)
- 600 DPI PNG (`FigureX_*.png`)
- CMYK PNG clone if Pillow is installed (`FigureX_*_CMYK.png`)
- `figure_legends.txt`

## Notes
- Figure 1 panel A uses PyMOL when available on PATH (`pymol -cq`).
- If PyMOL is unavailable, the script renders a placeholder panel so automation can still run.
- Customize fonts/colors/layout through `scripts/nature_publication.mplstyle`.
