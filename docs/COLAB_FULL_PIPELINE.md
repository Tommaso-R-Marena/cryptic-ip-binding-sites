# Run the full pipeline on Google Colab

## One-click open

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Tommaso-R-Marena/cryptic-ip-binding-sites/blob/main/notebooks/Colab_Full_Pipeline_Run.ipynb)

## What “entire pipeline” means

The notebook runs **every automated stage** in one pass via `scripts/colab_run_all.py`:

| # | Stage | Output |
|---|--------|--------|
| 1 | Environment + `fpocket` / APBS | verified tools |
| 2 | Download tier-1 PDBs (1ZY7, 1MAI, 1BWN) | `data/validation/` |
| 3 | Tier-1 gate (ADAR2 vs PLCδ1) | `validation/` — **must pass** |
| 4 | ML classifier benchmark | `ml_training/` + `models/*.pkl` |
| 5 | Yeast AlphaFold pilot screen | `yeast_pilot/` hits + summary |
| 6 | Publication package + **figures** | `publication/figures/*.pdf` |
| 7 | Optional MD pilot (`full` preset) | `md_validation/` |
| 8 | Zip archive | `colab_pipeline_results.zip` |

## Quick start

1. Open the notebook link above (or upload `notebooks/Colab_Full_Pipeline_Run.ipynb`).
2. **Runtime → Change runtime type → CPU** (High-RAM for `pilot`).
3. Run all cells.
4. Start with **`PRESET = 'quick'`** (25 proteins, ~15–30 min).
5. When satisfied, set **`PRESET = 'pilot'`** (500 proteins, ~2–4 h).

## Presets

| Preset | Proteins | ML | Figures | MD | Typical time |
|--------|----------|----|---------|----|--------------|
| `quick` | 25 | yes | yes | no | 15–30 min |
| `pilot` | 500 | yes | yes | no | 2–4 h |
| `full` | 500 | yes | yes | yes | 3–5 h |

## Google Drive (recommended for long runs)

In the notebook, set `USE_DRIVE = True` before running. Results persist under  
`MyDrive/cryptic_ip_results/`.

## Command-line equivalent (after clone)

```bash
bash scripts/colab_install.sh
python scripts/colab_run_all.py --preset quick
python scripts/colab_run_all.py --preset pilot --output-dir results/colab_run
```

### Useful flags

```bash
python scripts/colab_run_all.py --preset pilot \
  --score-threshold 0.75 \
  --skip-download \          # reuse AlphaFold cache
  --with-electrostatics \    # slower, includes APBS
  --skip-md
```

## Screening parameters (June 2026)

The yeast pilot uses **strict cryptic filters** (fixes the earlier 99% false-hit rate):

- Composite score ≥ **0.75**
- Pocket SASA ≤ **10 Ų**
- ≥ **4** basic residues (Arg/Lys/His) from fpocket lining residues
- Volume **300–800 Ų**
- Pocket pLDDT ≥ **70**

## Troubleshooting

| Issue | Fix |
|-------|-----|
| `fpocket not found` | Re-run `!bash scripts/colab_install.sh`; confirm **Pipeline Python** shows `/usr/local/miniforge3/bin/python` |
| Pipeline uses `/usr/bin/python3` | Re-run install cell; notebook must print `Pipeline Python: .../miniforge3/bin/python` |
| Tier-1 gate fails | Ensure 1ZY7/1MAI downloaded; set `WITH_ELECTROSTATICS = False` |
| Colab disconnects | Use Drive mount + `SKIP_DOWNLOAD = True` to resume yeast pilot |
| MD fails | Expected on some runtimes; use `--skip-md` or `PRESET = pilot` |
| Out of memory | Reduce `N_PROTEINS` or use High-RAM runtime |

## Resume yeast screen after disconnect

The yeast script writes checkpoints to `yeast_pilot/screen_checkpoint.json`. Re-run with  
`SKIP_DOWNLOAD = True` and the same `OUTPUT_DIR` — already-finished proteins are skipped.

## After Colab

Merge results into the repo publication folder or manuscript:

- `publication/RESULTS_SUMMARY.md`
- `publication/figures/Figure*.pdf`
- `yeast_pilot/yeast_pilot_hits.csv`
