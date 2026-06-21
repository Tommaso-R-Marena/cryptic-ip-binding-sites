# Run the full pipeline on Google Colab

## Quick start

1. Open [notebooks/Colab_Full_Pipeline_Run.ipynb](../notebooks/Colab_Full_Pipeline_Run.ipynb) in Colab:
   - Upload the notebook, or
   - Use **File → Open notebook → GitHub** and paste:
     `Tommaso-R-Marena/cryptic-ip-binding-sites`
2. Select a **CPU runtime** (GPU is not required).
3. Run all cells. Start with `N_PROTEINS = 25`, then set `N_PROTEINS = 500` for the full pilot.

## What runs

| Step | Script / module | Typical time (500 proteins) |
|------|-----------------|-----------------------------|
| Install | `scripts/colab_install.sh` | ~5 min |
| Tier-1 gate | `ValidationSuite` | ~3 min |
| ML retrain | `scripts/train_ml_classifier.py` | ~30–90 min |
| Yeast pilot | `scripts/run_yeast_pilot_screen.py` | ~1–3 hours |
| Publication | `scripts/run_publication_package.py` | ~10 min |

Add `--with-electrostatics` / `WITH_ELECTROSTATICS = True` for APBS — roughly **3–5× slower**.

## One-liner (after clone)

```bash
bash scripts/colab_install.sh
python scripts/run_yeast_pilot_screen.py --n-proteins 500 --workers 2
python scripts/train_ml_classifier.py --skip-build-dataset
python scripts/run_publication_package.py --skip-dataset-build --skip-figures
```

## Saving results

The Colab notebook zips `results/colab_run/` for download. For long runs, mount Google Drive:

```python
from google.colab import drive
drive.mount('/content/drive')
OUTPUT = Path('/content/drive/MyDrive/cryptic_ip_results')
```

## MD validation (optional)

```bash
pip install openmm mdtraj
python scripts/run_md_pilot_validation.py --top-n 5 --production-ns 1.0
```

OpenMM is not available on all Colab runtimes; use a local machine or HPC for production MD.
