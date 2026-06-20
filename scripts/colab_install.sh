#!/usr/bin/bash
# Bootstrap cryptic-ip-binding-sites on Google Colab (CPU or GPU runtime).
# Usage in Colab: !bash scripts/colab_install.sh

set -euo pipefail

echo "==> Installing conda/mamba if needed..."
if ! command -v conda &>/dev/null; then
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh
  bash /tmp/miniforge.sh -b -p /usr/local/miniforge3
  export PATH="/usr/local/miniforge3/bin:$PATH"
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "==> Installing structural biology tools..."
conda install -y -c conda-forge fpocket freesasa apbs pdb2pqr propka

echo "==> Installing Python dependencies..."
pip install -q -r requirements.txt
pip install -q -e .

echo "==> Verifying fpocket..."
fpocket -h >/dev/null && echo "fpocket OK"

echo "==> Done. Run notebooks/Colab_Full_Pipeline_Run.ipynb or:"
echo "    python scripts/run_yeast_pilot_screen.py --n-proteins 500 --workers 2"
echo "    python scripts/train_ml_classifier.py --skip-build-dataset"
echo "    python scripts/run_publication_package.py --skip-dataset-build --with-electrostatics"
