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
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "==> Installing structural biology tools..."
conda install -y -c conda-forge fpocket freesasa apbs pdb2pqr propka pip

echo "==> Installing Python dependencies..."
python -m pip install --upgrade pip setuptools wheel
python -m pip install -q -r requirements.txt
python -m pip install -q -e .
python -m pip install -q ipywidgets tqdm

echo "==> Verifying tools..."
fpocket -h >/dev/null && echo "  fpocket OK"
freesasa -h >/dev/null 2>&1 && echo "  freesasa OK" || echo "  freesasa: check manually"
apbs --version >/dev/null 2>&1 && echo "  apbs OK" || echo "  apbs: check manually"

echo "==> Verifying Python CLI..."
python -c "import cryptic_ip; from cryptic_ip.validation import ValidationSuite; print('  cryptic_ip OK')"

python scripts/colab_env.py

if command -v cryptic-ip &>/dev/null; then
  cryptic-ip check-dependencies || echo "  (check-dependencies reported missing optional tools — continue if fpocket works)"
fi

echo ""
echo "==> Done. Run the full pipeline:"
echo "    python scripts/colab_run_all.py --preset quick"
echo "    python scripts/colab_run_all.py --preset pilot   # ~2-4 h on Colab"
echo ""
echo "Or open notebooks/Colab_Full_Pipeline_Run.ipynb and Run All."
