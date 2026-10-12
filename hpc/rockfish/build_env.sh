#!/usr/bin/env bash
# Build the study environment on a Rockfish login node, and record what it got.
#
#   bash hpc/rockfish/build_env.sh            # create or update `cryptic-ip`
#   ENV_NAME=cryptic-ip-gpu bash hpc/rockfish/build_env.sh --gpu
#
# The recipe pins what it can. Three binaries (mmseqs2, foldseek, apbs) have no
# usable pin, so this script writes the resolved versions to
# hpc/rockfish/env_manifest.json: a run can then state what it used, which is
# the difference between an unpinned dependency and an unknown one.
set -euo pipefail

ENV_NAME="${ENV_NAME:-cryptic-ip}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WANT_GPU=0
for arg in "$@"; do
  case "$arg" in
    --gpu) WANT_GPU=1 ;;
    *) echo "unknown argument: $arg" >&2; exit 2 ;;
  esac
done

if ! command -v conda >/dev/null 2>&1; then
  echo "conda is not on PATH. On Rockfish: module load anaconda3" >&2
  exit 1
fi
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"

SOLVER=conda
if command -v mamba >/dev/null 2>&1; then
  SOLVER=mamba
fi
echo "== building $ENV_NAME with $SOLVER"

if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "== $ENV_NAME exists; updating it in place"
  "$SOLVER" env update -n "$ENV_NAME" -f "$REPO_ROOT/hpc/rockfish/environment.yml" --prune
else
  "$SOLVER" env create -n "$ENV_NAME" -f "$REPO_ROOT/hpc/rockfish/environment.yml"
fi

conda activate "$ENV_NAME"
python -m pip install -e "$REPO_ROOT"

if [ "$WANT_GPU" -eq 1 ]; then
  # Co-folding and MD only. Kept out of the default environment so a CPU study
  # does not drag a CUDA stack behind it, and left unpinned here on purpose:
  # which build works depends on the driver on the GPU nodes, which this
  # repository has not measured. Record what you install.
  echo "== GPU extras: install against the driver the GPU nodes actually run"
  echo "   openmm:  $SOLVER install -n $ENV_NAME -c conda-forge openmm cudatoolkit"
  echo "   boltz:   python -m pip install boltz"
  echo "   Then re-run: python scripts/rockfish.py preflight inside a GPU job."
fi

echo "== recording what was resolved"
python - "$REPO_ROOT/hpc/rockfish/env_manifest.json" <<'PY'
import json
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone

out = {
    "built_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    "python": sys.version.split()[0],
    "platform": platform.platform(),
    "packages": {},
    "binaries": {},
}
for module in ("numpy", "scipy", "pandas", "sklearn", "rdkit", "vina", "meeko", "gemmi", "openmm"):
    try:
        mod = __import__(module)
    except Exception as exc:  # noqa: BLE001
        out["packages"][module] = f"absent ({type(exc).__name__})"
    else:
        out["packages"][module] = getattr(mod, "__version__", "unknown")

# The unpinned three, plus the rest, by whatever each one answers to.
for binary, args in (
    ("fpocket", ["--help"]),
    ("apbs", ["--version"]),
    ("mmseqs", ["version"]),
    ("foldseek", ["version"]),
    ("mafft", ["--version"]),
    ("pdb2pqr", ["--version"]),
):
    path = shutil.which(binary)
    if path is None:
        out["binaries"][binary] = "absent"
        continue
    try:
        result = subprocess.run([path] + args, capture_output=True, text=True, timeout=30, check=False)
        text = (result.stdout or result.stderr).strip().splitlines()
        out["binaries"][binary] = {"path": path, "version": text[0][:120] if text else "unknown"}
    except (OSError, subprocess.SubprocessError) as exc:
        out["binaries"][binary] = {"path": path, "version": f"unreadable ({type(exc).__name__})"}

with open(sys.argv[1], "w") as handle:
    json.dump(out, handle, indent=2, sort_keys=True)
    handle.write("\n")
print(f"wrote {sys.argv[1]}")
PY

echo
echo "== done. Next:"
echo "   python scripts/rockfish.py preflight --site hpc/sites/rockfish.yml --out preflight.json"
