# AGENTS.md

## Cursor Cloud specific instructions

This is a Python (>=3.8; CI targets 3.10/3.11, VM uses 3.12) bioinformatics pipeline for
detecting buried inositol-phosphate binding sites. It ships two products:

- **`cryptic-ip` CLI** (entry point in `cryptic_ip/cli.py`) — e.g. `cryptic-ip analyze <pdb>`,
  `cryptic-ip check-dependencies`, `cryptic-ip validate`.
- **Streamlit web app** (`streamlit_app.py`) — single-protein analysis UI.

Standard commands are documented in `README.md`, `TESTING.md`, `CONTRIBUTING.md`, and
`.github/workflows/ci.yml`; prefer those. Notes below are the non-obvious, environment-specific bits.

### Python environment
- Dependencies are installed into a virtualenv at `/workspace/.venv` (the system Python is
  externally managed). Activate with `source /workspace/.venv/bin/activate` before running
  anything, or call binaries directly as `/workspace/.venv/bin/<cmd>`.
- The startup update script refreshes this venv (`pip install -r requirements.txt`,
  `pip install -e .`, plus dev tools). It does NOT install the external tools below.

### External bioinformatics tools (NOT pip-installable)
- `fpocket` is **required** for pocket detection (`cryptic-ip analyze/screen` and the Streamlit
  pipeline raise `RuntimeError` without it). `apbs` + `pdb2pqr` (electrostatics) are optional;
  electrostatics degrades gracefully to a warning when absent.
- These are installed once via Miniforge at `~/miniforge3` (from conda-forge) and persist in the
  VM snapshot. `~/miniforge3/bin` is added to `PATH` via `~/.bashrc`. If `fpocket` is not found,
  run `export PATH="$HOME/miniforge3/bin:$PATH"` (non-interactive shells may not source `.bashrc`).
- `freesasa` CLI is intentionally absent; SASA is computed with BioPython's `ShrakeRupley`, so
  `check-dependencies` listing `freesasa` as an optional missing tool is expected and harmless.

### Lint / test / build / run
- Lint (CI gate): `flake8 cryptic_ip --count --select=E9,F63,F7,F82 --show-source --statistics`.
  `black`/`isort`/`mypy` run non-blocking (`continue-on-error`) in CI.
- Tests: `pytest tests/ -m "not requires_network and not slow" --ignore=tests/integration -n auto`.
  Tests tagged `requires_network` hit RCSB/AlphaFold/UniProt; some tests skip when `fpocket` or a
  cached ADAR2 structure is unavailable. Benchmarks: `pytest tests/benchmarks`.
- Run CLI (needs a PDB file): download one first, e.g.
  `curl -fsSL -o data/1ZY7.pdb https://files.rcsb.org/download/1ZY7.pdb` then
  `cryptic-ip analyze data/1ZY7.pdb --output results/out.csv` (add `--use-ml-model` for the
  pre-trained classifier in `models/`).
- Run web app: `streamlit run streamlit_app.py` (defaults to PDB `1ZY7`). Requires `fpocket` on PATH.

### Known non-blocking warnings
- The bundled ML model (`models/cryptic_ip_classifier_v1.pkl`) was pickled with scikit-learn
  1.9.0 while `requirements.txt` pins 1.7.2, so `--use-ml-model` prints `InconsistentVersionWarning`.
  It still loads and runs. Do not "fix" this by re-pinning unless asked.
