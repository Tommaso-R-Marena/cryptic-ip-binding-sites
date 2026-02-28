# Integration test suite

This directory validates cross-module workflows rather than isolated units.

## Test files

- `test_full_pipeline.py`: ADAR2 download + ML scoring + batch output + mocked MD validation + comparative figure generation.
- `test_ml_integration.py`: ML training is wired into batch-style outputs and cache persistence.
- `test_comparative_workflow.py`: comparative outputs feed statistical validation and publication figure generation.

## Runtime expectations

- Local developer runs: ~30-90 seconds total.
- CI runs with optional dependencies missing: tests requiring scikit-learn are skipped.
- No GPU/OpenMM runtime is required (MD simulation is mocked in integration tests).

## Test data strategy

- Attempts to download the real ADAR2 AlphaFold structure once per test session.
- If network is unavailable, a tiny PDB stub is used so pipeline orchestration still validates.
- Remaining inputs are deterministic synthetic fixtures to keep CI reproducible.

## Common failures and troubleshooting

1. **`ModuleNotFoundError: sklearn`**
   - Install `scikit-learn` and `joblib`, or run only non-ML tests.
2. **Figure generation assertions fail**
   - Confirm matplotlib/seaborn are installed and writable temp directories are available.
3. **Unexpected empty comparative tables**
   - Check fixture hit labels still include both positive and negative outcomes.
4. **Network instability for ADAR2 fetch**
   - This is non-fatal in tests; fallback stub should be used automatically.
