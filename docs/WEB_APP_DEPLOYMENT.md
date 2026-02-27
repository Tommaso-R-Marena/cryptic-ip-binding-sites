# Streamlit Web App Deployment

This repository now includes a Streamlit interface (`streamlit_app.py`) for single-protein cryptic IP-binding site analysis.

## Run locally

```bash
pip install -r requirements.txt
streamlit run streamlit_app.py
```

## Streamlit Community Cloud

1. Push repository to GitHub.
2. In Streamlit Cloud, create a new app and select:
   - **Main file path**: `streamlit_app.py`
   - **Python version**: 3.11+
3. Add required system tools in your deployment image if you need full electrostatic pipeline support (`fpocket`, `pdb2pqr`, `apbs`).
4. Deploy.

## Heroku

The repository includes a `Procfile`:

```text
web: streamlit run streamlit_app.py --server.port=$PORT --server.address=0.0.0.0
```

Deploy with Heroku Container Registry or standard Python buildpack:

```bash
heroku create <app-name>
git push heroku <branch>:main
```

## Docker

Build and run with Docker:

```bash
docker build -t cryptic-ip-web .
docker run --rm -p 8501:8501 cryptic-ip-web
```

## Basic rate limiting

The Streamlit app includes per-session rate limiting (30-second cooldown between analysis runs), implemented in `streamlit_app.py`.

## Notes

- Structure inputs supported:
  - UniProt ID (AlphaFold fetch)
  - PDB ID (RCSB fetch)
  - custom uploaded structure file
- Results include:
  - 3D viewer with highlighted candidate pockets
  - electrostatic summary
  - CSV/JSON export
  - PyMOL `.pml` session script export
