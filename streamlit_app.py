"""Streamlit interface for cryptic IP-binding site prediction."""

from __future__ import annotations

import json
import tempfile
import time
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from cryptic_ip.analysis.analyzer import ProteinAnalyzer
from cryptic_ip.analysis.scorer import PocketScorer
from cryptic_ip.database.alphafold_client import AlphaFoldClient
from cryptic_ip.database.pdb_client import PDBClient
from cryptic_ip.errors import UserFacingError
from cryptic_ip.utils.input_validation import validate_structure_file


APP_TITLE = "Cryptic IP-Binding Site Explorer"
RATE_LIMIT_SECONDS = 30


@st.cache_resource
def get_clients() -> Tuple[AlphaFoldClient, PDBClient]:
    return AlphaFoldClient(), PDBClient()


def _set_progress(progress_bar, status, value: int, message: str) -> None:
    progress_bar.progress(value)
    status.info(message)


def _check_rate_limit() -> bool:
    last_run = st.session_state.get("last_run_ts")
    if last_run is None:
        return True

    elapsed = time.time() - last_run
    if elapsed >= RATE_LIMIT_SECONDS:
        return True

    st.warning(
        f"Rate limit active. Please wait {RATE_LIMIT_SECONDS - int(elapsed)}s "
        "before starting another job."
    )
    return False


def _prepare_structure(
    source_type: str,
    uniprot_id: str,
    pdb_id: str,
    uploaded_file,
) -> Tuple[Path, Dict[str, str]]:
    alphafold_client, pdb_client = get_clients()

    if source_type == "UniProt ID (AlphaFold)":
        path = alphafold_client.fetch_structure(uniprot_id.strip().upper())
        return path, {"source": "alphafold", "id": uniprot_id.strip().upper()}

    if source_type == "PDB ID (RCSB)":
        path = pdb_client.fetch_structure(pdb_id.strip().upper())
        return path, {"source": "pdb", "id": pdb_id.strip().upper()}

    if uploaded_file is None:
        raise ValueError("Please upload a structure file.")

    suffix = Path(uploaded_file.name).suffix or ".pdb"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as handle:
        handle.write(uploaded_file.getvalue())
        tmp_path = Path(handle.name)

    validate_structure_file(tmp_path)
    return tmp_path, {"source": "upload", "id": uploaded_file.name}


def _run_pipeline(pdb_path: Path, progress_bar, status) -> Tuple[pd.DataFrame, Optional[float], ProteinAnalyzer]:
    analyzer = ProteinAnalyzer(str(pdb_path))

    _set_progress(progress_bar, status, 15, "Validating structure and initializing analysis…")
    time.sleep(0.2)

    _set_progress(progress_bar, status, 40, "Running fpocket detection…")
    analyzer.detect_pockets()

    _set_progress(progress_bar, status, 65, "Calculating solvent accessibility metrics…")
    analyzer.calculate_sasa()

    _set_progress(progress_bar, status, 82, "Estimating electrostatic properties…")
    potential = analyzer.calculate_electrostatics()

    _set_progress(progress_bar, status, 95, "Scoring cryptic IP-binding candidates…")
    results = analyzer.score_all_pockets()
    scorer = PocketScorer()
    results["classification"] = results["composite_score"].apply(scorer.classify_site)

    _set_progress(progress_bar, status, 100, "Analysis complete.")
    return results, potential, analyzer


def _build_pymol_script(df: pd.DataFrame, pdb_name: str) -> str:
    top_rows = df.head(5)
    lines = [
        f"load {pdb_name}, protein",
        "hide everything",
        "show cartoon, protein",
        "color gray70, protein",
    ]
    for i, row in enumerate(top_rows.itertuples(index=False), start=1):
        color = "tv_red" if i == 1 else "orange"
        x, y, z = row.center
        lines.extend(
            [
                f"pseudoatom center_{row.pocket_id}, pos=[{x:.3f},{y:.3f},{z:.3f}]",
                f"select pocket_{row.pocket_id}, byres (protein within 5 of center_{row.pocket_id})",
                f"show sticks, pocket_{row.pocket_id}",
                f"color {color}, pocket_{row.pocket_id}",
            ]
        )
    lines.append("set ray_opaque_background, off")
    lines.append("bg_color white")
    return "\n".join(lines)


def _render_3d_viewer(pdb_text: str, candidates: pd.DataFrame) -> None:
    view = py3Dmol.view(width=920, height=560)
    view.addModel(pdb_text, "pdb")
    view.setStyle({"cartoon": {"color": "lightgray"}})

    for rank, row in enumerate(candidates.head(5).itertuples(index=False), start=1):
        color = "red" if rank == 1 else "orange"
        view.addSphere(
            {
                "center": {
                    "x": float(row.center[0]),
                    "y": float(row.center[1]),
                    "z": float(row.center[2]),
                },
                "radius": 2.3,
                "color": color,
                "alpha": 0.85,
            }
        )

    view.addSurface(
        py3Dmol.VDW,
        {"opacity": 0.65, "colorscheme": "RWB"},
        {"hetflag": False},
    )
    view.zoomTo()
    components.html(view._make_html(), height=580, width=940)


def main() -> None:
    st.set_page_config(page_title=APP_TITLE, layout="wide")
    st.title(APP_TITLE)
    st.caption("Single-protein cryptic inositol phosphate (IP) binding site prediction in ~2-3 minutes.")

    with st.expander("System health", expanded=False):
        health = {
            "status": "ok",
            "timestamp": datetime.utcnow().isoformat(),
            "last_run_ts": st.session_state.get("last_run_ts"),
            "python_tmp_dir": tempfile.gettempdir(),
        }
        st.json(health)

    with st.sidebar:
        st.header("Input")
        source_type = st.radio(
            "Structure source",
            ["UniProt ID (AlphaFold)", "PDB ID (RCSB)", "Upload structure"],
        )

        uniprot_id = st.text_input("UniProt ID", value="P78563", disabled=source_type != "UniProt ID (AlphaFold)")
        pdb_id = st.text_input("PDB ID", value="1ZY7", disabled=source_type != "PDB ID (RCSB)")
        uploaded_file = st.file_uploader(
            "Upload structure (.pdb/.cif)",
            type=["pdb", "cif", "ent"],
            disabled=source_type != "Upload structure",
        )

        st.markdown("---")
        st.subheader("Example gallery")
        st.markdown(
            "- **ADAR2** (UniProt: `P78563`)\n"
            "- **Pds5B** (UniProt: `Q9NTI5`)\n"
            "- **User validations**: share your validated IDs via contact form below."
        )

        run_analysis = st.button("Run cryptic site analysis", type="primary")

    docs_tab, results_tab, contact_tab = st.tabs(["Tutorial & docs", "Results", "Collaborate"])

    with docs_tab:
        st.subheader("Tutorial")
        st.markdown(
            "1. Choose UniProt, PDB, or upload a structure.\n"
            "2. Click **Run cryptic site analysis**.\n"
            "3. Inspect ranked pockets and interactive structure view.\n"
            "4. Export JSON/CSV and PyMOL script for downstream validation."
        )
        st.markdown("Documentation: `docs/TUTORIAL.md`, `docs/API.md`, and `docs/INSTALLATION.md` in the repository.")
        st.subheader("Citation")
        st.markdown(
            "If this interface helps your work, please cite the project and associated manuscript metadata from `CITATION.cff`."
        )

    with contact_tab:
        st.subheader("Contact for collaborations")
        with st.form("contact_form"):
            name = st.text_input("Name")
            email = st.text_input("Email")
            org = st.text_input("Institution")
            message = st.text_area("Project details / validation data")
            submitted = st.form_submit_button("Submit")

        if submitted:
            if not email or not message:
                st.error("Email and project details are required.")
            else:
                payload = {
                    "timestamp": datetime.utcnow().isoformat(),
                    "name": name,
                    "email": email,
                    "institution": org,
                    "message": message,
                }
                st.success("Thanks! Please copy/send this payload to the maintainers.")
                st.code(json.dumps(payload, indent=2), language="json")

    if run_analysis:
        if not _check_rate_limit():
            return

        with results_tab:
            progress_bar = st.progress(0)
            status = st.empty()

            try:
                structure_path, source_meta = _prepare_structure(source_type, uniprot_id, pdb_id, uploaded_file)
                results_df, potential, analyzer = _run_pipeline(structure_path, progress_bar, status)
                st.session_state["last_run_ts"] = time.time()
                st.session_state["results_df"] = results_df
                st.session_state["structure_path"] = str(structure_path)
                st.session_state["source_meta"] = source_meta
                st.session_state["electrostatic_potential"] = potential
                st.session_state["work_dir"] = str(analyzer.work_dir)
            except UserFacingError as exc:
                st.error(str(exc))
                return
            except Exception as exc:
                st.error(
                    "Analysis failed while processing your structure. "
                    "Please verify the file format/identifier and retry. "
                    f"Details: {exc}. See docs/TROUBLESHOOTING.md"
                )
                return

    with results_tab:
        results_df = st.session_state.get("results_df")
        structure_path = st.session_state.get("structure_path")

        if results_df is None:
            st.info("Run an analysis to see candidate cryptic IP-binding sites.")
            return

        st.success("Prediction complete.")
        st.metric("Candidate pockets", value=len(results_df))

        top_hits = results_df.head(10).copy()
        st.subheader("Score breakdown table")
        st.dataframe(
            top_hits[["pocket_id", "composite_score", "classification", "volume", "depth", "sasa", "basic_residues"]],
            use_container_width=True,
        )

        if structure_path and Path(structure_path).exists():
            st.subheader("Interactive 3D structure and site highlights")
            pdb_text = Path(structure_path).read_text(errors="ignore")
            _render_3d_viewer(pdb_text, top_hits)

        st.subheader("Electrostatic surface summary")
        potential = st.session_state.get("electrostatic_potential")
        if potential is None:
            st.warning("APBS electrostatics unavailable in this environment. Viewer still shows an electrostatic-style surface coloring.")
        else:
            st.write(f"Estimated electrostatic potential: **{potential:.2f} kT/e**")

        source_meta = st.session_state.get("source_meta", {})
        json_payload = {
            "source": source_meta,
            "generated_at": datetime.utcnow().isoformat(),
            "results": top_hits.to_dict(orient="records"),
            "electrostatic_potential": potential,
        }

        st.subheader("Exports")
        col1, col2, col3 = st.columns(3)
        with col1:
            st.download_button(
                "Download CSV",
                data=top_hits.to_csv(index=False),
                file_name="cryptic_ip_sites.csv",
                mime="text/csv",
            )
        with col2:
            st.download_button(
                "Download JSON",
                data=json.dumps(json_payload, indent=2),
                file_name="cryptic_ip_sites.json",
                mime="application/json",
            )
        with col3:
            pymol_script = _build_pymol_script(top_hits, Path(structure_path).name)
            st.download_button(
                "Download PyMOL session script",
                data=StringIO(pymol_script).getvalue(),
                file_name="cryptic_ip_session.pml",
                mime="text/plain",
            )


if __name__ == "__main__":
    main()
