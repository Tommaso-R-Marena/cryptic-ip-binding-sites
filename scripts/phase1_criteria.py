#!/usr/bin/env python3
"""Measure the project plan's Phase 1 criteria on the full control panel.

For each control: pick the inositol phosphate copy with the most protein
contacts, reduce the entry to that chain with every non-polymer atom removed,
detect and score pockets as the proteome screen does, and report where the
true site ranks. Then fetch the control's AlphaFold model, superpose the
binding region to measure AlphaFold-versus-crystal RMSD, carry the ligand into
the model's frame, and repeat the site measurements on the model - the input
the screen actually uses.

Criteria that cannot be measured are reported as unmeasured, never as met.
The report is evidence, not a gate: it exits zero whichever way the criteria
fall, so that the result is recorded rather than tuned until it passes.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cryptic_ip.utils.json_io import write_json_strict  # noqa: E402
from cryptic_ip.validation.phase1_criteria import (  # noqa: E402
    CONTROLS,
    binding_region_superposition,
    critical_test,
    evaluate_site,
    fetch_alphafold_model,
    judge,
    select_ligand_copy,
    site_residue_measurements,
    uniprot_for_entry,
    write_apo_chain,
    write_holo_site,
)

LOGGER = logging.getLogger("phase1_criteria")


def _download_pdb(pdb_id: str, directory: Path) -> Optional[Path]:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{pdb_id}.pdb"
    if path.exists() and path.stat().st_size > 0:
        return path
    try:
        with urllib.request.urlopen(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=120) as r:
            path.write_bytes(r.read())
        return path
    except Exception as exc:
        LOGGER.warning("could not download %s: %s", pdb_id, exc)
        return None


def measure(control, structures: Path, work: Path, use_apbs: bool) -> Dict[str, Any]:
    result: Dict[str, Any] = {"control": control.__dict__}
    path = _download_pdb(control.pdb_id, structures)
    if path is None:
        result["error"] = "structure unavailable"
        return result
    ligand = select_ligand_copy(path)
    if ligand is None:
        result["error"] = "no inositol phosphate in this entry"
        return result
    chain = ligand["chain"]
    result["ligand"] = {k: v for k, v in ligand.items() if k != "coords"}

    apo = write_apo_chain(path, chain, work / f"{control.pdb_id}_{chain}_apo.pdb")
    holo = write_holo_site(path, chain, ligand["residue_key"], work / f"{control.pdb_id}_{chain}_site.pdb")
    crystal: Dict[str, Any] = {"chain": chain}
    crystal["site"] = evaluate_site(apo, ligand["coords"], use_apbs=use_apbs)
    crystal["residues"] = site_residue_measurements(holo, apo, ligand["coords"], chain)
    result["crystal"] = crystal

    uniprot = uniprot_for_entry(control.pdb_id, chain)
    result["uniprot_id"] = uniprot
    model_result: Optional[Dict[str, Any]] = None
    if uniprot:
        model_path = fetch_alphafold_model(uniprot, structures / "alphafold")
        if model_path is not None:
            superposition = binding_region_superposition(path, chain, model_path, ligand["coords"])
            model_result = {"model": model_path.name, "superposition": {
                k: v for k, v in superposition.items() if k != "ligand_in_model_frame"
            }}
            if superposition.get("ok"):
                model_result["site"] = evaluate_site(
                    model_path, superposition["ligand_in_model_frame"], use_apbs=use_apbs
                )
    result["alphafold"] = model_result
    result["criteria"] = judge(control, crystal, model_result)
    return result


def _fmt(value: Any, digits: int = 2) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, float):
        return f"{value:.{digits}f}" if np.isfinite(value) else "-"
    return str(value)


def digest(results: List[Dict[str, Any]], tests: List[Dict[str, Any]]) -> str:
    lines = ["## Phase 1 criteria (project plan, sections 5-6)\n"]
    lines.append(
        "Each control reduced to the ligand's chain with all non-polymer atoms removed, "
        "scored as the proteome screen scores. `rank` is the true site's rank among all "
        "pockets by composite score.\n"
    )
    lines.append(
        "| control | pdb | role | ligand | pockets | site rank | site score | top score | "
        "APBS kT/e | Coulomb kT/e | basic <=5 A | coord. SASA holo / apo (basic) | "
        "AF model | AF RMSD (n) | AF site rank | AF site score |"
    )
    lines.append("|" + " --- |" * 16)
    for r in results:
        c = r["control"]
        if "error" in r:
            lines.append(f"| {c['name']} | {c['pdb_id']} | {c['role']} | {r['error']} |" + " |" * 12)
            continue
        s = r["crystal"]["site"]
        res = r["crystal"]["residues"]
        af = r.get("alphafold") or {}
        sup = af.get("superposition") or {}
        afs = af.get("site") or {}
        lines.append(
            f"| {c['name']} | {c['pdb_id']} | {c['role']} | {r['ligand']['comp_id']} "
            f"| {_fmt(s.get('n_pockets'))} | {_fmt(s.get('site_rank'))} | {_fmt(s.get('site_score'), 3)} "
            f"| {_fmt(s.get('top_score'), 3)} | {_fmt(s.get('apbs_potential_kT'))} "
            f"| {_fmt(s.get('coulomb_potential_kT'))} | {_fmt(res.get('n_basic_within_5A'))} "
            f"| {_fmt(res.get('basic_coordinating_sasa_holo_mean'), 1)} / "
            f"{_fmt(res.get('basic_coordinating_sasa_apo_mean'), 1)} "
            f"| {r.get('uniprot_id') or '-'} "
            f"| {_fmt(sup.get('rmsd'))} ({_fmt(sup.get('n_region_residues'))}) "
            f"| {_fmt(afs.get('site_rank'))}/{_fmt(afs.get('n_pockets'))} | {_fmt(afs.get('site_score'), 3)} |"
        )
    lines.append("\n### What the site pocket measures (apo), crystal and AlphaFold\n")
    lines.append(
        "| control | source | volume | depth | enclosure | lining SASA | basic (pocket) | "
        "Coulomb kT/e | APBS kT/e | pLDDT | score |"
    )
    lines.append("|" + " --- |" * 11)
    for r in results:
        for source in ("crystal", "alphafold"):
            site = ((r.get(source) or {}).get("site")) or {}
            if not site.get("ok"):
                continue
            lines.append(
                f"| {r['control']['name']} | {source} | {_fmt(site.get('pocket_volume'), 0)} "
                f"| {_fmt(site.get('burial_depth'))} | {_fmt(site.get('enclosure'))} "
                f"| {_fmt(site.get('sasa_mean'), 1)} | {_fmt(site.get('n_basic_residues_pocket'))} "
                f"| {_fmt(site.get('coulomb_potential_kT'))} | {_fmt(site.get('apbs_potential_kT'))} "
                f"| {_fmt(site.get('plddt_mean'), 1)} | {_fmt(site.get('site_score'), 3)} |"
            )
    lines.append("\n### Criteria\n")
    lines.append("| control | criterion | target | value | met |")
    lines.append("| --- | --- | --- | --- | --- |")
    for r in results:
        for name, item in (r.get("criteria") or {}).items():
            lines.append(
                f"| {r['control']['name']} | {name} | {item['target']} | {_fmt(item['value'], 3)} "
                f"| {_fmt(item['passes']) if item['passes'] is not None else 'unmeasured'} |"
            )
    lines.append("\n### The critical test\n")
    for test in tests:
        lines.append(f"**{test['source']}**\n")
        for key, value in test.items():
            if key != "source":
                lines.append(f"- {key}: {_fmt(value, 3)}")
        lines.append("")
    return "\n".join(lines) + "\n"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structures-dir", default="data/validation")
    parser.add_argument("--output-dir", default="results/phase1")
    parser.add_argument("--no-apbs", action="store_true")
    parser.add_argument("--only", nargs="*", help="Control names to run")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")

    out = Path(args.output_dir)
    work = out / "structures"
    results = []
    for control in CONTROLS:
        if args.only and control.name not in args.only:
            continue
        LOGGER.info("measuring %s (%s)", control.name, control.pdb_id)
        try:
            results.append(measure(control, Path(args.structures_dir), work, not args.no_apbs))
        except Exception as exc:  # one control failing must not hide the rest
            LOGGER.exception("control %s failed", control.name)
            results.append({"control": control.__dict__, "error": f"{type(exc).__name__}: {exc}"})

    tests = [critical_test(results, "crystal")]
    model_results = [
        {"control": r["control"], "alphafold": r.get("alphafold") or {}} for r in results
    ]
    tests.append(critical_test(model_results, "alphafold"))
    write_json_strict(out / "phase1_criteria.json", {"controls": results, "critical_test": tests}, indent=2)
    text = digest(results, tests)
    (out / "PHASE1.md").write_text(text, encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
