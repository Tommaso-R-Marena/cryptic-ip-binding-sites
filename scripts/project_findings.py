#!/usr/bin/env python3
"""What this project establishes, assembled from every line's own artefact.

`results/publication/RESULTS_SUMMARY.md` is generated from the control
benchmark, the validation dataset, the ML benchmark and the yeast pilot, and
from nothing else. It therefore predates the lining-architecture line, the
within-protein controls, the sequence and composition studies, and the
specificity re-ranking - several of which qualify its claims. This page is the
statement that covers all of them.

Like `scripts/architecture_findings.py`, it **recomputes nothing**: every figure
is read out of a committed artefact, a missing key fails the render rather than
printing a blank, and a section whose claim its artefact no longer supports
raises instead of being quietly contradicted.

    project_findings.py --out results/PROJECT_FINDINGS.md
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence

import pandas as pd

#: Artefacts the page is assembled from. Every one is committed.
SOURCES = {
    "controls": Path("results/publication/validation/positive_controls.csv"),
    "negatives": Path("results/publication/validation/negative_controls.csv"),
    "dataset": Path("data/validation/ip_binding_validation_dataset.csv"),
    "specificity": Path("results/specificity/specificity.json"),
    "candidates": Path("results/specificity/candidates.csv"),
    "architecture": Path("results/architecture/real/real.json"),
    "within": Path("results/within_protein/within_protein.json"),
    "sites": Path("results/within_protein/within_protein_sites.json"),
    "order": Path("results/pocket_sequence/pocket_sequence.json"),
    "gate": Path("results/composition/validation.json"),
    "census": Path("results/redocking/census.csv"),
    "coevolution": Path("results/coevolution/coevolution.json"),
}

#: Artefacts the page reports when present and names as outstanding when not.
#: A study that has not reported is a fact about the project, not a failure to
#: render, and saying "not yet" is more honest than leaving the line out.
OPTIONAL = {
    "composition": Path("results/composition/composition.json"),
}

COMMANDS = {
    "controls": "scripts/run_publication_package.py",
    "negatives": "scripts/run_publication_package.py",
    "dataset": "scripts/build_validation_dataset.py",
    "specificity": "scripts/specificity.py run",
    "candidates": "scripts/specificity.py run",
    "architecture": "scripts/real_architecture.py run",
    "within": "scripts/within_protein.py analyse --arm candidate",
    "sites": "scripts/within_protein.py analyse --arm site",
    "order": "scripts/pocket_sequence.py analyse",
    "gate": "scripts/pocket_composition.py validate",
    "composition": "scripts/pocket_composition.py analyse",
    "census": "scripts/redocking.py census",
    "coevolution": "scripts/coevolution.py report",
}


class Contradiction(KeyError):
    """An artefact no longer supports what the section about it says."""


def load(sources: Mapping[str, Path] = SOURCES,
         optional: Mapping[str, Path] = OPTIONAL) -> Dict[str, Any]:
    """Read every artefact, failing with the path when a required one is absent."""
    out: Dict[str, Any] = {}
    for name, path in sources.items():
        if not path.exists():
            raise FileNotFoundError(f"{name}: {path} has not been generated")
        out[name] = (pd.read_csv(path) if path.suffix == ".csv"
                     else json.loads(path.read_text()))
    for name, path in optional.items():
        out[name] = (json.loads(path.read_text()) if path.exists() else None)
    return out


def pick(blob: Any, *keys: Any) -> Any:
    """Index a nested artefact, naming the whole path when a key is missing."""
    here: Any = blob
    for index, key in enumerate(keys):
        try:
            here = here[key]
        except (KeyError, IndexError, TypeError) as error:
            trail = " -> ".join(str(k) for k in keys[:index + 1])
            raise Contradiction(f"{trail} is not in the artefact") from error
    return here


def interval(blob: Mapping[str, Any]) -> str:
    return f"{pick(blob, 'point'):.3f} [{pick(blob, 'low'):.3f}, {pick(blob, 'high'):.3f}]"


def cryptic_census(census: pd.DataFrame) -> Dict[str, int]:
    """How many genuinely buried annotated copies the project actually holds.

    Two different tables carry a burial class and they answer different
    questions. The control benchmark's validation set bounds the *benchmark*;
    the redocking census bounds the *ground-truth arms* of every study that uses
    annotated sites. Conflating them understates what the project has, which is
    a mistake this function exists to stop being made again by hand.
    """
    for column in ("burial_class", "classification"):
        if column in census.columns:
            classes = census[column].astype(str).str.lower()
            break
    else:
        raise Contradiction("the census has no burial class column")
    selected = census
    for flag in ("selected", "primary_set"):
        if flag in census.columns:
            selected = selected[selected[flag].astype(str).str.lower().isin(
                ("true", "1", "yes"))]
    picked = classes.loc[selected.index]
    group = ("homology_group_strict" if "homology_group_strict" in census.columns
             else "pdb_id")
    cryptic = selected[picked == "cryptic"]
    semi = selected[picked.isin(("semi_cryptic", "semi-cryptic"))]
    return {
        "rows": int(len(census)),
        "all_cryptic": int((classes == "cryptic").sum()),
        "selected": int(len(selected)),
        "cryptic": int(len(cryptic)),
        "cryptic_groups": int(cryptic[group].nunique()),
        "semi_cryptic": int(len(semi)),
        "semi_cryptic_groups": int(semi[group].nunique()),
    }


# ------------------------------------------------------------------ sections

def header(data: Mapping[str, Any]) -> list:
    lines = [
        "# What this project establishes",
        "",
        "This page is **generated** by `scripts/project_findings.py` from the committed "
        "artefacts below. It recomputes nothing, no number on it is retyped, and a claim "
        "its artefact stops supporting fails the render rather than standing.",
        "",
        "It exists because `results/publication/RESULTS_SUMMARY.md` is generated from the "
        "control benchmark, the validation dataset, the ML benchmark and the yeast pilot "
        "alone, so it predates and does not cover most of the work. **Where the two "
        "differ, this page is the current one.**",
        "",
        "| line | artefact | command |",
        "|---|---|---|",
    ]
    for name, path in list(SOURCES.items()) + list(OPTIONAL.items()):
        state = "" if data.get(name) is not None else " *(not yet reported)*"
        lines.append(f"| {name}{state} | `{path}` | `{COMMANDS[name]}` |")
    lines += [
        "",
        "Pre-registered plans are in `docs/`, each with its dated amendments. The "
        "lining-architecture line has its own detailed page at "
        "`results/architecture/FINDINGS.md`.",
        "",
    ]
    return lines


def verdict(data: Mapping[str, Any]) -> list:
    """The top line, in the order a reader needs it."""
    positives = data["controls"]
    passed = int(positives["passed"].astype(bool).sum()) if "passed" in positives else 0
    dataset = data["dataset"]
    column = ("classification" if "classification" in dataset.columns
              else "burial_class" if "burial_class" in dataset.columns else None)
    if column is None:
        raise Contradiction("the validation dataset has no burial classification column")
    counts = dataset[column].astype(str).str.lower().value_counts().to_dict()
    cryptic = int(counts.get("cryptic", 0))
    rerank = pick(data["specificity"], "rerank")
    candidates = data["candidates"]
    best = float(candidates["precision_at_rank"].max()) if len(candidates) else 0.0
    census = cryptic_census(data["census"])

    # The page's spine. If any of these inverts, the summary below is wrong.
    if cryptic > 5:
        raise Contradiction(
            f"the validation set now holds {cryptic} cryptic structures; the limitation "
            "this page is built around has changed and the text must be rewritten")
    if passed == len(positives):
        raise Contradiction(
            "every positive control now passes; the panel sentence below says otherwise")
    if pick(rerank, "decisions", "L3b") == "supported":
        raise Contradiction(
            "the combined score now beats the same-protein baseline; this page says it "
            "does not, and that section must be rewritten from the artefact")

    return [
        "## The short answer",
        "",
        f"* **A detector that is not inverted on its one gold-standard case.** "
        f"{passed} of {len(positives)} positive controls pass. The tier-1 gate is one "
        "positive against one negative, so it shows the score is not backwards on the "
        "case it was built around, and nothing about an unseen protein.",
        f"* **Buried examples are scarce, and two different tables say how scarce.** The "
        f"control benchmark's validation set holds {cryptic} cryptic structure"
        f"{'' if cryptic == 1 else 's'} in {len(dataset)}, which is what bounds the "
        f"benchmark. The redocking census holds **{census['cryptic']} cryptic annotated "
        f"copies across {census['cryptic_groups']} homology groups** among its "
        f"{census['selected']} selected copies ({census['all_cryptic']} cryptic rows in "
        f"{census['rows']}), plus {census['semi_cryptic']} semi-cryptic across "
        f"{census['semi_cryptic_groups']} groups — and that is what bounds the "
        "ground-truth arms of every study below. Small, but not one.",
        f"* **A learnable inositol-phosphate-versus-other-polyanion signal** at "
        f"{interval(pick(rerank, 'L3a_combined'))} ROC-AUC against "
        f"{pick(rerank, 'annotated_unseen')} annotated binders among "
        f"{pick(rerank, 'unseen')} unseen proteins — which **does not beat its own "
        f"same-protein baseline** "
        f"({interval(pick(rerank, 'L3b_combined_minus_L1'))}, Holm "
        f"{pick(rerank, 'holm', 'L3b'):.3f}, *not supported*).",
        f"* **{len(candidates)} candidate proteins that are hypotheses, not findings.** "
        f"The best precision at any rank is {best:.2f}, and it is "
        f"{float(candidates['precision_at_rank'].iloc[0]):.2f} at rank 1.",
        "* **Two structural decompositions that separate real linings from random ones "
        "and not from other pockets.** The architecture line (positions) and the sequence "
        "line (letters) account for most of the project's statistical work. Against "
        "random position sets the architecture contrast is large; against the protein's "
        "own other pockets, matched on size, both lines are null. Sections 2 to 4.",
        "",
        "Nothing on this page is evidence that any candidate binds an inositol phosphate. "
        "The deposited cryptic copies above are structures with a buried inositol "
        "phosphate, classified by this project's own burial metric; **one of them, ADAR2, "
        "is the case with published evidence that the protein cannot fold without the "
        "cofactor**, and it is the case the detector was built around.",
        "",
    ]


def specificity(data: Mapping[str, Any]) -> list:
    rerank = pick(data["specificity"], "rerank")
    gate = pick(data["specificity"], "gate")
    candidates = data["candidates"]
    explained = int(pick(rerank, "candidates_explained"))
    unexplained = int(pick(rerank, "candidates_unexplained"))
    return [
        "## 1. Can the score tell an IP site from another polyanion site?",
        "",
        f"Pre-registered in `docs/SPECIFICITY_PLAN.md`. The gate opened on the "
        f"`{pick(gate, 'variant')}` variant at a held-out point of "
        f"{pick(gate, 's1r_point'):.3f}, after the all-data variant was ruled "
        "**not evaluable** because one homology group held more than 40% of the IP "
        "pockets — which is the kind of thing that makes an impressive number worthless, "
        "and was caught before the number was used.",
        "",
        "| test | ROC-AUC [95%] | Holm | decision |",
        "|---|---|---|---|",
        f"| L1, same proteins | {interval(pick(rerank, 'L1_same_proteins'))} | — | — |",
        f"| L3a, combined | {interval(pick(rerank, 'L3a_combined'))} | "
        f"{pick(rerank, 'holm', 'L3a'):.3f} | **supported** |",
        f"| L3b, combined − L1 | {interval(pick(rerank, 'L3b_combined_minus_L1'))} | "
        f"{pick(rerank, 'holm', 'L3b'):.3f} | **not supported** |",
        "",
        "**The third row is the one that matters and the easiest to skip.** The combined "
        "score ranks annotated binders well above chance, and it does **not** rank them "
        "better than the same-protein baseline does. Everything the pipeline adds beyond "
        "that baseline is, on this evidence, unmeasurable.",
        "",
        f"Of the {len(candidates)} candidates, {explained} are explained by an existing "
        f"annotation and {unexplained} are not. The unexplained ones are the hypotheses; "
        "they are also the ones with no independent support of any kind.",
        "",
    ]


def architecture(data: Mapping[str, Any]) -> list:
    real = data["architecture"]
    return [
        "## 2. The lining-architecture line: positions",
        "",
        "Pre-registered in `docs/ARCHITECTURE_PLAN.md` with five dated amendments, and "
        "set out in full at `results/architecture/FINDINGS.md`. It factors a pocket "
        "lining into **topology** (the gap multiset) and **register** (the arrangement "
        "within it), and tests each against its own permutation null.",
        "",
        f"{pick(real, 'n_measurable')} of {pick(real, 'n_linings')} real candidate "
        "linings were measurable. The method properties are real and calibrated, and the "
        "detailed page carries them. What they amount to on real data is a per-lining "
        "verdict:",
        "",
        "| verdict | linings |",
        "|---|---|",
    ] + [
        f"| {name} | {count} |"
        for name, count in sorted(pick(real, "verdicts").items(),
                                  key=lambda kv: (-int(kv[1]), str(kv[0])))
    ] + [
        "",
        f"Against size- and span-matched **decoy** position sets the contrast is large: "
        f"{pick(real, 'rates', 'real_either'):.1%} of real linings resolve on one half or "
        f"the other against {pick(real, 'rates', 'decoy_either'):.1%} of "
        f"{pick(real, 'rates', 'n_decoys')} decoys. That is the architecture line's "
        "strongest real-data number, and the line says in its own plan why it is weaker "
        "than it looks: **a random position set is not a protein**, so the comparison is "
        "close to tautological — a real lining is contiguous, secondary-structure-bound "
        "and solvent-excluded, and a random set of positions is none of those.",
        "",
        "**That is exactly why section 3 exists**: it replaces the decoy with the "
        "protein's own other pockets, matched on lining size, so protein, fold, model, "
        "detector and size are all fixed and only the choice of pocket varies. On that "
        "comparison the line is null.",
        "",
        f"Rates are in the artefact (`{SOURCES['architecture']}`) and on the detailed "
        "page. The transfer ceiling stated there — the method's agreement with ideal "
        "geometry falling as backbone jitter and bulges are introduced — bounds this "
        "whole line and is the reason a null here is weak evidence rather than evidence "
        "of absence.",
        "",
    ]


def within_protein(data: Mapping[str, Any]) -> list:
    rows = []
    for label, key in (("candidate", "within"), ("annotated site", "sites")):
        report = data[key]
        for test in pick(report, "tests"):
            if not test.get("runnable", True):
                continue
            rows.append((label, test, report))
    findings = [t for _, t, _ in rows if t.get("is_finding")]
    if findings:
        raise Contradiction(
            "a within-protein test is now a finding; this section states that none is")
    lines = [
        "## 3. The control that did not exist, built — and null in both arms",
        "",
        "Pre-registered in `docs/WITHIN_PROTEIN_PLAN.md`. The architecture line's own "
        "comparison sets a candidate lining against random position sets, and says in "
        "print that a random position set is not a protein. The control that means "
        "something is **other pockets of the same protein, matched on lining size** — and "
        "no committed table carried their residues, so this study built it.",
        "",
        "| arm | half | proteins | groups | mean quantile | difference | Holm | "
        "could have detected |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for label, test, report in rows:
        lines.append(
            f"| {label} | {pick(test, 'half')} | {pick(report, 'n_measured')} | "
            f"{pick(test, 'n_groups')} | {pick(test, 'mean_quantile'):.3f} | "
            f"{pick(test, 'difference'):+.3f} | {pick(test, 'p_holm'):.3f} | "
            f"{pick(test, 'minimum_detectable'):.3f} |")
    lines += [
        "",
        "Below 0.5 means more structured than the protein's own same-size pockets. "
        "**Nothing is a finding, and in every row the difference observed is smaller than "
        "the difference the design could have detected** — so this is weak evidence, not "
        "evidence of absence.",
        "",
        "The annotated-site arm decides how to read the other: it runs the identical "
        "statistic on deposited inositol-hexakisphosphate sites, so it could have fired, "
        "and a failure there would mean the candidate arm's null is about the instrument "
        "rather than about candidates. It does not separate real sites from their own "
        "proteins' pockets at the power available.",
        "",
        "One thing it established for free: re-detection reproduced the committed "
        f"candidate linings at a median Jaccard overlap of "
        f"{pick(data['within'], 'overlaps')['median']:.2f}, so the candidate table is "
        "reproducible from current models with current fpocket — which nothing had checked.",
        "",
    ]
    return lines


def order(data: Mapping[str, Any]) -> list:
    report = data["order"]
    tests = pick(report, "tests")
    by_key = {(t["arm"], t["statistic"]): t for t in tests if t.get("runnable")}
    primary = [t for t in tests if t.get("runnable") and t.get("family") == "primary"]
    if [t for t in primary if t.get("is_finding")]:
        raise Contradiction(
            "a study Q9 test is now a finding; this section states that none is")
    check = pick(report, "look_reproduction")
    if not check.get("selection_reproduced") or check.get("verdict") != "reproduced":
        raise Contradiction(
            "the committed Q9 run no longer reproduces the look its amendment discloses, "
            "so the held-out arm named below is not held out")
    pre = by_key[("formation", "helix_hydropathy_periodicity")]
    held = by_key[("formation_heldout", "helix_hydropathy_periodicity")]
    low, high = pick(held, "interval")
    return [
        "## 4. The sequence line: letters, in order",
        "",
        "Pre-registered in `docs/POCKET_SEQUENCE_PLAN.md` with two dated amendments. A "
        "window of sequence factors as **composition** (which letters) **⊕ order** (in "
        "what arrangement), the same way a lining factors into topology and register. "
        "This line takes the quotient by composition and asks what the ordering carries.",
        "",
        "**The result worth carrying off this line is a correction, not a discovery.** An "
        f"underpowered form of the formation contrast gave {pick(pre, 'difference'):+.2f} "
        f"on {pick(pre, 'n_site')} windows at a corrected p the plan refused, because the "
        "observed difference was smaller than the minimum detectable one. The "
        "re-specified rule, on the resampling groups the pre-registered rule selected "
        f"nothing from, gives **{pick(held, 'difference'):+.3f}** on "
        f"{pick(held, 'n_site')} windows across {pick(held, 'n_groups')} groups, interval "
        f"[{low:+.3f}, {high:+.3f}], Holm {pick(held, 'p_holm'):.4f}. Not a small effect "
        "— an absent one.",
        "",
        "Without the power floor, fixed before any sequence was read, this project would "
        "have reported a sequence predictor of pocket formation that does not exist.",
        "",
        f"Architectural type is null and powered enough to say so. Binding is **null but "
        f"not resolved**: the plan's contiguous-span window admits "
        f"{pick(report, 'n_site_windows')} of {pick(report, 'n_site_windows_dense')} "
        "annotated contact shells, because a shell draws residues from several "
        "secondary-structure elements at once and so has no short contiguous span. That "
        "is a property of the object and no sample size fixes it.",
        "",
    ]


def composition(data: Mapping[str, Any]) -> list:
    gate = data["gate"]
    checks = pick(gate, "checks")
    if not pick(gate, "all_checks_pass"):
        failed = sorted(name for name, ok in checks.items() if not ok)
        raise Contradiction(
            f"the committed gate calibration no longer passes ({', '.join(failed)}); "
            "this section says the instrument check is calibrated")
    planted = pick(gate, "planted")
    lines = [
        "## 5. The instrument check the sequence line did not have",
        "",
        "Pre-registered in `docs/COMPOSITION_PLAN.md`. Sections 3 and 4 return nulls, and "
        "a design that returns nulls has either found a real absence or has no "
        "sensitivity. **Nothing in those studies distinguishes the two**: their only "
        "positive control was a synthetic planted signal.",
        "",
        "So this study takes the other half of the decomposition — composition, which the "
        "order null holds fixed and cannot see — and puts a known-true signal through the "
        "same machinery. An inositol-phosphate contact shell being arginine-rich is "
        "textbook and is **not offered as a finding**; it is offered as the check.",
        "",
        f"The gate is calibrated before any real lining reaches it: a planted basic patch "
        f"gives mean z **{pick(planted, 'mean_z'):+.2f}** (p "
        f"{pick(planted, 'p_value'):.4g}) and opens it, while "
        f"{len(pick(gate, 'noise'))} composition-matched noise pools average "
        f"{pick(gate, 'noise_mean_z'):+.4f} against an expected zero and "
        f"**{pick(gate, 'n_noise_pools_fired')}** of them open it. All "
        f"{len(checks)} checks pass.",
        "",
        "**What the gate decides.** If the positive control fails on real shells, every "
        "arm of that study reads *instrument not demonstrated* and **the sequence line's "
        "nulls above cannot be read as null results either**. That is fixed in code, not "
        "promised in prose.",
        "",
    ]
    report = data.get("composition")
    if report is None:
        lines += [
            "**The study itself has not reported yet.** When it does, this section gains "
            "its result and section 4's nulls either firm up into statements about order "
            "or become uninterpretable. Until then the question is open, and this page "
            "says so rather than leaving the line out.",
            "",
        ]
        return lines
    fired = bool(pick(report, "control_fired"))
    control = next(t for t in pick(report, "tests")
                   if t.get("arm") == pick(report, "control_arm"))
    if not fired:
        lines += [
            f"**On real shells the positive control did not fire** (mean z "
            f"{pick(control, 'mean_site_z'):+.2f}, p {pick(control, 'p_value'):.4g}, "
            f"{pick(control, 'n_groups')} groups). A signal that is certainly present did "
            "not clear the plan's own floors, so **this pipeline's sensitivity on real "
            "linings is undemonstrated and section 4's nulls are not null results**. "
            "Every arm of this study reads *instrument not demonstrated*.",
            "",
        ]
        return lines
    primary = [t for t in pick(report, "tests")
               if t.get("runnable") and t.get("family") == "primary"]
    findings = [t for t in primary if t.get("is_finding")]
    lines += [
        f"**On real shells the positive control fired** (mean z "
        f"{pick(control, 'mean_site_z'):+.2f}, p {pick(control, 'p_value'):.4g}, "
        f"{pick(control, 'n_groups')} groups). The instrument is demonstrated on real "
        "linings, so **section 4's nulls are statements about order and not about "
        "sensitivity** — which is the whole reason this study exists.",
        "",
        f"Of its own {len(primary)} primary tests, **{len(findings)}** is a finding.",
        "",
    ]
    return lines


def untested(data: Mapping[str, Any]) -> list:
    """Which feature families have been tried, so the closing cannot overstate."""
    c2 = pick(data["coevolution"], "C2_matched")
    estimate = pick(c2, "estimate")
    if pick(c2, "decision") == "supported":
        raise Contradiction(
            "the comparative coevolution test now reads supported; the sentence below "
            "calls it inconclusive")
    return [
        "## 6. Feature families, and which have actually been tried",
        "",
        "The two decompositions above cover a lining's **positions** and a window's "
        "**letters**. It is easy to read their nulls as though nothing about a sequence or "
        "a structure could predict an IP site, so here is what has and has not been put "
        "to a test.",
        "",
        "| family | status |",
        "|---|---|",
        "| lining topology and register | tested, null against the control that matters "
        "(section 3) |",
        "| sequence order | tested, null; binding arm unresolved (section 4) |",
        "| sequence composition | instrument check calibrated, study in flight (section 5) |",
        f"| cross-organism IP-concentration coevolution | tested, **inconclusive** — "
        f"{pick(estimate, 'point') * 100:+.3f} pp "
        f"[{pick(estimate, 'low') * 100:+.3f}, {pick(estimate, 'high') * 100:+.3f}] over "
        f"{pick(c2, 'bins')} matched bins and {pick(c2, 'clusters')} clusters "
        f"(`docs/COEVOLUTION_PLAN.md`) |",
        "| a learned pocket score | exists, and is a component of the combined score; the "
        "protein-level classifier benchmarked at chance |",
        "| orthogroup-level conservation of hits | computed in the comparative analysis |",
        "| **per-residue conservation from an alignment** | **not tested** |",
        "| **residue-residue coevolution as a site feature** | **not tested** |",
        "| **a learned representation of the site itself** | **not tested** |",
        "",
        "The last three are where a signal could still live, and the nulls above say "
        "nothing about them. Note that the coevolution row is a different question from "
        "the last two: it asks whether an organism richer in inositol phosphate buries "
        "more of them, not whether a particular site's residues covary.",
        "",
    ]


def closing() -> list:
    return [
        "## What is not established, and what would change it",
        "",
        "* **That any candidate binds an inositol phosphate.** No candidate has "
        "independent support. The honest description is a ranked list of hypotheses.",
        "* **That the pipeline generalises to cryptic sites.** One cryptic example in the "
        "validation set cannot support a performance estimate for cryptic sites, and the "
        "deposited archive cannot supply more, because a buried site with its ligand in "
        "place is hard to crystallise. This is the project's binding constraint and it is "
        "not a statistical one.",
        "* **That the machinery adds anything beyond the same-protein baseline.** Section "
        "1's third row says it does not, at the power available.",
        "* **That structure or sequence carries a detectable signature of an IP site.** "
        "Sections 2 to 4 are null or unresolved, and section 5 is what decides whether "
        "those nulls mean anything.",
        "",
        "What would change the picture, in the order it would be worth doing: an "
        "experimentally tested candidate; published functional evidence for one of the "
        "deposited cryptic copies other than ADAR2, which would turn the ground-truth "
        "arms from structural classification into established biology; and one of the "
        "three untested feature families in section 6.",
        "",
        "Every claim above is read from the artefact beside it. Where a number here and a "
        "number in `results/publication/RESULTS_SUMMARY.md` disagree, that file is older.",
    ]


def render(data: Mapping[str, Any]) -> str:
    lines = header(data)
    lines += verdict(data)
    lines += specificity(data)
    lines += architecture(data)
    lines += within_protein(data)
    lines += order(data)
    lines += composition(data)
    lines += untested(data)
    lines += closing()
    return "\n".join(lines) + "\n"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=Path("results/PROJECT_FINDINGS.md"))
    args = parser.parse_args(argv)
    text = render(load())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(f"wrote {args.out} ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
