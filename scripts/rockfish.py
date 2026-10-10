#!/usr/bin/env python3
"""Run the pre-registered studies on a SLURM cluster.

    python scripts/rockfish.py validate   --site hpc/sites/rockfish.yml
    python scripts/rockfish.py preflight  --site hpc/sites/rockfish.yml
    python scripts/rockfish.py plan       --study hpc/studies/flexible.yml
    python scripts/rockfish.py submit     --study hpc/studies/flexible.yml
    python scripts/rockfish.py status     --study hpc/studies/flexible.yml
    python scripts/rockfish.py collect    --study hpc/studies/flexible.yml --stage dock
    python scripts/rockfish.py resume     --study hpc/studies/flexible.yml --stage dock

The science lives in `scripts/<study>.py` and is identical to what the GitHub
workflows run. This script only decides where it runs, how it is sharded, and
whether the results add up.

Three refusals are deliberate and are the reason to go through this script
rather than writing sbatch files by hand:

* a study whose pre-registered plan documents do not exist is not submitted;
* a site whose claims have never been checked against the cluster is not
  submitted against unless that is stated on the command line;
* a stage whose shards do not reconcile does not produce a merged result.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from cryptic_ip.hpc.collect import (  # noqa: E402
    ACCOUNTED, CollectError, collect_stage, provenance, take_census,
)
from cryptic_ip.hpc.manifest import ManifestError, StudyManifest, load_manifest, total_tasks  # noqa: E402
from cryptic_ip.hpc.preflight import markdown, run_preflight  # noqa: E402
from cryptic_ip.hpc.render import render_study, write_scripts  # noqa: E402
from cryptic_ip.hpc.site import SiteConfig, SiteError, load_site  # noqa: E402

DEFAULT_SITE = "hpc/sites/rockfish.yml"
SUBMISSION_FILE = "submission.json"


# ------------------------------------------------------------------- utilities

def _load(site_path: str, study_path: Optional[str]):
    site = load_site(ROOT / site_path if not Path(site_path).is_absolute() else Path(site_path))
    manifest = None
    if study_path:
        path = Path(study_path)
        manifest = load_manifest(path if path.is_absolute() else ROOT / path)
    return site, manifest


def _work_dir(site: SiteConfig, manifest: StudyManifest, override: Optional[str]) -> Path:
    if override:
        return Path(override).expanduser()
    return Path(site.scratch_for(manifest.study.lower()))


def _echo(text: str = "") -> None:
    print(text)


def _plan_gate(manifest: StudyManifest) -> List[str]:
    """Pre-registration documents that must exist before anything is submitted."""
    return manifest.missing_plan_docs(ROOT)


# ------------------------------------------------------------------- validate

def cmd_validate(args: argparse.Namespace) -> int:
    try:
        site, _ = _load(args.site, None)
    except (SiteError, OSError) as exc:
        _echo(f"site file: {exc}")
        return 1
    _echo(f"site {site.name}: {len(site.partitions)} partition(s), "
          f"max_array_size {site.max_array_size}, throttle {site.max_array_throttle}")
    if not site.is_verified:
        _echo("  warning: this site file has never been checked against the live cluster "
              "(run `preflight` on a login node)")

    paths = [Path(p) for p in args.study] if args.study else sorted((ROOT / "hpc" / "studies").glob("*.yml"))
    if not paths:
        _echo("no study manifests found under hpc/studies/")
        return 1
    failures = 0
    for path in paths:
        full = path if path.is_absolute() else ROOT / path
        try:
            manifest = load_manifest(full)
        except (ManifestError, OSError) as exc:
            _echo(f"{path}: {exc}")
            failures += 1
            continue
        problems = manifest.problems_against(site)
        missing = _plan_gate(manifest)
        order = " -> ".join(stage.name for stage in manifest.order())
        _echo(f"study {manifest.study} ({path.name}): {total_tasks(manifest)} task(s); {order}")
        for stage in manifest.stages:
            shards = f"{stage.task_count} shards" if stage.is_array else "1 task"
            res = stage.resources
            gpu = f", {res.gpus} gpu" if res.gpus else ""
            _echo(f"    {stage.name}: {shards}, {res.cpus} cpu, {res.mem_gb:g} GB, "
                  f"{res.walltime_text}, partition {res.partition}{gpu}")
        for problem in problems:
            _echo(f"    does not fit: {problem}")
            failures += 1
        for doc in missing:
            _echo(f"    missing pre-registered plan: {doc}")
            failures += 1
    _echo()
    _echo("ok" if not failures else f"{failures} problem(s)")
    return 0 if not failures else 1


# ------------------------------------------------------------------- preflight

def cmd_preflight(args: argparse.Namespace) -> int:
    site, _ = _load(args.site, None)
    manifests = []
    for path in sorted((ROOT / "hpc" / "studies").glob("*.yml")):
        try:
            manifests.append(load_manifest(path))
        except ManifestError:
            continue
    report = run_preflight(site, manifests=manifests)
    text = markdown(report)
    _echo(text)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(report.as_dict(), indent=2, sort_keys=True) + "\n")
        _echo(f"wrote {args.out}")
    if args.markdown:
        Path(args.markdown).parent.mkdir(parents=True, exist_ok=True)
        Path(args.markdown).write_text(text)
        _echo(f"wrote {args.markdown}")
    if report.verdict == "blocked":
        return 1
    if report.verdict == "indeterminate":
        _echo("Run this on a Rockfish login node to actually check the site file's claims.")
    return 0


# ------------------------------------------------------------------------ plan

def cmd_plan(args: argparse.Namespace) -> int:
    site, manifest = _load(args.site, args.study)
    assert manifest is not None
    work_dir = _work_dir(site, manifest, args.work_dir)
    rendered = render_study(
        manifest, site, work_dir=str(work_dir), repo_root=str(args.repo_root or ROOT),
    )
    _echo(f"study {manifest.study}: {manifest.title}")
    _echo(f"work dir {work_dir}")
    _echo(f"repo     {args.repo_root or ROOT}")
    _echo()
    for item in rendered:
        stage = manifest.stage(item.stage)
        dep = f" after {item.after}:{','.join(item.depends_on)}" if item.depends_on else ""
        _echo(f"  {item.script_name}: {len(item.shards)} task(s), "
              f"shards {item.shards[0]}..{item.shards[-1]}, {stage.resources.walltime_text}{dep}")
    missing = _plan_gate(manifest)
    for doc in missing:
        _echo(f"  BLOCKED: missing pre-registered plan {doc}")
    if args.out_dir:
        paths = write_scripts(rendered, Path(args.out_dir))
        _echo()
        _echo(f"wrote {len(paths)} script(s) to {args.out_dir}")
    return 1 if missing else 0


# ---------------------------------------------------------------------- submit

def _sbatch(script: Path, dependency: Optional[str], dry_run: bool) -> str:
    command = ["sbatch", "--parsable"]
    if dependency:
        command.append(f"--dependency={dependency}")
    command.append(str(script))
    if dry_run:
        _echo("  would run: " + " ".join(shlex.quote(part) for part in command))
        return f"DRY-{script.stem}"
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"sbatch refused {script.name}: {result.stderr.strip() or result.stdout.strip()}")
    return result.stdout.strip().split(";")[0]


def cmd_submit(args: argparse.Namespace) -> int:
    site, manifest = _load(args.site, args.study)
    assert manifest is not None

    missing = _plan_gate(manifest)
    if missing:
        _echo(f"Refusing to submit study {manifest.study}: its pre-registered plan is not in the repository.")
        for doc in missing:
            _echo(f"  missing: {doc}")
        _echo("Write and commit the plan first. That order is the point of pre-registration.")
        return 2

    if not site.is_verified and not args.allow_unverified_site:
        _echo(f"Refusing to submit against site {site.name!r}: its limits have never been checked "
              "against the live cluster.")
        _echo(f"  Run: python scripts/rockfish.py preflight --site {args.site}")
        _echo("  Then set verified.date and verified.by in the site file, or pass --allow-unverified-site "
              "to submit anyway.")
        return 2

    if not args.skip_preflight and not args.dry_run:
        report = run_preflight(site, manifests=[manifest])
        blocking = report.blocking
        if blocking:
            _echo(f"Refusing to submit: {len(blocking)} preflight check(s) failed.")
            for check in blocking:
                _echo(f"  {check.name}: expected {check.expected}, got {check.observed}")
                if check.remedy:
                    _echo(f"    {check.remedy}")
            _echo("Pass --skip-preflight only if you have a reason these do not apply.")
            return 2

    work_dir = _work_dir(site, manifest, args.work_dir)
    repo_root = Path(args.repo_root or ROOT)
    only: Optional[Dict[str, Sequence[int]]] = None
    if args.stage:
        only = {
            stage.name: (range(stage.task_count) if stage.name in args.stage else ())
            for stage in manifest.stages
        }
    rendered = render_study(
        manifest, site, work_dir=str(work_dir), repo_root=str(repo_root), only_shards=only,
    )
    if not rendered:
        _echo("nothing to submit")
        return 1

    # SLURM writes a task's output file before the script runs, so the log
    # directory has to exist at submit time rather than be made by the script.
    script_dir = work_dir / "scripts"
    if not args.dry_run:
        (work_dir / "logs").mkdir(parents=True, exist_ok=True)
    paths = {item.script_name: path for item, path in zip(rendered, write_scripts(rendered, script_dir))} \
        if not args.dry_run else {item.script_name: script_dir / item.script_name for item in rendered}

    _echo(f"study {manifest.study} -> {work_dir}")
    job_ids: Dict[str, List[str]] = {}
    submitted: List[Dict[str, object]] = []
    for item in rendered:
        deps: List[str] = []
        for dep in item.depends_on:
            deps.extend(job_ids.get(dep, []))
        dependency = f"after{item.after}:" + ":".join(deps) if deps else None
        try:
            job_id = _sbatch(paths[item.script_name], dependency, args.dry_run)
        except RuntimeError as exc:
            _echo(f"  {exc}")
            _echo("  Earlier jobs in this study are already queued; cancel them with scancel if you stop here.")
            return 1
        job_ids.setdefault(item.stage, []).append(job_id)
        submitted.append({
            "stage": item.stage,
            "chunk": item.chunk_index,
            "job_id": job_id,
            "tasks": len(item.shards),
            "shards": list(item.shards),
            "script": str(paths[item.script_name]),
            "dependency": dependency,
        })
        _echo(f"  {item.stage:<10} {job_id:>12}  {len(item.shards):>4} task(s)"
              + (f"  depends on {dependency}" if dependency else ""))

    record = {
        "study": manifest.study,
        "site": site.name,
        "work_dir": str(work_dir),
        "repo_root": str(repo_root),
        "manifest": str(manifest.source_path),
        "dry_run": bool(args.dry_run),
        "jobs": submitted,
    }
    if not args.dry_run:
        (work_dir / SUBMISSION_FILE).write_text(json.dumps(record, indent=2) + "\n")
        _echo(f"\nwrote {work_dir / SUBMISSION_FILE}")
    _echo("\nFollow it with: python scripts/rockfish.py status --study "
          f"{args.study}" + (f" --work-dir {args.work_dir}" if args.work_dir else ""))
    return 0


# ---------------------------------------------------------------------- status

def cmd_status(args: argparse.Namespace) -> int:
    site, manifest = _load(args.site, args.study)
    assert manifest is not None
    work_dir = _work_dir(site, manifest, args.work_dir)
    _echo(f"study {manifest.study} in {work_dir}")
    if not work_dir.exists():
        _echo("  nothing there yet")
        return 0
    for stage in manifest.order():
        stage_dir = work_dir / stage.name
        if not stage.is_array:
            output = stage_dir / stage.output if stage.output else None
            done = "done" if output and output.is_file() else "not finished"
            _echo(f"  {stage.name:<10} single task: {done}")
            continue
        census = take_census(manifest, stage, stage_dir)
        counts = census.counts()
        _echo(f"  {stage.name:<10} " + ", ".join(f"{k}={v}" for k, v in counts.items() if v))
        if census.missing and args.verbose:
            _echo(f"      missing: {', '.join(str(s) for s in census.missing[:20])}")
    submission = work_dir / SUBMISSION_FILE
    if submission.is_file():
        record = json.loads(submission.read_text())
        ids = [str(job["job_id"]) for job in record.get("jobs", [])]
        _echo(f"  jobs: {', '.join(ids)}")
    return 0


# --------------------------------------------------------------------- collect

def cmd_collect(args: argparse.Namespace) -> int:
    site, manifest = _load(args.site, args.study if args.study else args.manifest)
    assert manifest is not None
    work_dir = Path(args.work_dir) if args.work_dir else _work_dir(site, manifest, None)
    try:
        census = collect_stage(
            manifest, args.stage, work_dir,
            records_out=Path(args.records) if args.records else None,
            census_out=Path(args.census_out) if args.census_out else None,
            allow_missing=args.allow_missing,
        )
    except CollectError as exc:
        _echo(str(exc))
        return 1
    _echo(census.summary())
    if args.provenance:
        Path(args.provenance).parent.mkdir(parents=True, exist_ok=True)
        Path(args.provenance).write_text(provenance(manifest, census, site_name=site.name))
        _echo(f"wrote {args.provenance}")
    if args.then:
        command = shlex.split(args.then)
        _echo(f"running: {args.then}")
        result = subprocess.run(command, check=False)
        return result.returncode
    return 0


# ---------------------------------------------------------------------- resume

def cmd_resume(args: argparse.Namespace) -> int:
    site, manifest = _load(args.site, args.study)
    assert manifest is not None
    work_dir = _work_dir(site, manifest, args.work_dir)
    stage = manifest.stage(args.stage)
    if not stage.is_array:
        _echo(f"stage {args.stage!r} is not sharded; re-submit it with `submit --stage {args.stage}`")
        return 1
    census = take_census(manifest, stage, work_dir / stage.name)
    wanted = sorted(census.missing + [s for state in ACCOUNTED for s in census.shards_in(state)]) \
        if args.include_failed else sorted(census.missing + census.shards_in("killed"))
    if not wanted:
        _echo(f"nothing to resume: {census.summary()}")
        return 0
    _echo(f"{len(wanted)} shard(s) to re-run: {', '.join(str(s) for s in wanted[:20])}"
          + (" ..." if len(wanted) > 20 else ""))
    rendered = render_study(
        manifest, site,
        work_dir=str(work_dir), repo_root=str(args.repo_root or ROOT),
        only_shards={s.name: (wanted if s.name == stage.name else ()) for s in manifest.stages},
    )
    script_dir = work_dir / "scripts" / "resume"
    paths = write_scripts(rendered, script_dir) if not args.dry_run else [
        script_dir / item.script_name for item in rendered
    ]
    for item, path in zip(rendered, paths):
        job_id = _sbatch(path, None, args.dry_run)
        _echo(f"  {item.stage} {job_id} {len(item.shards)} task(s)")
    return 0


# ---------------------------------------------------------------- check-inputs

def cmd_check_inputs(args: argparse.Namespace) -> int:
    """Fail before an array queues, not after every task has queued and died."""
    work_dir = Path(args.work_dir)
    found: Dict[str, object] = {}
    problems: List[str] = []
    for relative in args.require:
        path = work_dir / relative
        ok = path.is_file()
        found[relative] = {"path": str(path), "exists": ok, "bytes": path.stat().st_size if ok else 0}
        if not ok:
            problems.append(f"missing file {path}")
    for relative in args.require_dir:
        path = work_dir / relative
        ok = path.is_dir()
        count = sum(1 for _ in path.iterdir()) if ok else 0
        found[relative] = {"path": str(path), "exists": ok, "entries": count}
        if not ok:
            problems.append(f"missing directory {path}")
        elif count == 0:
            problems.append(f"empty directory {path}")
    payload = {"work_dir": str(work_dir), "ok": not problems, "problems": problems, "inputs": found}
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for problem in problems:
        _echo(problem)
    if problems:
        _echo("Stage the inputs into the work directory before submitting the array.")
        return 1
    _echo(f"inputs ok in {work_dir}")
    return 0


# ------------------------------------------------------------------------ main

def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--site", default=os.environ.get("ROCKFISH_SITE", DEFAULT_SITE),
                        help="site YAML (default: %(default)s)")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("validate", help="check the site file and every study manifest")
    p.add_argument("--study", action="append", help="one manifest (default: all under hpc/studies)")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("preflight", help="check the site file's claims against the cluster")
    p.add_argument("--out", help="write the report as JSON")
    p.add_argument("--markdown", help="write the report as markdown")
    p.set_defaults(func=cmd_preflight)

    p = sub.add_parser("plan", help="show, and optionally write, the sbatch scripts")
    p.add_argument("--study", required=True)
    p.add_argument("--work-dir")
    p.add_argument("--repo-root")
    p.add_argument("--out-dir", help="write the scripts here for inspection")
    p.set_defaults(func=cmd_plan)

    p = sub.add_parser("submit", help="submit a study, stages wired by dependency")
    p.add_argument("--study", required=True)
    p.add_argument("--work-dir")
    p.add_argument("--repo-root")
    p.add_argument("--stage", action="append", help="submit only these stages")
    p.add_argument("--dry-run", action="store_true", help="print the sbatch calls instead of running them")
    p.add_argument("--skip-preflight", action="store_true")
    p.add_argument("--allow-unverified-site", action="store_true",
                   help="submit although nobody has checked the site file against the cluster")
    p.set_defaults(func=cmd_submit)

    p = sub.add_parser("status", help="per-stage shard census")
    p.add_argument("--study", required=True)
    p.add_argument("--work-dir")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(func=cmd_status)

    p = sub.add_parser("collect", help="reconcile and merge a sharded stage")
    p.add_argument("--study", help="manifest path")
    p.add_argument("--manifest", help="manifest path (alias, for use inside a stage command)")
    p.add_argument("--stage", required=True)
    p.add_argument("--work-dir", required=True)
    p.add_argument("--records", help="write the merged records here")
    p.add_argument("--census-out", help="write the shard census here")
    p.add_argument("--provenance", help="write a provenance stub here")
    p.add_argument("--allow-missing", type=int, default=0,
                   help="how many shards may have left nothing at all (default 0)")
    p.add_argument("--then", help="command to run once the records are merged")
    p.set_defaults(func=cmd_collect)

    p = sub.add_parser("resume", help="re-run the shards that did not land")
    p.add_argument("--study", required=True)
    p.add_argument("--stage", required=True)
    p.add_argument("--work-dir")
    p.add_argument("--repo-root")
    p.add_argument("--include-failed", action="store_true",
                   help="also re-run shards that failed in the code, not only the ones the cluster took")
    p.add_argument("--dry-run", action="store_true")
    p.set_defaults(func=cmd_resume)

    p = sub.add_parser("check-inputs", help="verify staged inputs before an array queues")
    p.add_argument("--work-dir", required=True)
    p.add_argument("--require", action="append", default=[])
    p.add_argument("--require-dir", action="append", default=[])
    p.add_argument("--out")
    p.set_defaults(func=cmd_check_inputs)

    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except (SiteError, ManifestError) as exc:
        _echo(f"error: {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
