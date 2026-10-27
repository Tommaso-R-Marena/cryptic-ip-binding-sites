"""Check a site file's claims against the cluster it describes.

The site file is a hypothesis. Every value in it - partition names, walltime
caps, how many GPUs a node has and how much memory each one holds, SLURM's
`MaxArraySize` - was written without access to the machine, and a wrong value
does not announce itself: it is a job that sits in the queue against a
partition that does not exist, or an array that is refused at submit time, or a
co-folding job that runs out of VRAM four hours in.

So the parsers here are pure functions over scheduler output, testable without
a cluster, and the runner is a thin shell around them. `submit` will not run a
study against an unverified site unless told to in so many words.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from cryptic_ip.hpc.manifest import StudyManifest
from cryptic_ip.hpc.site import SiteConfig, format_walltime, parse_walltime

PASS = "pass"
FAIL = "fail"
WARN = "warn"
UNKNOWN = "unknown"

#: Checks whose failure means do not submit. A warning is for something a
#: person should see; a failure is for something that will waste an allocation.
BLOCKING = (FAIL,)


@dataclass
class Check:
    """One claim, what was expected, and what the cluster actually said."""

    name: str
    status: str
    expected: str = ""
    observed: str = ""
    detail: str = ""
    remedy: str = ""

    def as_dict(self) -> Dict[str, str]:
        return {
            "name": self.name,
            "status": self.status,
            "expected": self.expected,
            "observed": self.observed,
            "detail": self.detail,
            "remedy": self.remedy,
        }


@dataclass
class Report:
    """Every check, and whether the site may be submitted against."""

    site: str
    checks: List[Check] = field(default_factory=list)
    site_path: Optional[str] = None

    def add(self, check: Check) -> Check:
        self.checks.append(check)
        return check

    @property
    def blocking(self) -> List[Check]:
        return [check for check in self.checks if check.status in BLOCKING]

    @property
    def unknown(self) -> List[Check]:
        return [check for check in self.checks if check.status == UNKNOWN]

    @property
    def scheduler_reachable(self) -> bool:
        """Whether this machine can see the scheduler at all."""
        return any(c.name == "scheduler" and c.status == PASS for c in self.checks)

    @property
    def verdict(self) -> str:
        """``blocked``, ``indeterminate`` or ``ok``.

        ``indeterminate`` is the honest answer from a machine with no scheduler:
        nothing failed, but nothing about the cluster was actually checked, and
        reporting that as "ok to submit" would be a green light nobody earned.
        """
        if self.blocking:
            return "blocked"
        if not self.scheduler_reachable:
            return "indeterminate"
        return "ok"

    @property
    def ok(self) -> bool:
        """Whether a study may be submitted on the strength of this report."""
        return self.verdict == "ok"

    def counts(self) -> Dict[str, int]:
        out = {PASS: 0, FAIL: 0, WARN: 0, UNKNOWN: 0}
        for check in self.checks:
            out[check.status] = out.get(check.status, 0) + 1
        return out

    def as_dict(self) -> Dict[str, object]:
        return {
            "site": self.site,
            "site_path": self.site_path,
            "ok": self.ok,
            "verdict": self.verdict,
            "counts": self.counts(),
            "checks": [check.as_dict() for check in self.checks],
        }


# --------------------------------------------------------------------- parsers

def parse_sinfo_partitions(text: str) -> Dict[str, Dict[str, object]]:
    """Parse `sinfo -h -O 'Partition:|,Time:|,CPUs:|,Memory:|,Gres:'`-style rows.

    Expected one row per partition, pipe separated:
    ``name|timelimit|cpus|memory_mb|gres``. A trailing ``*`` marks the default
    partition and is not part of the name.
    """
    out: Dict[str, Dict[str, object]] = {}
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = [part.strip() for part in line.split("|")]
        if len(fields) < 2:
            continue
        name = fields[0].rstrip("*")
        if not name:
            continue
        row: Dict[str, object] = {}
        try:
            row["max_walltime"] = parse_walltime(fields[1])
        except ValueError:
            row["max_walltime"] = None
        if len(fields) > 2 and fields[2].isdigit():
            row["cpus"] = int(fields[2])
        if len(fields) > 3:
            mem = fields[3].rstrip("+")
            if mem.isdigit():
                row["mem_gb"] = int(mem) / 1024.0
        if len(fields) > 4:
            row["gres"] = fields[4]
        # Several rows can describe one partition, one per node type. Keep the
        # most generous, because that is what a job may actually ask for.
        if name in out:
            row = _merge_partition_rows(out[name], row)
        out[name] = row
    return out


def _merge_partition_rows(a: Mapping[str, object], b: Mapping[str, object]) -> Dict[str, object]:
    merged: Dict[str, object] = dict(a)
    for key, value in b.items():
        current = merged.get(key)
        if value is None:
            continue
        if current is None:
            merged[key] = value
        elif key in {"max_walltime", "cpus", "mem_gb"}:
            if current == -1 or value == -1:
                merged[key] = -1
            else:
                merged[key] = max(current, value)  # type: ignore[call-overload]
        elif key == "gres" and str(value) not in str(current):
            merged[key] = f"{current},{value}"
    return merged


def parse_scontrol_config(text: str) -> Dict[str, str]:
    """Parse `scontrol show config` into its ``key = value`` pairs."""
    out: Dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        if key and " " not in key:
            out[key] = value.strip()
    return out


def parse_gres_models(gres: str) -> Dict[str, int]:
    """Models and counts in a SLURM GRES string such as ``gpu:a100:4(S:0-1)``."""
    out: Dict[str, int] = {}
    for item in str(gres).split(","):
        item = item.strip()
        if not item or item in {"(null)", "null"}:
            continue
        item = item.split("(")[0]
        parts = item.split(":")
        if parts[0] != "gpu":
            continue
        if len(parts) >= 3:
            model, count = parts[1], parts[2]
        elif len(parts) == 2:
            model, count = "", parts[1]
        else:
            continue
        try:
            out[model.lower()] = out.get(model.lower(), 0) + int(count)
        except ValueError:
            continue
    return out


def parse_nvidia_smi(text: str) -> List[Dict[str, object]]:
    """Parse `nvidia-smi --query-gpu=name,memory.total --format=csv,noheader`."""
    out: List[Dict[str, object]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = [part.strip() for part in line.split(",")]
        if len(fields) < 2:
            continue
        name = fields[0]
        memory = fields[1].replace("MiB", "").replace("MB", "").strip()
        try:
            vram_gb = int(memory) / 1024.0
        except ValueError:
            continue
        out.append({"name": name, "vram_gb": vram_gb})
    return out


# ---------------------------------------------------------------------- runner

def _run(command: Sequence[str], timeout: int = 60) -> Optional[str]:
    """Run a read-only query. ``None`` when the command is absent or fails."""
    if shutil.which(command[0]) is None:
        return None
    try:
        result = subprocess.run(
            list(command), capture_output=True, text=True, timeout=timeout, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if result.returncode != 0:
        return None
    return result.stdout


SINFO_QUERY = ["sinfo", "-h", "-a", "-o", "%R|%l|%c|%m|%G"]


def run_preflight(
    site: SiteConfig,
    *,
    manifests: Sequence[StudyManifest] = (),
    runner=_run,
) -> Report:
    """Check ``site`` against the live cluster, as far as this node can see it."""
    report = Report(site=site.name, site_path=str(site.source_path) if site.source_path else None)

    report.add(_check_verified(site))
    report.add(_check_scheduler(runner))

    sinfo_text = runner(SINFO_QUERY)
    partitions = parse_sinfo_partitions(sinfo_text) if sinfo_text else {}
    if sinfo_text is None:
        report.add(Check(
            name="partitions",
            status=UNKNOWN,
            expected=", ".join(sorted(p.name for p in site.partitions.values())),
            observed="sinfo is not available on this machine",
            detail="Partition names, walltime caps and node sizes could not be checked.",
            remedy="Run this on a Rockfish login node.",
        ))
    else:
        for check in _check_partitions(site, partitions):
            report.add(check)

    config_text = runner(["scontrol", "show", "config"])
    report.add(_check_max_array(site, config_text))
    report.add(_check_account(site, runner))
    report.add(_check_scratch(site))
    for check in _check_gpu(site, partitions, runner):
        report.add(check)
    for check in _check_environment():
        report.add(check)
    for manifest in manifests:
        report.add(_check_manifest_fits(manifest, site))
    return report


def _check_verified(site: SiteConfig) -> Check:
    if site.is_verified:
        return Check(
            name="site file verified",
            status=PASS,
            expected="a dated verification",
            observed=f"verified {site.verified_on}" + (f" by {site.verified_by}" if site.verified_by else ""),
        )
    return Check(
        name="site file verified",
        status=WARN,
        expected="a dated verification",
        observed="never verified against the live cluster",
        detail="Every limit in the site file is still an assumption.",
        remedy="Once this report passes, set verified.date and verified.by in the site file.",
    )


def _check_scheduler(runner) -> Check:
    version = runner(["sbatch", "--version"])
    if version:
        return Check(name="scheduler", status=PASS, expected="sbatch on PATH", observed=version.strip())
    return Check(
        name="scheduler",
        status=UNKNOWN,
        expected="sbatch on PATH",
        observed="sbatch not found",
        detail="Nothing can be submitted from this machine.",
        remedy="Run on a cluster login node.",
    )


def _check_partitions(site: SiteConfig, observed: Mapping[str, Dict[str, object]]) -> List[Check]:
    out: List[Check] = []
    for role, partition in sorted(site.partitions.items()):
        row = observed.get(partition.name)
        if row is None:
            out.append(Check(
                name=f"partition {role} ({partition.name})",
                status=FAIL,
                expected=f"partition {partition.name!r} exists",
                observed="not in sinfo output",
                detail=f"Known partitions: {', '.join(sorted(observed)) or 'none reported'}.",
                remedy=f"Correct partitions.{role}.name in the site file.",
            ))
            continue
        problems: List[str] = []
        actual_walltime = row.get("max_walltime")
        if isinstance(actual_walltime, int) and partition.max_walltime != -1:
            if actual_walltime != -1 and partition.max_walltime > actual_walltime:
                problems.append(
                    f"site claims {format_walltime(partition.max_walltime)} but the partition caps at "
                    f"{format_walltime(actual_walltime)}"
                )
        actual_cpus = row.get("cpus")
        if isinstance(actual_cpus, int) and partition.max_cpus_per_node and partition.max_cpus_per_node > actual_cpus:
            problems.append(f"site claims {partition.max_cpus_per_node} cpus per node, sinfo reports {actual_cpus}")
        actual_mem = row.get("mem_gb")
        if isinstance(actual_mem, float) and partition.max_mem_gb_per_node:
            if partition.max_mem_gb_per_node > actual_mem + 1:
                problems.append(
                    f"site claims {partition.max_mem_gb_per_node:g} GB per node, sinfo reports {actual_mem:.0f} GB"
                )
        out.append(Check(
            name=f"partition {role} ({partition.name})",
            status=FAIL if problems else PASS,
            expected=(
                f"<= {format_walltime(partition.max_walltime)}, {partition.max_cpus_per_node or '?'} cpus, "
                f"{partition.max_mem_gb_per_node or '?'} GB"
            ),
            observed="; ".join(
                f"{key}={row[key]}" for key in ("max_walltime", "cpus", "mem_gb", "gres") if key in row
            ),
            detail="; ".join(problems),
            remedy=f"Lower partitions.{role} in the site file to what the cluster allows." if problems else "",
        ))
    return out


def _check_max_array(site: SiteConfig, config_text: Optional[str]) -> Check:
    if not config_text:
        return Check(
            name="MaxArraySize",
            status=UNKNOWN,
            expected=str(site.max_array_size),
            observed="scontrol is not available on this machine",
            detail=(
                "This is the limit study J's GitHub equivalent exceeded silently. "
                "Chunking is sized from it, so it is worth confirming."
            ),
            remedy="Run on a login node, or read `scontrol show config | grep MaxArraySize`.",
        )
    config = parse_scontrol_config(config_text)
    raw = config.get("MaxArraySize")
    if raw is None:
        return Check(
            name="MaxArraySize",
            status=UNKNOWN,
            expected=str(site.max_array_size),
            observed="not reported by scontrol",
        )
    try:
        actual = int(raw)
    except ValueError:
        return Check(name="MaxArraySize", status=UNKNOWN, expected=str(site.max_array_size), observed=raw)
    if site.max_array_size > actual:
        return Check(
            name="MaxArraySize",
            status=FAIL,
            expected=str(site.max_array_size),
            observed=str(actual),
            detail="Chunks sized from the site file would be refused at submit time.",
            remedy=f"Set max_array_size: {actual} in the site file.",
        )
    return Check(name="MaxArraySize", status=PASS, expected=str(site.max_array_size), observed=str(actual))


def _check_account(site: SiteConfig, runner) -> Check:
    if not site.account:
        return Check(
            name="account",
            status=WARN,
            expected="an allocation to charge",
            observed="no account set in the site file",
            detail="Clusters that require an account will reject every job.",
            remedy="Run `sacctmgr show assoc user=$USER format=Account` and set `account` in the site file.",
        )
    text = runner(["sacctmgr", "-nP", "show", "assoc", f"user={os.environ.get('USER', '')}", "format=Account"])
    if text is None:
        return Check(
            name="account",
            status=UNKNOWN,
            expected=site.account,
            observed="sacctmgr is not available on this machine",
        )
    accounts = {line.strip() for line in text.splitlines() if line.strip()}
    if site.account in accounts:
        return Check(name="account", status=PASS, expected=site.account, observed=", ".join(sorted(accounts)))
    return Check(
        name="account",
        status=FAIL,
        expected=site.account,
        observed=", ".join(sorted(accounts)) or "none",
        detail="Jobs would be rejected, or charged somewhere unintended.",
        remedy="Set `account` to one of the accounts listed above.",
    )


def _check_scratch(site: SiteConfig) -> Check:
    root = Path(os.path.expanduser(os.path.expandvars(site.scratch_root)))
    parent = root if root.exists() else root.parent
    if not parent.exists():
        return Check(
            name="scratch",
            status=FAIL,
            expected=f"{root} is creatable",
            observed=f"{parent} does not exist",
            detail="Study outputs have nowhere to go.",
            remedy="Correct `scratch_root`, or create the parent directory.",
        )
    if not os.access(parent, os.W_OK):
        return Check(
            name="scratch",
            status=FAIL,
            expected=f"{root} is writable",
            observed=f"{parent} is not writable",
            remedy="Correct `scratch_root` to a filesystem you can write to.",
        )
    try:
        stat = os.statvfs(str(parent))
        free_gb = stat.f_bavail * stat.f_frsize / 1e9
    except OSError:
        free_gb = float("nan")
    status = PASS
    detail = ""
    if free_gb == free_gb and free_gb < 50:
        status = WARN
        detail = "Under 50 GB free. Docking shards write receptors and poses per copy."
    return Check(
        name="scratch",
        status=status,
        expected=f"{root} writable with room",
        observed=f"{parent} writable, {free_gb:.0f} GB free" if free_gb == free_gb else f"{parent} writable",
        detail=detail,
    )


def _check_gpu(site: SiteConfig, partitions: Mapping[str, Dict[str, object]], runner) -> List[Check]:
    gpu_roles = {role: part for role, part in site.partitions.items() if part.has_gpu}
    if not gpu_roles:
        return []
    out: List[Check] = []
    for role, partition in sorted(gpu_roles.items()):
        row = partitions.get(partition.name) or {}
        models = parse_gres_models(str(row.get("gres", "")))
        if not models:
            out.append(Check(
                name=f"gpu model ({partition.name})",
                status=UNKNOWN,
                expected=f"{partition.gpus_per_node} x {partition.gpu_model}",
                observed=str(row.get("gres", "no GRES reported")),
                remedy="Check `sinfo -p <partition> -o %G` on a login node.",
            ))
        else:
            claimed = (partition.gpu_model or "").lower()
            matched = [name for name in models if claimed and claimed in name]
            out.append(Check(
                name=f"gpu model ({partition.name})",
                status=PASS if matched else FAIL,
                expected=f"{partition.gpus_per_node} x {partition.gpu_model}",
                observed=", ".join(f"{name or 'gpu'}:{count}" for name, count in sorted(models.items())),
                detail="" if matched else "The site file names a GPU model this partition does not advertise.",
                remedy="" if matched else f"Correct partitions.{role}.gpu_model.",
            ))
        # VRAM is the binding constraint for co-folding and SLURM does not
        # report it, so it can only be confirmed from a node that has a GPU.
        smi = runner(["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"])
        if smi is None:
            out.append(Check(
                name=f"gpu VRAM ({partition.name})",
                status=UNKNOWN,
                expected=f"{partition.gpu_vram_gb:g} GB" if partition.gpu_vram_gb else "stated VRAM",
                observed="nvidia-smi not present on this machine",
                detail=(
                    "VRAM decides whether the 19 candidates >= 600 residues and the 5 >= 1000 can be "
                    "co-folded whole or must be truncated, and truncation has to be pre-registered."
                ),
                remedy="Run this preflight inside a short GPU job: srun -p <gpu> --gres=gpu:1 --pty bash.",
            ))
            continue
        found = parse_nvidia_smi(smi)
        if not found:
            out.append(Check(
                name=f"gpu VRAM ({partition.name})",
                status=UNKNOWN,
                expected=f"{partition.gpu_vram_gb:g} GB" if partition.gpu_vram_gb else "stated VRAM",
                observed=smi.strip()[:120],
            ))
            continue
        smallest = min(item["vram_gb"] for item in found)  # type: ignore[type-var]
        claimed_vram = partition.gpu_vram_gb or 0.0
        out.append(Check(
            name=f"gpu VRAM ({partition.name})",
            status=PASS if smallest + 2 >= claimed_vram else FAIL,
            expected=f"{claimed_vram:g} GB",
            observed=f"{smallest:.0f} GB on {found[0]['name']}",
            detail="" if smallest + 2 >= claimed_vram else "Less VRAM than the site file claims.",
            remedy="" if smallest + 2 >= claimed_vram else f"Set partitions.*.gpu_vram_gb to {smallest:.0f}.",
        ))
    return out


#: Imports a study needs on the cluster, and what breaks without each.
REQUIRED_IMPORTS = (
    ("numpy", "every study"),
    ("pandas", "per-copy tables and bootstraps"),
    ("yaml", "reading these manifests on the cluster"),
)
OPTIONAL_IMPORTS = (
    ("vina", "docking (studies A, F, G, J)"),
    ("rdkit", "ligand construction and symmetry-corrected RMSD"),
    ("meeko", "flexible receptors (study J)"),
    ("openmm", "molecular dynamics"),
)
OPTIONAL_BINARIES = (
    ("fpocket", "pocket detection"),
    ("pdb2pqr", "protonation for docking"),
    ("mmseqs", "sequence homology groups"),
    ("foldseek", "structural homology groups"),
)


def _check_environment() -> List[Check]:
    out: List[Check] = []
    import importlib

    for module, why in REQUIRED_IMPORTS:
        try:
            importlib.import_module(module)
        except Exception as exc:  # noqa: BLE001 - any import failure is the same answer
            out.append(Check(
                name=f"import {module}",
                status=FAIL,
                expected="importable",
                observed=type(exc).__name__,
                detail=f"Needed for {why}.",
                remedy="Build the environment: bash hpc/rockfish/build_env.sh",
            ))
        else:
            out.append(Check(name=f"import {module}", status=PASS, expected="importable", observed="ok"))
    for module, why in OPTIONAL_IMPORTS:
        try:
            importlib.import_module(module)
        except Exception as exc:  # noqa: BLE001
            out.append(Check(
                name=f"import {module}",
                status=WARN,
                expected="importable",
                observed=type(exc).__name__,
                detail=f"Only {why} needs it.",
                remedy="bash hpc/rockfish/build_env.sh",
            ))
        else:
            out.append(Check(name=f"import {module}", status=PASS, expected="importable", observed="ok"))
    for binary, why in OPTIONAL_BINARIES:
        found = shutil.which(binary)
        out.append(Check(
            name=f"binary {binary}",
            status=PASS if found else WARN,
            expected="on PATH",
            observed=found or "not found",
            detail="" if found else f"Only {why} needs it.",
            remedy="" if found else "bash hpc/rockfish/build_env.sh",
        ))
    return out


def _check_manifest_fits(manifest: StudyManifest, site: SiteConfig) -> Check:
    problems = manifest.problems_against(site)
    return Check(
        name=f"study {manifest.study} fits the site",
        status=FAIL if problems else PASS,
        expected="every stage inside its partition's limits",
        observed=f"{len(problems)} problem(s)" if problems else "ok",
        detail="; ".join(problems),
        remedy="Lower the stage's resources, or move it to another partition." if problems else "",
    )


def markdown(report: Report) -> str:
    """The report as a table, for pasting into an issue or a provenance file."""
    counts = report.counts()
    headline = {
        "ok": "OK to submit",
        "blocked": "DO NOT SUBMIT",
        "indeterminate": "NOT VERIFIED - no scheduler on this machine, so nothing about the "
                         "cluster was actually checked",
    }[report.verdict]
    lines = [
        f"# Preflight: {report.site}",
        "",
        f"**{headline}** - "
        f"{counts[PASS]} passed, {counts[FAIL]} failed, {counts[WARN]} warnings, "
        f"{counts[UNKNOWN]} unknown.",
        "",
        "| check | status | expected | observed | what to do |",
        "|---|---|---|---|---|",
    ]
    for check in report.checks:
        remedy = check.remedy or check.detail or ""
        lines.append(
            f"| {check.name} | {check.status} | {check.expected or '-'} | "
            f"{check.observed or '-'} | {remedy or '-'} |"
        )
    if report.unknown:
        lines += [
            "",
            "`unknown` means this machine could not answer, not that the claim is wrong. "
            "A login node answers most of them; the GPU checks need a GPU job.",
        ]
    return "\n".join(lines) + "\n"
