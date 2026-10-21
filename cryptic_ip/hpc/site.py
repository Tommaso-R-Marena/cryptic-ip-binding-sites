"""What a cluster is, as data rather than as strings scattered through scripts.

Every value here is a claim about a machine this code cannot see. A wrong claim
is not a crash - it is a job that queues forever against a partition that does
not exist, or an array that is silently truncated. So each site file carries a
`verified` block saying when its claims were last checked against the live
cluster, and `cryptic_ip.hpc.preflight` is what checks them.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

import yaml

#: SLURM accepts several walltime spellings; these are the ones sbatch documents.
_WALLTIME = re.compile(
    r"^(?:(?P<days>\d+)-)?(?P<a>\d+)(?::(?P<b>\d+))?(?::(?P<c>\d+))?$"
)

UNLIMITED = -1


class SiteError(ValueError):
    """A site file that cannot be used as written."""


def parse_walltime(text: str) -> int:
    """Seconds in a SLURM time string, or ``UNLIMITED``.

    SLURM's own grammar is ambiguous without the day part: ``30:00`` is thirty
    minutes, ``30:00:00`` is thirty hours, and ``30`` is thirty minutes. We
    follow sbatch: with days present the next fields are H, M, S; without days
    the fields are M, S then H, M, S.
    """
    value = str(text).strip()
    if value.upper() in {"INFINITE", "UNLIMITED", "NONE"}:
        return UNLIMITED
    match = _WALLTIME.match(value)
    if not match:
        raise SiteError(f"not a SLURM time: {text!r}")
    days = int(match.group("days") or 0)
    parts = [match.group(k) for k in ("a", "b", "c")]
    present = [int(p) for p in parts if p is not None]
    if match.group("days") is not None:
        # days-H, days-H:M, days-H:M:S
        hours, minutes, seconds = (present + [0, 0, 0])[:3]
    elif len(present) == 1:
        hours, minutes, seconds = 0, present[0], 0
    elif len(present) == 2:
        hours, minutes, seconds = 0, present[0], present[1]
    else:
        hours, minutes, seconds = present[0], present[1], present[2]
    return days * 86400 + hours * 3600 + minutes * 60 + seconds


def format_walltime(seconds: int) -> str:
    """A SLURM time string for ``seconds``, in the unambiguous D-HH:MM:SS form."""
    if seconds == UNLIMITED:
        return "UNLIMITED"
    if seconds < 0:
        raise SiteError(f"negative walltime: {seconds}")
    days, rest = divmod(int(seconds), 86400)
    hours, rest = divmod(rest, 3600)
    minutes, secs = divmod(rest, 60)
    if days:
        return f"{days}-{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def _expand(value: str) -> str:
    """Expand ``$USER`` and ``~`` without requiring the paths to exist.

    An unset variable is an error rather than a literal: `os.path.expandvars`
    leaves `$USER` alone when USER is not in the environment, which would send
    a study's output into a directory actually named `$USER`.
    """
    expanded = os.path.expanduser(os.path.expandvars(str(value)))
    if "$" in expanded:
        raise SiteError(
            f"{value!r} references a variable that is not set in this environment "
            f"(expanded to {expanded!r}); set it, or write the path out in full"
        )
    return expanded


@dataclass(frozen=True)
class Partition:
    """One SLURM partition, with the limits a job must fit inside.

    ``max_mem_gb_per_node`` is **GiB**, the same unit as sbatch's `--mem=<n>G`
    and as `sinfo`'s megabytes divided by 1024. A node sold as "192 GB" reports
    192000 MB, which is 187.5 GiB, so writing 192 here asks for more memory
    than the node can give and the job waits forever. Preflight compares
    against what sinfo reports, in GiB.
    """

    role: str
    name: str
    max_walltime: int
    max_cpus_per_node: Optional[int] = None
    max_mem_gb_per_node: Optional[float] = None
    gpus_per_node: int = 0
    gpu_model: Optional[str] = None
    gpu_vram_gb: Optional[float] = None
    gres_name: str = "gpu"

    @property
    def has_gpu(self) -> bool:
        return self.gpus_per_node > 0

    @classmethod
    def from_dict(cls, role: str, raw: Mapping[str, object]) -> "Partition":
        missing = {"name", "max_walltime"} - set(raw)
        if missing:
            raise SiteError(f"partition {role!r} is missing {sorted(missing)}")
        unknown = set(raw) - {
            "name", "max_walltime", "max_cpus_per_node", "max_mem_gb_per_node",
            "gpus_per_node", "gpu_model", "gpu_vram_gb", "gres_name", "note",
        }
        if unknown:
            raise SiteError(f"partition {role!r} has unknown keys {sorted(unknown)}")
        gpus = int(raw.get("gpus_per_node") or 0)
        vram = raw.get("gpu_vram_gb")
        if gpus and vram is None:
            raise SiteError(
                f"partition {role!r} offers GPUs but states no gpu_vram_gb; "
                "VRAM is the binding constraint for co-folding, so it may not be left blank"
            )
        return cls(
            role=role,
            name=str(raw["name"]),
            max_walltime=parse_walltime(str(raw["max_walltime"])),
            max_cpus_per_node=_opt_int(raw.get("max_cpus_per_node")),
            max_mem_gb_per_node=_opt_float(raw.get("max_mem_gb_per_node")),
            gpus_per_node=gpus,
            gpu_model=_opt_str(raw.get("gpu_model")),
            gpu_vram_gb=_opt_float(vram),
            gres_name=str(raw.get("gres_name") or "gpu"),
        )


def _opt_int(value: object) -> Optional[int]:
    return None if value is None else int(value)  # type: ignore[arg-type]


def _opt_float(value: object) -> Optional[float]:
    return None if value is None else float(value)  # type: ignore[arg-type]


def _opt_str(value: object) -> Optional[str]:
    return None if value is None else str(value)


@dataclass(frozen=True)
class SiteConfig:
    """A cluster, its limits, and how to get a working environment on it."""

    name: str
    scheduler: str
    partitions: Dict[str, Partition]
    scratch_root: str
    account: Optional[str] = None
    project_root: Optional[str] = None
    env_activate: Sequence[str] = field(default_factory=tuple)
    module_loads: Sequence[str] = field(default_factory=tuple)
    max_array_size: int = 1001
    max_array_throttle: int = 50
    max_submitted_jobs: Optional[int] = None
    verified_on: Optional[str] = None
    verified_by: Optional[str] = None
    source_path: Optional[Path] = None

    # ----------------------------------------------------------------- lookups
    def partition(self, role: str) -> Partition:
        try:
            return self.partitions[role]
        except KeyError:
            known = ", ".join(sorted(self.partitions)) or "none"
            raise SiteError(f"site {self.name!r} has no {role!r} partition (has: {known})") from None

    @property
    def is_verified(self) -> bool:
        """Whether anyone has checked these claims against the live cluster."""
        return bool(self.verified_on)

    def scratch_for(self, study: str) -> str:
        """Where a study's working files live. Never the home filesystem."""
        return str(Path(_expand(self.scratch_root)) / study)

    @classmethod
    def from_dict(cls, raw: Mapping[str, object], *, source_path: Optional[Path] = None) -> "SiteConfig":
        unknown = set(raw) - {
            "name", "scheduler", "partitions", "scratch_root", "account", "project_root",
            "env_activate", "module_loads", "max_array_size", "max_array_throttle",
            "max_submitted_jobs", "verified", "note",
        }
        if unknown:
            raise SiteError(f"site file has unknown keys {sorted(unknown)}")
        for required in ("name", "partitions", "scratch_root"):
            if required not in raw:
                raise SiteError(f"site file is missing {required!r}")
        partitions_raw = raw["partitions"]
        if not isinstance(partitions_raw, Mapping) or not partitions_raw:
            raise SiteError("site file needs at least one partition")
        partitions = {
            str(role): Partition.from_dict(str(role), part)
            for role, part in partitions_raw.items()
            if isinstance(part, Mapping)
        }
        if len(partitions) != len(partitions_raw):
            raise SiteError("every partition must be a mapping")
        verified = raw.get("verified") or {}
        if not isinstance(verified, Mapping):
            raise SiteError("'verified' must be a mapping with 'date' and 'by'")
        # YAML 1.1 reads the bare key `on` as the boolean true, so a
        # `verified: {on: 2026-10-10}` block parses as `{True: ...}` and the
        # date disappears, leaving the site looking unverified. Reject the
        # spelling rather than quietly mis-read it.
        if True in verified or "on" in verified:
            raise SiteError(
                "use 'verified.date', not 'verified.on': YAML reads the key 'on' as the boolean true, "
                "so the date would be dropped and the site would look unverified"
            )
        # `or` would turn an explicit 0 into the default and hide the mistake.
        max_array = 1001 if raw.get("max_array_size") is None else int(raw["max_array_size"])
        if max_array < 1:
            raise SiteError(f"max_array_size must be positive, got {max_array}")
        throttle = 50 if raw.get("max_array_throttle") is None else int(raw["max_array_throttle"])
        if throttle < 1:
            raise SiteError(f"max_array_throttle must be positive, got {throttle}")
        return cls(
            name=str(raw["name"]),
            scheduler=str(raw.get("scheduler") or "slurm"),
            partitions=partitions,
            scratch_root=str(raw["scratch_root"]),
            account=_opt_str(raw.get("account")),
            project_root=_opt_str(raw.get("project_root")),
            env_activate=tuple(str(line) for line in (raw.get("env_activate") or ())),
            module_loads=tuple(str(line) for line in (raw.get("module_loads") or ())),
            max_array_size=max_array,
            max_array_throttle=throttle,
            max_submitted_jobs=_opt_int(raw.get("max_submitted_jobs")),
            verified_on=_opt_str(verified.get("date")),
            verified_by=_opt_str(verified.get("by")),
            source_path=source_path,
        )


def load_site(path: Path) -> SiteConfig:
    """Read a site YAML file."""
    path = Path(path)
    try:
        raw = yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError as exc:
        raise SiteError(f"{path}: {exc}") from exc
    if not isinstance(raw, Mapping):
        raise SiteError(f"{path}: expected a mapping at the top level")
    return SiteConfig.from_dict(raw, source_path=path)


def available_sites(root: Path) -> List[Path]:
    """Site files under ``root``, sorted."""
    return sorted(Path(root).glob("*.yml"))
