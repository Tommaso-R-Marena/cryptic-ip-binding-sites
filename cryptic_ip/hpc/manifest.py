"""What a study is, as data: stages, resources, dependencies, and its plan.

A manifest is the cluster counterpart of a workflow file in
`.github/workflows/`. It says nothing about science - the science stays in
`scripts/<study>.py`, so a study run here and a study run on shared runners
execute the same code with the same arguments.

One rule is enforced here rather than left to discipline: a manifest must name
the pre-registered plan documents its study depends on, and
`missing_plan_docs` reports the ones that do not exist. The CLI refuses to
submit when any are missing, so the infrastructure cannot be used to run a
study that has not been pre-registered.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

import yaml

from cryptic_ip.hpc.site import Partition, SiteConfig, format_walltime, parse_walltime

#: Placeholders a stage command may use. Anything else is a typo, and a typo in
#: a command line is a job that burns its walltime producing nothing.
PLACEHOLDERS = frozenset({"shard", "shards", "stage_dir", "work_dir", "study", "repo"})

#: The one definition of a placeholder, shared with `cryptic_ip.hpc.render`.
PLACEHOLDER_PATTERN = re.compile(r"\{([a-z_]+)\}")
_PLACEHOLDER = PLACEHOLDER_PATTERN

_STAGE_NAME = re.compile(r"^[a-z][a-z0-9_]*$")

#: How a stage waits on the stages it depends on. `ok` is the default; `any`
#: is what a reporting stage wants, because the report matters most when some
#: shards failed - the same reason study J's report job used `always()` rather
#: than `!cancelled()`.
AFTER_KINDS = frozenset({"ok", "any", "notok"})


class ManifestError(ValueError):
    """A manifest that cannot be used as written."""


@dataclass(frozen=True)
class Resources:
    """What one task of a stage asks the scheduler for."""

    partition: str
    cpus: int = 1
    mem_gb: float = 4.0
    walltime: int = 3600
    gpus: int = 0

    @classmethod
    def from_dict(cls, raw: Mapping[str, object]) -> "Resources":
        unknown = set(raw) - {"partition", "cpus", "mem_gb", "walltime", "gpus"}
        if unknown:
            raise ManifestError(f"resources have unknown keys {sorted(unknown)}")
        if "partition" not in raw:
            raise ManifestError("resources must name a partition role")
        # Explicit zeros must be rejected, not replaced by the default.
        cpus = 1 if raw.get("cpus") is None else int(raw["cpus"])
        gpus = 0 if raw.get("gpus") is None else int(raw["gpus"])
        mem = 4.0 if raw.get("mem_gb") is None else float(raw["mem_gb"])
        if cpus < 1:
            raise ManifestError(f"cpus must be at least 1, got {cpus}")
        if gpus < 0:
            raise ManifestError(f"gpus cannot be negative, got {gpus}")
        if mem <= 0:
            raise ManifestError(f"mem_gb must be positive, got {mem}")
        walltime = parse_walltime(str(raw.get("walltime") or "01:00:00"))
        if walltime <= 0:
            raise ManifestError("walltime must be a positive duration")
        return cls(partition=str(raw["partition"]), cpus=cpus, mem_gb=mem, walltime=walltime, gpus=gpus)

    @property
    def walltime_text(self) -> str:
        return format_walltime(self.walltime)

    def problems_against(self, partition: Partition) -> List[str]:
        """Why this request cannot run in ``partition``, as plain sentences."""
        out: List[str] = []
        if partition.max_walltime != -1 and self.walltime > partition.max_walltime:
            out.append(
                f"asks {self.walltime_text} but partition {partition.name!r} caps at "
                f"{format_walltime(partition.max_walltime)}"
            )
        if partition.max_cpus_per_node is not None and self.cpus > partition.max_cpus_per_node:
            out.append(
                f"asks {self.cpus} cpus but partition {partition.name!r} has "
                f"{partition.max_cpus_per_node} per node"
            )
        if partition.max_mem_gb_per_node is not None and self.mem_gb > partition.max_mem_gb_per_node:
            out.append(
                f"asks {self.mem_gb:g} GB but partition {partition.name!r} has "
                f"{partition.max_mem_gb_per_node:g} GB per node"
            )
        if self.gpus and not partition.has_gpu:
            out.append(f"asks {self.gpus} gpu(s) but partition {partition.name!r} has none")
        if self.gpus and partition.gpus_per_node and self.gpus > partition.gpus_per_node:
            out.append(
                f"asks {self.gpus} gpus but partition {partition.name!r} has "
                f"{partition.gpus_per_node} per node"
            )
        return out


@dataclass(frozen=True)
class Stage:
    """One step of a study: a command, run once or once per shard."""

    name: str
    command: Sequence[str]
    resources: Resources
    description: str = ""
    shards: Optional[int] = None
    output: Optional[str] = None
    depends_on: Sequence[str] = field(default_factory=tuple)
    after: str = "ok"
    throttle: Optional[int] = None
    env: Mapping[str, str] = field(default_factory=dict)

    @property
    def is_array(self) -> bool:
        return self.shards is not None

    @property
    def task_count(self) -> int:
        return self.shards if self.shards is not None else 1

    def output_for(self, shard: Optional[int]) -> Optional[str]:
        """The output path this stage promises for ``shard``."""
        if self.output is None:
            return None
        if self.is_array:
            if shard is None:
                raise ManifestError(f"stage {self.name!r} is sharded; name a shard")
            return self.output.format(shard=shard)
        return self.output

    @classmethod
    def from_dict(cls, raw: Mapping[str, object]) -> "Stage":
        unknown = set(raw) - {
            "name", "command", "resources", "description", "shards", "output",
            "depends_on", "after", "throttle", "env", "note",
        }
        if unknown:
            raise ManifestError(f"stage has unknown keys {sorted(unknown)}")
        for required in ("name", "command", "resources"):
            if required not in raw:
                raise ManifestError(f"stage is missing {required!r}")
        name = str(raw["name"])
        if not _STAGE_NAME.match(name):
            raise ManifestError(f"stage name {name!r} must be lowercase letters, digits and underscores")
        command = raw["command"]
        if isinstance(command, str) or not isinstance(command, Sequence):
            raise ManifestError(f"stage {name!r}: command must be a list of tokens, not one string")
        tokens = [str(token) for token in command]
        if not tokens:
            raise ManifestError(f"stage {name!r}: command is empty")
        resources_raw = raw["resources"]
        if not isinstance(resources_raw, Mapping):
            raise ManifestError(f"stage {name!r}: resources must be a mapping")
        shards = raw.get("shards")
        shard_count = None if shards is None else int(shards)  # type: ignore[arg-type]
        if shard_count is not None and shard_count < 1:
            raise ManifestError(f"stage {name!r}: shards must be at least 1, got {shard_count}")
        throttle = raw.get("throttle")
        throttle_count = None if throttle is None else int(throttle)  # type: ignore[arg-type]
        if throttle_count is not None and throttle_count < 1:
            raise ManifestError(f"stage {name!r}: throttle must be at least 1, got {throttle_count}")
        output = raw.get("output")
        output_text = None if output is None else str(output)
        env_raw = raw.get("env") or {}
        if not isinstance(env_raw, Mapping):
            raise ManifestError(f"stage {name!r}: env must be a mapping")
        after = str(raw.get("after") or "ok")
        if after not in AFTER_KINDS:
            raise ManifestError(
                f"stage {name!r}: after must be one of {sorted(AFTER_KINDS)}, got {after!r}"
            )
        if after != "ok" and not raw.get("depends_on"):
            raise ManifestError(f"stage {name!r}: after={after!r} means nothing without depends_on")

        stage = cls(
            name=name,
            command=tuple(tokens),
            resources=Resources.from_dict(resources_raw),
            description=str(raw.get("description") or ""),
            shards=shard_count,
            output=output_text,
            depends_on=tuple(str(dep) for dep in (raw.get("depends_on") or ())),
            after=after,
            throttle=throttle_count,
            env={str(k): str(v) for k, v in env_raw.items()},
        )
        stage.validate_placeholders()
        return stage

    def validate_placeholders(self) -> None:
        """Reject unknown placeholders, and `{shard}` where there are no shards."""
        used: Set[str] = set()
        for token in self.command:
            used.update(_PLACEHOLDER.findall(token))
        if self.output:
            used.update(_PLACEHOLDER.findall(self.output))
        unknown = used - PLACEHOLDERS
        if unknown:
            raise ManifestError(
                f"stage {self.name!r} uses unknown placeholder(s) {sorted(unknown)}; "
                f"known: {sorted(PLACEHOLDERS)}"
            )
        if not self.is_array and "shard" in used:
            raise ManifestError(f"stage {self.name!r} uses {{shard}} but declares no shards")
        if self.is_array and self.output and "shard" not in _PLACEHOLDER.findall(self.output):
            raise ManifestError(
                f"stage {self.name!r} is sharded but its output path does not vary by shard, "
                "so every task would write over the others"
            )
        if self.output is not None:
            # An output is stated relative to the stage directory so `collect`
            # can find it without re-deriving anyone's scratch layout.
            if self.output.startswith("/") or self.output.startswith("${"):
                raise ManifestError(
                    f"stage {self.name!r}: output must be relative to the stage directory, got {self.output!r}"
                )
            if ".." in Path(self.output).parts:
                raise ManifestError(f"stage {self.name!r}: output may not climb out of the stage directory")
            stray = set(_PLACEHOLDER.findall(self.output)) - {"shard"}
            if stray:
                raise ManifestError(
                    f"stage {self.name!r}: output may only use {{shard}}, not {sorted(stray)}"
                )


@dataclass(frozen=True)
class StudyManifest:
    """A study as the cluster sees it."""

    study: str
    title: str
    stages: Sequence[Stage]
    plan: Sequence[str] = field(default_factory=tuple)
    description: str = ""
    source_path: Optional[Path] = None

    def stage(self, name: str) -> Stage:
        for stage in self.stages:
            if stage.name == name:
                return stage
        known = ", ".join(s.name for s in self.stages)
        raise ManifestError(f"study {self.study!r} has no stage {name!r} (has: {known})")

    def order(self) -> List[Stage]:
        """Stages in an order where every dependency comes first."""
        return _topological(self.stages)

    def missing_plan_docs(self, repo_root: Path) -> List[str]:
        """Pre-registration documents this study names that do not exist."""
        root = Path(repo_root)
        return [doc for doc in self.plan if not (root / doc).is_file()]

    def problems_against(self, site: SiteConfig) -> List[str]:
        """Why this study cannot run on ``site``, as plain sentences."""
        out: List[str] = []
        for stage in self.stages:
            try:
                partition = site.partition(stage.resources.partition)
            except Exception as exc:  # SiteError, kept loose so one bad stage lists the rest
                out.append(f"stage {stage.name!r}: {exc}")
                continue
            out.extend(f"stage {stage.name!r} {problem}" for problem in stage.resources.problems_against(partition))
        return out

    @classmethod
    def from_dict(cls, raw: Mapping[str, object], *, source_path: Optional[Path] = None) -> "StudyManifest":
        unknown = set(raw) - {"study", "title", "stages", "plan", "description", "note"}
        if unknown:
            raise ManifestError(f"manifest has unknown keys {sorted(unknown)}")
        for required in ("study", "title", "stages"):
            if required not in raw:
                raise ManifestError(f"manifest is missing {required!r}")
        stages_raw = raw["stages"]
        if isinstance(stages_raw, str) or not isinstance(stages_raw, Sequence) or not stages_raw:
            raise ManifestError("manifest needs a non-empty list of stages")
        stages = [Stage.from_dict(item) for item in stages_raw if isinstance(item, Mapping)]
        if len(stages) != len(stages_raw):
            raise ManifestError("every stage must be a mapping")
        names = [stage.name for stage in stages]
        duplicates = sorted({name for name in names if names.count(name) > 1})
        if duplicates:
            raise ManifestError(f"duplicate stage name(s) {duplicates}")
        known = set(names)
        for stage in stages:
            for dep in stage.depends_on:
                if dep not in known:
                    raise ManifestError(f"stage {stage.name!r} depends on unknown stage {dep!r}")
                if dep == stage.name:
                    raise ManifestError(f"stage {stage.name!r} depends on itself")
        plan = raw.get("plan") or ()
        if isinstance(plan, str) or not isinstance(plan, Sequence):
            raise ManifestError("plan must be a list of document paths")
        if not plan:
            raise ManifestError(
                "manifest must name the pre-registered plan document(s) for its study; "
                "an empty plan list is how an un-pre-registered study gets run by accident"
            )
        manifest = cls(
            study=str(raw["study"]),
            title=str(raw["title"]),
            stages=tuple(stages),
            plan=tuple(str(doc) for doc in plan),
            description=str(raw.get("description") or ""),
            source_path=source_path,
        )
        manifest.order()  # raises on a dependency cycle
        return manifest


def _topological(stages: Sequence[Stage]) -> List[Stage]:
    by_name = {stage.name: stage for stage in stages}
    remaining: Dict[str, Set[str]] = {s.name: set(s.depends_on) for s in stages}
    out: List[Stage] = []
    while remaining:
        ready = sorted(name for name, deps in remaining.items() if not deps)
        if not ready:
            stuck = ", ".join(sorted(remaining))
            raise ManifestError(f"stage dependencies form a cycle among: {stuck}")
        for name in ready:
            out.append(by_name[name])
            del remaining[name]
        for deps in remaining.values():
            deps.difference_update(ready)
    return out


def load_manifest(path: Path) -> StudyManifest:
    """Read a study manifest YAML file."""
    path = Path(path)
    try:
        raw = yaml.safe_load(path.read_text()) or {}
    except yaml.YAMLError as exc:
        raise ManifestError(f"{path}: {exc}") from exc
    if not isinstance(raw, Mapping):
        raise ManifestError(f"{path}: expected a mapping at the top level")
    return StudyManifest.from_dict(raw, source_path=path)


def load_manifests(root: Path) -> List[StudyManifest]:
    """Every manifest under ``root``, sorted by file name."""
    return [load_manifest(path) for path in sorted(Path(root).glob("*.yml"))]


def total_tasks(manifest: StudyManifest) -> int:
    """How many tasks the whole study submits, array tasks included."""
    return sum(stage.task_count for stage in manifest.stages)


def iter_shards(stage: Stage) -> Iterable[Tuple[int, Optional[str]]]:
    """``(shard, output)`` for every task of ``stage``."""
    if not stage.is_array:
        yield 0, stage.output
        return
    for shard in range(stage.task_count):
        yield shard, stage.output_for(shard)
