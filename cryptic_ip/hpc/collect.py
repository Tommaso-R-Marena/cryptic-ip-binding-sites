"""Merge shard outputs, and refuse to pretend a missing shard never existed.

This module is the study J lesson in code. Two of its 272 docking copies died
when GitHub reclaimed their runners. Nothing was written for them, so they did
not appear among the 63 failures - they were simply absent, and the only thing
that revealed them was a merged record count of 270 against an expected 272.
The report at the time would have gone out on 207 copies either way.

So collection is a reconciliation, not a concatenation. Every shard is in
exactly one of five states, the merged file is only written when they add up,
and the census is published next to the result so a reader can check the
arithmetic:

    produced  an output file exists                      -> merged
    failed    the task recorded a non-zero exit          -> named
    killed    the scheduler took the task                -> named
    running   a record was written but never completed   -> named
    missing   no output and no record at all             -> named, and loudest
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from cryptic_ip.hpc.manifest import Stage, StudyManifest
from cryptic_ip.hpc.render import STATUS_DIR

PRODUCED = "produced"
FAILED = "failed"
KILLED = "killed"
RUNNING = "running"
MISSING = "missing"

#: States that mean the shard contributed nothing and we know why.
ACCOUNTED = (FAILED, KILLED, RUNNING)


class CollectError(RuntimeError):
    """A collection that would have published an unreconciled result."""


@dataclass
class Census:
    """Which shards landed, which did not, and what the scheduler said."""

    study: str
    stage: str
    expected: int
    states: Dict[int, str] = field(default_factory=dict)
    status_records: Dict[int, Mapping[str, object]] = field(default_factory=dict)
    outputs: Dict[int, str] = field(default_factory=dict)
    records_merged: int = 0

    def shards_in(self, state: str) -> List[int]:
        return sorted(shard for shard, value in self.states.items() if value == state)

    @property
    def produced(self) -> List[int]:
        return self.shards_in(PRODUCED)

    @property
    def missing(self) -> List[int]:
        return self.shards_in(MISSING)

    @property
    def reconciles(self) -> bool:
        """Whether every expected shard is in exactly one known state."""
        return len(self.states) == self.expected and not self.missing

    def counts(self) -> Dict[str, int]:
        out = {PRODUCED: 0, FAILED: 0, KILLED: 0, RUNNING: 0, MISSING: 0}
        for state in self.states.values():
            out[state] = out.get(state, 0) + 1
        return out

    def as_dict(self) -> Dict[str, object]:
        return {
            "study": self.study,
            "stage": self.stage,
            "expected_shards": self.expected,
            "counts": self.counts(),
            "records_merged": self.records_merged,
            "reconciles": self.reconciles,
            "produced": self.produced,
            "failed": self.shards_in(FAILED),
            "killed": self.shards_in(KILLED),
            "running": self.shards_in(RUNNING),
            "missing": self.missing,
            "status_records": {str(k): v for k, v in sorted(self.status_records.items())},
        }

    def summary(self) -> str:
        counts = self.counts()
        parts = [f"{counts[PRODUCED]} of {self.expected} shards produced output"]
        for state in (FAILED, KILLED, RUNNING, MISSING):
            if counts[state]:
                parts.append(f"{counts[state]} {state}")
        line = ", ".join(parts)
        if self.missing:
            shown = ", ".join(str(s) for s in self.missing[:10])
            more = f" and {len(self.missing) - 10} more" if len(self.missing) > 10 else ""
            line += (
                f"\nMissing shards wrote neither output nor a status record: {shown}{more}. "
                "A shard in this state is not censoring by the protocol - it is a shard that "
                "never ran, or whose node was taken, and it must not shrink the denominator quietly."
            )
        return line


def read_status(stage_dir: Path, shard: int) -> Optional[Mapping[str, object]]:
    """The status sentinel a task wrote for ``shard``, if any."""
    path = Path(stage_dir) / STATUS_DIR / f"shard-{shard}.json"
    if not path.is_file():
        return None
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {"status": "unreadable", "shard": shard}
    return value if isinstance(value, Mapping) else {"status": "unreadable", "shard": shard}


def take_census(manifest: StudyManifest, stage: Stage, stage_dir: Path) -> Census:
    """Classify every shard of ``stage`` by what it left behind."""
    stage_dir = Path(stage_dir)
    census = Census(study=manifest.study, stage=stage.name, expected=stage.task_count)
    for shard in range(stage.task_count):
        status = read_status(stage_dir, shard)
        if status is not None:
            census.status_records[shard] = status
        output_name = stage.output_for(shard) if stage.output else None
        output_path = stage_dir / output_name if output_name else None
        if output_path is not None and output_path.is_file():
            census.states[shard] = PRODUCED
            census.outputs[shard] = str(output_path)
            continue
        recorded = str((status or {}).get("status") or "")
        if recorded == "killed":
            census.states[shard] = KILLED
        elif recorded == "running":
            census.states[shard] = RUNNING
        elif recorded in {"failed", "unreadable"}:
            census.states[shard] = FAILED
        elif recorded == "ok":
            # The task said it succeeded and there is no output: that is a bug
            # in the stage, not an absent shard, and it is worth saying so.
            census.states[shard] = FAILED
        else:
            census.states[shard] = MISSING
    return census


def merge_records(census: Census) -> List[object]:
    """Concatenate the JSON lists the produced shards wrote.

    A shard whose file is not a list is kept as one record rather than
    discarded, because a stage may legitimately write one object per shard.
    """
    records: List[object] = []
    for shard in census.produced:
        path = Path(census.outputs[shard])
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise CollectError(f"shard {shard} wrote {path}, which is not readable JSON: {exc}") from exc
        if isinstance(payload, list):
            records.extend(payload)
        else:
            records.append(payload)
    census.records_merged = len(records)
    return records


def collect_stage(
    manifest: StudyManifest,
    stage_name: str,
    work_dir: Path,
    *,
    records_out: Optional[Path] = None,
    census_out: Optional[Path] = None,
    allow_missing: int = 0,
) -> Census:
    """Reconcile and merge one stage's shards.

    ``allow_missing`` is the number of shards permitted to have left nothing at
    all. It defaults to zero and has to be raised deliberately, because the
    whole point of this function is that a vanished shard cannot pass unnoticed.
    """
    stage = manifest.stage(stage_name)
    if not stage.is_array:
        raise CollectError(f"stage {stage_name!r} is not sharded; there is nothing to reconcile")
    if stage.output is None:
        raise CollectError(f"stage {stage_name!r} declares no output, so its shards cannot be collected")
    stage_dir = Path(work_dir) / stage.name
    census = take_census(manifest, stage, stage_dir)
    records = merge_records(census)

    if census_out is not None:
        Path(census_out).parent.mkdir(parents=True, exist_ok=True)
        Path(census_out).write_text(json.dumps(census.as_dict(), indent=2, sort_keys=True) + "\n")

    missing = len(census.missing)
    if missing > allow_missing:
        raise CollectError(
            f"{missing} shard(s) of stage {stage_name!r} left neither output nor a status record"
            f"{f' (allowing {allow_missing})' if allow_missing else ''}.\n"
            f"{census.summary()}\n"
            "Re-run them with `rockfish.py resume`, or, if they are genuinely unrecoverable, "
            "pass --allow-missing with that count so the number is a stated decision rather than "
            "a silently shorter dataset."
        )
    if not census.produced:
        raise CollectError(
            f"no shard of stage {stage_name!r} produced an output file.\n{census.summary()}\n"
            "A stage where every task ran and every copy failed still writes files, so this is a "
            "structural failure - the array did not run - rather than a degenerate result."
        )
    if records_out is not None:
        Path(records_out).parent.mkdir(parents=True, exist_ok=True)
        Path(records_out).write_text(json.dumps(records))
    return census


def provenance(
    manifest: StudyManifest,
    census: Census,
    *,
    site_name: str,
    job_ids: Sequence[str] = (),
    commit: Optional[str] = None,
) -> str:
    """A provenance stub naming the run, the census and the plan."""
    counts = census.counts()
    lines = [
        f"# Study {manifest.study}: {manifest.title}",
        "",
        f"**Run.** Site `{site_name}`, stage `{census.stage}`"
        + (f", SLURM job(s) {', '.join(job_ids)}" if job_ids else "")
        + (f", commit `{commit}`" if commit else "")
        + ".",
        "",
        "**Plan.** " + ", ".join(f"`{doc}`" for doc in manifest.plan) + ".",
        "",
        "## Shard census",
        "",
        "| state | shards | what it means |",
        "|---|---|---|",
        f"| produced | {counts[PRODUCED]} | an output file exists and was merged |",
        f"| failed | {counts[FAILED]} | the task recorded a non-zero exit |",
        f"| killed | {counts[KILLED]} | the scheduler took the task away |",
        f"| running | {counts[RUNNING]} | a record was written but never completed |",
        f"| missing | {counts[MISSING]} | nothing at all was written |",
        "",
        f"{census.expected} shards expected, {census.records_merged} records merged.",
        "",
    ]
    if census.missing:
        lines += [
            "The missing shards are named in `census.json`. They are not censoring by the "
            "protocol and must not be folded into any denominator.",
            "",
        ]
    return "\n".join(lines)
