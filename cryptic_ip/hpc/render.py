"""Turn a stage into an sbatch script.

Two things here are not cosmetic.

**The shard map.** An array task learns which shard it is from a bash array
written into the script, not from arithmetic on `$SLURM_ARRAY_TASK_ID`. That
keeps every index below the cluster's `MaxArraySize` and makes a resume
submission - an arbitrary, sparse set of shards - render exactly like a first
submission.

**The status sentinel.** Each task writes a record before it starts and
rewrites it when it ends, through a bash trap that fires on normal exit and on
the signals SLURM uses to kill a job. Study J lost two copies to reclaimed
runners: no record was written, so the copies did not appear among the failures,
they simply were not in the census at all, and only a count of 270 against an
expected 272 revealed it. A task that dies here leaves a file saying so.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from cryptic_ip.hpc.manifest import PLACEHOLDER_PATTERN, Stage, StudyManifest
from cryptic_ip.hpc.shards import ArrayChunk, chunk_shards
from cryptic_ip.hpc.site import SiteConfig, format_walltime

#: Signals a scheduler sends when it takes a node or a job away. Trapping them
#: is what turns "the task vanished" into "the task recorded that it was killed".
TRAPPED_SIGNALS = ("TERM", "INT", "HUP", "QUIT", "USR1", "XCPU")

STATUS_DIR = "_status"


@dataclass(frozen=True)
class RenderedStage:
    """One sbatch script, with what the submitter needs to know about it."""

    stage: str
    chunk_index: int
    script: str
    shards: Sequence[int]
    depends_on: Sequence[str]
    after: str
    script_name: str


def _escape_literal(text: str) -> str:
    """Escape a literal for inclusion inside a bash double-quoted word."""
    out = []
    for char in text:
        if char in '\\"$`':
            out.append("\\" + char)
        else:
            out.append(char)
    return "".join(out)


def expand_token(token: str, mapping: Mapping[str, str]) -> str:
    """Render one command token as a single bash double-quoted word.

    Literal text is escaped; a `{placeholder}` becomes the mapped text, which is
    either a shell variable reference such as `${SHARD}` or an escaped literal.

    A brace that is not a placeholder - the `{"shard": 3}` inside a JSON literal
    in a command, say - is literal text. The pattern is the same
    `PLACEHOLDER_PATTERN` the manifest validates against, so a command that
    loads cannot then fail to render. A hand-rolled scanner here disagreed with
    the validator, and the disagreement only showed up at submit time.
    """
    out: List[str] = []
    cursor = 0
    for match in PLACEHOLDER_PATTERN.finditer(token):
        name = match.group(1)
        if name not in mapping:
            raise KeyError(f"no value for placeholder {{{name}}}")
        out.append(_escape_literal(token[cursor:match.start()]))
        value = mapping[name]
        out.append(value if value.startswith("${") else _escape_literal(value))
        cursor = match.end()
    out.append(_escape_literal(token[cursor:]))
    return chr(34) + "".join(out) + chr(34)


def _directives(
    manifest: StudyManifest,
    stage: Stage,
    site: SiteConfig,
    chunk: Optional[ArrayChunk],
    log_dir: str,
) -> List[str]:
    partition = site.partition(stage.resources.partition)
    name = f"{manifest.study.lower()}-{stage.name}"
    if chunk is not None and chunk.chunk_index:
        name = f"{name}-c{chunk.chunk_index}"
    lines = [
        f"#SBATCH --job-name={name}",
        f"#SBATCH --partition={partition.name}",
        f"#SBATCH --time={format_walltime(stage.resources.walltime)}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={stage.resources.cpus}",
        f"#SBATCH --mem={int(round(stage.resources.mem_gb))}G",
    ]
    if site.account:
        lines.insert(1, f"#SBATCH --account={site.account}")
    if stage.resources.gpus:
        lines.append(f"#SBATCH --gres={partition.gres_name}:{stage.resources.gpus}")
    if chunk is not None:
        lines.append(f"#SBATCH --array={chunk.array_spec}")
        lines.append(f"#SBATCH --output={log_dir}/{stage.name}-%A_%a.out")
    else:
        lines.append(f"#SBATCH --output={log_dir}/{stage.name}-%j.out")
    # A requeued task re-runs its own shard, and the sentinel records the retry
    # rather than leaving two half-written outputs behind.
    lines.append("#SBATCH --requeue")
    # SLURM sends USR1 before the walltime kill, which is what lets the trap
    # write a record while the task is still alive.
    lines.append("#SBATCH --signal=B:USR1@120")
    return lines


def render_stage(
    manifest: StudyManifest,
    stage: Stage,
    site: SiteConfig,
    *,
    work_dir: str,
    repo_root: str,
    shards: Optional[Sequence[int]] = None,
    chunk: Optional[ArrayChunk] = None,
) -> str:
    """The sbatch script for one stage, or for one chunk of a sharded stage."""
    if stage.is_array and chunk is None:
        raise ValueError(f"stage {stage.name!r} is sharded; render one chunk at a time")
    if not stage.is_array and chunk is not None:
        raise ValueError(f"stage {stage.name!r} is not sharded; it takes no chunk")
    del shards  # the chunk carries them; kept in the signature for callers' clarity

    stage_dir = f"{work_dir}/{stage.name}"
    log_dir = f"{work_dir}/logs"
    mapping: Dict[str, str] = {
        "shard": "${SHARD}",
        "shards": str(stage.task_count),
        "stage_dir": "${STAGE_DIR}",
        "work_dir": "${WORK_DIR}",
        "study": manifest.study,
        "repo": "${REPO_ROOT}",
    }
    command = " ".join(expand_token(token, mapping) for token in stage.command)

    lines: List[str] = ["#!/bin/bash"]
    lines.extend(_directives(manifest, stage, site, chunk, log_dir))
    lines.extend(
        [
            "",
            f"# study {manifest.study}: {manifest.title}",
            f"# stage {stage.name}"
            + (f": {stage.description}" if stage.description else ""),
            "# Generated by cryptic_ip.hpc.render - edit the manifest, not this file.",
            "",
            "set -euo pipefail",
            "",
            f'REPO_ROOT="{_escape_literal(repo_root)}"',
            f'WORK_DIR="{_escape_literal(work_dir)}"',
            f'STAGE_DIR="{_escape_literal(stage_dir)}"',
            f'STATUS_DIR="$STAGE_DIR/{STATUS_DIR}"',
            'mkdir -p "$STAGE_DIR" "$STATUS_DIR" ' + f'"{_escape_literal(log_dir)}"',
            "",
        ]
    )

    if chunk is not None:
        lines.extend(
            [
                'if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then',
                '  echo "This is an array script: set SLURM_ARRAY_TASK_ID, or submit it with sbatch." >&2',
                "  exit 2",
                "fi",
                "",
                "# Array index to shard. Indices always start at zero so none can reach the",
                "# cluster's MaxArraySize, and a resume of a sparse shard set renders the same way.",
            ]
        )
        if _is_contiguous(chunk.shards):
            lines.append(f"SHARD=$(( SLURM_ARRAY_TASK_ID + {chunk.shards[0]} ))")
        else:
            lines.append(f"SHARD_MAP=({chunk.shard_map})")
            lines.append('SHARD="${SHARD_MAP[$SLURM_ARRAY_TASK_ID]}"')
        lines.extend(
            [
                'JOB_ID="${SLURM_ARRAY_JOB_ID:-local}_${SLURM_ARRAY_TASK_ID}"',
                "",
            ]
        )
    else:
        lines.extend(
            [
                "SHARD=0",
                'JOB_ID="${SLURM_JOB_ID:-local}"',
                "",
            ]
        )

    lines.extend(
        [
            'STATUS_FILE="$STATUS_DIR/shard-${SHARD}.json"',
            "",
            "write_status() {",
            '  printf \'{"study": "%s", "stage": "%s", "shard": %s, "status": "%s", '
            '"exit_code": %s, "job_id": "%s", "host": "%s", "attempt": %s, "at": %s}\\n\' \\',
            f'    "{_escape_literal(manifest.study)}" "{_escape_literal(stage.name)}" '
            '"$SHARD" "$1" "$2" "$JOB_ID" "$(hostname)" "${SLURM_RESTART_COUNT:-0}" "$(date +%s)" \\',
            '    > "$STATUS_FILE.tmp" && mv "$STATUS_FILE.tmp" "$STATUS_FILE"',
            "}",
            "",
            "# A signal trap's `exit` re-enters the EXIT trap, so without this flag the",
            "# \"killed\" record would be overwritten by a plain \"failed\" and the one fact",
            "# worth keeping - that the scheduler took the task, not that the code broke -",
            "# would be lost.",
            "KILLED=0",
            "",
            "on_exit() {",
            "  code=$?",
            '  if [ "$KILLED" -eq 1 ]; then exit "$code"; fi',
            '  if [ "$code" -eq 0 ]; then write_status ok 0; else write_status failed "$code"; fi',
            "  exit $code",
            "}",
            "",
            "on_signal() {",
            "  # The scheduler is taking this task away: record it before it goes, so a",
            "  # killed task is a known failure rather than a gap in the census.",
            "  KILLED=1",
            '  if [ -n "${CHILD:-}" ]; then kill -TERM "$CHILD" 2>/dev/null || true; fi',
            '  write_status killed "${1:-143}"',
            '  exit "${1:-143}"',
            "}",
            "",
            "trap on_exit EXIT",
        ]
    )
    for signal in TRAPPED_SIGNALS:
        lines.append(f"trap 'on_signal {128 + _signal_number(signal)}' {signal}")
    lines.extend(["", "write_status running -1", ""])

    if site.module_loads:
        lines.append("module purge 2>/dev/null || true")
        lines.extend(f"module load {name}" for name in site.module_loads)
        lines.append("")
    if site.env_activate:
        lines.extend(site.env_activate)
        lines.append("")

    lines.append('cd "$REPO_ROOT"')
    # One task owns its cores; threaded BLAS inside an array task otherwise
    # oversubscribes the node and every task slows down together. A stage that
    # states a thread count means it - study J pins OMP_NUM_THREADS=1 - so the
    # default goes first and the stage's own value overrides it.
    if "OMP_NUM_THREADS" not in stage.env:
        lines.append('export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"')
    for key, value in sorted(stage.env.items()):
        lines.append(f'export {key}="{_escape_literal(value)}"')
    # The command runs in the background and the script waits on it. A bash
    # trap does not run until the current foreground command returns, so with
    # the command in the foreground the USR1 that SLURM sends before the
    # walltime kill - the whole point of --signal=B:USR1@120 - would not be
    # handled until the task had already finished. A docking copy that takes
    # five hours would never get its warning, and the record would be written
    # too late to be written at all.
    lines.extend([
        "",
        command + " &",
        "CHILD=$!",
        "# `set -e` propagates the child's exit code, and the EXIT trap records it.",
        'wait "$CHILD"',
        "",
    ])
    return "\n".join(lines)


def _is_contiguous(shards: Sequence[int]) -> bool:
    """Whether ``shards`` is an unbroken run, so arithmetic can replace a lookup table."""
    return all(b - a == 1 for a, b in zip(shards, shards[1:]))


_SIGNAL_NUMBERS = {"HUP": 1, "INT": 2, "QUIT": 3, "USR1": 10, "TERM": 15, "XCPU": 24}


def _signal_number(name: str) -> int:
    return _SIGNAL_NUMBERS[name]


def render_study(
    manifest: StudyManifest,
    site: SiteConfig,
    *,
    work_dir: str,
    repo_root: str,
    only_shards: Optional[Mapping[str, Sequence[int]]] = None,
) -> List[RenderedStage]:
    """Every sbatch script this study needs, in dependency order.

    ``only_shards`` restricts a stage to a subset of its shards, which is how a
    resume re-runs the tasks that did not land without re-running the rest.
    """
    rendered: List[RenderedStage] = []
    for stage in manifest.order():
        if stage.is_array:
            wanted = list(only_shards.get(stage.name, range(stage.task_count))) if only_shards \
                else list(range(stage.task_count))
            if not wanted:
                continue
            throttle = stage.throttle or site.max_array_throttle
            for chunk in chunk_shards(wanted, max_array_size=site.max_array_size, throttle=throttle):
                suffix = f"-c{chunk.chunk_index}" if chunk.chunk_index else ""
                rendered.append(
                    RenderedStage(
                        stage=stage.name,
                        chunk_index=chunk.chunk_index,
                        script=render_stage(
                            manifest, stage, site,
                            work_dir=work_dir, repo_root=repo_root, chunk=chunk,
                        ),
                        shards=chunk.shards,
                        depends_on=tuple(stage.depends_on),
                        after=stage.after,
                        script_name=f"{manifest.study.lower()}-{stage.name}{suffix}.sbatch",
                    )
                )
        else:
            if only_shards is not None and stage.name in only_shards and not list(only_shards[stage.name]):
                continue
            rendered.append(
                RenderedStage(
                    stage=stage.name,
                    chunk_index=0,
                    script=render_stage(
                        manifest, stage, site, work_dir=work_dir, repo_root=repo_root,
                    ),
                    shards=(0,),
                    depends_on=tuple(stage.depends_on),
                    after=stage.after,
                    script_name=f"{manifest.study.lower()}-{stage.name}.sbatch",
                )
            )
    return rendered


def write_scripts(rendered: Sequence[RenderedStage], out_dir: Path) -> List[Path]:
    """Write the scripts to ``out_dir`` and return the paths, executable."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    for item in rendered:
        path = out_dir / item.script_name
        path.write_text(item.script)
        path.chmod(0o755)
        paths.append(path)
    return paths
