"""How a shard range maps onto SLURM array indices.

Study J is the reason this module exists rather than being one f-string at the
submission site. Its docking stage wanted 272 shards; GitHub generates at most
256 jobs from one matrix, and past that limit it does not fail the run - it
drops the job, so the report ran on zero shards and the run went green
(run 36883632819). SLURM's `MaxArraySize` is a different mechanism with the
same shape of consequence: an index at or above the limit is rejected, so a
naive `--array=0-1499` either fails at submit time or, if someone "fixes" it by
submitting the tail as `--array=1000-1499`, fails again for the same reason.

So shards are chunked, every chunk indexes from zero, and the chunk carries an
explicit map from array index to shard number. Two properties are pinned by
tests: the chunks tile the requested shards exactly - no gap, no overlap, no
duplicate - and no chunk ever emits an index the cluster will refuse.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


class ShardError(ValueError):
    """A chunking request that cannot be honoured."""


@dataclass(frozen=True)
class ArrayChunk:
    """One `sbatch --array` submission, and the shards its indices stand for."""

    chunk_index: int
    shards: Tuple[int, ...]
    throttle: int

    def __post_init__(self) -> None:
        if not self.shards:
            raise ShardError("an array chunk cannot be empty")
        if self.throttle < 1:
            raise ShardError(f"throttle must be at least 1, got {self.throttle}")

    @property
    def count(self) -> int:
        return len(self.shards)

    @property
    def array_spec(self) -> str:
        """The value for `sbatch --array`, throttled so one study cannot fill the queue."""
        last = self.count - 1
        throttle = min(self.throttle, self.count)
        span = "0" if last == 0 else f"0-{last}"
        return f"{span}%{throttle}"

    @property
    def shard_map(self) -> str:
        """The chunk's shards as a bash array literal, indexed by array task id."""
        return " ".join(str(shard) for shard in self.shards)

    def shard_for(self, task_id: int) -> int:
        """Which shard array index ``task_id`` runs."""
        if not 0 <= task_id < self.count:
            raise ShardError(f"array index {task_id} is outside this chunk of {self.count}")
        return self.shards[task_id]


def chunk_shards(
    shards: Iterable[int],
    *,
    max_array_size: int,
    throttle: int = 50,
) -> List[ArrayChunk]:
    """Split ``shards`` into array submissions no cluster will refuse.

    ``max_array_size`` is SLURM's own `MaxArraySize`, whose meaning is that the
    largest permitted index is one less than the value. Chunks are therefore at
    most ``max_array_size - 1`` tasks, which also leaves index headroom rather
    than sitting exactly on the limit.
    """
    if max_array_size < 2:
        raise ShardError(f"max_array_size must be at least 2, got {max_array_size}")
    wanted = list(shards)
    if not wanted:
        raise ShardError("no shards to run")
    if any(shard < 0 for shard in wanted):
        raise ShardError("shard numbers cannot be negative")
    duplicates = sorted({shard for shard in wanted if wanted.count(shard) > 1})
    if duplicates:
        raise ShardError(f"shard(s) {duplicates} requested more than once")
    ordered = sorted(wanted)
    per_chunk = max_array_size - 1
    chunks: List[ArrayChunk] = []
    for index, start in enumerate(range(0, len(ordered), per_chunk)):
        chunks.append(
            ArrayChunk(
                chunk_index=index,
                shards=tuple(ordered[start:start + per_chunk]),
                throttle=throttle,
            )
        )
    _assert_tiles(chunks, ordered)
    return chunks


def _assert_tiles(chunks: Sequence[ArrayChunk], wanted: Sequence[int]) -> None:
    """Fail loudly if the chunks do not cover exactly the requested shards.

    This is the check study J did not have. It is cheap, and the failure it
    catches is one that otherwise looks like a successful run on a short
    dataset.
    """
    seen: List[int] = []
    for chunk in chunks:
        seen.extend(chunk.shards)
    if sorted(seen) != sorted(wanted):
        missing = sorted(set(wanted) - set(seen))
        extra = sorted(set(seen) - set(wanted))
        raise ShardError(
            f"chunking did not tile the shard range: {len(seen)} tasks for {len(wanted)} shards"
            f"{f', missing {missing}' if missing else ''}"
            f"{f', unexpected {extra}' if extra else ''}"
        )
    if len(seen) != len(set(seen)):
        raise ShardError("chunking produced the same shard twice")


def index_map(chunks: Sequence[ArrayChunk]) -> Dict[Tuple[int, int], int]:
    """``(chunk_index, array_task_id) -> shard`` for the whole submission."""
    return {
        (chunk.chunk_index, task_id): shard
        for chunk in chunks
        for task_id, shard in enumerate(chunk.shards)
    }
