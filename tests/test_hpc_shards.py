"""Array chunking: the SLURM analogue of the matrix cap that lost study J a run.

GitHub drops the job when a matrix exceeds 256 entries, so study J's 272-shard
docking stage ran nothing and the run went green (run 36883632819). SLURM
refuses an index at or above MaxArraySize instead, which is louder but has the
same cause. These tests pin the two properties that make that impossible here:
the chunks cover exactly the shards asked for, and no chunk emits an index the
cluster would refuse.
"""

from __future__ import annotations

import pytest

from cryptic_ip.hpc.shards import ArrayChunk, ShardError, chunk_shards, index_map


@pytest.mark.parametrize("n_shards", [1, 2, 255, 256, 257, 272, 1000, 1001, 1002, 5000])
@pytest.mark.parametrize("max_array_size", [2, 11, 256, 1001])
def test_chunks_tile_the_shard_range_exactly(n_shards, max_array_size):
    chunks = chunk_shards(range(n_shards), max_array_size=max_array_size, throttle=10)
    covered = [shard for chunk in chunks for shard in chunk.shards]
    assert sorted(covered) == list(range(n_shards))      # nothing missing
    assert len(covered) == len(set(covered))             # nothing twice


@pytest.mark.parametrize("n_shards", [1, 272, 1001, 5000])
@pytest.mark.parametrize("max_array_size", [2, 256, 1001])
def test_no_chunk_emits_an_index_the_cluster_would_refuse(n_shards, max_array_size):
    for chunk in chunk_shards(range(n_shards), max_array_size=max_array_size, throttle=5):
        assert chunk.count <= max_array_size - 1
        highest_index = chunk.count - 1
        assert highest_index < max_array_size


def test_study_js_272_shards_fit_one_chunk_under_a_1001_cap():
    chunks = chunk_shards(range(272), max_array_size=1001, throttle=40)
    assert len(chunks) == 1
    assert chunks[0].array_spec == "0-271%40"


def test_the_same_272_shards_split_under_githubs_256_cap():
    # The exact case that was dropped silently. Here it becomes two submissions.
    chunks = chunk_shards(range(272), max_array_size=256, throttle=20)
    assert len(chunks) == 2
    assert [chunk.count for chunk in chunks] == [255, 17]
    assert sorted(s for c in chunks for s in c.shards) == list(range(272))


def test_a_sparse_resume_set_keeps_its_own_shard_numbers():
    chunk, = chunk_shards([3, 117, 245], max_array_size=1001, throttle=3)
    assert chunk.shards == (3, 117, 245)
    assert chunk.array_spec == "0-2%3"
    assert chunk.shard_for(0) == 3
    assert chunk.shard_for(2) == 245


def test_shards_are_ordered_even_when_requested_out_of_order():
    chunk, = chunk_shards([245, 3, 117], max_array_size=1001)
    assert chunk.shards == (3, 117, 245)


def test_the_throttle_never_exceeds_the_task_count():
    # "%50" on a two-task array is not wrong, but it reads as a bigger job than
    # it is, and the spec should say what will actually run.
    chunk, = chunk_shards([1, 2], max_array_size=1001, throttle=50)
    assert chunk.array_spec == "0-1%2"


def test_a_single_task_array_renders_a_single_index():
    chunk, = chunk_shards([5], max_array_size=1001, throttle=10)
    assert chunk.array_spec == "0%1"
    assert chunk.shard_map == "5"


def test_the_shard_map_is_a_bash_array_literal():
    chunk, = chunk_shards([2, 4, 6], max_array_size=1001)
    assert chunk.shard_map == "2 4 6"


def test_index_map_covers_every_task_of_every_chunk():
    chunks = chunk_shards(range(10), max_array_size=5, throttle=2)
    mapping = index_map(chunks)
    assert len(mapping) == 10
    assert sorted(mapping.values()) == list(range(10))
    assert mapping[(0, 0)] == 0


@pytest.mark.parametrize("bad_kwargs", [{"max_array_size": 1}, {"max_array_size": 0}])
def test_an_unusable_cap_is_refused(bad_kwargs):
    with pytest.raises(ShardError, match="max_array_size"):
        chunk_shards(range(4), **bad_kwargs)


def test_no_shards_is_refused_rather_than_submitting_an_empty_array():
    with pytest.raises(ShardError, match="no shards"):
        chunk_shards([], max_array_size=1001)


def test_a_duplicated_shard_is_refused():
    with pytest.raises(ShardError, match="more than once"):
        chunk_shards([1, 2, 2], max_array_size=1001)


def test_a_negative_shard_is_refused():
    with pytest.raises(ShardError, match="negative"):
        chunk_shards([-1, 0], max_array_size=1001)


def test_an_empty_chunk_cannot_be_constructed():
    with pytest.raises(ShardError, match="cannot be empty"):
        ArrayChunk(chunk_index=0, shards=(), throttle=1)


def test_a_nonpositive_throttle_is_refused():
    with pytest.raises(ShardError, match="throttle"):
        ArrayChunk(chunk_index=0, shards=(1,), throttle=0)


def test_shard_for_rejects_an_index_outside_the_chunk():
    chunk, = chunk_shards([1, 2], max_array_size=1001)
    with pytest.raises(ShardError, match="outside this chunk"):
        chunk.shard_for(2)
