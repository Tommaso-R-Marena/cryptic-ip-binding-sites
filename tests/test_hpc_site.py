"""The site layer: SLURM time strings, partition limits, and the verified marker."""

from __future__ import annotations

import pytest

from cryptic_ip.hpc.site import (
    UNLIMITED, Partition, SiteConfig, SiteError, format_walltime, load_site, parse_walltime,
)

ROCKFISH = "hpc/sites/rockfish.yml"


@pytest.mark.parametrize(
    "text,seconds",
    [
        ("30", 1800),               # bare field is minutes, per sbatch
        ("30:00", 1800),            # minutes:seconds
        ("05:55:00", 21300),        # hours:minutes:seconds
        ("72:00:00", 259200),
        ("1-00:00:00", 86400),      # days-hours:minutes:seconds
        ("2-12", 216000),           # days-hours
        ("1-06:30", 109800),        # days-hours:minutes
        ("INFINITE", UNLIMITED),
        ("unlimited", UNLIMITED),
    ],
)
def test_parse_walltime_covers_slurms_spellings(text, seconds):
    assert parse_walltime(text) == seconds


def test_the_ambiguous_case_follows_sbatch_not_intuition():
    # "30:00" is half an hour, not thirty hours. Reading it the other way would
    # ask for 60x the walltime and the job would never be scheduled.
    assert parse_walltime("30:00") == 1800
    assert parse_walltime("30:00:00") == 108000


@pytest.mark.parametrize("bad", ["", "soon", "1:2:3:4", "-5", "01:xx"])
def test_parse_walltime_rejects_nonsense(bad):
    with pytest.raises(SiteError):
        parse_walltime(bad)


@pytest.mark.parametrize("seconds", [0, 59, 60, 3599, 21300, 86399, 86400, 259200])
def test_format_walltime_round_trips(seconds):
    assert parse_walltime(format_walltime(seconds)) == seconds


def test_format_walltime_uses_the_day_form_past_a_day():
    assert format_walltime(86400) == "1-00:00:00"
    assert format_walltime(21300) == "05:55:00"


def test_unlimited_formats_back_to_a_word_slurm_accepts():
    assert format_walltime(UNLIMITED) == "UNLIMITED"


def test_partition_needs_a_name_and_a_cap():
    with pytest.raises(SiteError, match="missing"):
        Partition.from_dict("cpu", {"name": "defq"})


def test_partition_rejects_unknown_keys():
    with pytest.raises(SiteError, match="unknown keys"):
        Partition.from_dict("cpu", {"name": "defq", "max_walltime": "1:00:00", "cores": 48})


def test_a_gpu_partition_must_state_its_vram():
    # VRAM is the binding constraint for co-folding (docs/GPU_REQUIREMENTS.md),
    # so a GPU partition that leaves it blank is a planning error, not a default.
    with pytest.raises(SiteError, match="gpu_vram_gb"):
        Partition.from_dict("gpu", {"name": "a100", "max_walltime": "1:00:00", "gpus_per_node": 4})


def _site(**overrides):
    raw = {
        "name": "t",
        "scratch_root": "/scratch/$USER",
        "partitions": {"cpu": {"name": "defq", "max_walltime": "72:00:00"}},
    }
    raw.update(overrides)
    return raw


def test_site_requires_name_partitions_and_scratch():
    for missing in ("name", "partitions", "scratch_root"):
        raw = _site()
        del raw[missing]
        with pytest.raises(SiteError, match=missing):
            SiteConfig.from_dict(raw)


def test_site_rejects_unknown_keys_rather_than_ignoring_them():
    with pytest.raises(SiteError, match="unknown keys"):
        SiteConfig.from_dict(_site(queue="defq"))


@pytest.mark.parametrize("field,value", [("max_array_size", 0), ("max_array_throttle", 0)])
def test_site_rejects_nonpositive_limits(field, value):
    with pytest.raises(SiteError, match="positive"):
        SiteConfig.from_dict(_site(**{field: value}))


def test_verified_on_is_refused_because_yaml_reads_it_as_a_boolean():
    # In YAML 1.1 the bare key `on` is the boolean true, so `verified: {on: ...}`
    # parses as {True: ...}: the date disappears and the site silently reads as
    # unverified. Refusing the spelling is the only safe reading.
    with pytest.raises(SiteError, match="verified.date"):
        SiteConfig.from_dict(_site(verified={True: "2026-10-10", "by": "me"}))
    with pytest.raises(SiteError, match="verified.date"):
        SiteConfig.from_dict(_site(verified={"on": "2026-10-10"}))


def test_verified_date_is_what_marks_a_site_checked():
    site = SiteConfig.from_dict(_site(verified={"date": "2026-10-10", "by": "someone"}))
    assert site.is_verified
    assert site.verified_on == "2026-10-10"
    assert not SiteConfig.from_dict(_site()).is_verified


def test_unknown_partition_role_names_the_ones_that_exist():
    site = SiteConfig.from_dict(_site())
    with pytest.raises(SiteError, match="no 'gpu' partition"):
        site.partition("gpu")


def test_scratch_expands_the_user(monkeypatch):
    monkeypatch.setenv("USER", "tmarena")
    site = SiteConfig.from_dict(_site(scratch_root="/scratch/$USER"))
    assert site.scratch_for("flexible") == "/scratch/tmarena/flexible"


def test_scratch_refuses_an_unset_variable_rather_than_making_a_dollar_directory(monkeypatch):
    # Left to expandvars, "$NOPE" survives and the study writes into a
    # directory literally named "$NOPE".
    monkeypatch.delenv("NOPE", raising=False)
    site = SiteConfig.from_dict(_site(scratch_root="/scratch/$NOPE"))
    with pytest.raises(SiteError, match="not set in this environment"):
        site.scratch_for("flexible")


def test_the_shipped_rockfish_file_loads():
    site = load_site(ROCKFISH)
    assert site.name == "rockfish"
    assert {"cpu", "shared", "bigmem", "gpu"} <= set(site.partitions)
    assert site.partition("gpu").has_gpu
    assert site.partition("gpu").gpu_vram_gb


def test_the_shipped_rockfish_file_is_still_marked_unverified():
    # This test is the honesty marker. Nobody has checked these limits against
    # the cluster, and `submit` refuses an unverified site. When someone runs
    # preflight on a login node and fills in verified.date, this test should be
    # updated in the same commit - deliberately, not incidentally.
    assert not load_site(ROCKFISH).is_verified
