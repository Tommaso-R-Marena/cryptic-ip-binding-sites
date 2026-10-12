"""Preflight: the site file's claims against what the scheduler actually says."""

from __future__ import annotations

from cryptic_ip.hpc.manifest import StudyManifest
from cryptic_ip.hpc.preflight import (
    FAIL, PASS, UNKNOWN, WARN, markdown, parse_gres_models, parse_nvidia_smi,
    parse_scontrol_config, parse_sinfo_partitions, run_preflight,
)
from cryptic_ip.hpc.site import UNLIMITED, SiteConfig

SINFO = """\
defq*|72:00:00|48|192000|(null)
shared|72:00:00|48|192000|(null)
bigmem|72:00:00|48|1536000|(null)
a100|72:00:00|48|192000|gpu:a100:4(S:0-1)
"""


def _site(**overrides):
    raw = {
        "name": "rockfish",
        "scratch_root": "/tmp",
        "partitions": {
            # 187 GiB, not 192 GB: see Partition's docstring on the unit.
            "cpu": {"name": "defq", "max_walltime": "72:00:00", "max_cpus_per_node": 48,
                    "max_mem_gb_per_node": 187},
        },
    }
    raw.update(overrides)
    return SiteConfig.from_dict(raw)


def _runner(responses):
    def run(command, timeout=60):
        key = command[0]
        return responses.get(key)
    return run


# ----------------------------------------------------------------- the parsers

def test_sinfo_rows_become_partition_limits():
    rows = parse_sinfo_partitions(SINFO)
    assert set(rows) == {"defq", "shared", "bigmem", "a100"}
    assert rows["defq"]["max_walltime"] == 259200
    assert rows["defq"]["cpus"] == 48
    assert round(rows["defq"]["mem_gb"]) == 188
    assert "a100" in rows["a100"]["gres"]


def test_the_default_partition_marker_is_not_part_of_the_name():
    # sinfo prints "defq*" for the default partition. Keeping the star would
    # make every lookup of "defq" miss.
    assert "defq" in parse_sinfo_partitions("defq*|1:00:00|4|8000|(null)")


def test_an_infinite_time_limit_is_understood():
    rows = parse_sinfo_partitions("long|infinite|4|8000|(null)")
    assert rows["long"]["max_walltime"] == UNLIMITED


def test_several_node_types_in_one_partition_merge_to_the_most_generous():
    rows = parse_sinfo_partitions("mix|1:00:00|24|96000|(null)\nmix|8:00:00|48|192000|gpu:a100:2")
    assert rows["mix"]["cpus"] == 48
    assert rows["mix"]["max_walltime"] == 28800
    assert "a100" in rows["mix"]["gres"]


def test_memory_is_read_as_gibibytes_to_match_sbatchs_own_suffix():
    # sinfo reports MB; a "192 GB" node says 192000, which is 187.5 GiB. Reading
    # it as 192 would let a site file ask for more than the node has.
    rows = parse_sinfo_partitions("p|1:00:00|48|192000|(null)")
    assert 187.0 < rows["p"]["mem_gb"] < 188.0


def test_a_memory_plus_suffix_is_handled():
    assert round(parse_sinfo_partitions("p|1:00:00|4|8000+|(null)")["p"]["mem_gb"]) == 8


def test_blank_and_short_rows_are_skipped():
    assert parse_sinfo_partitions("\n\njunk\n") == {}


def test_scontrol_config_yields_its_pairs():
    config = parse_scontrol_config("MaxArraySize            = 1001\nMaxJobCount             = 10000\nnoise\n")
    assert config["MaxArraySize"] == "1001"
    assert config["MaxJobCount"] == "10000"


def test_gres_models_are_read_with_their_counts():
    assert parse_gres_models("gpu:a100:4(S:0-1)") == {"a100": 4}
    assert parse_gres_models("gpu:4") == {"": 4}
    assert parse_gres_models("gpu:a100:2,gpu:v100:2") == {"a100": 2, "v100": 2}


def test_a_null_gres_is_no_gpu_rather_than_an_error():
    assert parse_gres_models("(null)") == {}


def test_non_gpu_gres_is_ignored():
    assert parse_gres_models("nvme:1600") == {}


def test_nvidia_smi_output_becomes_names_and_vram():
    found = parse_nvidia_smi("NVIDIA A100-SXM4-40GB, 40960 MiB\nNVIDIA A100-SXM4-40GB, 40960 MiB\n")
    assert len(found) == 2
    assert round(found[0]["vram_gb"]) == 40


def test_unparseable_nvidia_smi_output_is_dropped_not_guessed():
    assert parse_nvidia_smi("no GPUs found") == []


# ------------------------------------------------------------------ the checker

def test_a_partition_the_cluster_does_not_have_blocks_submission():
    site = _site(partitions={"cpu": {"name": "nosuch", "max_walltime": "72:00:00"}})
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm 23.02"}))
    failing = [c for c in report.checks if c.status == FAIL]
    assert any("nosuch" in c.expected for c in failing)
    assert not report.ok


def test_claiming_more_walltime_than_the_partition_allows_blocks_submission():
    site = _site(partitions={"cpu": {"name": "defq", "max_walltime": "96:00:00"}})
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    assert any(c.status == FAIL and "defq" in c.name for c in report.checks)


def test_claiming_more_cores_or_memory_than_a_node_has_blocks_submission():
    site = _site(partitions={"cpu": {"name": "defq", "max_walltime": "72:00:00",
                                     "max_cpus_per_node": 96, "max_mem_gb_per_node": 512}})
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    detail = next(c.detail for c in report.checks if c.name.startswith("partition cpu"))
    assert "96 cpus" in detail and "512 GB" in detail


def test_a_matching_partition_passes():
    report = run_preflight(_site(), runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    assert next(c for c in report.checks if c.name.startswith("partition cpu")).status == PASS


def test_a_max_array_size_larger_than_the_cluster_allows_blocks_submission():
    site = _site(max_array_size=5000)
    report = run_preflight(site, runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm", "scontrol": "MaxArraySize = 1001\n",
    }))
    check = next(c for c in report.checks if c.name == "MaxArraySize")
    assert check.status == FAIL
    assert "max_array_size: 1001" in check.remedy


def test_a_max_array_size_within_the_limit_passes():
    report = run_preflight(_site(max_array_size=1001), runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm", "scontrol": "MaxArraySize = 1001\n",
    }))
    assert next(c for c in report.checks if c.name == "MaxArraySize").status == PASS


def test_an_account_the_user_does_not_have_blocks_submission():
    site = _site(account="wrong-alloc")
    report = run_preflight(site, runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm", "sacctmgr": "real-alloc\nother-alloc\n",
    }))
    check = next(c for c in report.checks if c.name == "account")
    assert check.status == FAIL and "real-alloc" in check.observed


def test_an_account_the_user_has_passes():
    site = _site(account="real-alloc")
    report = run_preflight(site, runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm", "sacctmgr": "real-alloc\n",
    }))
    assert next(c for c in report.checks if c.name == "account").status == PASS


def test_no_account_is_a_warning_because_some_clusters_require_one():
    report = run_preflight(_site(), runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    assert next(c for c in report.checks if c.name == "account").status == WARN


def test_a_gpu_model_the_partition_does_not_advertise_blocks_submission():
    site = _site(partitions={
        "cpu": {"name": "defq", "max_walltime": "72:00:00"},
        "gpu": {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4,
                "gpu_model": "H100", "gpu_vram_gb": 80},
    })
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    check = next(c for c in report.checks if c.name.startswith("gpu model"))
    assert check.status == FAIL and "a100" in check.observed


def test_less_vram_than_claimed_blocks_submission():
    # The case that decides whether the long-receptor tail can be co-folded at
    # all, or has to be truncated - and truncation has to be pre-registered.
    site = _site(partitions={
        "cpu": {"name": "defq", "max_walltime": "72:00:00"},
        "gpu": {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4,
                "gpu_model": "A100", "gpu_vram_gb": 80},
    })
    report = run_preflight(site, runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm",
        "nvidia-smi": "NVIDIA A100-SXM4-40GB, 40960 MiB\n",
    }))
    check = next(c for c in report.checks if c.name.startswith("gpu VRAM"))
    assert check.status == FAIL
    assert "gpu_vram_gb to 40" in check.remedy


def test_vram_is_unknown_rather_than_assumed_when_there_is_no_gpu_here():
    site = _site(partitions={
        "cpu": {"name": "defq", "max_walltime": "72:00:00"},
        "gpu": {"name": "a100", "max_walltime": "72:00:00", "gpus_per_node": 4,
                "gpu_model": "A100", "gpu_vram_gb": 40},
    })
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    check = next(c for c in report.checks if c.name.startswith("gpu VRAM"))
    assert check.status == UNKNOWN
    assert "GPU job" in check.remedy


def test_a_study_that_does_not_fit_the_site_is_reported_as_one_check():
    manifest = StudyManifest.from_dict({
        "study": "X", "title": "t", "plan": ["README.md"],
        "stages": [{"name": "big", "command": ["echo", "x"],
                    "resources": {"partition": "cpu", "cpus": 200, "walltime": "01:00:00"}}],
    })
    report = run_preflight(_site(), manifests=[manifest], runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    check = next(c for c in report.checks if c.name.startswith("study X"))
    assert check.status == FAIL and "200 cpus" in check.detail


# ---------------------------------------------------------------- the verdict

def test_no_scheduler_means_indeterminate_not_ok():
    # A machine with no sbatch has checked nothing about the cluster. Reporting
    # that as "ok to submit" would be a green light nobody earned.
    report = run_preflight(_site(), runner=_runner({}))
    assert report.verdict == "indeterminate"
    assert not report.ok
    assert "NOT VERIFIED" in markdown(report)


def test_a_failure_blocks_whatever_else_passed():
    site = _site(partitions={"cpu": {"name": "nosuch", "max_walltime": "72:00:00"}})
    report = run_preflight(site, runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    assert report.verdict == "blocked"
    assert "DO NOT SUBMIT" in markdown(report)


def test_everything_matching_on_a_login_node_is_ok_to_submit():
    site = _site(account="real-alloc", max_array_size=1001,
                 verified={"date": "2026-10-10", "by": "test"})
    report = run_preflight(site, runner=_runner({
        "sinfo": SINFO, "sbatch": "slurm", "scontrol": "MaxArraySize = 1001\n",
        "sacctmgr": "real-alloc\n",
    }))
    assert report.verdict == "ok"
    assert "OK to submit" in markdown(report)


def test_an_unverified_site_is_a_warning_in_its_own_report():
    report = run_preflight(_site(), runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    assert next(c for c in report.checks if c.name == "site file verified").status == WARN


def test_the_report_round_trips_through_its_dict_form():
    report = run_preflight(_site(), runner=_runner({"sinfo": SINFO, "sbatch": "slurm"}))
    payload = report.as_dict()
    assert payload["site"] == "rockfish"
    assert payload["verdict"] == report.verdict
    assert len(payload["checks"]) == len(report.checks)
