"""Every workflow parses, and every job id is a string.

YAML reads bare ``null``, ``true``, ``on`` and the like as non-strings; a job
named ``null`` makes GitHub reject the whole file at push time.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

WORKFLOWS = sorted((Path(__file__).resolve().parents[1] / ".github" / "workflows").glob("*.yml"))


@pytest.mark.parametrize("path", WORKFLOWS, ids=lambda p: p.name)
def test_workflow_jobs_are_named_by_strings(path):
    data = yaml.safe_load(path.read_text())
    jobs = data["jobs"]
    assert jobs and all(isinstance(name, str) for name in jobs), list(jobs)
    for name, job in jobs.items():
        needs = job.get("needs", [])
        for dependency in [needs] if isinstance(needs, str) else needs:
            assert dependency in jobs, f"{path.name}: {name} needs unknown job {dependency!r}"
