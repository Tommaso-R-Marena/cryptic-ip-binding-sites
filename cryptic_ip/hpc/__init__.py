"""Cluster execution for the pre-registered studies.

The GitHub Actions workflows in `.github/workflows/` are the reference
implementation of every study: they define the stages, the sharding and the
marker-block reporting. This package runs the same stages under SLURM, where
the compute is large enough that free shared runners cannot reach it.

The modules are deliberately separate so each is testable without a cluster:

* `site`      - what a cluster is: partitions, limits, modules, filesystems.
* `manifest`  - what a study is: stages, resources, dependencies, plan docs.
* `shards`    - how a shard range maps onto SLURM array indices.
* `render`    - how a stage becomes an sbatch script.
* `preflight` - whether the cluster matches what the site file claims.
* `collect`   - how shard outputs merge, and how a missing shard is caught.
"""

from cryptic_ip.hpc.manifest import Stage, StudyManifest, load_manifest
from cryptic_ip.hpc.shards import ArrayChunk, chunk_shards
from cryptic_ip.hpc.site import Partition, SiteConfig, format_walltime, load_site, parse_walltime

__all__ = [
    "ArrayChunk",
    "Partition",
    "SiteConfig",
    "Stage",
    "StudyManifest",
    "chunk_shards",
    "format_walltime",
    "load_manifest",
    "load_site",
    "parse_walltime",
]
