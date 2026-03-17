# Performance and Scaling Guide

## Expected Runtime Targets

Assuming 16 CPU cores, SSD storage, and electrostatics skipped for first pass:

| Proteins | Expected runtime | Throughput |
|---|---:|---:|
| 100 | 6-8 min | 750-1000 proteins/hour |
| 1,000 | 55-75 min | 800-1100 proteins/hour |
| 10,000 | 9-12 hours | 830-1100 proteins/hour |

## Profiling

Run:

```bash
python scripts/profile_pipeline.py --pdb tests/data/structures/1ZY7.pdb --report profiling_report.json --skip-electrostatics
```

The report includes:
- single protein runtime
- memory peak per structure
- function-level bottlenecks
- optimization recommendations

## Benchmark Tests

```bash
pytest tests/benchmarks -v
```

## HPC Deployment

- SLURM: `scripts/hpc/slurm_submit.sh`
- AWS Batch template: `scripts/hpc/aws_batch_job.json`
- GCP Batch template: `scripts/hpc/gcp_batch_job.json`
- Work queue runner: `scripts/hpc/work_queue_runner.py`
