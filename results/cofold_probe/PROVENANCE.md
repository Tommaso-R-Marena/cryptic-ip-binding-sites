# Boltz-2 CPU feasibility probe: co-folding is not affordable on these runners

**Run.** Cofold probe run 36729986766, Probe job 109936721893, commit `006d8b8`,
branch `claude/ip-binding-studies-26p6u4`. Started 2026-09-30 14:33 UTC, concluded
16:05 UTC, conclusion `success` (the probe is written to report a timeout as data
rather than to fail the job).

**Machine.** 4 CPUs, 15 GB RAM, ubuntu-latest. No GPU.

**Result.** One record, and it is decisive:

| length | ligand | affinity head | completed | wall time |
|---|---|---|---|---|
| 100 residues | IHP (IP6) as a CCD ligand | requested | **no** | **timed out at 5,400 s** |

The probe stops at the first length that does not complete, so 300 and 600 residues
were never attempted. Boltz-2 itself installed cleanly in 1 m 39 s, so this is a
compute limit, not a packaging problem. A partial model (`probe_100_model_0.pdb`) was
written, so the prediction was making progress — just nowhere near the budget.

**The settings were already the cheapest possible**: `--accelerator cpu`,
`--diffusion_samples 1`, `--recycling_steps 0`, `msa: empty`. There is no cheaper
configuration to retreat to.

**Decision.** Study M (Boltz-2 co-folding of IP6 into candidate pockets) is **not
affordable on GitHub-hosted runners** and is not run. A 100-residue protein is smaller
than every candidate in the shortlist; at 90 minutes without completing one of them,
even a token subset would cost more than the budget and would not be a result worth
having. Per the pre-agreed rule, no reduced version is run and presented as a finding.

The geometric question co-folding would have answered is taken up instead by
**study N** (`docs/TEMPLATE_PLAN.md`), which transplants real crystal IP6 sites onto
candidate pockets by rigid superposition and a clash test, and costs minutes rather
than GPU-days.

Study M remains viable on a GPU. It is reserved, not abandoned.

**Extraction.** `scripts/extract_log_block.py read <log> COFOLD_PROBE_JSON`. The probe
streams each record as it completes *and* re-emits them in the marker block; the two
independent copies were compared programmatically and agree, so no figure here was
retyped. The job-log artifact itself is not reachable from the analysis sandbox
(blob storage returns 403 through the proxy), which is why the block is extracted from
the log text rather than downloaded.
