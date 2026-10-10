# Running the studies on Rockfish

This describes the cluster harness in `cryptic_ip/hpc/`, `scripts/rockfish.py`,
`hpc/sites/` and `hpc/studies/`. It is engineering, not science: the studies
themselves stay in `scripts/<study>.py` and run byte-identically to the way the
GitHub workflows run them. What this adds is where they run, how they are
sharded, and whether the results add up afterwards.

**Nothing here has run on Rockfish yet.** Every limit in
`hpc/sites/rockfish.yml` is an assumption written without access to the
machine, and the harness knows it: `submit` refuses an unverified site. The
first real step is `preflight` on a login node.

---

## 1. What this is for

docs/GPU_REQUIREMENTS.md sets out the two things the project cannot do on free
shared runners, and one it can:

| | needs | why |
|---|---|---|
| **Study M**, co-folding candidates, controls, binders and ligand decoys | GPU, VRAM-bound | the measured CPU probe did not finish one 100-residue protein in 5,400 s, and the candidates have median length 436 |
| **Study P**, molecular dynamics on the shortlist | GPU, hours-bound | every candidate is a single apo conformer; this is the project's biggest standing weakness |
| **Study J**, flexible-receptor redocking | CPU only | Vina is CPU-only; this is a core request, not a GPU one |

Study J is the one that has already run to completion elsewhere, which makes it
the harness's test case: a port that reproduces 209 of 272 copies scored over 29
strict homology groups, J1 −0.078 [−0.140, −0.031] per group, is a working
harness. One that does not is a broken harness, and we would know which.

---

## 2. First run, in order

```bash
# On a Rockfish login node, in a clone of this repository.
module load anaconda3
bash hpc/rockfish/build_env.sh                 # writes hpc/rockfish/env_manifest.json
conda activate cryptic-ip

python scripts/rockfish.py preflight --out preflight.json --markdown preflight.md
```

Preflight compares every claim in the site file against what the scheduler
actually reports, and tells you which line to correct. Expect to correct
several. Then, for the GPU numbers specifically, which `sinfo` cannot answer:

```bash
srun -p a100 --gres=gpu:1 -t 00:10:00 --pty bash
python scripts/rockfish.py preflight            # now the VRAM check can answer
```

When the report comes back `OK to submit`, record that in the site file:

```yaml
verified:
  date: "2026-10-20"
  by: "tmarena, preflight on login02"
```

Note `date`, not `on`: in YAML the bare key `on` is the boolean true, so a
`verified: {on: ...}` block parses as `{True: ...}` and the date silently
disappears. The loader refuses that spelling rather than reading it as
unverified.

Then stage study J's inputs and go:

```bash
WORK=$HOME/scratch4-cryptic-ip/j
mkdir -p "$WORK"
# census.csv, census/ccd and the structures, from the redocking run that
# selected the 272 copies (run 36020459246 on GitHub).
python scripts/rockfish.py validate
python scripts/rockfish.py plan   --study hpc/studies/flexible.yml --out-dir /tmp/j-scripts
python scripts/rockfish.py submit --study hpc/studies/flexible.yml
python scripts/rockfish.py status --study hpc/studies/flexible.yml
```

---

## 3. What refuses to run, and why

Three refusals are deliberate. They are the reason to go through this script
rather than writing sbatch files by hand.

**A study whose pre-registered plan is not in the repository is not submitted.**
Every manifest names the plan documents its study depends on, and `submit`
checks they exist. Study M and study P are fully engineered in
`hpc/studies/cofold.yml` and `hpc/studies/md.yml` and both are refused today,
because `docs/COFOLD_PLAN.md` and `docs/MD_PLAN.md` have not been written.
That is the project's own rule — pre-registration before data — enforced by the
tool instead of by memory.

**A site nobody has checked is not submitted against.** Until `verified.date` is
filled in, `submit` stops and points at `preflight`. `--allow-unverified-site`
overrides it, and says so in the output.

**A stage whose shards do not reconcile does not produce a merged result.** See
§5.

---

## 4. How a study is described

A manifest is a list of stages. Each stage is a command, run once or once per
shard, with the resources one task needs:

```yaml
- name: dock
  shards: 272          # one array task per copy
  throttle: 40         # at most 40 running at once
  output: flexible-{shard}.json
  resources:
    partition: shared
    cpus: 2
    mem_gb: 8
    walltime: "05:55:00"
  env:
    OMP_NUM_THREADS: "1"
  command: [python, scripts/flexible.py, dock, --shard, "{shard}", ...]
```

Placeholders are `{shard}`, `{shards}`, `{stage_dir}`, `{work_dir}`, `{study}`
and `{repo}`; anything else is refused at load time, because a typo in a command
line is a task that burns its walltime producing nothing. Outputs are relative
to the stage directory, and a sharded stage whose output does not vary by shard
is refused — every task would write over the others.

`depends_on` wires the stages, and `after: any` is how a reporting stage says it
must run even when shards failed. That is the same decision the GitHub workflow
made with `always()` rather than `!cancelled()`: the report matters most when
something went wrong, and a report you only see on a clean sweep is a report you
never see.

Three units and limits worth knowing, each of which has already cost somebody a
submission somewhere:

* `mem_gb` is **GiB**, like sbatch's own `G` suffix. A node sold as 192 GB
  reports 192000 MB, which is 187.5 GiB; writing 192 asks for more than the node
  has and the job is never scheduled.
* `walltime` follows sbatch's grammar, where `30:00` is thirty **minutes** and
  `30:00:00` is thirty hours.
* `walltime` must exceed the stage's own internal budget. Study J's docking
  shards record their own censoring when they run out of `--shard-budget`; if
  the scheduler kills the task first, that record is never written, and a
  reported timeout becomes a vanished shard.

---

## 5. Why collection is a reconciliation

Study J merged 270 records from a 272-shard run. The two absent copies had lost
their runners mid-job, so nothing was written for either — they did not appear
among the 63 failures, they were simply not in the census. The only thing that
revealed them was the count.

So every task writes a status record before it starts and rewrites it when it
ends, through a trap that fires on the signals SLURM uses to take a job away.
Each shard then lands in exactly one of five states:

| state | meaning |
|---|---|
| `produced` | an output file exists, and is merged |
| `failed` | the task recorded a non-zero exit |
| `killed` | the scheduler took the task away |
| `running` | a record was written but never completed |
| `missing` | nothing at all was written |

`collect` refuses to write a merged record file while any shard is `missing`.
To proceed anyway you pass `--allow-missing N`, which makes the shortfall a
stated decision rather than a quietly shorter dataset. `resume` re-runs exactly
the shards that did not land:

```bash
python scripts/rockfish.py resume --study hpc/studies/flexible.yml --stage dock
python scripts/rockfish.py resume --study hpc/studies/flexible.yml --stage dock --include-failed
```

The default re-runs what the cluster lost; `--include-failed` also re-runs what
failed in the code, which is usually not what you want until the code changed.

A second refusal sits beside it: a stage where **no** shard produced output is a
structural failure, not a degenerate result. A run where every task ran and
every copy failed still writes files. Study J's report once ran on zero shards,
because a 272-entry matrix exceeded GitHub's 256-job cap and the job was dropped
rather than failed — it reported "not evaluable" and the run went green.

---

## 6. Array sizes

SLURM refuses an array index at or above `MaxArraySize`. The chunker splits a
shard range into submissions no cluster will refuse, every chunk indexing from
zero, with the chunk carrying the map back to shard numbers. Two properties are
pinned by tests: the chunks cover exactly the shards requested — no gap, no
overlap, no duplicate — and no chunk emits an index that would be rejected.

This is the same failure mode as the GitHub matrix cap, which is why it is a
module with tests rather than an f-string at the submission site.

---

## 7. What is not done yet

* **No part of this has run on Rockfish.** The site file is unverified; the
  numbers in it are to be corrected by `preflight`, not trusted.
* **`account` is null.** Fill it from `sacctmgr show assoc user=$USER
  format=Account`. A guessed account is a rejected job, which is the better
  failure, but still a failure.
* **The GPU environment is deliberately unpinned.** Which OpenMM and co-folding
  builds work depends on the driver on the GPU nodes, which this repository has
  not measured. `build_env.sh --gpu` prints the commands rather than running
  them, and records what ends up installed.
* **Study M and study P need their plans written** before either can run. The
  harder questions are in the plans, not the manifests: for study M, whether the
  long-receptor tail is truncated and whether real MSAs are used, both of which
  change what the result means; for study P, how long, how many replicates, and
  what result would count *against* the candidate list.
* **The input staging for study J is manual.** The workflow pulled the census and
  structures from a previous GitHub run's artifacts; on the cluster they have to
  be placed in the work directory, and `check-inputs` only verifies they are
  there.
