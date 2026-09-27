# Study I provenance

- Plan: `docs/ORTHOGONAL_PLAN.md`
- **Analysis run:** orthogonal 36259814742 (head caa92ac). Disorder job succeeded (184 proteins,
  184 sequences, no errors); **all 25 SHARK-capture shards succeeded** and contributed to I3.
- **Report run:** orthogonal 36279873788 (Report job 108509740644, head a7fa5a1), dispatched
  with `analysis_run=36259814742` so it reused the artifacts above and recomputed nothing.
- `orthogonal.json` extracted from `BEGIN_ORTHOGONAL_JSON`, `ORTHOGONAL.md` from
  `BEGIN_ORTHOGONAL_MD_B64` (gzip, CRC-checked) with `scripts/extract_log_block.py read`.
  Every figure in the JSON that also appears in the markdown was cross-checked against the
  CRC-verified markdown by script.
- `orthogonal_candidates.csv` is **not** committed here: the per-candidate table is in the
  `orthogonal-report` artifact of run 36279873788. Its base64 block is large enough that the
  three files together exceed what one log read returns, and nothing in the decisions depends
  on it.

## Two faults in the first attempt, and what they cost

1. **shark-dive did not finish within its job budget.** It ran 340 minutes against the 19
   candidates study H could not evaluate and was killed by `timeout-minutes`. **I4 therefore
   reports 0 queries** — not "no remote homologues found", but "the search did not complete".
   It is a pre-registered descriptive arm and no decision rests on it.
2. **The run concluded `cancelled`, so the report was skipped.** GitHub marks a job killed by
   `timeout-minutes` as cancelled, and the report job was gated on `if: !cancelled()`. The
   report was recovered by dispatching the report-only path; the gate is now `always()` and
   the emitted blocks are ordered largest-first (commit a7fa5a1). Pushing that fix
   auto-triggered full run 36279860843, which was cancelled at once; its outputs are not kept.

An earlier report-only run, 36279284598, produced the same numbers but emitted the blocks in
the old order, which made them unreadable from the log. It is superseded by 36279873788.
