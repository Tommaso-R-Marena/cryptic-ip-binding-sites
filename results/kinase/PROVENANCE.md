# Study K provenance

- Plan: `docs/KINASE_PLAN.md`
- Run: kinase 36264528692 (Report job 108471178863), head 8f69348
- Apo arm: rerank run 36055694750 `rerank-arms-*`, not re-docked
- Census/inputs: redocking run 36020459246
- `kinase.json` extracted from `BEGIN_KINASE_JSON`, `KINASE.md` from `BEGIN_KINASE_MD_B64`
  (gzip, CRC-checked) with `scripts/extract_log_block.py read`. Every figure in the JSON
  that also appears in the markdown was cross-checked against the CRC-verified markdown.

A second run (36265022079) was triggered by the plan-text correction touching
`docs/KINASE_PLAN.md`, which is one of the workflow's path triggers. It recomputes the
same arms and its outputs are not kept.
