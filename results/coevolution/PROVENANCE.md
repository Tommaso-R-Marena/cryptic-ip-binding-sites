# Study L provenance

- Plan: `docs/COEVOLUTION_PLAN.md`
- Run: coevolution 36264528735 (Report job 108469945312), head 8f69348
- Proteins: specificity run 36020458866 `specificity-rerank` (`proteins_combined.csv.gz`)
- Pocket shards for the rule-only arm: proteome screen run 35935291031
- `coevolution.json` extracted from `BEGIN_COEVOLUTION_JSON`, `COEVOLUTION.md` from
  `BEGIN_COEVOLUTION_MD_B64` (gzip, CRC-checked) with `scripts/extract_log_block.py read`.
  Every figure in the JSON that also appears in the markdown was cross-checked against the
  CRC-verified markdown.
