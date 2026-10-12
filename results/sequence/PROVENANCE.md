# Study Q6 method validation: provenance

**Run.** `scripts/sequence_signature.py validate --windows 40 --length 24 --shuffles 400`, seed 20261012, on the session container (single core), at commit `8f817d1`.

**Plan.** `docs/SEQUENCE_PLAN.md`, which pre-registers study Q6 and states that this validation is a property of the method and not a study result.

Every figure below is generated from `validation.json` by `scripts/sequence_signature.py provenance`. Nothing is transcribed, and a test checks that this file is what that subcommand produces.

## What this is

Three pools of synthetic sequence windows, matched on composition, one of which has a planted ordering signature. It answers one question: does the machinery find an ordering signature when one is there, and stay quiet when composition alone would explain the difference?

**It says nothing about whether annotated inositol-phosphate sites have a sequence signature.** That is study Q6, asked of real sequences, and it has not run.

## The pools (40 windows of 24 residues each)

| pool | mean raw | mean entropy (bits) | basic fraction | distinct letters | mean z | windows firing |
|---|---|---|---|---|---|---|
| planted helical face | 0.586 | 3.211 | 0.270 | 11 | +4.05 | 100% |
| same letters, shuffled | 0.172 | 3.211 | 0.270 | 11 | -0.24 | 0% |
| basic-rich, no ordering | 0.183 | 3.253 | 0.270 | 11 | -0.13 | 2% |

Statistic: `helix_charge_periodicity`. Null: 400 composition-preserving shuffles per window, so each window is compared only with reorderings of itself.

## Why the comparison is clean

The shuffled control holds the planted pool's letters exactly, so its entropy differs by 0.0e+00 bits and its basic fraction not at all. The basic-rich control is built from the planted pool's own per-window basic counts rather than from their mean fraction, which closed a residual of 0.021 in basic fraction - enough to fail the match at twelve windows. Across all three pools the entropy spread is 0.042 bits and the basic-fraction spread is 0.0e+00.

So the z column cannot be reporting composition: there is no composition difference left for it to report.

## Checks

| check | result |
|---|---|
| all pools have the same basic fraction | pass |
| all pools within 0.1 bits | pass |
| composition control has identical basic fraction | pass |
| composition control stays quiet | pass |
| planted pool fires | pass |
| shuffle control has identical basic fraction | pass |
| shuffle control has identical entropy | pass |
| shuffle control stays quiet | pass |

All 8 checks pass: **True**. A failed check exits non-zero and fails the workflow run.

## What this validation may not be used for

* It is not evidence that any real site has a sequence signature.
* It is not a benchmark of the detector: no structure is involved.
* It does not license adding a sequence term to the score. `signature.ESTABLISHED_SIGNATURES` stays empty until study Q6 reports, and `sequence_contribution` returns an abstention until then.
* The planted signal is a charge periodicity at exactly 3.6 residues, which is the cleanest possible case. A real helical face is noisier, so the firing rate here is an upper bound on sensitivity, not an estimate of it.
