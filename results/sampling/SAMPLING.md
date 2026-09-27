## The sampling ceiling (docs/SAMPLING_PLAN.md)

218 copies in 29 strict groups compared across arms; the re-ranking weight is frozen at w = 0.1.

| arm | copies | seed runs | mean poses | median s/copy |
|---|---|---|---|---|
| e32 | 218 | 654 | 39.9 | 448 |
| e128 | 218 | 654 | 39.0 | 1987 |
| e512 | 61 | 61 | 36.8 | 2253 |

### Decisions (Holm across G1-G3)

| | question | difference (group) | Holm p | decision |
|---|---|---|---|---|
| G1 | sampling ceiling, E128 - E32 | 0.019 [-0.001, 0.047] | 0.152 | **search-saturated** |
| G2 | top-pose success, E128 - E32 | 0.002 [-0.002, 0.009] | 0.715 | **no gain** |
| G3 | re-ranked (w = 0.1) top-pose success at E128 - Vina top-pose at E32 | 0.041 [0.004, 0.083] | 0.081 | **inconclusive** |

### Levels (group estimand)

| arm | ceiling | Vina top pose | re-ranked top pose |
|---|---|---|---|
| e32 | 0.385 [0.256, 0.516] | 0.092 [0.032, 0.172] | 0.140 [0.066, 0.224] |
| e128 | 0.404 [0.277, 0.534] | 0.094 [0.035, 0.173] | 0.133 [0.064, 0.216] |

**G4 (descriptive).** 420 seed runs had no near-native pose at E32.
- e128: 0.075 [0.032, 0.133] of 420 such runs find one.
- e512: 0.077 [0.000, 0.231] of 43 such runs find one.
- ladder e32 (seed 1, 68 copies): ceiling 0.348 [0.176, 0.535], top pose 0.078 [0.000, 0.203]
- ladder e128 (seed 1, 61 copies): ceiling 0.404 [0.221, 0.603], top pose 0.096 [0.004, 0.214]
- ladder e512 (seed 1, 68 copies): ceiling 0.373 [0.187, 0.567], top pose 0.077 [0.000, 0.202]

### Strata (ceiling and success)

| stratum | copies | ceiling E32 | ceiling E128 | top pose E128 |
|---|---|---|---|---|
| surface | 164 | 0.279 [0.162, 0.402] | 0.316 [0.196, 0.444] | 0.043 [0.010, 0.088] |
| semi_cryptic | 43 | 0.563 [0.372, 0.745] | 0.572 [0.382, 0.750] | 0.151 [0.069, 0.237] |
| cryptic (not evidence) | 11 | 0.810 [0.500, 1.000] | 0.839 [0.625, 1.000] | 0.345 [0.000, 0.750] |

Accounting: {"e32": {"records": 269, "errors": 19, "not_reached": 0}, "e128": {"records": 269, "errors": 51, "not_reached": 32}, "e512": {"records": 79, "errors": 11, "not_reached": 1}}

- G4 conditions on an E32 outcome, so it describes where a gain lands rather than evidencing it.
- Docking scores for a -9 polyanion from functions without explicit electrostatics are weak evidence whatever the search budget.
