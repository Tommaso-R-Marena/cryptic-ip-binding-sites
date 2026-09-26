## Does buried IP-site frequency track IP concentration? (docs/COEVOLUTION_PLAN.md)

Threshold: `combined` >= 0.519604 (the 1.0 % pooled false-positive anchor), pLDDT floor 70.0.

**C2 (primary, matched):** inconclusive - -0.017 [-0.341, 0.269] pp over 722 shared bins and 20080 clusters.

| organism | proteins | hits | hit rate |
|---|---|---|---|
| dictyostelium | 11679 | 106 | 0.908 % [0.702, 1.122] |
| human | 19899 | 243 | 1.221 % [1.018, 1.424] |
| yeast | 5806 | 85 | 1.464 % [1.113, 1.842] |

**C1 (raw, not the test):** -0.435 [-0.714, -0.150] pp.

**C4 (rule-only):** inconclusive - 0.183 [-0.081, 0.433] pp. Directions agree: False.

**C5 (deepest quartile):** inconclusive - 0.209 [-0.439, 0.906] pp.

| pLDDT decile | proteins | matched difference | decision |
|---|---|---|---|
| 70.0-75.0 | 3739 | 0.028 [-0.425, 0.532] pp | inconclusive |
| 75.0-79.3 | 3738 | 0.073 [-0.484, 0.678] pp | inconclusive |
| 79.3-83.2 | 3738 | 0.007 [-0.612, 0.616] pp | inconclusive |
| 83.2-86.3 | 3739 | 0.027 [-0.936, 1.078] pp | inconclusive |
| 86.3-88.9 | 3739 | 0.215 [-0.845, 1.262] pp | inconclusive |
| 88.9-91.3 | 3739 | 0.260 [-0.792, 1.589] pp | inconclusive |
| 91.3-93.3 | 3739 | 0.315 [-0.913, 1.325] pp | inconclusive |
| 93.3-95.2 | 3738 | 0.229 [-0.682, 1.185] pp | inconclusive |
| 95.2-97.0 | 3738 | -0.638 [-1.724, 0.439] pp | inconclusive |
| 97.0-98.9 | 3739 | -0.928 [-1.384, -0.526] pp | lower in dictyostelium |

- The primary test is C2, the matched comparison. C1 is reported only to show how much of the raw difference is model confidence, length and composition.
- Three organisms give two contrasts and no replication, so even a clean positive is a correlation across three points, not a demonstration of co-evolution.
- If C2 and C4 disagree in direction the primary is not robust to the learned model's training distribution and no co-evolution claim is made.
