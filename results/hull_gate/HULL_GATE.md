## The screen's hull-depth gate (docs/HULL_GATE_PLAN.md)

35158 unseen proteins, 35 annotated IP binders, 20073 sequence clusters; K = 546 (hits under the 10 Å gate).

| arm | ROC-AUC [95 %] | recall@K [95 %] | hits | annotated hits |
|---|---|---|---|---|
| gate10 | 0.701 [0.613, 0.789] | 0.114 [0.023, 0.233] | 546 | 4 |
| gate5 | 0.738 [0.645, 0.817] | 0.114 [0.023, 0.233] | 554 | 4 |
| none | 0.747 [0.663, 0.821] | 0.114 [0.023, 0.233] | 554 | 4 |

- Δ ROC-AUC, no gate − 10 Å: 0.046 [-0.015, 0.101]
- Δ ROC-AUC, 5 Å − 10 Å: 0.037 [-0.018, 0.091]
- Δ recall@K, no gate − 10 Å: 0.000 [0.000, 0.000]
- Δ recall@K, 5 Å − 10 Å: 0.000 [0.000, 0.000]

**Decision: keep** (neither removing nor relaxing the gate is shown non-inferior).
