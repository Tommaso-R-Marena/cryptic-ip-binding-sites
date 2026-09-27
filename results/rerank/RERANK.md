## Electrostatic re-ranking of docked IP poses (docs/RERANK_PLAN.md)

250 copies in 31 strict groups, 750 seed runs.

**F1:** no detectable difference.

| estimand | Vina | re-ranked (out of fold) | difference |
|---|---|---|---|
| group | 0.093 [0.034, 0.176] | 0.147 [0.090, 0.213] | 0.054 [-0.001, 0.114] |
| copy | 0.091 [0.050, 0.129] | 0.125 [0.091, 0.190] | 0.035 [0.000, 0.104] |

Chosen w by fold: [0.5, 0.1, 0.1, 0.1, 0.1].

**F2:** per-run shares {"reranked": {"ok": 0.12533333333333332, "scoring_failure": 0.228, "sampling_failure": 0.6466666666666666}, "vina": {"ok": 0.09066666666666667, "scoring_failure": 0.26266666666666666, "sampling_failure": 0.6466666666666666}}; sampling ceiling 0.399 [0.286, 0.521].
**F3:** 100 permutations, mean 0.002, max 0.023, fraction ≥ observed 0.00.
- surface: 186 copies, 25 groups; Vina 0.030, re-ranked 0.099
- semi_cryptic: 51 copies, 14 groups; Vina 0.154, re-ranked 0.216
- cryptic: 13 copies, 4 groups; Vina 0.343, re-ranked 0.241 (fewer than 5 groups: not evidence)
- classic arrestins: 7 copies, 1 groups; Vina 0.000, re-ranked 0.048 (fewer than 5 groups: not evidence)

Accounting: {"records": 269, "errors": 19, "not_reached": 0, "error_examples": ["ReceptorError: pdb2pqr failed:   File \"/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/pdb2pqr/main", "ReceptorError: PDB2PQR outputs disagree: 137887 atoms, 137507 charges", "ReceptorError: pdb2pqr failed:     coords = [bondatom.coords, nextatom.coords] |                                ^^^^^^^^", "ReceptorError: pdb2pqr failed:   File \"/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/pdb2pqr/main", "ReceptorError: pdb2pqr failed:     coords = [bondatom.coords, nextatom.coords] |                                ^^^^^^^^", "ReceptorError: PDB2PQR outputs disagree: 141887 atoms, 141367 charges", "ReceptorError: pdb2pqr failed:   File \"/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/pdb2pqr/main", "ReceptorError: pdb2pqr failed:     newcoords = hatom.coords |                 ^^^^^^^^^^^^ | AttributeError: 'NoneType' ", "ReceptorError: PDB2PQR outputs disagree: 141708 atoms, 141148 charges", "ReceptorError: pdb2pqr failed:   File \"/opt/hostedtoolcache/Python/3.11.16/x64/lib/python3.11/site-packages/pdb2pqr/main"]}
Reproducibility vs the redocking primary arm: {"pairs": 750, "same_success": 0.9866666666666667, "within_0.5A": 0.892}
