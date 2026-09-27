## Triaging the candidates by pocket conservation (docs/TRIAGE_PLAN.md)

Roles: {"candidate": 75, "control": 75, "positive": 34}

**H1:** enriched - 0.646 [0.485, 0.808] over 38 matched pairs in 33 clusters.

| role | n | evaluable | conserved-basic rate |
|---|---|---|---|
| candidate | 75 | 48 | 0.925 [0.825, 1.000] |
| control | 75 | 39 | 0.365 [0.216, 0.514] |
| positive | 34 | 26 | 0.714 [0.524, 0.905] |

**H2 (filter calibration):** annotated binders conserved-basic 0.714 [0.524, 0.905]; informative = True. the filter keeps most annotated binders, so it can triage

**H3:** {"explained": 30, "conserved basic pocket": 24, "not evaluable: 2 homologues (fewer than 10)": 4, "not evaluable: 6 homologues (fewer than 10)": 3, "not evaluable: 0 homologues (fewer than 10)": 2, "not evaluable: 9 homologues (fewer than 10)": 2, "not evaluable: 7 homologues (fewer than 10)": 2, "not evaluable: 1 homologues (fewer than 10)": 2, "not conserved": 2, "not evaluable: 5 homologues (fewer than 10)": 1, "not evaluable: 4 homologues (fewer than 10)": 1, "not evaluable: 3 homologues (fewer than 10)": 1, "not evaluable: 8 homologues (fewer than 10)": 1}

| organism | rank | protein | homologues | basic positions | conserved |
|---|---|---|---|---|---|
| human | 1 | ARRDC2 (Q8TBH0) | 119 | 8 | 4 |
| human | 2 | RLBP1 (P12271) | 360 | 9 | 9 |
| human | 4 | GRTP1 (Q5TC63) | 22 | 10 | 8 |
| human | 8 | PROM2 (Q8N271) | 145 | 5 | 5 |
| human | 13 | CSTF1 (Q05048) | 112 | 8 | 6 |
| human | 15 | WDR46 (O15213) | 62 | 8 | 8 |
| human | 20 | SLC12A1 (Q13621) | 420 | 5 | 5 |
| human | 21 | RRP12 (Q5JTH9) | 500 | 5 | 5 |
| human | 22 | ELMOD2 (Q8IZ81) | 419 | 8 | 6 |
| human | 23 | PRLHR (P49683) | 114 | 4 | 4 |
| human | 24 | TBC1D3H (P0C7X1) | 34 | 6 | 5 |
| human | 25 | TBC1D3L (B9A6J9) | 34 | 5 | 5 |
| human | 26 | ATCAY (Q86WG3) | 112 | 5 | 5 |
| yeast | 5 | COP1 (P53622) | 500 | 17 | 13 |
| yeast | 6 | MUB1 (Q03162) | 19 | 6 | 6 |
| yeast | 7 | NMD5 (P46970) | 45 | 4 | 4 |
| yeast | 8 | IMG1 (P25626) | 32 | 9 | 6 |
| yeast | 9 | VSB1 (P53273) | 43 | 5 | 5 |
| yeast | 15 | MSN5 (P52918) | 46 | 4 | 4 |
| yeast | 17 | NFS1 (P25374) | 107 | 7 | 7 |
| yeast | 18 | UTP30 (P36144) | 39 | 5 | 4 |
| yeast | 20 | CSC1 (Q06538) | 53 | 6 | 6 |
| yeast | 21 | KAP95 (Q06142) | 49 | 5 | 5 |
| yeast | 24 | PEX11 (Q12462) | 37 | 9 | 9 |

- A conserved basic pocket is not evidence of binding: study C showed these descriptors do not separate IP from other polyanions well enough to change a ranking.
- The candidate list was read before this plan was written, so this is a filter, not a test of the screen.
