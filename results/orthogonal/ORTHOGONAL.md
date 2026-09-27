## Orthogonal sequence evidence on the candidates (docs/ORTHOGONAL_PLAN.md)

**I1:** candidates more ordered - -0.093 [-0.149, -0.040] over 75 matched pairs in 65 clusters.

| role | n | scored | disordered-pocket rate |
|---|---|---|---|
| candidate | 75 | 75 | 0.015 [0.000, 0.046] |
| control | 75 | 75 | 0.111 [0.042, 0.181] |
| positive | 34 | 34 | 0.000 [0.000, 0.000] |

**I2 (disorder QC):** informative = True. annotated binders are rarely flagged, so a flagged candidate is demoted Demoted: P53244.

**I3 (motif support):** -0.132 [-0.348, 0.053] over 16 pairs in 15 clusters; 74 proteins gave no result. descriptive: orthogonal to study H's MAFFT criterion, not combined with it

**I4 (remote homology):** 0 of 0 unevaluable candidates have a remote homologue in the sampled proteomes. a within-proteome remote hit is not an orthologue set: study H's conservation criterion is not re-run on these

- The controls are already pLDDT-matched, so I1 asks whether sequence disorder carries information beyond AlphaFold's confidence. A null is a real answer.
- No sequence tool speaks to ligand identity; study C measured that and found the descriptors too weak to change a ranking.
