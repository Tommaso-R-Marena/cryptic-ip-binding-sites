## The α-arrestin lead (docs/ARRESTIN_PLAN.md)

**B1 (family ranks high among unseen proteins):** supported - ROC-AUC 0.757, 95 % [0.683, 0.856], 5th percentile 0.694, 20 unseen α-arrestins in 14 clusters.

Protocol validity (crystal redocks of arrestin-IP sites): {'crystal_sites': 24, 'mean_top_pose_success': 0.0, 'valid': False}.

| protein | overlap | conservation | convergence | scores | verdict |
|---|---|---|---|---|---|
| ARRDC2 (Q8TBH0) | no | no | no | no | **not supported** |
| ART5 (P53244) | no | no | no | no | **not supported** |

## Dossiers

### ARRDC2 (Q8TBH0): **not supported**

Failing criteria: 1_overlap, 2_conservation, 3_convergence, 4_scores.

| criterion | pass | detail |
|---|---|---|
| 1_overlap | no | {"jaccard": 0.0, "tm_score": 0.64057, "reference": "5TV1_A.pdb", "detail": ""} |
| 2_conservation | no | {"decision": "not conserved", "conserved_positions": [], "basic_positions": []} |
| 3_convergence | no | {"decision": "not evaluable: protocol validity failed"} |
| 4_scores | no | {"decision": "not evaluable: protocol validity failed"} |

Mapped site (reference 5TV1_A.pdb, copy IHP_A_401; TM-score 0.64057): reference residue → ARRDC2 residue: 226→211, 227→212, 332→290.

Residues within 4 Å of the best docked IP6 pose: PRO214, VAL215, LEU216, ASP287, ILE288, PRO289, GLY290, THR291, LYS293.

What would test it: isothermal titration calorimetry or a fluorescence binding assay of the purified protein with IP6 (and ATP as a polyanion control); charge-reversal mutants (K/R to E) of the mapped basic residues, which should abolish binding if the site is real. Docking scores for a -9 polyanion from a scoring function without electrostatics are weak evidence whatever they say.

### ART5 (P53244): **not supported**

Failing criteria: 1_overlap, 2_conservation, 3_convergence, 4_scores.

| criterion | pass | detail |
|---|---|---|
| 1_overlap | no | {"jaccard": 0.0, "tm_score": 0.678, "reference": "7F1W_D.pdb", "detail": ""} |
| 2_conservation | no | {"decision": "not evaluable: 8 homologues (fewer than 10)", "conserved_positions": [], "basic_positions": []} |
| 3_convergence | no | {"decision": "not evaluable: protocol validity failed"} |
| 4_scores | no | {"decision": "not evaluable: protocol validity failed"} |

Mapped site (reference 7F1W_D.pdb, copy IHP_D_501; TM-score 0.678): reference residue → ART5 residue: 171→306.

Residues within 4 Å of the best docked IP6 pose: GLU15, GLU18, ARG19, ILE22, SER23, TYR24, PHE25, ILE303, LYS304, LYS305, PHE306, ASP486, LYS487.

What would test it: isothermal titration calorimetry or a fluorescence binding assay of the purified protein with IP6 (and ATP as a polyanion control); charge-reversal mutants (K/R to E) of the mapped basic residues, which should abolish binding if the site is real. Docking scores for a -9 polyanion from a scoring function without electrostatics are weak evidence whatever they say.

