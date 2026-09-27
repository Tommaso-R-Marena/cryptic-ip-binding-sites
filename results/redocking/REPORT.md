## Redocking benchmark (docs/REDOCKING_PLAN.md)

Copies found: 662 in 367 entries; eligible 462; selected 272; primary set 269; incomplete (flagged stratum) 3.

| status | copies |
|---|---|
| eligible | 462 |
| excluded: configuration differs from the CCD | 167 |
| excluded: crystal artefact | 32 |
| excluded: crystal copy does not match the CCD template: Explicit valence for atom # 26 O, 3, is greater than permitted | 1 |

| arm group | records | docked | failed or not reached |
|---|---|---|---|
| alphafold | 269 | 269 | {} |
| pockets | 269 | 249 | {'ReceptorError: pdb2pqr failed:   File "/usr/share/miniconda/': 10, 'ReceptorError: pdb2pqr failed:     coords = [bondatom.coords': 2, 'ReceptorError: pdb2pqr failed:     newcoords = hatom.coords ': 2, 'RuntimeError: fpocket exceeded 2700 s': 1, 'ReceptorError: PDB2PQR outputs disagree: 137887 atoms, 13750': 1, 'ReceptorError: PDB2PQR outputs disagree: 141887 atoms, 14136': 1, 'ReceptorError: PDB2PQR outputs disagree: 141708 atoms, 14114': 1, 'ReceptorError: PDB2PQR outputs disagree: 133912 atoms, 13353': 1, 'ReceptorError: PDB2PQR outputs disagree: 136188 atoms, 13580': 1} |
| primary | 272 | 253 | {'ReceptorError: pdb2pqr failed:   File "/usr/share/miniconda/': 10, 'ReceptorError: pdb2pqr failed:     coords = [bondatom.coords': 2, 'ReceptorError: pdb2pqr failed:     newcoords = hatom.coords ': 2, 'ReceptorError: PDB2PQR outputs disagree: 136188 atoms, 13580': 1, 'ReceptorError: PDB2PQR outputs disagree: 137887 atoms, 13750': 1, 'ReceptorError: PDB2PQR outputs disagree: 141887 atoms, 14136': 1, 'ReceptorError: PDB2PQR outputs disagree: 141708 atoms, 14114': 1, 'ReceptorError: PDB2PQR outputs disagree: 133912 atoms, 13353': 1} |
| secondary | 272 | 253 | {'ReceptorError: pdb2pqr failed:   File "/usr/share/miniconda/': 10, 'ReceptorError: pdb2pqr failed:     coords = [bondatom.coords': 2, 'ReceptorError: pdb2pqr failed:     newcoords = hatom.coords ': 2, 'ReceptorError: PDB2PQR outputs disagree: 136188 atoms, 13580': 1, 'ReceptorError: PDB2PQR outputs disagree: 137887 atoms, 13750': 1, 'ReceptorError: PDB2PQR outputs disagree: 141887 atoms, 14136': 1, 'ReceptorError: PDB2PQR outputs disagree: 141708 atoms, 14114': 1, 'ReceptorError: PDB2PQR outputs disagree: 133912 atoms, 13353': 1} |

### Decisions (Holm across R1-R4)

| | question | estimate (group) | Holm p | decision |
|---|---|---|---|---|
| R1 | protocol reliability (group estimand) | 0.110 [0.046, 0.193] | 0 | **unreliable** |
| R2 | success(cryptic) - success(surface) | 0.299 [-0.056, 0.942] | – | **not evaluable: [4, 25] strict groups (cryptic, surface); need 5** |
| R3 | AlphaFold cross-docking success (group estimand) | 0.045 [0.000, 0.122] | 0 | **not trustworthy** |
| R4 | Vina score separates true site from decoy (ROC-AUC) | 0.744 [0.658, 0.873] | 0 | **discriminates** |

### Outcomes (primary set)

| outcome | copies | groups | per copy | per group |
|---|---|---|---|---|
| success | 250 | 31 | 0.096 [0.056, 0.141] | 0.110 [0.046, 0.193] |
| best20_success | 250 | 31 | 0.295 [0.236, 0.419] | 0.324 [0.221, 0.441] |
| success_1A | 250 | 31 | 0.012 [0.002, 0.025] | 0.016 [0.001, 0.039] |
| success_3A | 250 | 31 | 0.199 [0.135, 0.317] | 0.202 [0.115, 0.302] |
| p_success | 250 | 31 | 0.100 [0.059, 0.143] | 0.111 [0.046, 0.194] |
| spearman | 250 | 31 | 0.072 [0.013, 0.134] | 0.058 [-0.025, 0.146] |

### Top-pose success by stratum

| stratum | level | copies | groups | per copy | per group | |
|---|---|---|---|---|---|---|
| burial class | cryptic | 13 | 4 | 0.333 [0.000, 0.600] | 0.343 [0.000, 0.750] | fewer than 5 groups: not evidence |
| burial class | semi_cryptic | 51 | 14 | 0.150 [0.101, 0.216] | 0.175 [0.085, 0.285] |  |
| burial class | surface | 186 | 25 | 0.065 [0.008, 0.088] | 0.044 [0.010, 0.084] |  |
| metal | metal within 3 Å | 42 | 10 | 0.151 [0.076, 0.225] | 0.202 [0.067, 0.402] |  |
| metal | no metal | 208 | 28 | 0.085 [0.034, 0.131] | 0.079 [0.029, 0.142] |  |
| interface | interface | 58 | 13 | 0.029 [0.000, 0.056] | 0.041 [0.000, 0.103] |  |
| interface | single chain | 192 | 27 | 0.116 [0.065, 0.164] | 0.124 [0.050, 0.223] |  |
| species | InsP3 | 39 | 10 | 0.085 [0.000, 0.148] | 0.049 [0.000, 0.116] |  |
| species | InsP4 | 23 | 5 | 0.101 [0.000, 0.126] | 0.052 [0.000, 0.106] |  |
| species | InsP5 | 10 | 4 | 0.000 [0.000, 0.000] | 0.000 [0.000, 0.000] | fewer than 5 groups: not evidence |
| species | InsP6 | 154 | 19 | 0.106 [0.059, 0.176] | 0.114 [0.047, 0.199] |  |
| species | other | 24 | 13 | 0.083 [0.000, 0.208] | 0.115 [0.000, 0.269] |  |
| resolution | ≤ 2.0 Å | 85 | 18 | 0.106 [0.028, 0.159] | 0.077 [0.016, 0.165] |  |
| resolution | 2.0-2.5 Å | 64 | 19 | 0.109 [0.067, 0.174] | 0.140 [0.040, 0.268] |  |
| resolution | 2.5-3.0 Å | 57 | 14 | 0.105 [0.021, 0.179] | 0.072 [0.016, 0.139] |  |
| resolution | > 3.0 Å | 44 | 6 | 0.045 [0.015, 0.111] | 0.077 [0.007, 0.188] |  |
| method | X-ray | 191 | 28 | 0.115 [0.061, 0.162] | 0.117 [0.045, 0.210] |  |
| method | cryo-EM | 59 | 6 | 0.034 [0.000, 0.058] | 0.035 [0.000, 0.091] |  |

### Controls

```json
{
 "failure_decomposition": {
  "failed_seed1": 242,
  "scoring": 230,
  "sampling": 12
 },
 "seed_noise": {
  "top_score_sd": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.2326982968879921,
    "low": 0.20207922907278567,
    "high": 0.25172413913142017,
    "p_value": 0.0
   },
   "per_group": {
    "point": 0.2218528641703246,
    "low": 0.1845231859516681,
    "high": 0.26478592261642464,
    "p_value": 0.0
   }
  },
  "seed_agreement": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.772,
    "low": 0.6540517961570593,
    "high": 0.8776978417266187,
    "p_value": 0.001
   },
   "per_group": {
    "point": 0.8184163347210297,
    "low": 0.699263661965275,
    "high": 0.9152335429395474,
    "p_value": 0.0
   }
  }
 },
 "vinardo_success": {
  "arm": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.056,
    "low": 0.02100402661064426,
    "high": 0.10257575757575751,
    "p_value": 0.0
   },
   "per_group": {
    "point": 0.05623849725715321,
    "low": 0.006505828137706696,
    "high": 0.13178513297250025,
    "p_value": 0.0
   }
  },
  "arm_minus_primary_seed1": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.024,
    "low": 0.004607242808745239,
    "high": 0.06349567099567097,
    "p_value": 0.036
   },
   "per_group": {
    "point": 0.014240380394775721,
    "low": 0.00031938677738741617,
    "high": 0.034052608296882815,
    "p_value": 0.036
   }
  }
 },
 "ad4_success": {
  "arm": {
   "copies": 248,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.0967741935483871,
    "low": 0.05142857142857143,
    "high": 0.13793670598911068,
    "p_value": 0.0
   },
   "per_group": {
    "point": 0.09984775924054862,
    "low": 0.025480352453786995,
    "high": 0.1953032630451985,
    "p_value": 0.0
   }
  },
  "arm_minus_primary_seed1": {
   "copies": 248,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.06451612903225806,
    "low": 0.029406397595534563,
    "high": 0.10714594661393967,
    "p_value": 0.006
   },
   "per_group": {
    "point": 0.05782383334565498,
    "low": 0.00884819323623878,
    "high": 0.13537285295349805,
    "p_value": 0.006
   }
  }
 },
 "deprotonated_success": {
  "arm": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.044,
    "low": 0.015870927318295737,
    "high": 0.10001322751322744,
    "p_value": 0.0
   },
   "per_group": {
    "point": 0.0552321456054524,
    "low": 0.006801075268817205,
    "high": 0.12820274911741644,
    "p_value": 0.0
   }
  },
  "arm_minus_primary_seed1": {
   "copies": 250,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.012,
    "low": -0.014862193071768523,
    "high": 0.06802721088435375,
    "p_value": 0.587
   },
   "per_group": {
    "point": 0.013234028743074904,
    "low": -0.00294292102021262,
    "high": 0.035190679206169466,
    "p_value": 0.155
   }
  }
 },
 "metals_success": {
  "arm": {
   "copies": 42,
   "groups": 10,
   "evidence": true,
   "per_copy": {
    "point": 0.14285714285714285,
    "low": 0.09523809523809523,
    "high": 0.19753086419753085,
    "p_value": 0.0
   },
   "per_group": {
    "point": 0.17592592592592587,
    "low": 0.07407407407407407,
    "high": 0.30555555555555547,
    "p_value": 0.0
   }
  },
  "arm_minus_primary_seed1": {
   "copies": 42,
   "groups": 10,
   "evidence": true,
   "per_copy": {
    "point": 0.09523809523809523,
    "low": 0.018472222222222213,
    "high": 0.13581649831649822,
    "p_value": 0.027
   },
   "per_group": {
    "point": 0.050925925925925895,
    "low": -0.06020833333333332,
    "high": 0.14912037037037013,
    "p_value": 0.334
   }
  }
 },
 "site_finding": {
  "lands_in_true_site": {
   "copies": 244,
   "groups": 31,
   "evidence": true,
   "per_copy": {
    "point": 0.3073770491803279,
    "low": 0.19013415938669995,
    "high": 0.553687651331719,
    "p_value": 0.137
   },
   "per_group": {
    "point": 0.329732538176561,
    "low": 0.19672370561127678,
    "high": 0.47233535887473277,
    "p_value": 0.019
   }
  },
  "any_positive_pocket_fraction": 0.38934426229508196
 }
}
```

- The cryptic class holds few strict groups; any stratum marked 'evidence: false' has fewer than 5 independent groups and its interval is not evidence.
- Docking scores for a -9 polyanion from a scoring function without electrostatics are weak evidence.
