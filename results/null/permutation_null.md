## Permutation-control null (diagnostic D1)

Run 35965940997 (commit ad7013f), on the table of full benchmark run 35949588200.
Each row is `benchmark.py run --permute --arm full --grouping sequence --n-draws 10`
with a different repeat (permutation and fold seed).

| task | permutations | pooled ROC-AUC mean ± sd [2.5, 97.5 pct] | > 0.60 | within-fold mean ± sd | verdict |
|---|---|---|---|---|---|
| burial | 30 | 0.497 ± 0.051 [0.406, 0.600] | 1 | 0.500 ± 0.067 | **chance** |
| cryptic_ip_site | 10 | 0.498 ± 0.049 [0.413, 0.567] | 0 | 0.507 ± 0.060 | **chance** |

Pooled ROC-AUC per permutation, repeats 0 upward:

- burial: 0.613, 0.447, 0.595, 0.47, 0.543, 0.511, 0.452, 0.437, 0.507, 0.554, 0.411, 0.5, 0.508, 0.567, 0.45, 0.485, 0.532, 0.562, 0.464, 0.483, 0.542, 0.47, 0.485, 0.504, 0.505, 0.463, 0.513, 0.393, 0.468, 0.476
- cryptic_ip_site: 0.506, 0.565, 0.568, 0.495, 0.494, 0.469, 0.526, 0.396, 0.484, 0.478

Repeat 0 is the full run's own permutation. It reproduces that run's 0.613
exactly, and it is the single most extreme of the 30.
