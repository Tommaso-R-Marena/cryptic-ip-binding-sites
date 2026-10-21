# Study Q6 method validation

Synthetic pools only. No annotated sequence is read here; this is the
evidence that the machinery finds an ordering signature when one is
planted and reports nothing when composition alone would explain the
difference. Statistic: `helix_charge_periodicity`, 400 shuffles per window.

| pool | n | mean raw | mean entropy | basic fraction | mean z | fires |
|---|---|---|---|---|---|---|
| planted helical face | 40 | 0.586 | 3.211 | 0.270 | +4.05 | 100% |
| same letters, shuffled | 40 | 0.172 | 3.211 | 0.270 | -0.24 | 0% |
| basic-rich, no ordering | 40 | 0.183 | 3.253 | 0.270 | -0.13 | 2% |

The three pools are matched on composition, so the separation in the z
column cannot be composition.

| check | result |
|---|---|
| shuffle control has identical entropy | pass |
| shuffle control has identical basic fraction | pass |
| composition control has identical basic fraction | pass |
| all pools within 0.1 bits | pass |
| all pools have the same basic fraction | pass |
| planted pool fires | pass |
| shuffle control stays quiet | pass |
| composition control stays quiet | pass |

All checks pass: **True**
