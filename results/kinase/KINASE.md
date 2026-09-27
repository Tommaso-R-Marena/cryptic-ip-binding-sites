## The IP kinases and the missing cosubstrate (docs/KINASE_PLAN.md)

33 of the primary-set copies carry a nucleotide cofactor within 6.0 Å of the ligand: {"ADP": 15, "ATP": 9, "ANP": 8, "ACP": 1}.

**K0 (descriptive):** the IPK superfamily is 0 of 18 on top-pose success — EhIP6KA, IPMK, ITPKA, ITPKC, PPIP5K2. Failure kinds: {"scoring": 18}.

**K1 (primary):** not evaluable: 4 strict groups (fewer than 5) — -0.012 [-0.068, 0.031] over 32 copies in 4 strict groups, Holm p 0.655.

**K2 (ceiling):** not evaluable: 4 strict groups (fewer than 5) — -0.010 [-0.031, 0.000].

| arm | ceiling | Vina top pose | re-ranked |
|---|---|---|---|
| apo | 0.542 [0.354, 0.854] | 0.273 [0.000, 0.751] | 0.363 [0.068, 0.796] |
| holo | 0.531 [0.354, 0.844] | 0.260 [0.000, 0.750] | 0.358 [0.102, 0.781] |

| ligand stratum | copies | top-pose success |
|---|---|---|
| PP-IP (InsP7 + InsP8) | 6 | 0.000 [0.000, 0.000] |
| InsP6 | 154 | 0.114 [0.043, 0.196] |

**K4 (electrostatics on the holo receptor):** 0.098 [0.030, 0.161]; 0 copies carried cofactor charges. no copy carried cofactor charges, so this repeats study G's re-ranking on a receptor whose cofactor is uncharged and is not evidence about electrostatics

**Metals only (quoted from study A):** mean 0.16666666666666666 over 12 copies. study A's metals-only arm on the same copies, quoted not re-docked, so the cofactor's contribution can be separated from the metals'

- The IPK slice of study A was read before this plan was written, so K0 and K3 are a re-analysis of data in hand, not a test. K1, K2 and K4 are blind.
- Vina's scoring function does not read partial charges, so in K1 and K2 the cofactor is a shaped, typed occluder: it restores the site's shape and hydrogen bonding, not its electrostatics.
- The holo arm adds the cofactor and the catalytic metals together. Study A's metals-only arm is quoted on the same copies so the two contributions are not confounded.
- A win here would narrow where the protocol may be trusted; it would not rescue the proteome ranking, whose candidate pockets are not catalytic sites.
