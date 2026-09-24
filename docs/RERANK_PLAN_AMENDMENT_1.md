# Electrostatic re-ranking: amendment 1

Dated 2026-09-24. Written before any pose for study F was generated or scored. No
redocking outcome had been read either. It amends `docs/RERANK_PLAN.md`; everything
not named here stands.

## What prompted it

The unit tests run cross-fitting on synthetic pose lists. In one of them Vina always
ranks a wrong pose first, and E_el is pure noise. There F1 still decided **improves**:
any perturbation of a systematically wrong ranking lifts success above zero. So F1 on
its own can credit electrostatics with a gain that any added noise would produce. In
the original plan the permutation control (F3) was only reported. It has to be part
of the decision.

## Changes

1. **F3 runs 100 permutations** (seeds 1–100) instead of 20. The permutation p-value is
   (1 + #{permutations with a difference ≥ the observed one}) / 101.
2. **F1 decision:**
   - **Improves**: the 95 % lower bound of the group difference is > 0 **and** the
     permutation p-value is < 0.05.
   - **Gain not specific to electrostatics**: the lower bound is > 0 but the permutation
     p-value is ≥ 0.05. The re-ranking helps, but shuffled energies help as much.
   - **Worsens**: the 95 % upper bound is < 0.
   - **No detectable difference**: otherwise.
   - **Not evaluable**: fewer than 5 groups.

The F1 point estimate is still compared with the permutation distribution, as F3
already specified. Holm across the study's primary tests stays trivial: F1 is one
decision with two conditions.
