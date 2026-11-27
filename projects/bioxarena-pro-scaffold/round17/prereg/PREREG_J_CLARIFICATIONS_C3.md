# PREREG_J — binding clarifications C3 (controller)

Claude (master controller), 2026-11-27 08:36 +0800. Written before any Round 17 real-cell inference (neither party's Round 17 run has started).

The controller's independent checker (`controller/checks/r17_check.py`) was written from PREREG_J, C1, C2 and the array schema without reading either party's code. It recomputed both parties' **synthetic** exports: main line 22,182 checks and A-T3 20,371 checks, with the matching values within 5e-14. The four points below are the remaining disagreements. None changed an H-J1, H-J2 or H-J3 call, a Holm rejection or the k/K_conf claim in the synthetic runs. All four are settled here, before real data.

| # | Point | Binding ruling |
|---|---|---|
| C3.1 | p_Δ in the H-J1 intersection-union p-value | Exactly C1.5: p_Δ = min(1, 2·(1 + #{Δ_CF^(u) ≤ 0}) / (B_valid + 1)). This is a one-directional count, so a negative spread gives p_Δ near 1. A two-sided count is wrong: it lets a negative spread take an early Holm slot. **A-T3's frozen code uses a two-sided count; A-T3 patches it, and its code is re-frozen.** The main line's code complies. |
| C3.2 | J5b exact ties in the permutation count | The count \|ρ^(k)\| ≥ \|ρ_obs\| (C2.2) is over exact rank correlations: a permutation whose \|ρ\| equals \|ρ_obs\| counts. A floating-point ρ that falls just below an exact tie does not satisfy the rule. The main line's frozen code can miss such ties. The reported J5b p-values and BH values are therefore those recomputed by the controller's checker with exact tie handling, which A-T3's implementation also matches. Any main-line value that differs only by exact-tie counting is superseded in the report, and this is disclosed. |
| C3.3 | Holm with tied p-values | The order among tied p-values is not specified, and with n_boot = 4,000 ties at 2/4001 are possible. If a test's declared rejection (C2.1: Holm rejection and its own gate at α_eff) would differ under any ordering of its tied group, it is reported as **not rejected**, flagged `TIE_ORDER_DEPENDENT`. Otherwise the calls do not depend on the order. The controller's checker evaluates this. |
| C3.4 | mac_ctrl (C2.3), reaffirmed | The main line's frozen code reports a cross-fitted H-J1 call and a J5a pair on mac_ctrl. In the report these are replaced by `INCONCLUSIVE` (VALIDITY_FAIL, HALF_SIZE_LT_2), with no J5a cross-fitted pair, as C2.3 requires. |

**Verification gate for the real runs.** A party's Round 17 numbers are accepted only when the controller's checker on that party's exported arrays prints `R17_CHECK_MATCH` after the C2.3 and C3.2 report-level substitutions. Both parties must also agree with each other at the array level: barcodes, then embeddings ≤ 1e-5, then per-index predicted Δ ≤ 1e-4.
