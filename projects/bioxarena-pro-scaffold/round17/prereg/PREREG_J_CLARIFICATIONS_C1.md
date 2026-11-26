# PREREG_J — binding clarifications C1 (controller)

Claude (master controller), 2026-11-27 07:44 +0800. Written before any Round 17 real-cell inference. These clarifications resolve the interpretation points the main line listed in `R17_MAIN_IMPL_READY` (07:41). They bind both parties: the main line's implementation, and A-T3's clean-room re-implementation. They change no hypothesis, estimator, threshold or line role in the hash-locked PREREG_J v2 (`PREREG_J_LOCK.md`, commit 9789b0b).

| # | Point | Binding ruling |
|---|---|---|
| C1.1 | Control-sampling seed | PREREG_J governs: primary control sampling, SE-600M pools, half-splits, bootstrap and permutations use seed **20261701**; the within-line resample uses **20261702**. The controller's dispatch text "seed 20261601" was a controller error and is withdrawn. Every Round 17 pool is sampled and encoded fresh, including the R16 lines used as mismatched sources. No R16 embedding is reused. |
| C1.2 | Gate files | Real-cell execution requires non-empty `round17/status/R17_CODE_FREEZE.txt`, `round17/status/R17_GO.txt` and `round17/prereg/PREREG_J_LOCK.md`. The lock record plays the role PREREG_J §5.2 gives `PREREG_J_COMMIT.txt`, and the controller also writes `round17/status/PREREG_J_COMMIT.txt` pointing to 9789b0b. |
| C1.3 | Treated TeloHAEC line key | `pfizer_telohaec_treated` (PREREG_J prose) and `pfizer_telohaec_cytokine` (locked reader, data directory) are the same line: Treated (IL-1β + TNF-α), a within-family secondary replication. |
| C1.4 | Stratified half-split, s = 0..19 | Sort the evaluated perturbations by canonical target symbol. With `rng = np.random.default_rng(20261701 + s)`, permute the half_valid = True indices, then the half_valid = False indices, and concatenate [perm_true, perm_false]. Positions 1, 3, 5, … go to P1 (⌊N/2⌋) and positions 0, 2, 4, … to P2 (⌈N/2⌉). Sort the indices within each half before within-half bootstrapping. Export the s = 0 index lists. |
| C1.5 | p-values for Holm (K_conf > 1) | H-J1 is a conjunction. Holm orders it by the intersection-union p-value p_H-J1 = max(p_var, p_Δ), where p_Δ = min(1, 2·(1 + #{Δ_CF^(u) ≤ 0}) / (B_valid + 1)). Rejection also requires the (1 − α_adj) percentile lower bound of Δ_CF,0 to be > 0. H-J3a/b use the two-sided bootstrap p = min(1, 2·(1 + min(#{d ≤ 0}, #{d ≥ 0})) / (B_valid + 1)), and EXACT_TIE_SELECTED_ZERO gives p = 1. Non-confirmatory lines use α = 0.05/3 per line and are not part of Holm. |
| C1.6 | V_noise and validity | Var_boot uses ddof = 1 over the joint-finite bootstrap rows (all 56 indices finite). p_var = (1 + #{S²_null ≥ S²_obs}) / (n_valid_rows + 1). Validity requires the per-index valid_frac ≥ 0.95 for every index and a joint valid_frac ≥ 0.95. |
| C1.7 | H-J1 "INCONCLUSIVE" outcomes | The validity-fail and insufficient-precision outcomes both report `INCONCLUSIVE`, with `subreason` set to `VALIDITY_FAIL` or `INSUFFICIENT_PRECISION`. `INCONCLUSIVE_LOW_SIGNAL` (J2) is separate. |
| C1.8 | H-J2 Bonferroni-56 companion | Per-index two-sided α = 0.025/56. Normal-calibrated: z = Φ⁻¹(1 − 0.025/112) (exact value; the prose's 3.5117 is that value rounded) × boot_se from the 4,000-draw bootstrap. Empirical: 20,000 draws from `np.random.default_rng(20261701).integers(0, N, size=(20000, N))`, quantiles 0.025/112 and 1 − 0.025/112. Both are reported. The primary H-J2 call uses the PREREG_I 97.5% rule only. |
| C1.9 | medoid_idx (secondary comparator) | Euclidean medoid of the 56 rows of the 4-fold mean `batch_encoder.weight`: argmin_b Σ_b′ ‖W̄_b − W̄_b′‖₂, ties to the lowest index. The cosine medoid is logged as a diagnostic only. |
| C1.10 | J5b permutation null | One shared (10,000 × 56) permutation table from `np.random.default_rng(20261701)` serves every line pair. a549_ko_cas9pos is excluded from the "all fresh lines with N ≥ 16" pair set (it is not a fresh line). It stays in the cross-model within-line comparison. |
| C1.11 | J5a auxiliary P2→P1 G4 | Run it, with no skip; its bootstrap seed is 20261701 + 10000. The primary J5a pair stays (b*₁, b†₁) selected on P1 and evaluated on P2. |
| C1.12 | Array layout and gene loading | Stacked per-model files (`preds_{MODEL}_sweep56.npz`, `preds_{MODEL}_g4_mismatched.npz`) as specified in `round17/src/R17_ARRAY_SCHEMA.md` are accepted. Both parties export in that schema, and the controller recomputes from it. Ensembl-only var names are canonicalized through the frozen HGNC bridge. `min_gene_overlap = 1` is a loading setting only. The PREREG_I floors that decide each line's role were applied at counts-only QC and are unchanged. |
| C1.13 | Run configuration | `--cohort all` (11 lines; wtc11_i3neuron omitted), `--n-boot 4000 --n-boot-bonf56 20000 --n-perm 10000 --n-splits-diag 20 --seed 20261701 --ctrl-seed 20261701 --ctrl-resample-seed 20261702 --device mps`, with the auxiliary G4 included. |

**Staged verification** (PREREG_J §5.2), as in Round 16:
- control barcodes exact;
- SE-600M embeddings ≤ 1e-5;
- per-index predicted Δ ≤ 1e-4;
- metrics ≤ 1e-5, with identical calls.

A-T3's arrays are compared with the main line's before any metric is compared.
