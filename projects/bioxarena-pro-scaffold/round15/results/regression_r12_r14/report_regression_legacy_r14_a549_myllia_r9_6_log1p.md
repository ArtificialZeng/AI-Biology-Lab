# Round 15 `PREREG_H` Report (`regression_legacy_r14_a549_myllia_r9_6_log1p`)

- **Wave ID**: `r14_a549_myllia` | **Status Label**: `interim, more lines pending`
- **Confirmatory Claim 1 (`H-E1`)**: **CONFIRMATORY_ABORTED_K0**
- **Confirmatory Claim 2 (`H-E2`)**: **CONFIRMATORY_ABORTED_K0**
- **Disclosure**: K_conf = 0 (CONFIRMATORY_ABORTED_K0): zero confirmatory claims tested; descriptive reporting only.

---

## Round 15 (`PREREG_H` v2) Governance, Multi-Line, & Wave Disclosures

- **Protocol**: `PREREG_H` v2 (`seed = 20261224`, `n_boot = 4000`, `n_perm = 100`)
- **Release Wave**: `r14_a549_myllia` (`interim, more lines pending`) | **Pending Candidate Lines**: `['cal33', 'hela', 'ac16', 'molm13', 'u87mg']` | **Excluded `jurkat`**: `True`
- **Mandatory Scope Notice (`§8.2`)**: K_conf = 0 (CONFIRMATORY_ABORTED_K0): zero confirmatory claims tested; descriptive reporting only.
- **Arc Challenge Curation Disclosure (`§4.2` / `§10.3`)**: Competition challenge panels filter targets by `>= 80%` median knockdown, select the stronger of 2 sgRNAs per gene, and uniformly downsample to `400` cells/perturbation.
- **Prior-Art Notice (`§10.6`)**: Inspect NeurIPS 2026 `OpenReview lnpzb1HvwX` and `CanPert/CASTER` boundaries (`PRIOR_ART_H.md` / `PRIOR_ART_H_CONTROLLER_FAN2026.md`) cleanly.

## `SA-H1` (Leave-`vcc2025_h1`-Out Source Sensitivity — `5` Remaining Sources, Strictly Descriptive)

| Target Line | `n_targets` | `FCT r_orth_cal` [97.5% CI] | `rho_half_ab_orth` | Paired `FCT - Mis_s` Ranges (Min/Max 97.5% CI lo) |
|---|---:|---|---:|---|
| `a549_myllia` | 28 | +0.1049 [+0.0714, +0.1371] | +0.3768 | [+0.0462, +0.0625] (`5` sources) |

## `SA-Panel` (`129` Shared Val-Test Targets vs Remainder — Strictly Descriptive)

| Target Line | Subgroup | `n_targets` | `rho_half_ab_orth` | `FCT r_orth_cal` [97.5% CI] |
|---|---|---:|---:|---|
| `a549_myllia` | `in_val_test_overlap_panel` | 4 | +0.5386 | nan [+0.0609, +0.1338] |
| `a549_myllia` | `panel_unique_remainder` | 28 | +0.3755 | nan [+0.0704, +0.1436] |
