# Round 12 Clean-Room Audit Code — Old-Data Consistency & Blinding Report

**Author:** Agent `A-T3` (Team A · Pre-Registration & Statistical Audit)  
**Locked Protocol:** [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) (`commit 4c606f4732ea080768e27757c689d40d7892cd3d`)  
**Clean-Room Invariant:** Implemented 100% independently under `/Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/` without reading, importing, or inspecting any file inside `round12/src/`. Tested strictly on historical Round 9–11 standardized caches and synthetic fixtures.

---

## 1. Old-Data Consistency Verification (`round12/audit/` vs `round12/results/` Outputs, `tol = 1e-10`)

Verified via [`verify_vs_r12_results.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/audit/verify_vs_r12_results.py) at identical `seed = 20261201`:

| Cohort (`Phase`, `Source Mode`) | `n_boot / n_perm` | Point Checks (`max abs diff`) | CI97.5 Checks (`max abs diff`) | Verdicts | `<= 1e-10` Status |
|---|---|---:|---:|:---:|---|
| `vcc2025_h1` (`dev_r9`, `r9_6`) | `400 / 20` | `56` (**`5.55e-17`**) | `82` (**`2.78e-17`**) | `0` mismatch | **`PASS`** (`dev_smoke_post_r9_vcc`) |
| `feng_ipsc` (`dev_r11`, `r9_6`) | `400 / 20` | `59` (**`1.11e-16`**) | `88` (**`9.71e-17`**) | `0` mismatch | **`PASS`** (`dev_smoke_post_r11_feng`) |
| `calu3` (`dev_r10`, `r9_r10_8`, `6,101` joint genes, `9` targets) | `400 / 20` | `65` (**`2.22e-16`**) | `100` (**`1.11e-16`**) | `0` mismatch | **`PASS`** (`dev_smoke_post_src8_calu3`) |
| `feng_ipsc` (`dev_r11`, `r9_6`, full spec) | `4,000 / 100` | `59` (**`1.11e-16`**) | `88` (**`1.11e-16`**) | `0` mismatch | **`PASS`** (`dev_smoke_r11_feng`) |
| `cd4t_rest` (`dev_r10`, `r9_6`, pre-refactor `results/`) | `4,000 / 100` | `59` (**`1.73e-17`**) | `88` primary/non-degenerate (**`2.78e-17`**); `2` `GlobalScale r_orth_cal` `0.00462` | `0` mismatch | **`PRIMARY/N5 PASS`** (`dev_smoke_r10_cd4t_rest`)* |

*\*Note on pre-refactor `dev_smoke_r10_cd4t_rest` (`22:39` timestamp prior to `round12/results/dev_smoke_post_*`):* Every confirmatory & non-degenerate metric (`FCT`, `E1_ConstMed`, `TMVS`, `Hybrid_ProjShared`, `Source_Grand_Mean`, all `6` `H-E1b` paired `FCT - Mis_s`, `H-E2a E1(split)`, `H-E2b E1(split)-E1(all)`, `H-D3`, `Perm100`, and `N_5` including `N_5 Transfer_Mean_GlobalScale = 0.0`) agrees to **`<= 2.78e-17`**. Only the degenerate `0/0` correlation row `Transfer_Mean_GlobalScale` under `r_orth_cal` (`valid_frac=0.338` in pre-refactor output vs frozen R11 `_moments_stats` zero-vector convention `0.0` in `round12/audit/`, while yielding `r_orth_cal = ""` in all post-refactor `dev_smoke_post_*` runs) differed on CI bounds.

---

## 2. Unit Tests, Blinding Gate, & Synthetic Cohort Verification
- **`python3 -I -m pyflakes round12/audit/*.py`**: `PYFLAKES_OK` (`0` warnings/errors).
- **`python3 -I round12/audit/test_r12_audit.py -v`**: **`4/4 OK`**:
  - `test_1`: Hard runtime gate blocks `--phase holdout_r12` (`RuntimeError: PREREG_E hard gate blocked`) while `CODE_FREEZE_R12.txt` is absent; `dev_*` phases refuse Round 12 targets and `round12/shared/` directories; verifies `CONFIRMATORY_LINES_R12.txt` parser across `K_conf = 0/1/2` and `EXCLUDED_FEASIBILITY_Q3`.
  - `test_2`: Synthetic `a549_ko` reader verifies unique `cell_id` vs duplicate `cell_barcode`, `(genes x cells)` tab-MTX orientation, synthetic `Cas9` exclusion, `Condition == "Nontargeting"`, `Half A = batch1+batch3` vs `Half B = batch2+batch4`, `Lane` centering, and `SA-Cas9` (`exclude_zero_cas9=True`).
  - `test_3`: Synthetic `hek293_xucao` reader verifies `E3` pre-metric exclusion (`EXCLUDED_FEASIBILITY_Q3` on float-only matrix), CSC/CSR `layers/counts` ingestion, default `apply_author_kd40=False` (retaining `<40%` KD singlets), Stage-1 vs deterministic Stage-2 `seed=42` split, **Ruling A1** AnnData 0.2.0 `attrs["_index"]` resolution (`cell_barcode` / `gene_symbol` without `_index` dataset, `nperts=0` + `"control"` NTC retention with `"nan"` guide, and strict `RuntimeError` `MISSING_ANNDATA_INDEX` when real index metadata/datasets are absent — zero synthetic fallback), and **Ruling B** (`max_pct_mt = None`, `>55%`-MT cells retained while `<300` UMI or `<200` genes are filtered out).
  - `test_4`: Confirms `adjudicate_round12` (`K_conf = 0/1/2` + modality boundary notice), exact `<= 1e-12` parity vs frozen Round 11 `exploratory_core` / `run_n5_audit`, and joint `align_source_and_target_r12` gene intersection when an extension source omits target genes.
- **Post-patch (`RULING_R12_PREFREEZE` A1 + B & `RULING_R12_HEK293_UNITS` §1–2) checks:** Verified `dev_verify_r9_vcc` (`vcc2025_h1`, `r9_6`, `400/20`) against `round12/results/dev_smoke_post_r9_vcc` (`56` pt `5.55e-17`, `82` CI `2.78e-17`, `PASS`), plus post-`RULING_R12_HEK293_UNITS` `pyflakes` (`PYFLAKES_OK`), `test_r12_audit.py` (`4/4 OK`, including Case G `/X` vs `nascent_counts` + odd-stratum `perm[:n//2]` Stage-2 `"__NTC__"` split verification), `test_exploratory_check.py` (`4/4 OK`), `FREEZE_MANIFEST_R10.txt` (`22/22 OK`), and `FREEZE_MANIFEST_R11.txt` (`41/41 OK`).
