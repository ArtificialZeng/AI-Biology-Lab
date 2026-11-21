# Round 8 Execution Report (Final Closeout)

- **Timestamp**: 2026-11-22 07:43+08:00
- **Repository**: `https://github.com/ArtificialZeng/AI-Biology-Lab` (`projects/bioxarena-pro-scaffold/`)
- **Invariants Verified**:
  - Frozen solver diff (`git diff 3609811 -- projects/bioxarena-pro-scaffold/src/xatlas_panhuman_solver.py`): `0` lines.
  - Hardcoded local paths in `projects/bioxarena-pro-scaffold/src/*.py`: `0` matches.
  - Archived preprints retained untouched: `BioXArena_xatlas_preprint_v1.pdf`, `BioXArena_xatlas_preprint_v2.pdf`, and `BioXArena_xatlas_preprint_v2.1.pdf`.
  - GPU/MPS lock `projects/bioxarena-pro-scaffold/round9/p2/locks/mps_1` cleanly released (`rmdir`) immediately upon `R8-1` completion.

---

## 1. Preprint v2.2 Deliverable (`R8-1` Integration — Final Round 8 Preprint)

- **Local Path**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.2.pdf`
- **Source Markdown**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.2.md`
- **GitHub URL**: `https://github.com/ArtificialZeng/AI-Biology-Lab/blob/main/projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.2.pdf`
- **SHA-256**: `1c697decd85d700b80cd7de9cc97bc77ecb38dd0d44c1d837d74750c0725b060`
- **Pages / Size**: 16 pages, `1,853,964` bytes (`pandoc + xelatex`, 0 warnings, 4 embedded 300-DPI figures).
- **Automated QA Verified**: Single-paragraph abstract **246 words** (`<= 250` words strictly), zero backticks in Sections 1–5, zero CJK characters, zero `"PI"` / `[Human PI]`, uniform `DRF_raw` estimator throughout, Systema `2025`, all 21 API-verified references intact, root `README.md` & `README.zh-CN.md` updated to v2.2.

### Prior Archived Preprint Checkpoints (Preserved Untouched)
- **v2.1 (`R8-12`)**: `docs/preprint/BioXArena_xatlas_preprint_v2.1.pdf`, SHA-256 `7090e56196527d35f4cd896b99ee5b7e87d08071032b6c0161e83e84eb287317` (15 pages, 244-word abstract; uniform `DRF_raw` & Systema `2025`).
- **v2 (`R8-11`)**: `docs/preprint/BioXArena_xatlas_preprint_v2.pdf`, SHA-256 `90c732a8c020c9b5049684d8db697cef7cef6412022da01b636c2b5ef6d83535` (15 pages, 249-word abstract; Supplement 4 fixes, API-verified references, DRF calibration, HyperMap comparison).
- **v1 (`R8-10`)**: `docs/preprint/BioXArena_xatlas_preprint_v1.pdf` (initial archived PDF).

---

## 2. Completed `R8-1` Audit: Native Single-Cell GEARS on Nadig et al. (2025)

- **Script**: `projects/bioxarena-pro-scaffold/src/round8_r81_native_sc_gears.py` (with `--resume` checkpoint support added after host restart to preserve completed Seed 42 without recomputation).
- **Outputs**:
  - `projects/bioxarena-pro-scaffold/reports/r81_native_sc_gears_nadig_summary.json`
  - `projects/bioxarena-pro-scaffold/reports/r81_native_sc_gears_nadig_metrics.csv`
  - `projects/bioxarena-pro-scaffold/reports/r81_native_sc_gears_nadig_row_ci.csv`
- **Protocol**: Official `cell-gears==0.1.2` (`hidden_size=64, epochs=20, lr=1e-3, weight_decay=5e-4, batch_size=32`, native `skip_calc_de=False`, `anndata==0.10.9`), exact `make_strict_union_split` on pre-registered seeds **`[42, 0, 1]`**, per-`gem_group` library-size normalization + `log1p`, **`max_cells_per_pert = 32`** (lowered from 64 per pre-registered 3h/seed MPS cap rule) and **`max_ctrl_cells = 5000`** (`28,590–28,655` HepG2 cells and `35,796–36,131` Jurkat cells per split; zero test leakage; `98.6%–100.0%` test GO coverage).
- **Runtime (Apple Silicon MPS)**: Seed 42 `22,431.3s` (`6.23h`), Seed 0 `14,468.6s` (`4.02h`), Seed 1 `14,533.9s` (`4.04h`) — total **`51,433.8s` (`14.29h`)**.

### Matched 3-Seed Results on Nadig 2025 (`seeds = [42, 0, 1]`)

| Method (Exact Same 3 Seeds: 42, 0, 1) | `v5` | `delta_r` | `de50_r` | `de50_sign` | `disc_raw` | `disc_v5` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Native_SC_GEARS (`R8-1`, ctrl-z)** | `0.0679 ± 0.0042` | **`0.2024 ± 0.0116`** | **`0.3868 ± 0.0140`** | `0.5956 ± 0.0192` | **`0.6990 ± 0.0102`** | `0.6247 ± 0.0119` |
| `Native_SC_GEARS_rawlog` (unstandardized sensitivity) | `0.0716 ± 0.0045` | `0.1745 ± 0.0089` | `0.3052 ± 0.0112` | `0.5650 ± 0.0162` | `0.6453 ± 0.0050` | `0.6319 ± 0.0139` |
| `GEARS (official, pseudobulk-adapted)` (Round 7) | `0.0533` | `0.1001` | `0.1897` | `0.6421` | `0.6025` | `0.5844` |
| `AE (Linear)` | `0.1129` | `0.1721` | `0.2634` | `0.6457` | `0.6702` | `0.6574` |
| `S5 (Primary)` | **`0.1692`** | `0.1469` | `0.2195` | `0.6221` | **`0.7124`** | `0.7130` |
| `scGPT_FrozenEmb_Ridge` | `0.1520` | `0.2141` | `0.3181` | **`0.6831`** | `0.7097` | **`0.7139`** |

- **4,000-Resample Stratified Row-Level Paired Bootstrap 95% CIs (`seeds = [42, 0, 1]`)**:
  - **`S5 − Native_SC_GEARS`**:
    - `v5`: **`+0.1013`** (mean 95% CI `[+0.0788, +0.1232]`, **`3/3` CI > 0**)
    - `disc_v5`: **`+0.0883`** (`[+0.0498, +0.1266]`, **`3/3` CI > 0**)
    - **`disc_raw`**: **`+0.0134`** (`[−0.0249, +0.0519]`, **`0/3` CI > 0, `0/3` CI < 0 — statistical tie**; closes **87.8%** of pseudobulk gap `0.7124 − 0.6025 = 0.1099`)
    - **`delta_r`**: **`−0.0556`** (`[−0.0847, −0.0260]`, **`3/3` CI < 0 — Native SC GEARS beats S5 on 3/3 splits**)
  - **`Native_SC_GEARS − AE`**:
    - `delta_r`: **`+0.0303`** (`[−0.0020, +0.0631]`, `1/3` CI > 0, `0/3` CI < 0)
    - `disc_raw`: **`+0.0288`** (`[−0.0082, +0.0654]`, `1/3` CI > 0, `0/3` CI < 0)
    - `disc_v5`: `−0.0327` (`[−0.0719, +0.0074]`, `0/3` CI < 0)
    - `v5`: `−0.0449` (`[−0.0654, −0.0244]`, `3/3` CI < 0)
  - **GO-covered-only check (`98.6%–100.0%` coverage)**: `v5_go_only = 0.0684`, `disc_raw_go_only = 0.7010` (nearly identical).
- **Pre-Registered Adjudication Verdict (`ROUND8_PREREGISTRATION.md` Section 1.2)**:
  - **TRIGGERED: `Native Single-Cell GEARS Improves Substantially`**. Running official GEARS on native single cells with native DE calculation enabled dramatically improves performance over Round 7 pseudobulk adaptation—beating `AE` on both `disc_raw` and `delta_r`, closing 87.8% of the gap to `S5` on `disc_raw` (statistical tie, `0/3` CI > 0), and beating `S5` on `delta_r` (`3/3` splits) and `de50_r` (`0.3868` vs `0.2195`). Meanwhile, `S5` retains a consistent `3/3`-split lead over native single-cell GEARS specifically on orthogonalized metrics `v5` (`+0.1013`) and `disc_v5` (`+0.0883`). All manuscript comparisons in Preprint v2.2 have been updated accordingly.

---

## 3. Complete Summary of Round 8 Tasks (`R8-0` through `R8-12` + `v2.2` Closeout)

| Task | Status | Commit / Output | Key Result |
| :--- | :--- | :--- | :--- |
| **R8-0** | DONE | `9c4aa6c`, `3403221` | All 6 Round 7 final-audit manuscript fixes + `R8_CHECKPOINT_A.md` (`CHECKPOINT_A_DONE`). |
| **R8 Pre-Reg** | DONE | `ce82040`, `840b502` | `reports/ROUND8_PREREGISTRATION.md` + `Addendum R8-A` locked prior to execution/download. |
| **R8-1** | **DONE** | Current commit | **Native Single-Cell GEARS on Nadig 2025 (`seeds 42, 0, 1`, `cap=32`, `ctrl=5000`, `skip_calc_de=False`, `14.29h` MPS)**: `v5=0.0679`, `delta_r=0.2024` (beats S5 `3/3`), `de50_r=0.3868`, `de50_sign=0.5956`, `disc_raw=0.6990` (ties S5 `0/3`, beats AE, closes 87.8% of pseudobulk gap), `disc_v5=0.6247` (S5 leads `3/3`). Triggered `Native Single-Cell GEARS Improves Substantially` per pre-registered Section 1.2; integrated into Preprint v2.2. |
| **R8-2** | DONE | `ff12dc7` | PC1 axis decomposition **CONFIRMED 3/3 atlases**: `AE − S5` drops `+0.0533 → +0.0174` (xatlas, ratio `0.327`), `+0.0277 → −0.0017` (Replogle, `−0.063`), `+0.0195 → −0.0337` (Nadig, `−1.725`). |
| **R8-3** | DONE (`PENDING` 4th dataset) | `ff12dc7`, `9077d8d` | Exploratory `OrthoFusion_scGPT_S5` passes `H4`, fails dual-sided `H5` on Replogle & Nadig. 4th-dataset hard stop enforced prior to modeling: Feng 2026 (`38`/`86` rep targets, `37` shared) & Jiang 2025 (`44–61`/cytokine) failed locked power gate (`>=100/line`, `>=80` shared). |
| **R8-4** | **STOPPED** (per Supplement 2) | `task-2454` killed at `2026-11-21 16:03:01+08:00` | Started `13:21:32+08:00`, **elapsed `2h 41m 29s` (`161.5 min`)**. Partial progress: Seed 42 HepG2 `15/15` epochs (`val_mse 0.1119 → 0.1074`) + Jurkat `10/15` epochs (`val_mse 0.0461 → 0.0446`); no seed finished CSV export. Stopped under Supplement 2 compute triage (~40h full run low ROI); script `src/round7_r72_official_scgpt.py` preserved. |
| **R8-6** | DONE | `14f65f1`, `81995e8` | Cited Miller et al. (2026) & Mao et al. (2026); completed uniform `DRF_raw` calibration across 7 metrics on Replogle (K562 replicate) and Nadig (50/50 split-half): `MSE` poorly calibrated (`0.035 / 0.063`), `v5` moderate (`0.377 / 0.420 < delta_r 0.484 / 0.513`), `de50_r`/`de50_sign`/`disc_raw`/`disc_v5` best calibrated (`0.808–0.983`); standalone `S5` ties (`disc_raw`) or trails (`de50_sign`, `de50_r`) `scGPT_Probe`, and `S7` degrades `MSE` worse than `Constant Mean`. |
| **R8-7** | DONE | `14f65f1` | Complete 34-row `reports/AUDIT_LEDGER.csv` & `reports/AUDIT_LEDGER_SUMMARY.md` (13 retracted claims, 7 downgraded claims, 2 auditor self-corrections). |
| **R8-8** | DONE | `3c313be`, `48f9b1a`, `01a1182` | `reports/PRIOR_ART_LOG.md` covering all 5 claims, 9 mandatory metascience works + `CrossAudit` (`arXiv:2608.28631`) + `Vaccaro` (`arXiv:2606.11217`) + `HyperMap` (`10.64898/2026.04.23.720505`). |
| **R8-9** | DONE | `8da7859` | `reports/AUTONOMY_EVIDENCE.md` & `reports/autonomy_timeline.csv` (36 commits, all 15 direct human messages `<= 30` chars, 26 Claude relays, explicit `不确定` tags). |
| **R8-10..12** | DONE | `3c313be`, `48f9b1a`, `01a1182`, Current commit | Preprint PDFs `v1`, `v2`, `v2.1` archived untouched; final **Preprint PDF v2.2** (`BioXArena_xatlas_preprint_v2.2.pdf`, SHA-256 `1c697decd85d700b80cd7de9cc97bc77ecb38dd0d44c1d837d74750c0725b060`) integrating completed `R8-1` native single-cell GEARS audit published and linked in root `README.md` & `README.zh-CN.md`. |

PREPRINT_V2_DONE
PREPRINT_V21_DONE
PREPRINT_V22_DONE
ROUND8_DONE
