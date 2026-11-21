# Round 8 Execution Report (Preprint v2 Checkpoint)

- **Timestamp**: 2026-11-21 16:48+08:00
- **Repository**: `https://github.com/ArtificialZeng/AI-Biology-Lab` (`projects/bioxarena-pro-scaffold/`)
- **Invariants Verified**:
  - Frozen solver diff (`git diff 3609811 -- projects/bioxarena-pro-scaffold/src/xatlas_panhuman_solver.py`): `0` lines.
  - Hardcoded local paths in `projects/bioxarena-pro-scaffold/src/*.py`: `0` matches.
  - Preprint PDF v1 retained untouched at `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v1.pdf`.

---

## 1. Preprint v2 Deliverable (`R8-11`)

- **Local Path**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.pdf`
- **Source Markdown**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.md`
- **GitHub URL**: `https://github.com/ArtificialZeng/AI-Biology-Lab/blob/main/projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.pdf`
- **SHA-256**: `90c732a8c020c9b5049684d8db697cef7cef6412022da01b636c2b5ef6d83535`
- **Pages / Size**: 15 pages, `1,844,079` bytes (`pandoc + xelatex`, 0 warnings, 4 embedded 300-DPI figures).
- **All 7 Mandatory Supplement 4 Fixes Verified**:
  1. **Abstract length**: Single paragraph, **249 words** (`<= 250` words strictly), covering problem → 3 contributions → Codabench `0.0998` (with noise threshold and comparison category) → H1–H3 → PC1 decomposition → boundary vs. frozen scGPT probe → DRF conclusion (including unfavorable `de50_sign`/`de50_r` finding).
  2. **Comparison category & AI governance wording**: Updated throughout to *"multi-round iteration by a dual-AI executor–auditor system under a human-set goal, not comparable to single-shot 2-hour autonomous agents"*; Section 2.4 explicitly states human goal-setting, manual forwarding of Claude's critiques in Rounds 1–3, automated system relay from Round 6 onward, and Codabench Submission 1062221 submitted by the auditor on the human's behalf under human authorization.
  3. **API-verified references (`reports/REFERENCE_CHECK.csv`)**: All 21 references verified directly via Crossref API, arXiv API, bioRxiv API, NCBI PubMed E-utilities, ICML 2026 official virtual page (`https://icml.cc/virtual/2026/73543`), and alphaXiv HTTP 200 (`https://www.alphaxiv.org/abs/2611.society-ai-scientists`). Separated **Miller et al. (2026)** (`doi:10.1038/s41587-026-03307-w`; `bioRxiv 10.1101/2025.10.20.683304`) from **Mao et al. (2026)** (`arXiv:2604.27646`); replaced all placeholder author entries with exact API-verified authors (`Dong & Chen 2026`, `Dhaka, Gao & Ideker 2026`, `Zhu & Li 2026`, `Luo et al. 2026`, `Ahmad, Anderson, Qu et al. 2026`); standardized `Viñas Torné et al. (2026)` single year, `Adduri et al. (2025)` bioRxiv DOI (removed unverified NeurIPS), and exact PubMed-verified titles/years for Nadig (2025), Feng (2026), and Jiang (2025); removed unverified authors/workshop for ICML 2026 PCF and unverified institutions for A Society of AI Scientists.
  4. **HyperMap direct comparison**: Positioned **HyperMap** (Dhaka, Gao & Ideker, *bioRxiv* `10.64898/2026.04.23.720505`) as the closest prior work to `S5` in Section 1 (Introduction), Section 4 (Discussion), and `reports/PRIOR_ART_LOG.md`.
  5. **Integrated completed `R8-6` Miller DRF calibration**: Added Section 3.3 + Table 3 reporting DRF and rankings across all 7 metrics on Replogle 2022 and Nadig 2025, prominently reporting unfavorable findings for standalone `S5` (`de50_sign`, `de50_r`, `MSE`) and `S7` (`MSE` worse than `Constant Mean`). Kept `R8-1` (native single-cell GEARS) accurately marked as **in progress** in Section 5.
  6. **Unified figure & table GEARS labels**: Regenerated all 4 300-DPI figures via `src/round7_r75_generate_figures.py` replacing `"Official GEARS"` with `"GEARS (official package, pseudobulk-adapted)"` in Figure 1 legend, Figure 2B y-axis, and Table 1.
  7. **Root `README.md` & `README.zh-CN.md`**: Updated BioXArena-Pro project links to Preprint v2 while preserving `BioXArena_xatlas_preprint_v1.pdf`.

---

## 2. Summary of Round 8 Tasks Status (`R8-0` through `R8-11`)

| Task | Status | Commit / Output | Key Result |
| :--- | :--- | :--- | :--- |
| **R8-0** | DONE | `9c4aa6c`, `3403221` | All 6 Round 7 final-audit manuscript fixes + `R8_CHECKPOINT_A.md` (`CHECKPOINT_A_DONE`). |
| **R8 Pre-Reg** | DONE | `ce82040`, `840b502` | `reports/ROUND8_PREREGISTRATION.md` + `Addendum R8-A` locked prior to execution/download. |
| **R8-1** | **RUNNING** | `task-2432` | Native single-cell GEARS on Nadig 2025 (`seeds 42, 0, 1`, `cap=32`, `ctrl=5000`, `skip_calc_de=False`, `anndata==0.10.9`). Currently Seed 42 HepG2 Epoch 18/20. |
| **R8-2** | DONE | `ff12dc7` | PC1 axis decomposition **CONFIRMED 3/3 atlases**: `AE − S5` drops `+0.0533 → +0.0174` (xatlas, ratio `0.327`), `+0.0277 → −0.0017` (Replogle, `−0.063`), `+0.0195 → −0.0337` (Nadig, `−1.725`). |
| **R8-3** | DONE (`PENDING` 4th dataset) | `ff12dc7`, `9077d8d` | Exploratory `OrthoFusion_scGPT_S5` passes `H4`, fails dual-sided `H5` on Replogle & Nadig. 4th-dataset hard stop enforced prior to modeling: Feng 2026 (`38`/`86` rep targets, `37` shared) & Jiang 2025 (`44–61`/cytokine) failed locked power gate (`>=100/line`, `>=80` shared). |
| **R8-4** | **STOPPED** (per Supplement 2) | `task-2454` killed at `2026-11-21 16:03:01+08:00` | Started `13:21:32+08:00`, **elapsed `2h 41m 29s` (`161.5 min`)**. Partial progress: Seed 42 HepG2 `15/15` epochs (`val_mse 0.1119 → 0.1074`) + Jurkat `10/15` epochs (`val_mse 0.0461 → 0.0446`); no seed finished CSV export. Stopped under Supplement 2 compute triage (~40h full run low ROI); script `src/round7_r72_official_scgpt.py` preserved. |
| **R8-6** | DONE | `14f65f1`, `81995e8` | Cited Miller et al. (2026) & Mao et al. (2026); completed DRF calibration across 7 metrics on Replogle (K562 replicate) and Nadig (50/50 split-half): `MSE` poorly calibrated (`DRF 0.0631 / 0.0354`), `disc_v5`/`disc_raw` (`0.973–0.989`) and `de50_sign`/`de50_r` (`0.703–0.902`) well calibrated; standalone `S5` loses to `scGPT_Probe` & `AE` on `de50_sign`. |
| **R8-7** | DONE | `14f65f1` | Complete 34-row `reports/AUDIT_LEDGER.csv` & `reports/AUDIT_LEDGER_SUMMARY.md` (13 retracted claims, 7 downgraded claims, 2 auditor self-corrections). |
| **R8-8** | DONE | `3c313be`, updated in v2 | `reports/PRIOR_ART_LOG.md` covering all 5 claims, 9 mandatory metascience works + `CrossAudit` (`arXiv:2608.28631`) + `Vaccaro` (`arXiv:2606.11217`) + `HyperMap` (`10.64898/2026.04.23.720505`). |
| **R8-9** | DONE | `8da7859` | `reports/AUTONOMY_EVIDENCE.md` & `reports/autonomy_timeline.csv` (36 commits, all 15 direct human messages `<= 30` chars, 26 Claude relays, explicit `不确定` tags). |
| **R8-10** | DONE | `3c313be` | Preprint PDF v1 archived (`BioXArena_xatlas_preprint_v1.pdf`). |
| **R8-11** | DONE | `48f9b1a` | Preprint PDF v2 (`BioXArena_xatlas_preprint_v2.pdf`, SHA-256 `90c732a8c020c9b5049684d8db697cef7cef6412022da01b636c2b5ef6d83535`) + `reports/REFERENCE_CHECK.csv` + regenerated figures. |
| **R8-12** | DONE | Current commit | **Preprint PDF v2.1** (`BioXArena_xatlas_preprint_v2.1.pdf`, SHA-256 `7090e56196527d35f4cd896b99ee5b7e87d08071032b6c0161e83e84eb287317`, 15 pages, 244-word single-paragraph abstract): uniform `DRF_raw` estimator throughout abstract/prose (removed mixed `DRF_interp` `0.703`), explicitly reports `v5` (`DRF_raw 0.38–0.42 < delta_r 0.48–0.51`) and states that S5's sole clear advantage over the frozen scGPT probe is confined to a moderately calibrated metric while tying (`disc_raw`) or trailing (`de50_sign`, `de50_r`) on the best-calibrated metrics (`0.81–0.98`); full 7-metric DRF table ordered by `DRF_raw` in Section 3.3; Systema year corrected to `2025` across manuscript, `REFERENCE_CHECK.csv`, and `PRIOR_ART_LOG.md`; root README updated to v2.1 (`v1` and `v2` preserved). |

---

## 3. Preprint v2.1 Deliverable (`R8-12`)

- **Local Path**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.1.pdf`
- **Source Markdown**: `projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.1.md`
- **GitHub URL**: `https://github.com/ArtificialZeng/AI-Biology-Lab/blob/main/projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.1.pdf`
- **SHA-256**: `7090e56196527d35f4cd896b99ee5b7e87d08071032b6c0161e83e84eb287317`

PREPRINT_V2_DONE
PREPRINT_V21_DONE
