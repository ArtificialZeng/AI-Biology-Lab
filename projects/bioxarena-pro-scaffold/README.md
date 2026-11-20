# BioXArena-Pro Autonomous Research Project — `PanPerturb-Bridge`

> **Project Lead**: Antigravity  
> **BioXArena-Pro Milestone**: **Level 2 Confirmed Breakthrough** (`xatlas-orion-perturb`: 20-split `+0.0481 ± 0.0068`, `20/20`, $p = 1.09 \times 10^{-17}$; pending live Codabench upload for L3)  
> **External Generalization**: **Zero-Shot Confirmed on Replogle et al. (*Cell* 2022) K562 & RPE1** (frozen hyperparameters: `S5 = 0.1401`, `S7 = 0.1359` vs `SVD+Ridge 0.0883`, Ahlmann-Eltze *Nature Methods* 2025 Linear `0.0967`, `GEARS-GNN 0.1233`)  
> **Disciplines**: Official grader only (D1), reproducible code for every number (D2), full progressive ablation with exact test-model parity (D3), preflight ≠ score (D4), zero leakage (D5), honest staging (D6), strict D7 comparability.

---

## 核心突破一览（4 个人类细胞系 · 2 大独立图谱 · 实时 Reactome 机理验证）

1. **BioXArena-Pro `xatlas-orion-perturb`（HCT116 + HEK293T）**：
   - **20 分割总审计**：`100-SVD + Ridge (0.0605)` → **`PanPerturb-Bridge S7 (0.1086 ± 0.0080)`**（净增量 **`+0.0481`**, `20/20` 全胜, $p = 1.09 \times 10^{-17}$）。
   - **正面对决主流文献与深度学习（同 10 分割）**：碾压 `PyTorch ResMLP (0.0482)`、`GEARS-GNN (0.0528)`、`Pan-Human GNN (0.0538)` 与 **Ahlmann-Eltze (*Nature Methods* 2025) 双线性线性模型 (`0.0672`)**（`+63.5%`, $p = 1.92 \times 10^{-7}$）。
   - **超严格双细胞系联合留出 + 重拟合轴审计（10 分割）**：`0.0658` → **`0.1115 ± 0.0089`**（`+0.0457`, `10/10`, $p = 3.77 \times 10^{-8}$）。
2. **外部独立图谱零调参迁移 — Replogle et al. (*Cell* 2022) K562 & RPE1**：
   - 零超参微调直接迁移，10 轮双细胞系严格联合留出：`S0 (0.0883)` / `Ahlmann-Eltze Linear (0.0967)` / `ResMLP (0.1091)` / `GEARS-GNN (0.1233)` → **`PanPerturb-Bridge S7 (0.1359)` / `S5 (0.1401)`**（`10/10` 显著优于全部基线）。
3. **实时 Reactome 生物机理验证**：
   - 在真实人源基因上高显著富集 **蛋白酶体降解（`29/37`, `FDR = 4.44e-16`）**、**核糖体质量控制（`70/106`, `FDR = 4.44e-16`）**、**mRNA 剪接（`39/213`, `FDR = 9.88e-15`）**；真实案例 **`ADAM10` (K562)** 由 `GATA1` 与 3 个面板外中介体亚基（`MED30*/21*/9*`）从 `-0.118` 精准还原至 `0.579`。

## 快速导航
- **[完整学术手稿 (docs/PAPER_DRAFT.md)](docs/PAPER_DRAFT.md)**：含完整方法论、4 张核心定量表与 3 套 300-DPI 论文主图。
- **[300-DPI 论文主图 (docs/figures/)](docs/figures/)**：`fig1_orion_ablation_and_sota.png/pdf`、`fig2_replogle2022_external_generalization.png/pdf`、`fig3_biological_mechanism_and_reactome.png/pdf`。
- **[成果验收登记表 (ACCEPTANCE.md)](ACCEPTANCE.md)** | **[任务看板 (TASK_DISPATCH.md)](TASK_DISPATCH.md)** | **[作战总纲 (MISSION_CHARTER.md)](MISSION_CHARTER.md)**。
- **[一键复现源码 (src/)](src/)**：
  - `src/xatlas_panhuman_solver.py`（主算法 + 10 分割消融 + 超严格重拟合审计 + 官方 455 行提交）
  - `src/xatlas_sota_dl_benchmark.py`（DL/SOTA 正面对决 + 10 个全新确认分割）
  - `src/replogle_external_validation.py`（Replogle *Cell* 2022 零调参外部泛化）
  - `src/xatlas_biological_mechanism_audit.py` & `src/replogle_mechanism_reactome.py`（OOF 机理归因 + 实时 Reactome 富集）
  - `src/generate_publication_figures.py`（300-DPI 论文主图生成）
