# BioXArena-Pro Autonomous Research Project — Cross-Line Perturbation Modeling

> **Current Stage**: **Level 2** (`xatlas-orion-perturb`: internal v5 gain confirmed; Codabench submission pack frozen at commit `2a8b2e0` awaiting single upload by human PI)  
> **External Validation**: Frozen hyperparameters evaluated zero-shot on **Replogle et al. (*Cell* 2022) K562 & RPE1** (10 strict simultaneous union-held-out splits)  
> **Strict Disciplines**: Frozen solver hyperparameters, zero gene-identity lookup on xatlas, unfavorable results reported in main tables, row-level paired bootstrap 95% CIs (replacing invalid split-level paired $t$-tests).

---

## 核心实测摘要（含不利指标与外部模块迁移性）

1. **BioXArena-Pro `xatlas-orion-perturb`（HCT116 & HEK293T，官方 v5 残差评分器）**：
   - **10 个开发分割**：`100-SVD + Ridge (S0: 0.0613 ± 0.0040)` / `Ahlmann-Eltze 双线性基线 (0.0672 ± 0.0052)` → **`S7 (0.1099 ± 0.0087)`**（行级 95% CI：`[T3-CI 待填]`）。
   - **10 个额外随机重切分（与开发分割共享基因池，Seeds 101–110）**：`S0 (0.0596)` → `S7 (0.1073 ± 0.0072)`；超严格双细胞系同时留出 + 重拟合轴：`S0 (0.0658)` → `S7 (0.1115 ± 0.0089)`（行级 95% CI：`[T3-CI 待填]`）。
   - 自实现、未调参的简化神经网络基线在 10 个开发分割上得分：`ResMLP 0.0482`、简化单系 `GNN 0.0528`、简化跨系 `GNN 0.0538`。
2. **外部独立数据集零调参验证（Replogle et al., *Cell* 2022，K562 & RPE1，10 严格联合分割）**：
   - **v5 正交残差指标**：`S5 (0.1401 ± 0.0057)` 与 `S7 (0.1359 ± 0.0076)` 明显优于 `S0 (0.0883)`、`Ahlmann-Eltze 双线性基线 (0.0967)`、`30-kNN (0.0971)` 与自实现未调参简化基线（`ResMLP 0.1091`、`GNN 0.1233`）（行级 95% CI：`[T3-CI 待填]`）。
   - **传统未投影指标（S7 落后于线性基线）**：标准 Δ 表达 Pearson `S7 0.1162` 低于 `30-kNN 0.1479`、`AE 0.1470`、`S0 0.1408`；Top-50 DE Pearson `S7 0.2013` 低于 `S0 0.2119`（`S5` 为 `0.2170`）；Top-50 方向准确率 `S7 0.5962` 低于 `S0 0.6256`、`AE 0.6255`。
   - **模块迁移性**：外部稳定正迁移主要来自 **S1 (`+0.0356`)** 与 **S5 跨细胞系桥接 (`+0.0134`)**；**S2 (`-0.0024`)** 与 **S6 互惠响应借用 (`-0.0109`)** 在 Replogle 上为**负贡献**。

---

## 一键复现与依赖锁定

- **依赖锁定**：[requirements-lock.txt](requirements-lock.txt)（Python 3.13.5）
- **环境变量**（已彻底移除全部硬编码绝对路径）：
  - `BIOXARENA_DIR`: 指向本地 `BioXArena-Pro` 目录
  - `REPLOGLE_DATA_DIR`: 指向存放 Replogle 2022 官方 `.h5ad` 的目录
- **复现命令**：
  - 复现 xatlas 主表与提交核验：
    ```bash
    BIOXARENA_DIR=/path/to/BioXArena-Pro python3 src/xatlas_panhuman_solver.py
    ```
  - 复现 Replogle 2022 缓存（含 SHA-256 校验）与外部主表：
    ```bash
    REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/prepare_replogle2022_cache.py --verify-sha256
    REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode pan
    REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode dl
    ```
  - 复现行级配对自助法 95% CI（T3）：
    ```bash
    BIOXARENA_DIR=/path/to/BioXArena-Pro REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/compute_row_bootstrap_ci.py --mode xatlas
    ```

## 目录导航
- [学术手稿草稿 (docs/PAPER_DRAFT.md)](docs/PAPER_DRAFT.md) | [内部进展备忘 (docs/PRESS_RELEASE_CN.md)](docs/PRESS_RELEASE_CN.md)
- [成果验收台账 (ACCEPTANCE.md)](ACCEPTANCE.md) | [任务看板 (TASK_DISPATCH.md)](TASK_DISPATCH.md) | [通信日志 (COMMUNICATION_LOG.md)](COMMUNICATION_LOG.md)
- [Codabench 冻结清单 (CODABENCH_FREEZE_MANIFEST.json)](outputs/perturbation_baselines/xatlas-orion-perturb/CODABENCH_FREEZE_MANIFEST.json)
