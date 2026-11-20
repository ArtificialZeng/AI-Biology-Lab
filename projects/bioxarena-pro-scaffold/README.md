# BioXArena-Pro Autonomous Research Project — Cross-Line Perturbation Modeling

> **Current Stage**: **Level 2** (`xatlas-orion-perturb` Official Codabench Private Test `Submission 1062221`: **`0.0998`**, leading within the `0.02` noise band; human-AI multi-round iteration, NOT the same comparison category as single-shot 2h autonomous agents)  
> **External Validation**: Frozen hyperparameters evaluated zero-shot on **Replogle et al. (*Cell* 2022) K562 & RPE1** (10 strict simultaneous union-held-out splits; S5 and S7 reported as co-primary results)  
> **Strict Disciplines**: Frozen solver hyperparameters, zero gene-identity lookup on xatlas, unfavorable metrics and negative modules reported in main tables, row-level paired bootstrap 95% CIs (4,000 resamples).

---

## 核心实测摘要（含官方私测、不利指标与行级 95% CI）

1. **BioXArena-Pro `xatlas-orion-perturb`（官方 v5 残差评分器）**：
   - **官方 Codabench 私有测试实测（Submission ID `1062221`，唯一正式提交）**：**`0.0998`**（相对 22 个 2 小时自主单次运行 agent 最佳 `0.0825` 高 `+0.0173`，位于 `0.02` 官方噪声阈值内，属「噪声范围内领先」；实测验证→测试偏移 `−0.0101`；人机协作多轮迭代 ~25 个手调常数，与 2 小时自主 agent 非同一比较类别）。
   - **10 个开发分割**：`S0 (0.0613)` / `AE 双线性 (0.0672)` → `S7 (0.1099 ± 0.0087)`（`S7 vs S0` 行级中位 95% CI `[+0.0318, +0.0677]`, 10/10 CI > 0；`S7 vs AE` `[+0.0256, +0.0606]`, 10/10 CI > 0）。
   - **10 个额外随机重切分（与开发分割共享基因池）**：`S0 (0.0596)` → `S7 (0.1073 ± 0.0072)`；**超严格联合留出 + 重拟合轴**：`S0 (0.0658)` → `S7 (0.1115 ± 0.0089)`（行级中位 95% CI `[+0.0276, +0.0635]`, 10/10 CI > 0）。
   - 自实现、未调参简化基线（10 开发分割）：`ResMLP 0.0482`、简化单系 `GNN 0.0528`、简化跨系 `GNN 0.0538`。
2. **外部独立数据集零调参验证（Replogle et al., *Cell* 2022，K562 & RPE1，10 严格联合分割）**：
   - **v5 正交残差指标**：**`S5 (0.1401 ± 0.0057)`**（`S5 vs AE` 行级中位 95% CI `[+0.0323, +0.0549]`, 10/10 CI > 0；`S5 vs GNN` `[+0.0038, +0.0302]`, 7/10 CI > 0）与 **`S7 (0.1359 ± 0.0076)`**（`S7 vs AE` `[+0.0240, +0.0537]`, 10/10 CI > 0；`S7 vs GNN` `[-0.0041, +0.0286]`, 3/10 CI > 0）优于线性基线（`S0 0.0883`、`AE 0.0967`、`30-kNN 0.0971`）。
   - **传统未投影指标（S7 落后于线性基线）**：标准 Δ 表达 Pearson `S7 0.1162` 低于 `30-kNN 0.1479`、`AE 0.1470`、`S0 0.1408`；Top-50 DE Pearson `S7 0.2013` 低于 `S0 0.2119`（`S5` 为 `0.2170`）；Top-50 方向准确率 `S7 0.5962` 低于 `S0 0.6256`、`AE 0.6255`。
   - **模块迁移性**：外部稳定正迁移来自 **S1 (`+0.0356`, `[+0.0256, +0.0450]`, 10/10)** 与 **S5 跨细胞系桥接 (`+0.0134`)**；**S2 (`−0.0024`)** 与 **S6 互惠响应借用 (`−0.0109`, `[-0.0169, -0.0050]`, 9/10 CI < 0)** 为**负贡献**。

---

## 一键复现与依赖锁定

- **依赖锁定**：[requirements-lock.txt](requirements-lock.txt)（Python 3.13.5）
- **复现命令**（零硬编码绝对路径）：
  ```bash
  BIOXARENA_DIR=/path/to/BioXArena-Pro python3 src/xatlas_panhuman_solver.py
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/prepare_replogle2022_cache.py --verify-sha256
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode pan
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode dl
  ```

## 目录导航
- [第六轮阶段 A 检查点报告 (reports/ROUND6_PHASE_A.md)](reports/ROUND6_PHASE_A.md) | [第六轮执行报告 (reports/ROUND6_EXECUTION_REPORT.md)](reports/ROUND6_EXECUTION_REPORT.md)
- [学术手稿草稿 (docs/PAPER_DRAFT.md)](docs/PAPER_DRAFT.md) | [内部进展备忘 (docs/PRESS_RELEASE_CN.md)](docs/PRESS_RELEASE_CN.md) | [验收台账 (ACCEPTANCE.md)](ACCEPTANCE.md)
- [Codabench 冻结清单 (CODABENCH_FREEZE_MANIFEST.json)](outputs/perturbation_baselines/xatlas-orion-perturb/CODABENCH_FREEZE_MANIFEST.json)
