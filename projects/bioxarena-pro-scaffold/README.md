# BioXArena-Pro Autonomous Research Project — Cross-Line Perturbation Modeling

> **Current Stage**: **Level 2** (`xatlas-orion-perturb` Official Codabench Private Test `Submission 1062221`: **`0.0998`**, leading within the `0.02` noise band; human-AI multi-round iteration, NOT comparable to single-shot 2h autonomous agents)  
> **External Validation**: Frozen hyperparameters evaluated zero-shot on **Replogle et al. (*Cell* 2022) K562 & RPE1** (10 strict simultaneous union-held-out splits; **S5 and S7 co-primary**)  
> **Round 6 Complete Report**: [reports/ROUND6_EXECUTION_REPORT.md](reports/ROUND6_EXECUTION_REPORT.md) (`ROUND6_DONE`)

---

## 核心实测摘要（全五指标 + 留一消融 + 诚实边界）

1. **官方 Codabench 私有测试集（唯一提交 `1062221`，冻结包 `2a8b2e0`）**：`perturbation/xatlas-orion-perturb = 0.0998`（较 22 个 2 小时自主 agent 最佳 `0.0825` 高 `+0.0173`，属 `0.02` 噪声范围内领先；实测验证→测试偏移 `−0.0101`；人机协作多轮迭代与 2 小时自主 agent 非同一比较类别）。
2. **残差 / 区分度 vs 传统未投影指标的权衡（T2 & T3）**：
   - **v5 残差指标**：xatlas 10 开发分割 `S7 = 0.1099`（vs `S0 0.0613` 行级中位 95% CI `[+0.0318, +0.0677]` 10/10；vs `AE 0.0672` `[+0.0256, +0.0606]` 10/10）；Replogle 10 严格联合分割 **`S5 = 0.1401`**（vs `AE 0.0967` `[+0.0323, +0.0549]` 10/10；vs 简化 GNN `[+0.0038, +0.0302]` 7/10）与 **`S7 = 0.1359`**（vs `AE` `[+0.0240, +0.0537]` 10/10；vs 简化 GNN `[-0.0041, +0.0286]` 仅 3/10）。
   - **独立于 v5 的原始 Δ 空间扰动区分度 (`Discrim(Delta)`)**：Replogle 上 **`S5 = 0.6890`**（vs `AE 0.6461`，`+0.0417`, `[+0.0198, +0.0721]`, 9/10 CI > 0），`S7 = 0.6762`（vs AE 仅 4/10 CI > 0）；xatlas 上 `S5 = 0.6262`、`S7 = 0.6286`（vs `AE 0.5856`）。
   - **传统未投影指标（落后于线性基线）**：标准 Δ 表达 Pearson（Replogle `S7 0.1162` / `S5 0.1193` vs `AE 0.1470`、`30-kNN 0.1479`，8/10 CI < 0）、Top-50 DE Pearson 与 Top-50 方向准确率上，S7 均低于线性基线。
3. **双数据集留一消融（T4）**：真正跨数据集稳定正贡献的是 **`S1_HybridTOM` (`+0.0109 / +0.0103`)**、**`S7_WienerSnap` (`+0.0072 / +0.0067` on v5，降低区分度)**、**`S5_CrossLineBridge` (`+0.0020 / +0.0061`)** 与 **`S2_Deflation` (`+0.0108 / +0.0035`)**；**`S6` 互惠响应借用**在 Replogle 上为强负贡献（`−0.0153`, 9/10 CI < 0，去掉 S6 后 Replogle v5 升至 `0.1512`），且在匿名 xatlas 上无法排除 CRISPRi 局部顺式抑制混杂（T6）。

---

## 一键复现与依赖锁定（T7）

- **依赖锁定**：[requirements-lock.txt](requirements-lock.txt)（Python 3.13.5；所有源码无硬编码绝对路径）
- **复现命令**：
  ```bash
  # xatlas 主表
  BIOXARENA_DIR=/path/to/BioXArena-Pro python3 src/xatlas_panhuman_solver.py
  # Replogle 2022 缓存与外部主表
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/prepare_replogle2022_cache.py --verify-sha256
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode pan
  REPLOGLE_DATA_DIR=/path/to/replogle2022 python3 src/replogle_external_validation.py --mode dl
  ```

## 目录导航
- [第六轮完整执行报告 (reports/ROUND6_EXECUTION_REPORT.md)](reports/ROUND6_EXECUTION_REPORT.md) | [阶段 A 检查点 (reports/ROUND6_PHASE_A.md)](reports/ROUND6_PHASE_A.md)
- [学术手稿草稿 (docs/PAPER_DRAFT.md)](docs/PAPER_DRAFT.md) | [内部备忘 (docs/PRESS_RELEASE_CN.md)](docs/PRESS_RELEASE_CN.md) | [验收台账 (ACCEPTANCE.md)](ACCEPTANCE.md)
