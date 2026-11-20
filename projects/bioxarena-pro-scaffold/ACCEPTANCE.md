# BioXArena-Pro 成果验收与评审登记表 (Acceptance Ledger)

本文档由研发团队维护，记录真实实验数据、官方口径复核结果、全流程消融与阶段定级。全部数字均由仓库内可执行脚本直接产出并落盘至 CSV/JSON（遵守 D1–D7 铁律）。

---

## 阶段定位与真实定级：L2（官方口径下单项显著突破官方基线，2.43× 噪声阈值）

按五级科研标准定位：
- L0: 代码跑通、产物合规
- L1: 复现了基准官方已公布的参考基线
- **L2: 官方口径下单项显著超过官方参考基线且超出 0.02 噪声阈值 ← 当前位置（已确认达成）**
- L3: 官方口径下单项超过排行榜 22 配置最佳且超出噪声（待 Codabench 私有 test 实测或内部 ≥0.1200 门槛确认）
- L4: Codabench 实测、Mean64 层面成立（顶刊/顶会）

> **D7 严格可比性声明**：
> 本文档所有增量比较均在**同一数据集（public train 10 个独立随机分割）、同一官方评分器（`tasks.utils.strict_group_mean_from_row_scores` + `grade.py` 正交投影逻辑）、且 100% 推理时可得输入（两个细胞系训练集中均零泄漏验证/测试靶基因）**下完成。内部 10 分割官方 SVD 基线为 `0.0613 ± 0.0040`，较官方文档公布参考值 `0.051` 系统性偏高约 `+0.010`；因此本结果**不直接宣称与私有 Codabench test 排行榜 SOTA（`0.0825`）数值等价**，L3 升级须以 Codabench 私有 test 实测为准。

---

## 一、主攻任务 `perturbation/xatlas-orion-perturb` 核心突破：`PanPerturb-Bridge`

可执行代码：[src/xatlas_panhuman_solver.py](src/xatlas_panhuman_solver.py)（与最终提交产物 `outputs/perturbation_baselines/xatlas-orion-perturb/solution.py` 完全一致）  
结构化结果：[reports/xatlas_10split_ablation.csv](reports/xatlas_10split_ablation.csv)、[reports/xatlas_strict_union_refit_audit.csv](reports/xatlas_strict_union_refit_audit.csv)、[outputs/perturbation_baselines/xatlas-orion-perturb/metrics.json](outputs/perturbation_baselines/xatlas-orion-perturb/metrics.json)

### 1. Benchmark A：10 分割递进消融实验（与 `xatlas_10split_bench.py` 完全同分割）

10 个独立随机分割分为两部分：
- **5 个开发分割（Dev-5）**：Seeds `[42, 0, 1, 2, 3]`
- **5 个未触碰留出确认分割（Untouched Held-5）**：Seeds `[7, 11, 99, 2024, 31337]`（超参冻结前完全未参与开发调优）

| 阶段 | 模型架构增量（细胞系无关统一超参） | 10 分割均值 ± SD | HCT116 均值 | HEK293T 均值 | Dev-5 均值 | **Untouched Held-5 均值** | 相对 S0 净增量 (Δ ± SD) | 正增益比例 | 配对 $t$ 检验 $p$ 值 |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **S0** | `100-SVD + Ridge (α=10)`（官方参考基线） | 0.0613 ± 0.0040 | 0.0767 | 0.0460 | 0.0606 | 0.0621 | +0.0000 ± 0.0000 | — | — |
| **S1** | `+ SelfMask + BaseSub + √n_cells + HybridTOM(p=1.5)` | 0.0677 ± 0.0068 | 0.0866 | 0.0487 | 0.0704 | 0.0649 | +0.0063 ± 0.0057 | 8/10 | $8.87 \times 10^{-3}$ |
| **S2** | `+ Multi-Scale Spectral Deflation Pyramid (d0..d10)` | 0.0722 ± 0.0061 | 0.0935 | 0.0509 | 0.0737 | 0.0707 | +0.0109 ± 0.0046 | **10/10** | $5.69 \times 10^{-5}$ |
| **S3** | `+ Pan-Human Cross-Line Gene Graph (4,040 共有基因)` | 0.0748 ± 0.0053 | 0.0941 | 0.0555 | 0.0767 | 0.0729 | +0.0135 ± 0.0048 | **10/10** | $1.46 \times 10^{-5}$ |
| **S4** | `+ Empirical-SNR Weighted Non-Reproducible Anchors` | 0.0811 ± 0.0038 | 0.1033 | 0.0589 | 0.0829 | 0.0793 | +0.0197 ± 0.0040 | **10/10** | $1.14 \times 10^{-7}$ |
| **S5** | `+ GeneSlope-CosKRR 跨系翻译器与对称双向桥接` | 0.0887 ± 0.0046 | 0.1138 | 0.0635 | 0.0906 | 0.0868 | +0.0273 ± 0.0040 | **10/10** | $7.03 \times 10^{-9}$ |
| **S6** | `+ 阈值互惠因果律 (ReLU(|z|-4)² 激活全部 6,018 可重复敲除)` | 0.1027 ± 0.0065 | 0.1371 | 0.0682 | 0.1039 | 0.1014 | +0.0413 ± 0.0067 | **10/10** | $1.85 \times 10^{-8}$ |
| **S7** | **`+ 输出 SVD Wiener 降噪 + 经验流形吸附 (FULL)`** | **0.1099 ± 0.0087** | **0.1462** | **0.0736** | **0.1134** | **0.1064** | **+0.0485 ± 0.0085** | **10/10** | **$3.63 \times 10^{-8}$** |

**单调性与无过拟合验证**：
- `S0 → S7` 连续 7 个算法模块在 **总体 10 分割**、**HCT116**、**HEK293T**、**Dev-5** 与 **Untouched Held-5** 五个维度上全部实现**严格单调正增长**。
- 在 5 个从未触碰的确认分割上，`S7` 取得 `0.1064 ± 0.0069`（相对 `S0` 的 `0.0621` 净增量 `+0.0443`，最小单分割增量 `+0.0333`，最大 `+0.0581`），每个分割均显著超过官方 `0.02` 噪声阈值。

### 2. Benchmark B：超严格双细胞系联合留出 + 重拟合 `context_baseline` & `residual_axis` 审计（10 分割）

为彻底排除任何潜在的跨细胞系间接信息或公开 `context_baseline / residual_axis` 估计泄漏，在每轮分割中：
1. 将 HCT116 与 HEK293T 验证靶基因并集（每轮约 355–365 个唯一基因）**同时从两个细胞系训练集彻底移除**；
2. 仅在剩余训练集可重复行上按官方文档公式严格重拟合 `context_baseline` 与 `residual_axis`（先经核验与官方文件数值一致至 `5.95e-10` / `cosine=1.00000000`）。

| 指标（超严格联合留出 + 重拟合轴，10 分割） | S0 (`100-SVD + Ridge`) | **S7 (`PanPerturb-Bridge`)** | 配对净增量 (Δ ± SD) | 正增益比例 | 配对 $t$ / $p$ 值 |
|:---|:---:|:---:|:---:|:---:|:---:|
| **10 分割总体均值 ± SD** | 0.0658 ± 0.0033 | **0.1115 ± 0.0089** | **+0.0457 ± 0.0081** | **10/10** | $t=17.01,\; p=3.77 \times 10^{-8}$ |
| Dev-5 均值 (Seeds 42, 0, 1, 2, 3) | 0.0653 | 0.1132 | +0.0479 | 5/5 | — |
| **Untouched Held-5 均值 (Seeds 7..31337)** | 0.0663 | **0.1097** | **+0.0434** | **5/5** | — |

### 3. 跨细胞系翻译器 5 折交叉验证独立审计（781 对共享可重复敲除）

脚本：[src/xatlas_translator_cv.py](src/xatlas_translator_cv.py) | 结果：[reports/xatlas_translator_5fold_cv.json](reports/xatlas_translator_5fold_cv.json)
- `HEK293T → HCT116`：原始共有基因复制 `0.2376 ± 0.1768` → `GeneSlope + Soft-Power(p=1.5) CosKRR` **`0.3783 ± 0.2635`**（+0.1407）
- `HCT116 → HEK293T`：原始共有基因复制 `0.2420 ± 0.1806` → `GeneSlope + Soft-Power(p=1.5) CosKRR` **`0.3551 ± 0.2431`**（+0.1131）
- （D7 注释：该指标度量“跨细胞系同靶基因翻译器精度”，不同于“未见靶基因零样本预测”指标，仅用于独立验证 S5 翻译器模块有效性，不与主任务排行榜数值混比。）

### 4. 官方 455 行测试集提交与 Preflight 核验

- 训练全量 100% 公开训练集，使用与消融阶段 `S7` 完全一致的 `PanPerturbBridgeModel` 类生成全部 455 个 `predictions/<id>.npz` 与 `submission.csv`（28.8s）。
- 官方 `training.submission_preflight.validate_public_submission` 校验 **100% 通过**。
- 独立产物检查：全部 455 个 `.npz` 均含有限 `float32` 向量 `(5000,)`，`max |y| = 1.1670`，最小正交残差比例 `0.08426`（远高于 `grade.py` 的 `1e-4` 退化门槛）。

---

## 二、历史 5 任务官方口径复核记录（维持 L1 基线复现）

1. **gutcrc-cross-cohort-detection**: 实测 `0.4453` vs 官方参考基线 `0.444`（持平，L1）
2. **cpjump-genetic-matching**: 实测 `0.6112` vs 官方参考基线 `0.613`（持平，L1）
3. **solubility**: 实测 `0.7150` vs 官方参考基线 `0.711`（持平，随机 KFold 未用骨架划分，L1）
4. **msk-met-tumor-origin**: 实测 `0.4443` vs 官方参考基线 `0.549`（低于基线，L1）
5. **gdsc-lineage-holdout-response**: 官方 `grade.py` 实测 `0.1199` vs 官方参考基线 `0.148`（低于基线，此前自写指标 `0.1942` 已彻底废止，L1）
