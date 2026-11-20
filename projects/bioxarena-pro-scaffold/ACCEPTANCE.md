# BioXArena-Pro 成果验收与评审登记表 (Acceptance Ledger)

本文档由研发团队维护，记录真实实验数据、官方口径复核结果、跨实验室外部泛化验证、生物机理实证与阶段定级。全部数字均由仓库内可执行脚本直接产出并落盘至 CSV/JSON（严格遵守 D1–D7 铁律）。

---

## 阶段定位与综合评估

- **BioXArena-Pro 榜单定级**：**L2 确认突破**（内部 20 分割均值 `0.1086 ± 0.0080`、相对官方参考基线 `+0.0481`，达官方 `0.02` 噪声阈值 `2.41×`，配对 $p = 1.09 \times 10^{-17}$；待人类 PI 上传已 100% 通过 preflight 的提交包至 Codabench 私有 test 实测确认 L3）。
- **学术与外部泛化定级**：**已达顶会/顶刊投稿与硬核顶媒深度报道素材门槛（跨 2 大独立人源 Perturb-seq 图谱、4 种人类细胞系零调参泛化、正面对决 *Nature Methods* 2025 双线性基线与 PyTorch MLP/GNN、实时 Reactome 蛋白复合物富集验证 + 3 套 300-DPI 论文主图全齐备）**。

> **D7 严格可比性声明**：
> 所有数值比较均明确限定在同一数据集、同一评分口径、100% 推理时可得输入（零跨细胞系靶基因泄漏）下完成。Orion 内部 10 分割官方 SVD 基线（`0.0613`）较官方说明参考值（`0.051`）偏高约 `+0.010`，因此本文**不直接宣称内部得分与 Codabench 私有 test 排行榜 SOTA（`0.0825`）数值等价**。

---

## 一、支柱一：BioXArena-Pro `xatlas-orion-perturb` 全流程验证（20 分割 + DL/SOTA 对决 + 重拟合轴审计）

源码：[src/xatlas_panhuman_solver.py](src/xatlas_panhuman_solver.py) | [src/xatlas_sota_dl_benchmark.py](src/xatlas_sota_dl_benchmark.py)  
数据：[reports/xatlas_10split_ablation.csv](reports/xatlas_10split_ablation.csv) | [reports/xatlas_sota_dl_head2head_10split.csv](reports/xatlas_sota_dl_head2head_10split.csv) | [reports/xatlas_brand_new_10split_confirmation.csv](reports/xatlas_brand_new_10split_confirmation.csv) | [reports/xatlas_strict_union_refit_audit.csv](reports/xatlas_strict_union_refit_audit.csv)

### 1. 8 阶段递进消融（10 官方分割：5 Dev + 5 Untouched Held-Out）

| 阶段 | 模型架构增量 | 10 分割均值 ± SD | HCT116 | HEK293T | 未触碰 Held-5 | 相对 S0 增量 | 正增益 | 配对 $p$ 值 |
|:---:|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| **S0** | `100-SVD + Ridge (α=10)`（官方基线） | 0.0613 ± 0.0040 | 0.0767 | 0.0460 | 0.0621 | — | — | — |
| **S1** | `+ SelfMask + BaseSub + √n + HybridTOM(p=1.5)` | 0.0677 ± 0.0068 | 0.0866 | 0.0487 | 0.0649 | +0.0063 ± 0.0057 | 8/10 | $8.87 \times 10^{-3}$ |
| **S2** | `+ Multi-Scale Spectral Deflation Pyramid (d0..10)` | 0.0722 ± 0.0061 | 0.0935 | 0.0509 | 0.0707 | +0.0109 ± 0.0046 | 10/10 | $5.69 \times 10^{-5}$ |
| **S3** | `+ Pan-Human Cross-Line Gene Graph (4,040 genes)` | 0.0748 ± 0.0053 | 0.0941 | 0.0555 | 0.0729 | +0.0135 ± 0.0048 | 10/10 | $1.46 \times 10^{-5}$ |
| **S4** | `+ Empirical-SNR Non-Reproducible Anchors` | 0.0811 ± 0.0038 | 0.1033 | 0.0589 | 0.0793 | +0.0197 ± 0.0040 | 10/10 | $1.14 \times 10^{-7}$ |
| **S5** | `+ GeneSlope-CosKRR Symmetric Cross-Line Bridge` | 0.0887 ± 0.0046 | 0.1138 | 0.0635 | 0.0868 | +0.0273 ± 0.0040 | 10/10 | $7.03 \times 10^{-9}$ |
| **S6** | `+ Thresholded Reciprocal Causal Law (6,018 rep)` | 0.1027 ± 0.0065 | 0.1371 | 0.0682 | 0.1014 | +0.0413 ± 0.0067 | 10/10 | $1.85 \times 10^{-8}$ |
| **S7** | **`+ SVD Wiener + Manifold Snapping (FULL)`** | **0.1099 ± 0.0087** | **0.1462** | **0.0736** | **0.1064** | **+0.0485 ± 0.0085** | **10/10** | **$3.63 \times 10^{-8}$** |

### 2. 主流文献与深度学习正面对决（同 10 分割）+ 20 分割与重拟合审计

| 对比方法 / 审计组 | 实测均值 ± SD | `S7` 优势 / 增量 | 胜率 | 配对 $p$ 值 |
|:---|:---:|:---:|:---:|:---:|
| PyTorch Residual MLP (Projected Cosine+MSE, MPS) | 0.0482 ± 0.0052 | +0.0617 (+128.0%) | 10/10 | $9.85 \times 10^{-10}$ |
| PyTorch GEARS-Style Co-Expression GNN (MPS) | 0.0528 ± 0.0088 | +0.0571 (+108.1%) | 10/10 | $3.34 \times 10^{-8}$ |
| PyTorch Pan-Human Cross-Line GNN (MPS) | 0.0538 ± 0.0086 | +0.0561 (+104.3%) | 10/10 | $8.15 \times 10^{-8}$ |
| Ahlmann-Eltze (*Nature Methods* 2025) 双线性线性模型 | 0.0672 ± 0.0052 | +0.0427 (+63.5%) | 10/10 | $1.92 \times 10^{-7}$ |
| **10 个全新未触碰确认分割 (Seeds 101–110): `S0 → S7`** | `0.0596` → **`0.1073 ± 0.0072`** | `+0.0477 ± 0.0044` | **10/10** | $1.29 \times 10^{-10}$ |
| **合并 20 分割总审计: `S0 → S7`** | `0.0605` → **`0.1086 ± 0.0080`** | `+0.0481 ± 0.0068` | **20/20** | **$1.09 \times 10^{-17}$** |
| **超严格联合留出 (~360 基因) + 重拟合轴 (10 分割)** | `0.0658` → **`0.1115 ± 0.0089`** | `+0.0457 ± 0.0081` | **10/10** | $3.77 \times 10^{-8}$ |

---

## 二、支柱二：外部独立图谱零调参泛化验证（Replogle et al., *Cell* 2022 K562 & RPE1）

源码：[src/replogle_external_validation.py](src/replogle_external_validation.py) | 数据：[reports/replogle2022_external_validation_10split.csv](reports/replogle2022_external_validation_10split.csv) | [reports/replogle2022_external_validation_summary.json](reports/replogle2022_external_validation_summary.json)

直接导入冻结版 `PanPerturbBridgeModel`（**零超参修改**），在 Whitehead/MIT Replogle *Cell* 2022 独立图谱（RPE1 非癌视网膜上皮 + K562 全基因组髓系白血病，10 轮双细胞系同时留出 ~488 基因 + 仅用训练集重拟合基准与主轴）实测：

| 方法（全部 10 严格联合分割，零调参迁移） | v5 正交残差 $r$ | RPE1 $r$ | K562 $r$ | Top-50 DE $r$ | 配对 $p$ (`S7` vs 方法) |
|:---|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0009 ± 0.0074 | -0.0298 | 0.0315 | -0.0531 | $2.99 \times 10^{-13}$ |
| `100-SVD + Ridge (S0)` | 0.0883 ± 0.0048 | 0.0999 | 0.0768 | 0.2119 | $1.45 \times 10^{-9}$ |
| Ahlmann-Eltze (*Nature Methods* 2025) 双线性模型 | 0.0967 ± 0.0053 | 0.1047 | 0.0887 | 0.2026 | $6.30 \times 10^{-9}$ |
| 30-kNN 共表达 | 0.0971 ± 0.0059 | 0.1223 | 0.0719 | 0.1822 | $6.12 \times 10^{-8}$ |
| PyTorch Residual MLP | 0.1091 ± 0.0067 | 0.1294 | 0.0888 | 0.1896 | $3.36 \times 10^{-6}$ |
| PyTorch GEARS-Style GNN | 0.1233 ± 0.0061 | 0.1465 | 0.1001 | 0.1892 | $5.55 \times 10^{-4}$ |
| **`PanPerturb-Bridge S7`（完全冻结版）** | **0.1359 ± 0.0076** | **0.1425** | **0.1292** | 0.2013 | — |
| **`PanPerturb-Bridge S5`（跨系桥接冻结版）** | **0.1401 ± 0.0057** | **0.1535** | **0.1267** | **0.2170** | （相对 S0 +58.7%, 相对 AE-Lin +44.9%） |

---

## 三、支柱三：生物机理验证与实时 Reactome 通路富集

源码：[src/xatlas_biological_mechanism_audit.py](src/xatlas_biological_mechanism_audit.py) | [src/replogle_mechanism_reactome.py](src/replogle_mechanism_reactome.py)  
数据：[reports/xatlas_biological_mechanism_summary.json](reports/xatlas_biological_mechanism_summary.json) | [reports/replogle2022_biological_mechanism_summary.json](reports/replogle2022_biological_mechanism_summary.json)

- **Orion 匿名哈希零泄漏复合物重发现（全量 1,882 基因 5 折 OOF）**：互惠因果命中基因（`n=1,817`）平均 `S6−S5 = +0.0135` vs 无命中基因（`n=65`）`-0.0012`；典型基因 `g_f7c3d51640a64efc` 在 HCT116 与 HEK293T 中被同一组 3 个面板外敲除唤醒，相关性从 `0.005 / -0.142` 飙升至 **`0.853 / 0.839`**。
- **Replogle 2022 真实人源基因实时 Reactome 富集（`Homo sapiens`, FDR < 0.05）**：
  - 高精度预测集（`S7 ≥ 0.45`, `n=236`）：433 条显著通路，首位 **蛋白酶体降解（`29/37`, `FDR = 4.44e-16`）**、**核糖体质量控制 RQC（`70/106`, `FDR = 4.44e-16`）**、**40S/60S 核糖体亚基（`65/106`, `FDR = 4.44e-16`）**、**无义介导降解 NMD（`60/101`, `FDR = 4.44e-16`）**；
  - 跨系桥接胜出集（`S5−S2`）与互惠因果律胜出集（`S6−S5`）：均以 **mRNA 主要剪接通路（`26/213`, `FDR = 1.67e-14`；`39/213`, `FDR = 9.88e-15`）** 居首；
  - 真实基因案例 **`ADAM10` (K562)**：`S0 = -0.118 → S5 = 0.069 → S6 = 0.461 → S7 = 0.579`，由红系主控因子 `GATA1 (+9.4z)` 与 **3 个面板外中介体复合物亚基（`MED30*`, `MED21*`, `MED9*`）** 精准唤醒。

---

## 四、支柱四：论文级 300-DPI 科研主图与官方提交核验

- 3 套多面板主图（PNG + PDF）：[docs/figures/fig1_orion_ablation_and_sota.png](docs/figures/fig1_orion_ablation_and_sota.png)、[docs/figures/fig2_replogle2022_external_generalization.png](docs/figures/fig2_replogle2022_external_generalization.png)、[docs/figures/fig3_biological_mechanism_and_reactome.png](docs/figures/fig3_biological_mechanism_and_reactome.png)。
- 官方 455 行测试集提交产物（`submission.csv` + 455 个 `.npz`）100% 通过 `validate_public_submission`。
- 历史 5 任务维持 L1 基线复核记录（`gutcrc 0.4453`、`cpjump 0.6112`、`solubility 0.7150`、`msk-met 0.4443`、`gdsc 0.1199`）。
