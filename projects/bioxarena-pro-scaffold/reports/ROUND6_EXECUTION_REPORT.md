# 第六轮执行报告（Round 6 Execution Report — 全部任务 T0–T7 完成）

**日期**：2026-11-21  
**主控**：Claude | **执行**：Antigravity  
**仓库作用域**：`AI-Biology-Lab/projects/bioxarena-pro-scaffold`  
**核心求解器状态**：超参 100% 冻结（`src/xatlas_panhuman_solver.py` 仅在 commit `0e7d3a7` 修改第 35 行 `BIOXARENA_DIR` 路径解析，未改动任何算法常数或逻辑；Seed 42 复跑 `S0=0.070151`、`S7=0.126020` 与 455 测试 `.npz` 最大差值 `0.0`）。

---

## 一、任务完成总表（T0–T7）

| 任务 | 状态 | 核心 Commit | 执行脚本 | 结果文件 | 复现命令 |
|:---:|:---:|:---:|:---|:---|:---|
| **T0** | **完成** | `2a8b2e0` | `src/xatlas_panhuman_solver.py` | `outputs/.../CODABENCH_FREEZE_MANIFEST.json`, `submission.csv.sha256` | `shasum -a 256 -c outputs/perturbation_baselines/xatlas-orion-perturb/submission.csv.sha256` |
| **T1** | **完成** | `5d96040`, `fe4a265`, 本次提交 | 人工审计 + `grep` | `docs/PRESS_RELEASE_CN.md`, `docs/PAPER_DRAFT.md`, `README.md`, `ACCEPTANCE.md`, `COMMUNICATION_LOG.md` | `grep -nE "p ?=|untouched|从未触碰|never-seen|\+0\.010|Causal Law|因果律|Pan-Human|全人源|全面击败|全面领跑|线性壁垒" docs/PRESS_RELEASE_CN.md docs/PAPER_DRAFT.md README.md ACCEPTANCE.md COMMUNICATION_LOG.md` |
| **T2** | **完成** | 本次提交 | `src/round6_t2_multi_metric_bench.py`, `src/round6_t2_discrim_bootstrap_ci.py` | `reports/t2_all_5metric_combined.csv`, `reports/t2_all_5metric_summary.json`, `reports/t2_discrim_delta_row_ci.csv`, `reports/t2_discrim_delta_row_ci_summary.json` | `BIOXARENA_DIR=... REPLOGLE_DATA_DIR=... python3 src/round6_t2_multi_metric_bench.py --mode replogle_pan`（依次跑 `replogle_dl`, `xatlas_pan`, `xatlas_dl`, `summarize` 与 `round6_t2_discrim_bootstrap_ci.py`） |
| **T3** | **完成** | `fe4a265` | `src/compute_row_bootstrap_ci.py` | `reports/row_bootstrap_ci_xatlas.csv`, `reports/row_bootstrap_ci_replogle.csv`, `reports/row_bootstrap_ci_summary.json` | `BIOXARENA_DIR=... REPLOGLE_DATA_DIR=... python3 src/compute_row_bootstrap_ci.py --mode xatlas`（及 `replogle_gnn`, `replogle`） |
| **T4** | **完成** | 本次提交 | `src/round6_t4_leave_one_out_ablation.py` | `reports/t4_lomo_xatlas.csv`, `reports/t4_lomo_replogle2022.csv`, `reports/t4_lomo_summary.json` | `BIOXARENA_DIR=... REPLOGLE_DATA_DIR=... python3 src/round6_t4_leave_one_out_ablation.py` |
| **T5** | **完成** | `fe4a265` | `src/round6_t5_local_reactome_enrichment.py` | `reports/t5_local_reactome_enrichment.csv`, `reports/t5_local_reactome_enrichment_summary.json`, `reports/t5_enrichment_inputs.json` | `REACTOME_GMT_PATH=... python3 src/round6_t5_local_reactome_enrichment.py` |
| **T6** | **完成** | `fe4a265` | `src/round6_t6_crispri_proximity_test.py` | `reports/t6_replogle_crispri_proximity_genes.csv`, `reports/t6_replogle_crispri_proximity_summary.json` | `REPLOGLE_DATA_DIR=... GRCH38_GTF_PATH=... python3 src/round6_t6_crispri_proximity_test.py` |
| **T7** | **完成** | `0e7d3a7` | `src/prepare_replogle2022_cache.py` | `requirements-lock.txt`, `reports/ROUND6_PHASE_A.md` | `REPLOGLE_DATA_DIR=... python3 src/prepare_replogle2022_cache.py --verify-sha256` |

---

## 二、官方测试结果（Codabench 唯一正式提交）

- **提交信息**：人类 PI 授权下由主控使用 T0 冻结包（commit `2a8b2e0`，`submission.csv` SHA-256 `e2613a9b6616368ba5ca29edf132fe023761f0f72072b25176724207daa6a207`，455 个 `.npz` 聚合 SHA-256 `42532c39...`）单次提交至官方 Codabench。
- **Submission ID**：`1062221`（配置名 `panperturb-bridge__human-ai-collab`，状态 `Finished`）。
- **官方存档文件路径**：`gemini/Codabench官方结果_submission1062221.html`。
- **官方私有测试得分**：`perturbation/xatlas-orion-perturb = 0.0998`（metric `per_source_mean_projected_residual_effect_pearson`，`status ok`；其余 63 项任务未提交按 0 计）。
- **按预先声明标准的严格定性**：
  1. `0.0998` 为该任务当前最高分，较 22 个自主单次运行 agent 最佳（`0.0825`，`mlevolve__claude-opus-5-5`）高 **`+0.0173`**，但**小于官方 `0.02` 噪声阈值**，落在 `0.0825–0.1025` 区间，定性为**「在噪声范围内领先」**，不得写成「显著超越」或「SOTA 突破」；
  2. **实测验证→测试偏移**：内部 10 开发分割均值 `0.1099` vs 官方私有测试 `0.0998`，实测偏移 **`−0.0101`**；
  3. **不可混比类别声明**：本方案为人机协作、多轮迭代、含约 25 个手调常数的专用管线，**与 2 小时自主单次运行 agent 不是同一比较类别**；
  4. 本次提交为唯一一次正式测试评估，严禁据此改模型或再次提交任何变体（包括 S5）至 Codabench。

---

## 三、阶段 B 预先声明的判定标准（于 `5d96040` 实验运行前锁定）

1. **T2**：全 5 指标（`v5`、`Δ 表达 Pearson`、`Top-50 DE Pearson`、`Top-50 方向准确率`、独立于 v5 的原始 Δ 空间扰动区分度 `disc_raw` 及 v5 空间 `disc_v5`）全部上主表；在 `disc_raw` 上超过全部线性基线 `>+0.02` 且 `10/10` 分割领先方可宣称独立区分度占优；传统三指标落后必须正面报告并解释。
2. **T4**：以冻结 `Full_S7` 为基准逐一留一移除 10 个模块/子组件；在**两个数据集**上均值净贡献均 `> +0.002` 且 `≥ 7/10` 分割为正判定为「跨数据集可迁移贡献」；仅在 xatlas 为正判定为「数据集特异」；在 Replogle `< 0` 如实写为「外部负贡献」。
3. **T5**：使用本地官方 `ReactomePathways.gmt`（CC BY 4.0，SHA-256 `a9feaa3a...`）与评测基因全集背景（`M_bg = 1,784`，同时核验 GMT 交集 `1,452`）做单侧超几何检验 + BH 校正；若 `S7−S0` 与 `S6−S5` Top-150 均无 `FDR < 0.05` 通路则彻底删除通路富集表述。
4. **T6**：仅在 Replogle 2022 上使用公开 GENCODE v46 GRCh38 GTF（SHA-256 `d620d548...`）检验 `< 1 kb` / `< 10 kb` / `< 100 kb` 顺式邻近富集与反事实去除效果；若 `< 10 kb` 显著高于零分布（>5×），必须如实报告 CRISPRi 顺式混杂且明确声明 xatlas 无法排除该混杂。

---

## 四、各任务完整实测数据（含全部不利结果）

### 1. T2：全五指标主表与行级 95% CI（`reports/t2_all_5metric_summary.json`, `reports/t2_discrim_delta_row_ci_summary.json`）

#### (a) Replogle et al. (*Cell* 2022, 10 严格同时联合留出分割，零调参；S5 与 S7 并列为主结果)

| 方法 | v5 残差 $r$ | Δ 表达 $r$ | Top-50 DE $r$ | Top-50 方向 | `Discrim(Delta)`（原始 Δ 空间） | `Discrim(v5)` |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0009 | −0.0322 | −0.0531 | 0.4992 | 0.5000 | 0.4999 |
| 30-kNN | 0.0971 | **0.1479** | 0.1822 | 0.6095 | 0.6365 | 0.6329 |
| AE 双线性基线 | 0.0967 | 0.1470 | 0.2026 | 0.6255 | 0.6461 | 0.6375 |
| 自实现未调参简化 ResMLP | 0.1091 | 0.0876 | 0.1896 | 0.5901 | 0.6447 | 0.6529 |
| 自实现未调参简化 GNN | 0.1233 | 0.0971 | 0.1892 | 0.5999 | 0.6403 | 0.6499 |
| S0 (`100-SVD + Ridge`) | 0.0883 | 0.1408 | 0.2119 | **0.6256** | 0.6505 | 0.6404 |
| S1 (`+预处理+HybridTOM`) | 0.1239 | 0.1054 | 0.1917 | 0.5968 | 0.6728 | 0.6735 |
| S2 (`+谱剥离`) | 0.1215 | 0.1033 | 0.1918 | 0.5922 | 0.6781 | 0.6790 |
| S3 (`+跨系基因图`) | 0.1267 | 0.1079 | 0.2067 | 0.5996 | 0.6865 | 0.6871 |
| S4 (`+非可重复锚点`) | 0.1267 | 0.1077 | 0.2088 | 0.5995 | 0.6861 | 0.6866 |
| **S5 (`+跨系桥接`，并列主结果)** | **0.1401** | 0.1193 | **0.2170** | 0.6064 | 0.6890 | 0.6893 |
| S6 (`+互惠响应借用`) | 0.1292 | 0.1106 | 0.2155 | 0.6020 | **0.6993** | **0.7001** |
| **S7 (`+降噪吸附`，并列主结果)** | 0.1359 | 0.1162 | 0.2013 | 0.5962 | 0.6762 | 0.6771 |

**Replogle 行级配对自助法 95% CI（4,000 次分层重抽，与主控 `replogle_rowlevel.py` 同口径）**：
- `disc_raw`（独立于 v5 的原始 Δ 空间区分度）：
  - **`S5 − AE`**：均值 **`+0.0417`**，中位 95% CI **`[+0.0198, +0.0721]`**（**9/10 分割 CI > 0**，10/10 点估 > 0；通过预先声明的 T2 胜出标准！）；
  - **`S5 − S0`**：均值 **`+0.0373`**，中位 95% CI **`[+0.0125, +0.0649]`**（**9/10 分割 CI > 0**）；
  - **`S7 − AE`**：均值 `+0.0268`，中位 95% CI **`[-0.0011, +0.0546]`**（仅 **4/10 分割 CI > 0**，CI 触零）；
- `disc_v5`（v5 投影空间区分度）：`S5 − AE = +0.0498`（`[+0.0279, +0.0730]`, **10/10 CI > 0**）；`S7 − AE = +0.0357`（`[+0.0117, +0.0594]`, **9/10 CI > 0**）；
- `delta_r`（传统未投影 Δ 表达 Pearson，**不利结果**）：`S5 − AE = −0.0277`（`[-0.0472, -0.0092]`, **8/10 CI < 0**）；`S7 − AE = −0.0308`（`[-0.0522, -0.0095]`, **8/10 CI < 0**）。

#### (b) xatlas-orion-perturb（10 开发分割 `dev10` 与 10 额外随机重切分 `resplit10`）

| 方法 | `dev10` v5 | `dev10` Δ $r$ | `dev10` DE50 | `dev10` 方向 | `dev10` `Discrim(D)` | `dev10` `Discrim(v5)` | `resplit10` v5 / `Discrim(D)` |
|:---|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
| Constant Mean | 0.0090 | 0.0608 | 0.0941 | 0.5420 | 0.5000 | 0.4999 | 0.0073 / 0.5000 |
| 30-kNN | 0.0410 | 0.0832 | 0.1327 | 0.5669 | 0.5698 | 0.5664 | 0.0391 / 0.5639 |
| AE 双线性基线 | 0.0672 | **0.1367** | 0.2039 | **0.6048** | 0.5856 | 0.5952 | 0.0646 / 0.5829 |
| 自实现未调参简化 ResMLP | 0.0486 | 0.0390 | 0.0964 | 0.5425 | 0.5951 | 0.5962 | — |
| 自实现未调参简化 GNN | 0.0543 | 0.0432 | 0.1025 | 0.5464 | 0.5900 | 0.5922 | — |
| S0 (`100-SVD + Ridge`) | 0.0613 | 0.1327 | **0.2050** | 0.6046 | 0.5891 | 0.5993 | 0.0596 / 0.5870 |
| **S5 (`+跨系桥接`)** | 0.0887 | 0.0813 | 0.1576 | 0.5768 | 0.6262 | 0.6262 | 0.0878 / 0.6237 |
| S6 (`+互惠响应借用`) | 0.1027 | 0.0937 | 0.1885 | 0.5875 | **0.6453** | **0.6471** | 0.1006 / **0.6439** |
| **S7 (`+降噪吸附`)** | **0.1099** | 0.0995 | 0.1789 | 0.5857 | 0.6286 | 0.6313 | **0.1073** / 0.6263 |

- **xatlas `dev10` 行级 95% CI**：`disc_raw S5−AE = +0.0376`（`[-0.0058, +0.0627]`, 4/10 CI > 0）；`S7−AE = +0.0453`（`[+0.0027, +0.0708]`, 5/10 CI > 0）；`delta_r S5−AE = −0.0554`（8/10 CI < 0）；`S7−AE = −0.0372`（4/10 CI < 0）。
- **为何传统三指标落后但区分度领先**：因为 `S1–S7` 在训练目标上显式正交投影掉了通用应激主轴（`residual_axis`），放弃了占表达方差大头、同时高频出现在单基因 Top-50 DE 中的共享应激响应分量，导致单条未投影相关性和 Top-50 DE 方向准确率低于保留该主轴的线性模型（甚至 `Constant Mean` 在 xatlas 上仅凭共享方向也能拿到 `DE50 = 0.0941`、`Sign = 0.5420`）；而在无法靠共享应激主轴作弊的跨扰动检索指标 `Discrim(Delta)` 上，残差空间方法（尤其 `S5`/`S6`）全面反超线性基线。同时，`S7` 的流形吸附在两数据集上均降低了 `Discrim(Delta)`（Replogle `0.6993 → 0.6762`，xatlas `0.6453 → 0.6286`），揭示了“提升平均余弦 vs 保持扰动个体区分度”的真实权衡。

---

### 2. T3：行级配对自助法 95% CI 汇总（`reports/row_bootstrap_ci_summary.json`）

- **xatlas BenchA (10 开发分割)**：`S7−S0 = +0.0485`（中位 95% CI `[+0.0318, +0.0677]`, **10/10 CI > 0**）；`S7−AE = +0.0427`（`[+0.0256, +0.0606]`, **10/10 CI > 0**）；`S5−AE = +0.0215`（`[+0.0091, +0.0317]`, **10/10 CI > 0**）；`S6−S5 = +0.0140`（`[+0.0050, +0.0244]`, 9/10 CI > 0）；`S1−S0 = +0.0063`（`[-0.0015, +0.0174]`, 4/10 CI > 0）。
- **xatlas BenchB (超严格联合留出 + 重拟合轴)**：`S7−S0 = +0.0457`（`[+0.0276, +0.0635]`, **10/10 CI > 0**）；`S7−AE = +0.0398`（`[+0.0215, +0.0569]`, **10/10 CI > 0**）。
- **Replogle 2022 (10 严格联合分割)**：
  - `S1−S0 = +0.0356`（`[+0.0256, +0.0450]`, **10/10 CI > 0**）
  - `S5−AE = +0.0434`（`[+0.0323, +0.0549]`, **10/10 CI > 0**）
  - `S7−AE = +0.0392`（`[+0.0240, +0.0537]`, **10/10 CI > 0**）
  - `S7−S0 = +0.0476`（`[+0.0326, +0.0610]`, **10/10 CI > 0**）
  - `S5−简化GNN = +0.0168`（`[+0.0038, +0.0302]`, **7/10 CI > 0**）
  - **`S7−简化GNN = +0.0126`（`[-0.0041, +0.0286]`, 仅 3/10 CI > 0）**
  - **`S6−S5 = −0.0109`（`[-0.0169, -0.0050]`, 9/10 CI < 0）**

---

### 3. T4：双数据集逐模块留一消融 LOMO（`reports/t4_lomo_summary.json`）

净贡献定义为 `Δ = Full_S7 − LOO_variant`（v5 指标；正值代表去掉该模块后性能下降，即模块为正贡献；负值代表去掉该模块后性能上升，即模块为负贡献）：

| 留一移除模块 | xatlas 净贡献（中位 95% CI, 正/CI>0） | Replogle 净贡献（中位 95% CI, 正/CI>0/CI<0） | 按预先声明标准判定 |
|:---|:---:|:---:|:---|
| **`S1_HybridTOM` ($p=1.5$)** | **`+0.0109`** (`[+0.0027, +0.0193]`, 10/10, CI>0: 7/10) | **`+0.0103`** (`[+0.0038, +0.0155]`, 10/10, CI>0: **10/10**) | **跨数据集核心贡献（最稳健）** |
| **`S7_WienerSnap`** | **`+0.0072`** (`[+0.0008, +0.0142]`, 10/10, CI>0: 6/10) | **`+0.0067`** (`[+0.0006, +0.0114]`, 10/10, CI>0: **8/10**) | **跨数据集核心贡献（仅 v5；降低区分度）** |
| **`S2_Deflation` (`d0..d10`)** | **`+0.0108`** (`[+0.0004, +0.0198]`, 10/10, CI>0: 5/10) | **`+0.0035`** (`[-0.0051, +0.0113]`, 9/10, CI>0: 1/10) | **跨数据集正贡献（Replogle 单分割 CI 跨 0）** |
| **`S5_CrossLineBridge`** | **`+0.0020`** (`[-0.0058, +0.0098]`, 8/10, CI>0: 0/10) | **`+0.0061`** (`[+0.0022, +0.0103]`, 10/10, CI>0: **9/10**) | **跨数据集核心贡献（Replogle 极强；xatlas LOMO 与 S3/S6 重叠）** |
| `S1_SqrtNCells` | `+0.0011` (`[-0.0020, +0.0050]`, 8/10, CI>0: 0/10) | **`+0.0068`** (`[+0.0009, +0.0123]`, 10/10, CI>0: 6/10) | Replogle 显著有效，xatlas 微正 |
| `S1_SelfMask` | `+0.0004` (`[-0.0012, +0.0022]`, 9/10, CI>0: 0/10) | `+0.0006` (`[-0.0003, +0.0018]`, 9/10, CI>0: 1/10) | 极微正（LOMO 效应 < 0.001） |
| `S1_BaseSub` | `−0.0007` (`[-0.0032, +0.0035]`, 5/10, CI<0: 2/10) | `+0.0007` (`[-0.0018, +0.0031]`, 8/10, CI>0: 0/10) | 中性 / 无稳定贡献 |
| `S3_CrossLineGraph` | `−0.0002` (`[-0.0055, +0.0040]`, 4/10, CI>0: 0/10) | `+0.0016` (`[-0.0029, +0.0059]`, 10/10, CI>0: 0/10) | 在 S5 存在时近似冗余 |
| `S4_NonRepAnchors` | **`+0.0043`** (`[+0.0006, +0.0073]`, 10/10, CI>0: 6/10) | **`−0.0006`** (`[-0.0018, +0.0008]`, 1/10, CI<0: 2/10) | **数据集特异（仅 xatlas 有效）** |
| **`S6_Reciprocal`** | **`+0.0130`** (`[+0.0016, +0.0253]`, 10/10, CI>0: **8/10**) | **`−0.0153`** (`[-0.0240, -0.0080]`, 0/10, CI<0: **9/10**) | **数据集特异 / 外部强负贡献（去掉 S6 后 Replogle v5 升至 `0.1512`）** |

---

### 4. T5：本地 `ReactomePathways.gmt` 超几何富集（`reports/t5_local_reactome_enrichment_summary.json`, `reports/t5_enrichment_inputs.json`）

- **GMT 来源**：`ReactomePathways.gmt`（CC BY 4.0，SHA-256 `a9feaa3ae9bbb746422cba77a669813de627fb15b7d7478b67c8d57020ec26ae`，2,874 条通路）；评测基因全集背景 `M_bg = 1,784`（GMT 交集背景 `1,452`）。
- **`S7_top_164` vs `S0_top_164`**：重叠 `113/164`（`68.9%`），两者最显著富集均为同一组基础代谢/核糖体/蛋白酶体通路（`S0` 富集倍数 `7.47×` 甚至高于 `S7` 的 `6.15×`），确认第五轮“最准基因富集不能证明 S7 特有机理”的判断。
- **增益集合（严格以 1,784 评测全集为背景，且在 1,452 GMT 背景及主控 50 次随机 150 基因零对照 0/50 显著下均通过核验）**：
  - `S7 − S0` Top-150：6 条显著通路（BH-FDR < 0.05），首位 **`mRNA Splicing - Major Pathway`**（`k = 38/145`, `3.12×`, `FDR = 1.18e-08`）；
  - `S6 − S5` Top-150：21 条显著通路（GMT 背景 15 条），首位 **`mRNA Splicing`**（`k = 34/150`, `2.70×`, `FDR = 9.64e-06`）；
  - `S5 − S2` Top-150：5 条显著通路，首位 **`Mitochondrial translation`**（`k = 20/82`, `2.90×`, `FDR = 2.76e-03`）。

---

### 5. T6：Replogle CRISPRi 基因组邻近 / 双向启动子检验（`reports/t6_replogle_crispri_proximity_summary.json`）

- **参考坐标**：`gencode.v46.basic.annotation.gtf.gz`（GRCh38，SHA-256 `d620d548...`），成功映射 `10,530/10,539`（`99.9%`）Replogle Ensembl 基因 TSS。**严禁且未在 xatlas 上做任何基因身份或坐标分析。**
- **邻近对占比与零分布对比**：
  - 100,000 随机基因对零分布：`< 1 kb = 0.002%`，`< 10 kb = 0.003%`，`< 100 kb = 0.022%`；
  - 全部激活供体对（$|z|>4, n=17,213$）：`< 1 kb = 0.320%`（权重占 `0.252%`），`< 10 kb = 0.331%`（权重占 `0.255%`），`< 100 kb = 0.343%`（权重占 `0.256%`）；
  - **强负 $z$ 供体对（$z < -4, n=3,939$）**：`< 1 kb = 1.371%`（**685× 随机**），`< 10 kb = 1.422%`（**474× 随机**，权重占 **`3.32%`**）；而强正 $z$ 对（$z > +4, n=13,272$）`< 10 kb` 仅 `0.008%`（权重 `0.000%`）。
- **Replogle 反事实去除 `< 10 kb` / `< 100 kb` 供体**：
  - 全体 2,614 扰动：`delta_all = −0.01045` → `delta_no10kb = −0.01084`（`< 10 kb` 对总体仅贡献 `+0.00039`）；
  - 前 40 `ownRec` 正增益基因：去除 `< 10 kb`（及 `< 100 kb`）后平均增益 `0.2536 → 0.2471`，**保留 `97.41%`**（且前 40 `S6−S5` 胜出基因中 `0.0%` 首要供体位于 `< 100 kb` 内）。
- **对 xatlas 的硬性局限声明**：尽管 Replogle 上顺式占比极低，但实验证实强负 $z$ 对确实对 `< 10 kb` CRISPRi 双向启动子/局部顺式抑制存在 474× 富集；而 xatlas 上 S6 增益前 40 基因中 78% 首要供体为强负 $z$，且基准规则禁止查 xatlas 基因身份与坐标，因此**该局部实验效应在 xatlas 上无法排除**。

---

## 五、这一轮让结论变弱的发现（全部不利结果集中披露）

1. **Codabench 官方私测未跨过 0.02 噪声门槛**：单次正式提交 `Submission ID 1062221` 实测 **`0.0998`**，较 22 个 2 小时自主 agent 最佳 `0.0825` 仅高 `+0.0173`（`< 0.02` 噪声阈值），属于「噪声范围内领先」，且与单次 2 小时自主 agent 非同一比较类别（本方案为人机协作多轮迭代 + 约 25 个手调常数），实测验证→测试下降 `−0.0101`。当前仍维持 **L2**，不具备顶媒宣布“SOTA 突破”的条件。
2. **传统未投影指标上 S5 与 S7 均显著输给线性基线**：在标准 Δ 表达 Pearson 上，Replogle 上 `S7 − AE = −0.0308`（8/10 CI < 0）、`S5 − AE = −0.0277`（8/10 CI < 0），xatlas 上 `S5 − AE = −0.0554`（8/10 CI < 0）、`S7 − AE = −0.0372`（4/10 CI < 0）；在 Top-50 DE Pearson 与 Top-50 方向准确率上，S7 在两数据集上同样低于 `S0` 和 `AE`。
3. **行级 95% CI 推翻了先前“S7 显著优于简化 GNN”的假象**：Replogle 上 `S7 vs 简化 GNN` 的行级中位 95% CI 为 `[-0.0041, +0.0286]`，仅 **3/10 分割 CI > 0**（此前无效配对 $t$ 检验误报为高度显著）；仅 `S5 vs 简化 GNN` 在 7/10 分割上 CI > 0。
4. **`S7` 的流形吸附在两数据集上均伤害了独立扰动区分度**：在 `Discrim(Delta)` 上，Replogle 从 `S6 0.6993 / S5 0.6890` 跌至 `S7 0.6762`（导致 `S7 − AE` 仅 4/10 分割 CI > 0，而 `S5 − AE` 达 9/10 CI > 0）；xatlas 从 `S6 0.6453` 跌至 `S7 0.6286`。
5. **T4 留一消融证明 10 个模块中仅有 4 个跨数据集稳定为正，其余 6 个为数据集特异、冗余或负贡献**：`S6`（互惠响应借用）在 Replogle 上是**强负贡献（`−0.0153`, 9/10 CI < 0）**，去掉 `S6` 后 Replogle v5 反升至 `0.1512`；`S4` 仅在 xatlas 有效（Replogle `−0.0006`）；`S3`、`S1_BaseSub`、`S1_SelfMask` 在完整系统留一消融中接近冗余或中性。
6. **xatlas 原始 Δ 空间区分度单分割 CI 宽度较大**：xatlas 上 `disc_raw S5−AE` 仅 4/10 分割 CI > 0、`S7−AE` 仅 5/10 分割 CI > 0（虽 10/10 点估为正，但样本量小使单分割行级抽样噪声偏大）。
7. **T6 证实强负 $z$ 供体对存在 474× 的 `< 10 kb` CRISPRi 顺式邻近富集**：导致匿名化 xatlas 上 `S6` 的 `+0.0130` 增益（前 40 胜出基因 78% 为强负 $z$ 供体）无法排除局部双向启动子/顺式抑制混杂。
8. **Replogle 缓存多 guide 去重使用了全屏结果统计量**：`prepare_replogle2022_cache.py` 在同一基因对应多条 guide 时按 `(reproducible, ad_counts, num_cells_filtered)` 选优，已在方法部分如实披露。

ROUND6_DONE
