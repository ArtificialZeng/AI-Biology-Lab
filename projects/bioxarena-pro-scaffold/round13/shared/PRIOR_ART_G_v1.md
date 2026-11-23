# R13 文献检索与 Prior-Art 审计报告（PRIOR_ART_G）

- **智能体**：A-T5（文献与官方指标）
- **检索截止日**：`2026-11-23`（覆盖截至当日之 `bioRxiv`、`arXiv`、`EuropePMC`、`PubMed`、`OpenReview` 及 `NeurIPS / ICML / ICLR / RECOMB / ISMB` 公开文献与基准代码仓）
- **检索对象**：针对 `R12_VERDICT.md` 两条核心主张（`Claim 2` 为本次优先重点，`Claim 1` 为 `2026-11-15` 后新增检索及 `Kawsar et al. 2026` 精准边界核查）进行投稿前先占评估与诚实重叠度判定（`完全 / 部分 / 无`；不写任何对外宣传文字）。

---

## 一、核心结论摘要

1. **Claim 2（优先重点：按响应幅度筛选基因 — 含 split-half / cross-fitting 筛选 — 在 `r_orth_cal` 下持续抬高确定性、与扰动靶点无关常数基线 `E1_ConstMed`）**：
   - 在截至 `2026-11-23` 的全部公开文献与会议记录中，**未发现**任何已发表论文或预印本报告：在零共享控制抽样噪声的确定性设定下，仅因使用 `log1p(c_med * FC̄_{p,g}) - log1p(c_med)` 这类**与目标基因完全无关的确定性常数基线**（`E1_ConstMed`）配合非线性 `log1p` Jacobian 权重 `w(μ) = μ/(1+μ)` 与表达量/响应幅度正相关结构，会在**独立划分（split-half）按响应幅度选取评估基因集（`resp_top50_split`）** 时仍将 `r_orth_cal` 显著推高（我方实测：`A549 N=445` 上 `E1(resp_top50_split) = +0.2565` vs `E1(all_nontarget) = -0.0201`，`H-E2b = +0.2766 [+0.2412, +0.3164]`）。
   - 对关键对标论文 **Nicol, Shivakumar & Irizarry (2026)** (`doi:10.64898/2026.05.07.723486`, `PMC13192888`) 的全文与源码核查表明：其数学推导（Assumptions 1–3, Propositions 1–2）完全围绕**共享控制样本随机抽样噪声协方差 `σ² c_0`** 展开；且在其自身 `Split DEG + control`（Half A 选 top-20 DEG + Half B 拆分控制 `ctl-mean`）实验中，作者观察到零信息基线 `ctl-mean` 坍缩至 `[-0.044, +0.047]`（近零），**既未推导亦未发现**确定性常数基线在 `split-half` 选基因下仍维持显著正相关的机制。
   - **总体重叠评级**：现有文献与 Claim 2 最高仅为 **`部分`**（机制核心无 **`完全`** 先占；但撰写时必须精确区分共享控制噪声文献 [Nicol 2026 / Fukutani 2026 / SafeDelta]、克隆假重复 [Huang 2026]、信号稀释/DEG 加权文献 [Miller 2026 / SBB / scPertEval / AmbiModBench / VCC 2026] 及 `log1p` 变换经典文献）。

2. **Claim 1（`2026-11-15` 至 `2026-11-23` 增量检索 + `Kawsar et al. 2026` 边界核查）**：
   - `2026-11-15` 至 `2026-11-23` 新增检索（覆盖 `bioRxiv` `bioinformatics/genomics/systems_biology/genetics/cell_biology` 全部 458 篇条目、`EuropePMC`、`PubMed` 及 `arXiv` 新近预印本）中，**无任何新增研究**实现或测试跨细胞系 / 跨扰动模态（`CRISPRi -> CRISPRko`）的“纯未扰动基线生物变换（`FCT`）对全部错配源基线（`Mis_s`）的配对 Intersection-Union 检验（`r_orth_cal`）”。
   - 对 **Kawsar et al. (2026, `CanPert/CASTER`)** (`doi:10.64898/2026.11.15.766815`, bioRxiv `2026-11-17`, 代码仓 `v5.1.0`) 源码与全量表格逐项核查确认：`CASTER_wrongctrl` 仅用于**部分筛选场景 `t1` / `t1_grouped`** 下的神经网络控制集输入消融（固定目标基线 `ctrl + shift`、单次循环替换 `1` 个邻近系 `(idx+1)%4`），**从未**运行于零样本未见细胞系 `t2`（且在 `t2` 中 `CASTER` 本身以 `0.1709` 显著劣于 `PerturbMean` `0.1967`，Holm `p=1.0`），也未测试无监督“纯基线变换”或跨模态泛化（`CRISPRi -> CRISPRko`）。其与 Claim 1 为 **`部分`** 重叠，不构成完全先占。

---

## 二、Claim 2 专项审计：基因筛选偏差与确定性常数基线抬高

### 2.1 Nicol, Shivakumar & Irizarry (2026) 全文与官方代码范围精准核查

- **文献标识**：`doi:10.64898/2026.05.07.723486`（bioRxiv v1 `2026-05-12`，v2 `2026-08-11`；`PMC13192888`；GitHub 源码 `nicol_systema` / `spurious_correlation_paper` & `split_control`）。
- **作者实际做了什么（原文与代码可证）**：
  1. **共享控制随机抽样噪声理论（Assumptions 1–3, Propositions 1–2）**：假定 `d/n_0 -> c_0`、控制细胞独立同分布且协方差为 `σ²I_d`、真实均值变化 `μ̃` 为确定性向量。Proposition 1 证明当训练预测与评估真值相减共享同一控制组均值 `\hat{μ}_0` 时，即使零信号（`κ = 0`），`PearsonΔ` 极限仍收敛至非零正偏置 `σ² c_0 / sqrt((r_1² + σ²(c_0 + c_1))(r_2² + σ² c_0))`；Proposition 2 证明将控制细胞划分为互不重叠的两半可消除分子上的共享抽样噪声 `σ² c_0`（零信号下一致收敛于 `0`）。
  2. **“Top-k DEG 筛选 + 拆分控制”实证实验（`split_control/ctl-mean-topk.ipynb`, `eval_split_control_topk.ipynb`, `spurious_correlation_paper/figure_2_recompute.py`）**：作者在 10 个 Perturb-seq 数据集（含 `TianCRISPRa/i`, `ReplogleRPE1_v2`, `ReplogleK562_gwps`, `XuKinetics2024`, `Frangieh*`, `Norman2019`, `Adamson2016`）上比较 `Shared DEG + control` 与 `Split DEG + control`（在 Half A 用 Scanpy `get_DE_genes` 选 top-20 DEG，在 Half B 将控制组 25/25 拆分计算 `ctl-mean`）。结果显示：`Shared DEG + control` 下 `ctl-mean` 为正（TianCRISPRi `+0.408225`，TianCRISPRa `+0.322631`，XuKinetics `+0.154863`，FrangiehControl `+0.117597`，ReplogleK562_gwps `+0.065688`，ReplogleRPE1_v2 `+0.059049`），而切换至 **`Split DEG + control` 后，`ctl-mean` 在全部数据集中坍缩至零附近（`-0.0441` 至 `+0.0471`）**（TianCRISPRi `-0.030718`，TianCRISPRa `+0.005242`，XuKinetics `-0.006349`，FrangiehControl `+0.047121`，ReplogleK562_gwps `-0.004059`，ReplogleRPE1_v2 `-0.002293`）。
- **作者明确没有做什么**：
  1. 全部数学分析仅针对**有限样本共享控制抽样噪声（`σ² c_0`）**；**未涉及**零噪声极限下依然存在的确定性 `log1p` Jacobian 压缩权重 `w(μ) = μ/(1+μ)`。
  2. **未构造也未检验**任何基于表达基准均值和全局倍比变化的确定性、与扰动靶点无关常数预测器（如 `E1_ConstMed`: `log1p(c_med * FC̄_{p,g}) - log1p(c_med)`）。
  3. **未使用**去除了全细胞系平均偏差的正交校准相关系数 `r_orth_cal`。
  4. 由于作者测试的 `ctl-mean`（控制半组对另一半组差值）在 `Split DEG + control` 下直接归零，**Nicol et al. (2026) 事实上将“split-half 选基因 + split-half 拆控制”视作彻底消除了伪相关膨胀**，完全未察觉确定性常数基线在 `split-half` 响应幅度选基因下（`resp_top50_split`）仍因基线丰度–幅度耦合而维持高正偏置。
- **重叠度判定**：**`部分`**。

---

### 2.2 Claim 2 全领域对照文献清单（含截至 2026-11-23 最新会议与预印本）

| 文献 / 工作 | 精确标识符 | 作者做了什么（已核实） | 未做什么（与 Claim 2 差异） | 重叠度 |
|---|---|---|---|---|
| **Nicol, Shivakumar & Irizarry (2026)** | `doi:10.64898/2026.05.07.723486` (`PMC13192888`) | 推导共享控制抽样噪声导致的 `PearsonΔ` 偏置，验证 `split-control` 与 `Split DEG + control` 使 `ctl-mean` 归零（`[-0.044, +0.047]`）。 | 未涉及确定性常数基线 `E1_ConstMed`、`log1p` Jacobian、正交校准 `r_orth_cal`；且未发现 `split-half` 选基因下确定性基线仍显著抬高。 | **`部分`** |
| **Fukutani et al. (2026)** (`npj Syst Biol Appl`) | `doi:10.1038/s41540-026-00722-2` (`PMC13287743`) | 解析推导标准化中共享不稳定参照 `Z, W`（`X-Z`, `(X-Z)/W`）对独立随机变量诱导虚高相关性。 | 纯共享随机参照噪声效应；不涉及 `top-k` 响应幅度选基因、确定性常数基线或 `log1p` Jacobian。 | **`部分`** |
| **Liang et al. (2026, `SafeDelta`)** | `doi:10.64898/2026.10.05.756906` (`bioRxiv`, `2026-10-13`) | 针对 few-shot 细胞系微调，证明共享同一控制组在反比例测序深度下携带控制噪声并虚高内积/PCC；提出双半控制门控 `SafeDelta` 决策准则。 | 同样研究控制抽样噪声携带（代码中将常数向量 PCC 定义为 `NaN` 弃权），无 `top-k` / `resp_top50_split` 选基因偏差与 `E1_ConstMed` 现象。 | **`部分`** |
| **Huang (2026)** | `arXiv:2610.32962v2` (`2026-10-23` / v2 `10-26`) | 证明在混合 CRISPR 筛选中随机拆半（random halves）共享同一转导克隆而导致克隆假重复（`G σ_a² / K_eff`），使空扰动亦具“可重复性”并虚高 `cell-eval` 上限与迁移模型；用多组静态基因集及交叉验证检验迁移。 | 机理为**细胞克隆假重复**而非确定性非线性表达尺度假象；基因集为静态全局过滤（`top 2,000` 高控制均值等），未研究按单扰动响应幅度 `top-k` 筛选及 `E1_ConstMed` 在 `r_orth_cal` 下的抬高。 | **`部分`** |
| **Miller, Mejia et al. (2026)** (`Nat Biotechnol` / 先前预印本) | `doi:10.1038/s41587-026-03307-w` (`arXiv:2506.22641v1`) | 指出未校准 `MSE` 与控制偏移指标因“信号稀释（signal dilution）”与“控制偏置（control bias）”偏袒均值基线；提出 DEG 加权 `WMSE` 与 `$R^2_w(\Delta)$`，使均值基线得分跌至零/负区间。 | 其结论是 DEG 加权**压低**均值基线至 null，与 Claim 2 揭示的“按响应幅度筛选 `top-k` 基因配合 `r_orth_cal` 反而**系统性抬高**确定性常数基线 `E1_ConstMed`”作用方向与机制均相反。 | **`部分`** |
| **Vollenweider & Bühlmann (2026, `SBB`)** | `doi:10.64898/2026.04.20.719650` (`bioRxiv`, v2 `2026-05-27`) | 提出 Signal–Bounds–Baselines 原则、`GlobalShiftModel`（`X_0 + b`）、DEG 加权/PDS 过滤指标与 BDS/DRF 校准。 | 未考察 `split-half` 选基因（`resp_top50_split`）对确定性 `log1p(c_med * FC̄)` 基线在 `r_orth_cal` 下的抬高效应。 | **`部分`** |
| **Schäfer et al. (2026, `scPertEval`)** | `doi:10.64898/2026.07.23.740433` (`bioRxiv`, `2026-07-27`) | 构建协议分类法与 Python 包，系统比较 `pearson_pert_top_k`、`pearson_ctrl_top_k`、`heg_k` 等空间对 DRF/BDS 排序的影响。 | 其 `top_k` 直接取于同侧 GT 效应大小（无 `split-half` 去泄漏检验），亦未发现或报告 `E1_ConstMed` 在 `resp_top50_split` 下于 `r_orth_cal` 维持虚高正值。 | **`部分`** |
| **Viñas Torné et al. (2026, `Systema`)** | `doi:10.1038/s41587-025-02777-8` (`Nat Biotechnol`) | 识别扰动与控制间的系统性共享变异（systematic variation）并提出以 `all_perturbed_mean` 为参考的 `pearson_pert` 与非靶点参考。 | 关注全转录组去系统偏移；未研究按响应幅度选基因后确定性 `log1p` 变换基线 `E1_ConstMed` 在正交化相关 `r_orth_cal` 下的再度膨胀。 | **`部分`** |
| **Huang, Yang, Yu, Chen & Li (2026, `AmbiModBench`)** | `arXiv:2609.32527v1` (`2026-09-26`) | 批评绝对指标混淆共享背景响应，采用 Top-100 / Top-2% / Top-20（`|y_t|`）方向一致率并与同 split 训练集均值作差比较；用基因置换实验（Fig 3c）揭露部分幅度相关指标对基因身份不敏感。 | 将训练集均值作为**减项参考**（未审计其绝对或正交校准值是否被 `top-k` 选基因抬高），且未涉及 `split-half` 下 `E1_ConstMed` + `log1p` Jacobian 耦合。 | **`部分`** |
| **Panagopoulos, Omidi, Tong & Pang (2026)** | `arXiv:2610.26267v1` (`2026-10-19`; OpenReview withdrawn `1YbFu1Ad0l`) | 揭示 Varying top-N DEG 时 Precision 排名 `86.7%` 翻转，因将单一预测均值重复填充至细胞行会虚高 Mann-Whitney 检验统计量（`67.1%` 基因对过阈）。 | 揭示的是**样本量假重复膨胀 DE 显著集**，与 `resp_top50_split` 下单向量 `E1_ConstMed` 的 `log1p` Jacobian–丰度相关性无关。 | **`部分`** |
| **Gao (2026, `The Delta Gap`)** | `arXiv:2611.10658v1` (`2026-11-07`) | 量化全局 Pearson 与 delta Pearson 鸿沟（GEARS `0.977` vs `0.076`），并加测 top-20（`|δ_obs|`）Pearson；指出训练平均 delta 基线在 Replogle K562 上达 `0.267`。 | 未剥离共享响应，亦未检测 `split-half` 按幅度选基因对确定性 `E1_ConstMed` 在 `r_orth_cal` 下的抬高效应。 | **`部分`** |
| **Arc Virtual Cell Initiative & Goodarzi (2026, `cell-eval2` / VCC 2026)** | `doi:10.1016/j.cell.2026.08.004` (`Cell`, `2026-08-26`; `vcc2026-metrics-brief`) | 规定 DE 指标先由参照控制组过滤均值 `> 5 CPM`（参考单侧门控，独立于预测）后行 BH 检验与方向/LFC 评分。 | 仅基于对照均值设下限门控，未研究按扰动响应幅度 `top-k` 选基因诱导确定性常数基线虚高。 | **`部分`** |
| **Ahlmann-Eltze, Huber & Anders (2025) / Wong et al. (2025) / Csendes et al. (2025)** | `doi:10.1038/s41592-025-02772-6` (`Nat Methods`); `doi:10.1093/bioinformatics/btaf317` (`Bioinformatics`); `doi:10.1186/s12864-025-11600-2` (`BMC Genomics`) | 在千余高表达基因（`heg_k`）或基准集上展示简单均值/线性基线匹敌或超越深度学习/基础模型。 | 所用均值为训练组扰动均值（或 `top 1000 HEG` 固定集合），不涉及按单扰动幅度 `top-k` / `split-half` 选基因诱导 `E1_ConstMed` 膨胀。 | **`部分`** |
| **Sun, Zhu & Chen (2026) / Kalygina et al. (2026) / ICLR 2027 Anon.** | `OpenReview:EZ2KZwM1mR` (MLGenX@ICLR'26); `OpenReview:U5Gq071xcI` (GenBio'26); `OpenReview:DVGHFaScK9` (ICLR'27 sub.) | 分别揭示数据筛选对众数基线虚高（`0.62-0.80` 坍塌）、条件级基线饱和度分层（Baseline Saturation）、化学扰动伪散装/学习免费基线三轴诊断。 | 均不针对 `split-half` 按响应幅度选基因诱导 `E1_ConstMed` 在 `r_orth_cal` 下的确定性偏置。 | **`部分`/`无`** |
| **David & Peng (2026-11-18, `ReliFlow`)** | `arXiv:2611.23999v1` (`2026-11-18`) | 仅在**训练集单扰动内做 cross-fitting** 以估计残差方差不确定性，加权组合扰动流匹配传输场。 | 与评估阶段 `split-half` 基因筛选偏差及常数基线虚高完全无关。 | **`无`** |
| **Raza et al. (2026-11-20) / Sanchez et al. (2026)** | `doi:10.64898/2026.11.19.767931`; `doi:10.64898/2026.10.02.756269` | 分库水平假重复检验（单库 split-half 测试 `99-100%` 假阳性）；单细胞聚类后患者队列 DE 双重浸渍与批校正消除。 | 聚焦细胞分群/文库批次样本间检验，未涉及 Perturb-seq `E1_ConstMed` 或响应幅度选基因问题。 | **`无`** |
| **经典统计与表达尺度奠基工作**：Booeshaghi & Pachter (2021), Ahlmann-Eltze & Huber (2023), Townes et al. (2019), Kriegeskorte et al. (2009), Neufeld et al. (2022), Song et al. (2026, `ClusterDE`), Lawley (1943/1944) | `doi:10.1093/bioinformatics/btab085`; `doi:10.1038/s41592-023-01814-1`; `doi:10.1186/s13059-019-1861-6`; `doi:10.1038/nn.2303`; `doi:10.1093/biostatistics/kxac047`; `doi:10.1177/15578666251383562`; `doi:10.2307/3607974` | 分别奠定 `log1p` 依赖表达均值 `μ` 的方差压缩性质、双重浸渍（circular analysis）选择偏置理论、泊松计数拆分与选择后推断消除假阳性之方法论基础。 | **无一**将 `log1p` Jacobian 与“非随机（按响应幅度）选基因 + 目标无关常数预测器 + `split-half` 仍不可消除的正交校准偏置”关联起来。 | **`部分`**（理论基石） |

**Claim 2 综合判定**：现有关于统计“split-half / cross-fitting 消除选择泄漏”的直觉（含 Nicol et al. 2026 与 Neufeld et al. 2022）仅成立在“特征/靶标独立于确定性预测函数坐标权重”的随机噪声模型中；**Claim 2 发现的是由 `log1p` 非线性压制与真实生物响应幅度–基因基准丰度共变所致的结构性非零偏置（在 `split-half` 划分下仍顽固留存 `+0.2565`），截至 `2026-11-23` 无完全或本质等效的先占文献。**

---

## 三、Claim 1 审计：`Kawsar et al. 2026` 精准边界与 `2026-11-15` 后增量检索

### 3.1 Kawsar et al. (2026, `CanPert / CASTER`) 全源码与结果核查

- **文献标识**：Kawsar, Dip, Paul, Noman, Zhang & Moni (2026), *"Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line"*, `doi:10.64898/2026.11.15.766815`（bioRxiv v1 `2026-11-17`；代码仓 `CanPert` `CITATION.cff v5.1.0`, `date-released: 2026-11-15`）。
- **作者实际做了什么（核实 `experiments.py`, `train.py`, `results_v1.py`, `t1_fold_means.csv`, `t2_fold_means.csv`, `H2_tests.csv`, `calibrated_metrics.csv`）**：
  1. 覆盖 4 个人类细胞系 CRISPRi 筛选（`K562, RPE1, HepG2, Jurkat`，6,546 共享基因上取 top 2,000 HVGs）。
  2. 仅在**部分筛选适配任务 `t1`（含 `t1_grouped`，目标细胞系有 25% 扰动用于训练、5% 验证、70% 测试）** 加入了 `CASTER_wrongctrl` 对照组（学习曲线 `lc` 中明确关闭 `wrongctrl=False`；**零样本未见细胞系任务 `t2` 中完全不含 `CASTER_wrongctrl`**，已在 `t2_fold_means.csv` 与 `calibrated_metrics.csv` 核实）。
  3. `wrong_ctrl_predictor`（`train.py`）的实现为：用另一细胞系 `LINES[(idx+1)%4]` 的控制细胞输入已在部分筛选集上训练完毕的神经网络 `CASTER`，计算预测增量 `shift = base(p, wc).mean(0) - wc.mean(0)` 后**加回目标细胞系自身的正确控制细胞均值上**（`ctrl + shift`），原文注释写明 *"so only the model's use of the control set changes, not the baseline expression level"*。
  4. `t1` 下配对检验中 `CASTER` (`0.2759`) 优于单周期错配对照 `CASTER_wrongctrl` (`0.1516`)，同时未用任何单细胞控制信息的纯线性基线 `Linear` 在 `t1` 上反而略高（`0.2821`；在 `calibrated_metrics.csv` 为 `Linear 0.2350` vs `CASTER 0.2105`），且作者自己在 `results_v1.py` 指出 *"Linear sees no control cells at all and does as well. This suggests that the adaptation comes mainly from the screened knockdowns of the line"*，并在 `pearson_deltapert_top100` 下发现仅输入细胞系标签 `CASTER_label` 与真实控制细胞表现接近。
  5. **零样本未见细胞系任务 `t2`（`H2_tests.csv`）**：`CASTER` (`mean = 0.1709`) **显著败于简单训练集均值 `PerturbMean` (`0.1967`)**（Holm 校正单侧 `p = 1.0`, median diff `-0.0142`；`de_auprc`: `0.1074` vs `0.1359`）。
- **与 Claim 1 的本质边界差异**：
  1. **预测对象根本不同**：Kawsar et al. 测试的是**经过有监督微调的扰动特定神经网络**对“控制集特征条件化”的消融（保留目标系基线表达并依赖 25% 靶向扰动训练），而非**零靶向扰动训练、完全基于未扰动控制细胞群基线构成的纯生物变换基线（`FCT`）**；且其真实零样本预测 `CASTER` 在 `t2` 彻底失败于 `PerturbMean`。
  2. **实验设计与检验逻辑不同**：`CASTER_wrongctrl` 仅配对单个循环邻系，未执行针对全部 6 个源系错配均值基线（`Mis_s`）的 Intersection-Union Test（`min_s (FCT - Mis_s) > 0`），无 `r_orth_cal` 去偏移正交校准，亦未验证跨干扰模态迁移（如我方从 CRISPRi 跨至 `Norman A549` CRISPRko 的 `6/6` IUT 成立，最小边际 `+0.0796`）。
- **重叠度判定**：**`部分`**。

### 3.2 `2026-11-15` 至 `2026-11-23` 新增公开文献检索结果

| 文献 / 工作 | 精确标识符与发布日 | 核心内容（已核实） | 涉及 Claim 1 与否 / 重叠度 |
|---|---|---|---|
| **Sharma et al. (2026)** | `doi:10.64898/2026.10.20.760537` v2 (`2026-11-18`, bioRxiv) | 用 guide abundance 代理单细胞扰动强度，指出 `STATE` / `Cell2Sentence` 在 `HCT116/HEK293T` 下无法在检索中复现扰动身份差异。 | 无基线变换 / 错配控制检验；**`无`**。 |
| **Wei et al. (2026, `UPM`)** | `doi:10.64898/2026.10.25.761625` v2 (`2026-11-18`, bioRxiv) | 构建 110M 细胞、1B+ 参数统一 VAE 预测器，在 Tahoe-100M 与细胞因子基准上对比 `State`/`Stack`。 | 纯大规模端到端深度学习迁移，无纯基线变换 `FCT` vs `Mis_s`；**`无`**。 |
| **Cupaioli & Merelli (2026, `fmaudit`)** | `doi:10.64898/2026.11.03.763968` v3 (`2026-11-19`, bioRxiv) | 测试 5 个 scFM 与 `STATE` 在 K562、类器官与鼠脑皮层，发现基因敲除不优于表达匹配对照基因删除及训练均值。 | 用表达匹配对照基因作负对照，未研究跨细胞系基线变换 `FCT` vs 全部错配源基线；**`无`**（对 Claim 1）。 |
| **He et al. (2026)** | `doi:10.64898/2026.11.16.767041` v1 (`2026-11-17`, bioRxiv) | 用二维序列–功能模型条件化单一 IFN-β 单细胞状态预测表达响应（Spearman `0.616`），并检出比例统计量按较大项门控造假，及条件特异性 eQTL 失败。 | DNA 序列模型状态条件化与 eQTL；非 Perturb-seq 基线变换或 `Mis_s` 对照；**`无`**。 |
| **Bendidi et al. (2026, `Occam`)** | `arXiv:2611.21959v1` (`2026-11-16`) | 提出去噪 MLP 框架用于跨细胞系/双扰动预测，并采用独立 A/B 控制组避免 Nichol 噪声通胀；基线含 `Transfer` 与 `Mean perturbation`。 | 评估有监督模型及常规 delta 直接迁移，无基于未扰动控制的 `FCT` vs `Mis_s` IUT；**`无`**（纯基线变换机制方面；基准方法对照为 **`部分`**）。 |
| **David & Peng (2026, `ReliFlow`)** | `arXiv:2611.23999v1` (`2026-11-18`) | 训练集内 cross-fitting 估计单扰动可靠性用于组合扰动流匹配。 | 组合扰动方法，无跨细胞系/模态基线变换对比；**`无`**。 |
| **Tarikere (2026)** | `doi:10.21203/rs.3.rs-11768707/v1` (`2026-11-18`, Research Square) | Adamson 单数据集（82 扰动）上分析 OOD 预测器自适应选择（GO/GenePT MLP）与 Embedding Coverage 的相关性。 | 探索性单系 OOD 模型选择；**`无`**。 |
| **Liu et al. (`PerturbBridge`) / Hutchins et al. (`CIRCE`) / Raza et al.** | `doi:10.64898/2026.09.24.754259` v2 (`11-16`); `doi:10.64898/2026.11.19.767979` (`11-20`); `doi:10.64898/2026.11.19.767931` (`11-20`) | 分别聚焦潜在 Schrödinger 桥运输预测、组合信号强度的玻尔兹曼监督分级、单细胞耐药分库级噪音阈值框架。 | 均不涉及“正确未扰动控制 vs 错配源控制”基线变换；**`无`**。 |
| **更早 Claim 1 核心背景工作补充确认**：Huang et al. (`CGC`, `2026-10-01`)；Zhang et al. (`2026-11-04`, 低秩解码瓶颈) | `doi:10.64898/2026.09.30.755756`; `arXiv:2611.05849v1` | 分析跨上下文特异性压缩现象（CGC）；证明 STATE 解码器仅具 4–32 维低秩瓶颈且可靠性受训练最近邻距离限制。 | 解释复杂模型何以失败，均未建立或验证无参数 `FCT` vs `Mis_s` IUT 特异性；**`部分`/`无`**。 |

---

## 四、给主控决策的文献规范约束（无对外文字）

1. **引用 Nicol et al. (2026)** 必须精准限定为：“处理有限样本**共享控制抽样噪声**并证明 `Split DEG + control` 下 `ctl-mean` 坍缩至零，未触及确定性、靶点无关常数基线（`E1_ConstMed`）在 `split-half` 选基因（`resp_top50_split`）与 `log1p` Jacobian 下在 `r_orth_cal` 中保持虚高正偏置的机制。”
2. **引用 Kawsar et al. (2026, `CanPert`)** 必须精准限定为：“其 `CASTER_wrongctrl` 仅在部分筛选任务 `t1/t1_grouped` 下替换 1 个循环邻近细胞系的控制细胞输入以验证神经 FiLM 模块的控制条件化（且固定目标基线并依赖 25% 靶向扰动训练），未运行于零样本未见细胞系 `t2`（其 `t2` 上 `CASTER` 本身劣于 `PerturbMean`），亦未测试纯基线变换 `FCT` 对全源错配均值 `Mis_s` 的 IUT 或 `CRISPRi -> CRISPRko` 跨模态特异性。”
3. 截至 `2026-11-23`，`R12_VERDICT.md` 所保留的 **Claim 1** 与 **Claim 2** 均**不存在任何 `完全` 重叠的先占文献**；现有最高相关工作均属方法学/统计基础上的 **`部分`** 重叠。
