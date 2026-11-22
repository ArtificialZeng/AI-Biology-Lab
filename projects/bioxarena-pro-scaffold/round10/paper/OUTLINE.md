# 第九轮 + 第十轮合稿论文骨架（`round10/paper/OUTLINE.md`）

- **定位与范围**：第九轮（`PREREG_A v4`，6 条既有源线留一交叉验证：`K_conf = 3` 确认 co-fluctuation 与边际方差缩放）+ 第十轮（`PREREG_C`，4 条全新未见人类 CRISPRi 细胞系外推评测：`K_conf = 2` 前瞻检验 `FCT` 细胞系特异性与非劣效性）+ **R10-E 事后探索性基线消融预留槽位**（严格待 `A-T2` / `A-T3` 正式复核及 `A-T5` 查新完成后填入）。
- **版本与治理锚点**：
  - Round 9：`PREREG_A v4`（commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`），`CODE_FREEZE`（`1ca9f60f49adf1ac7cce8266e40533273411db30`）。
  - Round 10：`PREREG_C`（commit `dc4b4a5`，sha256 `0fb48aeb08aa6befe9583030f33ec668155e9e937616d291ea0c211e0d1bf4c6`），`QC_LOCK`（`86872ed741c3aa681a7a9091e60d03bd9922038e`），`CODE_FREEZE_R10`（`f8d2ebfc238d83a0c3a6c86696950ed292bf90e9`，22/22 文件哈希冻结），`CONFIRMATORY_LINES.txt`（`K_conf = 2`: `dld1`, `cd4t_rest`）。
- **AI 治理与作者披露声明（必保留）**：人类负责人 Zijian Zeng 无正规生物学背景；主控 `Claude`（Anthropic）负责预注册审核、结果核验及独立复算；执行智能体 `Antigravity`（Google Gemini 模型）负责数据标准化（`A-T1`）、主分析与基线消融（`A-T2`）、独立复核与 `C3` 交叉核对（`A-T3`）、VCC 2026 合规/前瞻管线与合稿骨架（`A-T4`）及 API 文献查新（`A-T5`）。所有阴性/未确认结果、预注册偏差与文献边界均照实披露。

---

## 1. 标题候选（Title Candidates）

> **措辞原则**（沿用 Round 9 `PAPER_REVIEW_R1.md` **B3**）：禁止使用 *"Rejects"*、*"Disproves"* 等过强证伪词；预注册阈值未通过统一表述为 **"Did Not Confirm"** / **"Fails Pre-Registered Specificity"**；且不把尚未由 `A-T2`/`A-T3` 复核完成的 `R10-E` 结论提前写死进主标题。

1. **保守写实型（主推，聚焦两轮预注册主判定）：**  
   *Pre-Registered Audits Across Ten Human Perturb-seq Screens Do Not Confirm Target-Line Specificity in Zero-Shot Control-Conditioned CRISPRi Predictors*
2. **方法学分解 + 实证并重型：**  
   *Separating Amplitude, Marginal Variance, and Context Specificity in Zero-Shot CRISPRi Perturbation Prediction: Two Pre-Registered Audits Across Ten Human Cell Lines*
3. **贴合 VCC 2026 瓶颈问题型：**  
   *Why Line-Specific Zero-Shot CRISPRi Generalization Remains Unconfirmed: Amplitude Orthogonalization and Matched-Control Audits Across Ten Perturb-seq Datasets*
4. **预留 `R10-E` 复核后可选后缀扩展（仅在 `E1–E7` 经 `A-T2`/`A-T3` 交叉复核成立后启用）：**  
   *... Despite Consistent Gains Over Unscaled Perturbation Means* `[OR if R10-E confirmed: ... Attributing Fold-Change Transfer Gains to Nonlinear Transformation Rather Than Target Baselines]`

---

## 2. 摘要骨架（Abstract Skeleton）

> **排版与篇幅要求**：成稿时主确认性摘要控制在 ~250–300 词内（若期刊允许结构化摘要可分成明确两段）；**必须严格把「预注册确认性结果」与「事后探索性 `R10-E` 结果」分成两段隔离**，且清晰标注 `K_conf = 3`（R9）与 `K_conf = 2 confirmatory cell lines`（R10）。

### 2.1 第一段：预注册确认性研究与主要判定（Confirmatory Study & Results）
1. **研究背景与核心困境**：单靠目标细胞系未扰动对照（zero-shot）预测基因型/细胞系特异的 CRISPRi 转录响应仍是核心难题（Arc VCC 2026 Wrap-Up 显示顶级模型跨细胞系特异性斜率仅 `0.02–0.13`，对照生物重复 `1.02`）。
2. **方法学混淆因素（`R18`）**：直接用单行标量缩放跨线平均（`Transfer_Mean_GlobalScale`）即可产生显著原始残差相关（R9 确认集 `r_resid = -0.126` 至 `+0.250`，全 6 线 `-0.316` 至 `+0.352`；R10 确认集 `dld1 +0.240`、`cd4t_rest -0.262`），而在单扰动内对共享跨线响应 $\hat{\Delta}^{\mathrm{TMean}}_p$ 做正交投影消幅（`r_orth_cal`）后精确归零（`0.000` / `NaN`）。
3. **第九轮（`PREREG_A v4`，6 线 `1,419,744` 细胞，`6,148` 基因，LOLO `K_conf = 3`）**：在发展线（`Replogle K562/RPE1`）锁定管线并排除低靶点 `KOLF2.1J`（`N = 18`）后，保留 `nadig_hepg2`（`N = 799`）、`nadig_jurkat`（`N = 934`）与 `vcc2025_h1`（`N = 80`）作盲态留出确认。在 Bonferroni `97.5%` 置信区间与 `3/3` 一致性规则下，**Group A（`Hybrid_ProjShared`，协方差共波动转移）与 Group B（`Transfer_Mean_VarScaled` / `TMVS`，边际方差缩放）均未确认（`0/4` 标准通过）**：二者均未达 `slope > 0.13`（`0/3`）、均未在 `H1` 上使 `r_orth_cal > 0`（`2/3`）、均未击败全部 5 条错配细胞系对照或 100 次基因置换零模型（`0/3`），且对 `Transfer_Mean` 全部呈现劣效（`3/3`）。在本数据集中，共波动结构在边际方差之外未显现可检测增量（`R14`，`2/3` 未一致性通过）。
4. **第十轮（`PREREG_C`，4 条全新人类单 sgRNA CRISPRi 数据集 `4,264,773` 质控细胞，**`K_conf = 2 confirmatory cell lines`**）**：针对 R9 事后观察到的二次基线 `Fold-Change Transfer`（`FCT`，匹配 VCC 2026 #1 `CellularAI` / #82 `AtlasShift` 缩放类别），预先锁定从零下载至评测的 6 步盲态管线；`H9` 神经元因高 MOI（`~4`）导致单 sgRNA 靶点稀疏且超出 `60 GB` 预算在下载前排除（`P2`），`calu3`（`N = 8`）与 `wtc11_i3neuron`（`N = 2`）归入次要报道，`dld1`（`N = 86`）与双供体 `cd4t_rest`（`N = 302`，`284` 个跨供体两半有效）通过确认准入门槛。在独立前瞻验证下：
   - **`H-C1` 确认（`2/2`）**：`FCT` 去幅度信号显著大于零（`dld1` `r_orth_cal = +0.1401 [+0.0870, +0.1917]`；`cd4t_rest` `+0.0560 [+0.0171, +0.0938]`），且在控制 `ln(n_cells)` 与敲低深度 `|k_p|` 后截距仍保持正值（`S7`）。
   - **`H-C2` 未确认（`0/2`）**：`FCT` `97.5%` CI 下界在两条线上均未同时高于 6 条错配源线对照均值的上界（`+0.1092`、`+0.0761`）及 100 次基因置换的 `97.5` 百分位（`+0.1179`、`+0.0554`）。**主要科学论断「`FCT` 捕捉到目标线特异的基因层面扰动信号」未确认（NOT CONFIRMED）。**
   - **`H-C3` 并行实用性非劣效通过（`12/12` 非劣）**：相对 `Transfer_Mean`，`FCT` 在全部 6 种评价口径上均无一例劣效（其中 `cd4t_rest` `cent_de50_r` 严格更优 `+0.0305 [+0.0177, +0.0435]`；仅 `S1` 纳入非可重复靶点后在 `cent_*` 出现 `4` 项劣效）。
   - 三方独立代码/复算（`A-T2`、`A-T3` 与主控 `Claude`）数值精确一致（R10 主分析 `349` 项跨组差异 `max <= 6.66e-16`）。

### 2.2 第二段：事后探索性基线消融（Explicitly Labeled Post-Hoc Exploratory Slot）
> **[Post-Hoc Exploratory Paragraph — PLACEHOLDER: DO NOT POPULATE CONCLUSIONS UNTIL `A-T2` (`R10E_ABLATION_DONE`), `A-T3` (`R10E_AUDIT_DONE`), AND `A-T5` (`R10E_LIT_DONE`) COMPLETE]**
- **结构占位文案骨架（数值与结论待填充）：**  
  *In a strictly post-hoc exploratory ablation executed after unblinding Round 10 confirmatory outcomes and prior to prospective replication, we decomposed why `FCT` preserves `H-C1` positive de-amplituded signal and `H-C3` non-inferiority despite failing target-line baseline specificity (`H-C2`). Across `[R9 LOLO lines + R10 4 lines]` evaluated under both `resp_top50` and selection-unbiased `all_nontarget` masks, removing gene-specific target control baseline information via `[E1 constant median baseline, E2 p10/p90 quantiles, E3 pooled source control mean, E5 rank-shuffled baseline, E6 log(FC_bar), E7 FC_bar - 1]` yielded `[INSERT VERIFIED r_orth_cal RANGE / PAIRWISE (FCT - E1) 97.5% CIs AND PERM100 COMPARISON HERE]`. `[INSERT VERIFIED CONCISE CONCLUSION ONCE R10E_AUDIT_DONE & R10E_LIT_DONE ARE COMPLETE — explicitly stating that this post-hoc mechanism remains a hypothesis requiring pre-registered prospective testing (planned Round 11) before any confirmatory claim.]`*

---

## 3. 论文总体结构与章节骨架（Manuscript Structure）

1. **Title & AI Governance / Role Disclosure**
2. **Abstract**（严格分段：Confirmatory Summary vs. Demarcated Post-Hoc Exploratory Summary）
3. **Section 1: Introduction**
   - 1.1 零样本跨细胞系 CRISPRi 预测的不对称性：守恒响应高相关 vs. 细胞系特异斜率瓶颈（VCC 2026 Wrap-Up `0.02–0.13` vs. `1.02`；`CASTER`、`Molina & Zhang 2026`、`Shoeibi & Yousefi 2026 arXiv:2608.00152v3`）。
   - 1.2 从 `CIPHER`（共波动理论）与边际方差缩放（`TMVS`）到 `CellularAI`/`AtlasShift` 类基线（`FCT`）：不依赖目标扰动样本的控制条件预测器家族梳理。
   - 1.3 零样本评测中的三大方法学陷阱：`R18` 全上下文标量振幅混杂、`R4/R14` 边际方差 vs. 非对角共波动、`R3/R12` 缺位错配对照与置换零模型及 `S6/arXiv:2608.00152` 选择/细胞数混杂。
   - 1.4 两轮递进预注册研究设计（Round 9 `PREREG_A v4` $\to$ Round 10 `PREREG_C` $\to$ `[R10-E post-hoc ablation rationale]`）。
4. **Section 2: Methods & Pre-Registration Governance**
   - 2.1 十套人类 Perturb-seq 数据集整合、单 sgRNA 严格筛选（`P4`）、`H9` 下载前预算/密度排除（`P2`）与跨轮 HGNC 符号桥接
   - 2.2 标准化表达空间、Addendum-B 可重复性筛选、`P5` 靶点准入及离散物理单元（通道 / 双供体 `P3`）Split-Half 设计
   - 2.3 两轮预注册锁链、运行时盲态闸门（`check_prereg_gate`）、机械准入判定（`P1`）与三方独立复算（`C3`）
   - 2.4 预测器与阴性对照数学定义（Group A/B/R14/R18 + Group C `FCT` + 错配均值/置换 + `[R10-E PLACEHOLDER: E1–E7 数学定义框架]`）
   - 2.5 去幅度评测指标（`r_orth_cal`）、联合重抽样 Bonferroni 置信区间、`H-C3` 并行非劣效、`R19` 三层分解与 `S7`（`P6`）混杂控制
   - **2.6 预注册偏差集中披露（Deviations from Pre-Registration）**：`R9 D1`（`span3` vs `rank1_q1`）、`R10 D-R10-1`（`H-C3` 命名映射）、`R10 D-R10-2`（`H-C1` 额外保护条件）、`FCT` $\mu_t$ 下限微小差异（`0` vs `1e-8`）及 `H9` P2 排除原因文字修正
5. **Section 3: Round 9 Results — Confirmatory Audit of Co-Fluctuation & Variance Transfer (`PREREG_A v4`, `K_conf = 3`)**
   - 3.1 `R18` 单标量幅度混杂证实（`Transfer_Mean_GlobalScale`）
   - 3.2 Group A（`Hybrid_ProjShared`）未确认（`H-A1..A4`）与 `CIPHER_Linear` 跨线零样本退化
   - 3.3 Group B（`TMVS`）未确认（`H-B1..B4`）与 `R14` 共波动增量检验阴性
   - 3.4 `R19` 三层分解、`S6` `resp_top50` 选择偏倚、全套敏感性分析（`S1, S2, S2b, S3, S4a, S4b, S5`）及二次探索线索（`FCT`）
6. **Section 4: Round 10 Results — Prospective Verification of Fold-Change Transfer (`PREREG_C`, `K_conf = 2`)**
   - 4.1 下载前排除 `H9`、机械准入计数（`dld1` `86` / `cd4t_rest` `302`）与跨轮 HGNC 符号挽回（`83->86`、`301->302`）
   - 4.2 主要终点裁决：`H-C1` 确认（`2/2`）vs. `H-C2` 未确认（`0/2`），`FCT` 细胞系特异性未通过
   - 4.3 `H-C3` 并行非劣效确认（`12/12` 非劣，`1` 项严格优）与 `S1` 含非可重复靶点时的稳健性退化（`4` 项劣效）
   - 4.4 `S7`（`P6`）细胞数与敲低深度混杂控制、`SA-1..SA-5`（含 `all_nontarget`、`R19` 及 `SA-3` R9 锚点迁移表现）
7. **Section 5: Post-Hoc Exploratory Analysis — Deconstructing `FCT` Without Target Baselines (`[R10-E PLACEHOLDER]`)**
   - *严格独立章节与全套免责声明，仅在 `A-T2` / `A-T3` / `A-T5` 全部签收后填充 `E1–E7` 正式表图、配对差 `FCT - E1` 与文献边界 `PRIOR_ART_D.md`。*
8. **Section 6: Discussion & Methodological Synthesis**
   - 6.1 两轮合卷视角下重估 VCC 2026 `0.02–0.13` 斜率瓶颈
   - 6.2 对 `CIPHER`、`CASTER`、`Molina & Zhang`、`Systema`、`CGC`/`AMM-SimWMag` 与 `arXiv:2608.00152v3` 的定位与互补
   - 6.3 **文献先验与方法归属核查汇总**（`PRIOR_ART_A–C` + `[PLACEHOLDER: PRIOR_ART_D]`）
   - 6.4 面向未来真值公布后的端到端前瞻基建（`A-T4` `cell-eval2` 本地复现、`vcc prep`/评分桥接及 VCC 2026 合规限制）
9. **Section 7: Limitations**（合并 R9 + R10 + `[R10-E Placeholder]`）
10. **Section 8: Data, Code, Reproducibility, and References**

---

## 4. 完整图表清单（Figure & Table Inventory）

### 4.1 主图清单（Main Figures）

| 编号 | 图题与核心展示内容 | 覆盖轮次 / 数据 | 关键信息与标注要求 |
| :--- | :--- | :--- | :--- |
| **Figure 1** | **Two-Round Pre-Registered Blinding Design & Geometric Decomposition** | R9 (`6` lines) + R10 (`4` lines + `H9` P2 excl.) | (A) 两轮 6 步锁定闸门与零泄漏信息流（R9 LOLO Dev $\to$ Held-Out；R10 6 源对 4 新线 + `P2` 预算/密度门控）；(B) 通道/双供体（`CE0006864` vs `CE0008162`）Split-Half 分割；(C) `r_orth_cal` 几何示意图（残差向 $\hat{\Delta}^{\mathrm{TMean}}_p$ 正交投影与 Spearman–Brown 去衰减）。 |
| **Figure 2** | **Scalar Amplitude Confounding (`R18`), Three-Tier Decomposition (`R19`), and `resp_top50` Ceiling Bias (`S6`)** | R9 (`6` 线) + R10 (`4` 线) | (A) `Transfer_Mean_GlobalScale` 原始 `r_cal`（R9 `[-0.330, +0.360]`、R10 `dld1 +0.240`/`cd4t -0.262`）对比 `r_orth_cal = 0.000`；(B) `R19` 堆叠条形图（`all_nontarget` 下 `L1` 占 `1.6%–51.3%`、两轮确认线 `11.3%–34.8%`，`resp_top50` 人为抬升 `L2` 并掩盖 `L1`）；(C) `S6` 噪声天花板 $\rho^{\perp}_{AB}$ 从 `resp_top50` 降至 `all_nontarget`（如 `dld1 0.828 -> 0.253`、`cd4t 0.648 -> 0.147`）。 |
| **Figure 3** | **Round 9 Confirmatory Results: Group A (`Hybrid_ProjShared`) and Group B (`TMVS`) Are Not Confirmed** | R9 (`K_conf = 3` 主 + 3 排除线) | (A) `r_orth_cal [97.5% CI]` 与 `slope_ols` 对标 5 条错配对照组 & `100` 次置换分布小提琴图（`H-A1..A3`、`H-B1..B3` 全标明 PASS/FAIL）；(B) `R14`（`Hybrid_on_TMVS - TMVS` `95% CI`）与 `SrcCorrTgtVar` 机制诊断；(C) `H-A4` / `H-B4` 6 口径全劣效森林图；(D) `S1–S4b` & `S2b`（秩 1 `q_1`）敏感性热图。 |
| **Figure 4** | **Round 10 Prospective Audit: `FCT` Passes Positive Signal (`H-C1`) and Non-Inferiority (`H-C3`) Yet Fails Cell-Line Specificity (`H-C2`)** | R10 (`dld1` `N=86`, `cd4t_rest` `N=302` + `calu3`, `i3Neuron`) | (A) `FCT` `r_orth_cal [97.5% CI]` vs. 全部 6 个错配源线均值（最高 `H1 +0.1092`、`Jurkat +0.0761`）与 `100` 次基因置换零分布（`p97.5 +0.1179`、`+0.0554`），清晰标示 `H-C1 = 2/2 PASS`、`H-C2 = 0/2 FAIL`；(B) `SA-3` R9 锚点方法在 R10 新线上的表现；(C) `H-C3` 6 项评价口径与 `S1` 全部靶点敏感性差值森林图（`FCT - Transfer_Mean`）。 |
| **Figure 5** | **Cell-Count / Knockdown-Depth Controls (`S7`) and Cross-Round Sensitivity Robustness** | R9 + R10 确认线 | (A) `dld1` 与 `cd4t_rest` 上单扰动 `r_orth_p` 对 $\ln(n_{\mathrm{cells}})$ 与敲低深度 $|k_p|$ 散点图 + OLS 回归线（展示 `Spearman` 无关/弱负相关及稳健正截距 `+0.1400`、`+0.0607`）；(B) 按 `n_cells` 三分位（`T1/T2/T3`）分层比较 `FCT` vs. 错配与置换阈值；(C) R10 `Primary` vs `S1` vs `S5` vs `all_nontarget` 一致性矩阵。 |
| **Figure 6** | **`[POST-HOC EXPLORATORY — PLACEHOLDER]` Baseline Ablation Deconstructing `FCT` (`R10-E`)** | **[待 R10-E 填充]**（R9 6 线 + R10 4 线；`resp_top50` & `all_nontarget`） | **必须标注「事后探索性（post hoc）」**。[预留 `(A)` `FCT` vs `E1 (median)`、`E2 (p10/p90)`、`E3 (src mean)`、`E5 (rank-shuffled)`、`E6 (log FC)`、`E7 (FC-1)` & `perm100` 在 10 条线上的 `r_orth_cal [97.5% CI]` 对比；`(B)` 配对差 `FCT - E1 [97.5% CI]` 森林图；`(C)` `resp_top50` vs `all_nontarget` 对照检验是否排除选择偏倚。] |

### 4.2 主表与附表清单（Main & Supplementary Tables）

| 编号 | 表题与核心列 | 状态 |
| :--- | :--- | :--- |
| **Table 1** | **Unified Dataset Inventory Across Rounds 9 & 10**（细胞系、研究出处/DOI/GEO/IGVF、角色、原始/质控单 sgRNA 细胞数、NTC 样本数、Split-Half 物理轴、`\|G_eval\|`、可重复与准入靶点数 `N`，含 `H9` P2 下载前排除披露） | 数值已完备（来自 R9 Table 1 + `R10_QC_COUNTS.md` + `DATA_SOURCES_R10.md`） |
| **Table 2** | **Round 9 Main Confirmatory Adjudication (`PREREG_A v4`, `K_conf = 3`)**（`Hybrid_ProjShared`、`TMVS`、`GlobalScale`、`R14`、`H-A1..A4` / `H-B1..B4` 与 3 条排除线对比） | 数值已完备（来自 R9 Table 2–4） |
| **Table 3** | **Round 10 Main Confirmatory & Secondary Adjudication (`PREREG_C`, `K_conf = 2`)**（`dld1`、`cd4t_rest` + `calu3`、`wtc11_i3neuron`：`FCT`、`Max Mismatch CI_hi`、`Perm100 med/p97.5`、`H-C1/H-C2/H-C3`、`slope_ols`、`S7` OLS 截距与 Spearman、`all_nontarget`、`SA-3` 锚点） | 数值已完备（来自 `R10_VERDICT.md` + `AUDIT_R10_REPORT.md` + `A-T2.txt`） |
| **Table 4** | **`H-C3` Parallel Non-Inferiority Across All Six Conventions (`FCT - Transfer_Mean`)**（`ctrl_disc_ownmasked` / `ctrl_disc_panel` / `ctrl_de50_r` / `cent_*` 及 `pds_cosine_panel`，注明 `D-R10-1` 与代码变量命名对照，附 `Primary` vs `S1` 对比） | 数值已完备 |
| **Table 5** | **Cross-Round `R19` Three-Tier Multiplicative Decomposition & `S6` Noise-Ceiling Comparison**（10 线 `resp_top50` vs `all_nontarget` 的 $\beta_t$、`L1/L2/L3` 及 $\rho^{\perp}_{AB}$） | R9 与 R10 `all_nontarget`/`top50` `L1/L2/L3` 与 $\rho^{\perp}_{AB}$ 数值已备（注：可由 `A-T2`/`A-T3` summary JSON 补全 R10 全 4 线 $\beta_t$ 与 `calu3`/`i3Neuron` `L1/L2/L3`） |
| **Table 6** | **`[POST-HOC EXPLORATORY — PLACEHOLDER]` `R10-E` Baseline Ablation (`E1–E7`) & Paired Difference (`FCT - E1`) Matrix** | **[严格留空占位]**；每张 E 表均须显著标注 `事后探索性（post hoc）`，待 `round10/exploratory/BASELINE_ABLATION.md` 与 `round10/audit/exploratory/` 交叉核对（点估计差 `< 1e-6`）完成后填入 |
| **Supp. Tables S1–S5** | `S1`: R9 全部敏感性表（`S1, S2, S2b, S3, S4a, S4b, S5`）；`S2`: R9 事后 `FCT` 表 & R10 `S1/S5` 完整表；`S3`: 两轮双组 `C3` 交叉核对汇总 ledger；`S4`: 文献查新矩阵（`PRIOR_ART_A–C` + `[PLACEHOLDER: PRIOR_ART_D]`）；`S5`: `A-T4` `cell-eval2 v0.22.0` VCC 2025 公开复现对照分与前瞻时间戳校验格式规格 | R9/R10 部分数值完备；`PRIOR_ART_D` 槽位留空待 `R10E_LIT_DONE` |

---

## 5. 方法与偏差小节细纲（Detailed Methods & Deviations Outline）

### 5.1 数据获取、准入闸门与基因空间统一（Data Provenance & Blinding Architecture）
1. **Round 9 数据集（6 线，`1,419,744` 细胞）**：`nadig_hepg2`（`GSE264667`，`145,473` 细胞，`N=799`）、`nadig_jurkat`（`GSE264667`，`262,956` 细胞，`N=934`）、`vcc2025_h1`（GCS VCC 2025 训练+验证+测试合并去重 `76,352` 重复 NTC 后保留 `414,694` 细胞，`N=80`）作为 `K_conf = 3`；`replogle_k562_essential`（`310,385` 细胞，`N=859`/`832` dev）与 `replogle_rpe1`（`247,914` 细胞，`N=1,070`/`832` dev）限定为开发集；`nourreddine_kolf21j`（`38,322` 细胞，`N=18 < 50`）作为低样本次要队列排除出确认裁决。
2. **Round 10 数据集（`34.080 GB <= 60.000 GB` 下载总额，全 6 条 R9 线固定作为外部源池 $\mathcal{S}$）**：
   - `dld1`（Yeung et al. 2026，`GSE337988`，仅取 Low-MOI `SAM24482910_1..24`，严格单 sgRNA 过滤后 `701,121` 细胞 [`685,008` targeting + `16,113` NTC]，`12` 奇对 `12` 偶通道分割，`4,934` 去重靶点 $\to$ `202` 可重复 $\to$ `93` 源线共享 $\to$ **`86`** `n_eval_targets`）。
   - `cd4t_rest`（Zhu et al. 2026 *Cell*，`GSE314342` 静息态，下载前按字典序锁定前 2 个人类供体 `CE0006864` [`D4`, `CD4i_R2`, 24 lanes] 与 `CE0008162` [`D1`, `CD4i_R1`, 23 lanes] —— 规程 `P3`，`3,502,875` 单 sgRNA+NTC 细胞 [`3,346,822` targeting + `156,053` NTC]，`11,903` 去重 $\to$ `3,675` 可重复 $\to$ `398` 共享 $\to$ **`302`** `n_eval_targets`，其中 **`284`** 个在两个供体中各具 $\ge 2$ 细胞）。
   - `wtc11_i3neuron`（Tian et al. 2021，`GSE152988`/Zenodo `13350497`，仅 CRISPRi、`31,407` 细胞、`N=2`）与 `calu3`（Sunshine et al. 2023，`GSE208240`/Zenodo `13350497`，去双靶向 `2,667` 细胞后 `29,370` 细胞、`N=8`）预锁定为仅次要报告。
   - **`h9_neuron`（IGVF `IGVFDS9825FOBP` / `PMC13629240`）下载前排除（`P2`）**：未下载任何表达矩阵。依据轻量元数据审计，`MOI ~ 4` 使单 HTO+单 sgRNA 细胞仅占 `5.53%`（`53,037 / 958,228`，均值约 `~19` 个质控后单 sgRNA 细胞/靶基因，不足 `>= 25` 门槛），且强制配对 `.h5`+`.pkl` 需 `50.999 GB`（与现有 `34.080 GB` 合计 `85.079 GB` 突破 `60.000 GB` 预算上限）。
   - **HGNC 跨轮统一桥接**：通过 `hgnc_complete_set.txt` 将 R9 `COMMON_GENES.txt` 中的 `164` 个旧 `prev_symbol` 与 R10 新注释映射到统一标准符号（`DLD-1`: `|G_eval| = 6,124`、`83 -> 86` 靶点，找回 `DARS->DARS1`、`NARS->NARS1`、`TAZ->TAFAZZIN`；`CD4+ T`: `|G_eval| = 6,123`、`301 -> 302` 靶点，找回 `H2AFZ->H2AZ1`）。
3. **机械确认判定门（`P1`）与单 sgRNA 锁定（`P4`）**：`A-T1` 在代码冻结前仅输出数量质控表 `R10_QC_COUNTS.md`；满足 `n_eval_targets >= 50` 且 `|G_eval| >= 3,000` 才进入确认集（`K_conf = 2`），且代码中强写 `K_conf < 2` 自动中止确认裁决分支。

### 5.2 预测模型与对照框架（Predictors & Nulls）
- **Round 9**：`Transfer_Mean`、`Source_Grand_Mean`、`CIPHER_Linear`、`Hybrid_Raw`、`Hybrid_ProjShared`（Group A）、`Transfer_Mean_VarScaled` (`TMVS`, Group B)、`Transfer_Mean_GlobalScale`（`R18`）、`Hybrid_on_TMVS`（`R14`）及 `MismatchedCov(s')` / `DiagCov` / `ShuffledGenes` / `SrcCorrTgtVar` / `MismatchedVar(s')` / `perm100`。
- **Round 10**：`FCT`（几何平均源线倍数变化映射至目标线未扰动均值的对数增量 `log1p(mu_t * FC_bar) - log1p(mu_t)`，明确无算法新颖性声明，归属 `CellularAI`/`AtlasShift` 同类基准）、`6` 条 R9 源线的 `NegCtrl_FCT_MismatchedMean(s')`、`100` 次置换 `NegCtrl_FCT_Perm(k)`（仅置换 $\mu_t$ 基因顺序，`default_rng(k)`）、`Transfer_Mean_GlobalScale` 与 `SA-3` R9 锚点对比。
- **`[R10-E PLACEHOLDER]` 探索性无目标基线消融家族（待 `A-T2` / `A-T3` 复核后填入定义检验细节与数值）**：
  - `E1`（目标对照均值中位数常数基线）、`E2`（10% 与 90% 分位数常数基线）、`E3`（源线对照均值的平均）、`E5`（目标均值分布按源线排序重排）、`E6`（$\ln \overline{\mathrm{FC}}$ 常数趋无穷极限）、`E7`（$\overline{\mathrm{FC}} - 1$ 常数趋零极限）及配对差 `FCT - E1` 联合 bootstrap（`n_boot = 4,000`, `seed = 20261122`, `resp_top50` 与 `all_nontarget` 双掩码）。

### 5.3 统计检验、混杂控制与三方复核（Endpoints, Confounder Controls, & Cross-Checks）
- **去幅度相关 `r_orth_cal`**：每扰动剔除自身靶基因 $g_p$，在掩码内（主 `resp_top50` / 敏感性 `all_nontarget` `S6` / 前 75% `S5`）对 $r^{\mathrm{obs}}_p$ 与 $r^{\mathrm{pred}}_p$ 正交投影去除 $\hat{\Delta}^{\mathrm{TMean}}_p$ 方向，池化计算 Pearson 相关并以联合 bootstrap（`4,000` 次，`default_rng(20261121)` for R9 / `default_rng(20261122)` for R10，每次重抽样同步重算 $\rho^{\perp}_{AB}$）的 Spearman–Brown 天花板校正；Bonferroni 双侧 `97.5%` 置信区间与 `K_conf / K_conf` 全票规则。
- **`H-C3` 并行实用性**：`FCT - Transfer_Mean` 在 6 种口径下检验非劣效（CI 上界 `>= 0`）。
- **`S7`（`P6`）混杂控制**：针对 `arXiv:2608.00152v3` 发现的细胞数纠缠，计算单扰动 `r_orth_p`、`Spearman(r_orth_p, ln n_cells)`、`Spearman(r_orth_p, |k_p|)`、中心化 OLS `[1, ln n_cells, |k_p|]` 截距 97.5% CI，以及 `n_cells` 三分位（`T1/T2/T3`）分层检验。
- **跨组复算（`C3`）**：R9 `A-T2` vs `A-T3`（`456` 项确定性零偏差 `<= 1.11e-14`，`S2b` `90` 项 `<= 1.11e-16`，`Claude` 独立复算 Group B `TMVS`）；R10 `A-T2` vs `A-T3`（主分析 `349` 项点估计最大差 `6.66e-16`、CI 最大差 `1.07e-05`；`S1` `364` 项差 `<= 3.89e-16`）并与主控自写脚本（`1,000` bootstrap，种子 `777`）判定完全一致。

### 5.4 与预注册的偏差（Deviations from Pre-Registration — 合并披露小节）
1. **`R9 D1`（`Hybrid_ProjShared` 投影子空间维度，`DEVIATIONS.md`）**：`PREREG_A v4 §4.2` 文字写为秩 1 $P_\perp = I_G - q_1 q_1^\top$，而冻结代码 `1ca9f60f`（`src/predictors.py::fit_shared_subspace` 及 `A-T3` `audit/independent_predictors.py`）实际去除的是 $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$ **三维子空间**。差异在冻结前即存在、审稿时（`PAPER_REVIEW_R1.md B1`）发现；R9 主判定以冻结代码为准，`A-T2`/`A-T3` 均补跑了字面秩 1 敏感性 **`S2b`**（点估计差 `1.11e-16`），结果与主运行完全一致（Group A 仍为 `未确认`，`R14` 仍为 `False`），且该点上 `A-T2`/`A-T3` 实现并不独立；Round 10 `PREREG_C P7` 已在预注册阶段明确锁定沿用 `1ca9f60f` 三维实现。
2. **`R10 D-R10-1`（`H-C3` 六种口径变量命名差异，`CODE_REVIEW_R10.md`）**：`PREREG_C §4.2` 文本使用 `ctrl_disc_ownmasked` / `ctrl_disc_panel` 与 `cent_disc_ownmasked` / `cent_disc_panel`，冻结代码沿用 R9 字段名 `ctrl_disc_raw` / `ctrl_disc_raw_panel` 与 `cent_disc_raw` / `cent_disc_raw_panel`。两者数学定义完全一致，论文正表统一采用预注册命名并注明一对一映射。
3. **`R10 D-R10-2`（`H-C1` 额外保守保护条件，`CODE_REVIEW_R10.md`）**：`round10/src/adjudicate_r10.py` 沿用 R9 要求除 `CI_lo > 0` 外还需 `valid_frac >= 0.95` 与 `rho_ab_orth > 0`。两条确认线均严格满足（`valid_frac = 1.0`，`rho_ab_orth = +0.8277 / +0.6479`），与审计线判定无任何分歧。
4. **`R10` 代码细微差异与前置规则修正记录**：`fct_from_fc_bar` 目标均值 $\mu_t$ 下限写成 `0` 而文本写为 `1e-8`（差异仅 `~1e-8*(FC-1)`，可忽略，已在代码审查记录在案）；`H9` 下载前排除原因由主控校订记录为 P2 双重理由（单 sgRNA 靶点过少 + 超出 60 GB 预算），详见 `CONFIRMATORY_LINES.txt` 与 `R10_QC_COUNTS.md`。

---

## 6. 局限性（Combined Limitations）

1. **执行智能体同源性及部分实现非独立性**：两轮核心执行 `A-T2` 与独立审计 `A-T3` 均由 `Antigravity`（Google Gemini）完成。虽有主控 `Claude`（Anthropic）另写代码独立复算 R9 Group B `TMVS` 及 R10 两条确认线的 `FCT` 与对照指标（R9 `HepG2/H1` `FCT` 与 R10 全确认集），但 R9 Group A（`Hybrid_ProjShared`）未由跨模型主控复算，且在 R9 `1ca9f60f` 中 `A-T3` 复刻了与 `A-T2` 相同的三维子空间投影（`D1`）。
2. **确认性细胞系数量有限（R9 `K_conf = 3`；R10 `K_conf = 2`）且谱系覆盖受限**：全票通过规则（`100%`）与 Bonferroni `97.5%` 区间以严格控制假阳性为代价降低了检出效能（如 R9 中 `vcc2025_h1` 仅 `N = 80` 个共享靶点导致方差大幅放宽）。在 R10 中，`H9` 神经元因高 MOI 与 `60 GB` 预算被预排除（`P2`），`wtc11_i3neuron`（`N = 2`）与 `calu3`（`N = 8`）受共享可重复靶点限制无法作为确认组，最终 R10 仅由一条癌细胞系（`DLD-1`）与一种原代静息免疫细胞（`CD4+ T`）构成，缺乏神经元或心肌系确认样本。
3. **`cd4t_rest` 跨供体+跨批次 Split-Half 放大噪声天花板分母与样本折损（`P3`）**：由于预锁定按两个不同人类供体（男性亚洲裔 `CE0006864` / `CD4i_R2` vs. 女性拉丁裔 `CE0008162` / `CD4i_R1`）划分两半，生物学供体差异及实验批次差异一并进入 $\rho^{\perp}_{AB}$（`0.6479`，低于通道分割的 `dld1` `0.8277`），且使 `302` 个评测靶点中仅 `284` 个同时在两供体各有 $\ge 2$ 细胞，增加了去衰减估计方差。
4. **共享基因空间与可重复滤波的外推边界**：R10 沿用 R9 确定的 `6,148` 交集基因（经 HGNC 桥接后 `|G_eval| = 6,124` / `6,123`），偏向高表达保守基因；且主评测限定经过 Addendum-B 筛选的单 sgRNA 强效应扰动（`n_cells >= 25`）。当 `S1` 放宽至全部靶点（`dld1` `N=500`，`cd4t_rest` `N=1,422`）时，`FCT` 信号减半（`+0.0717`、`+0.0314`）且 `H-C3` 在质心口径出现 `4` 项劣效。
5. **`resp_top50` 基因选择偏倚（`S6`）尚待在 `R10-E` 消融中彻底剥离**：按真实 $|\Delta^{\mathrm{obs}}|$ 排序筛选前 50 基因虽与主流差异基因评价相近，但会严重高估噪声音花板（R10 `dld1` $\rho^{\perp}_{AB}$ 从 `0.8277` 降至 `0.2533`、`cd4t_rest` 从 `0.6479` 降至 `0.1468`；R9 同样大幅下降）并放大 `R19` 逐基因增量 `L2`。尽管 R10 主 `FCT` 在 `all_nontarget` 下的点估计与 `resp_top50` 几乎重合（`0.1355` vs `0.1401`；`0.0556` vs `0.0560`），`[R10-E PLACEHOLDER: E1–E7 在 all_nontarget 下是否依然复现 FCT 是排除选择偏倚的关键，待 A-T2/A-T3 复核结果写入]`。
6. **VCC 2026 外部官方指标闭源与测试集真值未释出（`R6`, `A-T4`）**：`ArcInstitute/cell-eval2 v0.22.0` 公开源码未实现 Wrap-Up 所述 `0.02–0.13` 跨线斜率计算公式（参照 `10.64898/2026.11.17.767214` 尚未被公开文献库索引）；同时 `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` 截至 2026-11-22 仍为空且测试集 `controls-final` 受竞赛条款（`/terms`）限制，故无法在 `AC16`、`MOLM13`、`U-87MG` 上开展开源评测或公开生成基于受限对照的预测。
7. **`[POST-HOC R10-E & FUTURE R11 PLACEHOLDER]`**：任何来自 `R10-E` 的结论均属于看到第十轮结果之后发起的**事后探索性发现**，仅能作为假设生成，在独立新数据集上完成新的锁定预注册检验（计划 Round 11，`[Placeholder: ref R11_DATA_SCREEN.md once A-T1 finishes]`）之前，不可升级为确认性结论。
8. **文献查新边界（`[PLACEHOLDER: ref PRIOR_ART_D.md]`）**：所有“未见先行研究”判断均严格受限于检索日期与公开学术 API 覆盖范围（`PRIOR_ART_A–C` `[and PRIOR_ART_D.md once A-T5 completes R10E_LIT_DONE]`）。

---

## 7. 待完成与交接核对清单（Wait-Gate Before Drafting Combined Manuscript）

- [ ] **A-T2 (`R10E_ABLATION_DONE`)**：完成 `round10/exploratory/BASELINE_ABLATION.md` 与 CSV（`E1–E7`、配对差 `FCT - E1`、`resp_top50` & `all_nontarget`、`4,000` bootstrap）$\to$ 填充本文 §2.2、Fig. 6、Table 6。
- [ ] **A-T3 (`R10E_AUDIT_DONE`)**：完成 `round10/audit/exploratory/` 对 `E1, E6, E7, FCT - E1` 的独立复核（点估计差 `< 1e-6`）$\to$ 写入 C3 核对说明。
- [ ] **A-T5 (`R10E_LIT_DONE`)**：完成 `round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_D.md` $\to$ 填充 §1、§6.3 与 Supp. Table S4 中关于「去目标线基线 / 几何平均非线性压缩 / 无情境基线复现」的先行文献定位。
- [ ] **A-T1 (`R11_SCREEN_DONE`)**：完成 `round10/shared/R11_DATA_SCREEN.md` $\to$ 在 §6.4 / §7 补充后续第十一轮前瞻检验候选数据集规模估算。
- [ ] **Table 5 补全（可选快速核对）**：待确认需要时从 `round10/results/summary_holdout_r10_log1p_primary.json` 抽取 `calu3`、`wtc11_i3neuron` 与 `dld1`/`cd4t_rest` 完整 $\beta_t$ 与 `resp_top50`/`all_nontarget` `L1/L2/L3` 双组值填入 Table 5。

---

## 8. 核心参考文献目录（已验证引文，待合稿定稿按正文顺序重排编号）

- **Round 9 已有验证文献**：Adduri et al. (2026) *Cell* (`10.1016/j.cell.2026.08.004`); Arc Institute VCC 2026 Wrap-Up (`2026-11-19`) & `cell-eval2` (`v0.22.0`); Dhaka et al. (2026) `HyperMap` (`10.64898/2026.04.23.720505`); Huang et al. (2026) `CGC` (`10.64898/2026.09.30.755756`); Kawsar et al. (2026) `CASTER` (`10.64898/2026.11.15.766815`); Kuznets-Speck et al. (2025/2026) `CIPHER v2` (`10.1101/2025.06.27.661814`); Luan et al. (2026) (`10.64898/2026.10.20.760514`); Miller et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03307-w`, `PMID:42823510`); Molina & Zhang (2026) (`10.64898/2026.07.24.740459`); Nadig et al. (2025) *Nat. Genet.* (`10.1038/s41588-025-02169-3`, `GSE264667`); Nourreddine et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03199-w`); Replogle et al. (2022) *Cell* (`10.1016/j.cell.2022.05.013`); Shoeibi & Yousefi (2026) *arXiv*:`2608.00152v3`; Sun et al. (2026) `AMM-SimWMag` *Genes* (`10.3390/genes17070816`); Viñas Torné et al. (2025) `Systema` *Nat. Biotechnol.* (`10.1038/s41587-025-02777-8`).
- **Round 10 新增数据集与元数据出处**：Yeung J et al. (2026) `DLD-1` (*bioRxiv*, `doi:10.64898/2026.07.10.737863`, `GSE337988`); Zhu X et al. (2026) `CD4+ T` (*Cell*, `doi:10.1016/j.cell.2026.08.002`, `PMID:42664972`, `GSE314342`); Tian R, Kampmann M et al. (2021) `WTC11 i3Neurons` (*Nat. Neurosci.*, `doi:10.1038/s41593-021-00862-0`, `PMID:34031600`, `GSE152988`/Zenodo `13350497`); Sunshine S, Hein MY, Weissman JS et al. (2023) `Calu-3` (*Nat. Commun.*, `doi:10.1038/s41467-023-41788-4`, `PMID:37803001`, `GSE208240`/Zenodo `13350497`); `IGVFDS9825FOBP` / `PMC13629240` (`H9` 元数据预算/MOI 下载前审计依据).
- **`[PLACEHOLDER]`**：补充 `PRIOR_ART_D.md` 新查核文献（含竞争方案详情及基线批评文献如 `scGPT`/`GEARS` 批评类论文，待 `A-T5` 核实后写入）。
