# 第九轮 + 第十轮（含 R10-E）+ 第十一轮（描述性报告）合稿论文骨架（`round10/paper/OUTLINE.md`）

- **定位与核心叙事（Round 9 $\to$ Round 10 $\to$ R10-E Post-Hoc $\to$ Round 11 `K_conf = 1` Descriptive）**：
  1. **第九轮（`PREREG_A v4`，6 线 LOLO `K_conf = 3`）**：预注册的协方差共波动转移（`Hybrid_ProjShared`）与边际方差缩放（`TMVS`）方法**均未确认**（`0/4` 标准）；共波动在边际方差之外未提供可检测增量（`R14`）。
  2. **第十轮（`PREREG_C`，4 条全新单 sgRNA 人类线 `K_conf = 2`）**：与 VCC 2026 #1/#82 同类的基线倍数缩放模型 `Fold-Change Transfer`（`FCT`）在预注册主掩码 **`resp_top50`** 下通过去幅度正信号（`H-C1 2/2`）与并行非劣效（`H-C3 12/12`），但**未通过细胞系特异性检验（`H-C2 0/2 -> 未确认`）**。
  3. **第十轮事后消融与 `13:15` 主控更正（`R10-E`，严格标注事后探索性 *post hoc*）**：双团队正式复算证明（`510` 项最大差 `1.11e-16`），`resp_top50` 按目标观测 $|\Delta^{\mathrm{obs}}_t|$ 选基因造成了选择偏倚，使无目标基线信息的常数基线 `E1` 也能复现 `86.4%–115.5%` 的 `FCT` 分数并抬高置换阈值；而在无偏全非靶基因口径 **`all_nontarget` (`S6`)** 下，`E1/E2/E6/E7` 与置换对照塌缩至零附近或负值，真实 `FCT` 在全部 `8/8` 有效线（`N >= 18`）保持正值且配对差 `FCT - E1` 97.5% CI 全胜（`8/8`）。
  4. **第十一轮（`PREREG_D v2 + Addendum A`，`K_conf = 1` —— 严格按 `P1` 仅作描述性报告）**：在全新锁定候选池上检验全基因口径 `FCT` 与 `resp_top50` 选择陷阱。因 `hcasmc` 公开表达矩阵缺失细胞级别 guide 分配而在评估前依 `Q3` 可行性排除、`telohaec`（`N = 9 < 50`）受 GWAS 文库低敲低效应制约未过门槛、`feng_ipsc`（`N = 93`）预锁定为次要线，合格确认性线仅余 **`cd4t_stim48`** 一条（**`K_conf = 1`**）。依据 `P1` 规则，**第十一轮确认性检验因合格细胞系不足而没有进行，所有结果均为描述性结果**；经 `A-T5` 查新（`PRIOR_ART_E.md`）、`A-T3` 独立复核（`round11/audit/exploratory_check/REPORT.md`）及主控裁定（`EXPLORATORY_REVIEW_R11.md §D`），在大纲末尾设独立节**『事后探索』**仅收入 3 条事后发现与 1 条机制局限，而机制推导、解药类措辞及 `Q5a` 守恒论断均不写入正文。
- **锁定与版本锚点**：
  - **Round 9**：`PREREG_A v4`（commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`），`CODE_FREEZE`（`1ca9f60f49adf1ac7cce8266e40533273411db30`）。
  - **Round 10 & R10-E**：`PREREG_C`（commit `dc4b4a5`），`QC_LOCK`（`86872ed`），`CODE_FREEZE_R10`（`f8d2ebfc`，`22/22` 冻结文件 `OK`），`R10_VERDICT.md`（含 `13:15` 更正），`BASELINE_ABLATION.md`，`AUDIT_R10E_REPORT.md`，`PRIOR_ART_D.md`（Q1–Q5 完整）。
  - **Round 11**：`PREREG_D v2`（`commit 72fc0c937358b089e78664b464f9db5aafb790b1`，`sha256 1a054498...`）+ `PREREG_D_ADDENDUM_A.md`，`CODE_FREEZE_R11`（`cc834d863105c66874acb4c27cc214fd79704955`，`41/41` 冻结文件校验 `OK`），`DEVIATIONS_R11.md`（`D1–D4`），`CONFIRMATORY_LINES_R11.txt`（`K_conf = 1`），`R11_VERDICT.md`（`2026-11-22 21:05`），`EXPLORATORY_REVIEW_R11.md`（含 `21:55 §D` 补充），`round11/shared/PRIOR_ART_E.md`，`round11/audit/exploratory_check/REPORT.md`（`21:42`）。
- **AI 治理与作者披露声明（必保留；人类角色仅称 `Human project lead` / `人类项目负责人`）**：人类项目负责人 Zijian Zeng（Human project lead Zijian Zeng，无正规生物学背景）；主控 `Claude`（Anthropic）负责预注册审核、结果核验、偏差审查及独立复算；执行智能体 `Antigravity`（Google Gemini 模型）负责数据整合（`A-T1`）、主分析与消融（`A-T2`）、独立复核与 `C3` 交叉核对（`A-T3`）、VCC 2026 合规/前瞻管线与合稿骨架（`A-T4`）及 API 文献/源码审计（`A-T5`）。

---

## 1. 标题候选（Title Candidates）

> **措辞原则**（沿用 `PAPER_REVIEW_R1.md B3` 与 `R11_VERDICT.md`）：不使用过度绝对化的英文动词（如「驳倒/否定/证实」类词汇）；如实反映 R9/R10 未确认及 R11 因 `K_conf = 1` 仅作描述性报告与 `H-D2` 对 H1 基线仍未通过、总效应量微小的事实。

1. **三轮全覆盖保守写实型（主推）：**  
   *Multi-Round Pre-Registered Audits of Zero-Shot Control-Conditioned CRISPRi Predictors Reveal Modest Baseline Signal, Outcome-Selection Bias, and Unconfirmed Cell-Line Specificity*
2. **方法学 + 评估诊断聚焦型：**  
   *Separating Amplitude, Selection Bias, and Baseline Scaling in Zero-Shot CRISPRi Perturbation Prediction: Confirmatory and Descriptive Audits Across Human Perturb-seq Atlases*
3. **直面 VCC 2026 基准瓶颈型：**  
   *Why Zero-Shot Line-Specific CRISPRi Generalization Remains Unestablished: Co-Fluctuation Failures, `top-K` Mask Bias, and Fold-Change Baseline Audits*

---

## 2. 摘要骨架（Abstract Skeleton）

> **强制约束**：严格区分（1）R9 & R10 确认性结论；（2）R10-E 事后探索发现；（3）R11 锁定描述性结果（**必须写清 `PREREG_D Q7` 与 `R11_VERDICT.md §8` 全部必报项**）。摘要中**绝不出现** `exploratory/null_cone_theorem`（含 `N_5`、一阶泰勒启发式推导等）内容。

### 2.1 第一段：第九轮与第十轮已锁定确认性判定（Round 9 & Round 10 Confirmatory Outcomes）
1. **背景与 `R18` 标量混淆**：单靠目标细胞系未扰动对照零样本预测 CRISPRi 细胞系特异响应在 VCC 2026 中斜率偏低。每细胞系单一标量缩放（`Transfer_Mean_GlobalScale`）即可制造高原始残差相关（R9 确认集 `-0.126` 至 `+0.250`，全 6 线 `-0.316` 至 `+0.352`；R10 `dld1 +0.240`，`cd4t_rest -0.262`），而在单扰动内对共享跨线响应 $\hat{\Delta}^{\mathrm{TMean}}_p$ 正交消幅（`r_orth_cal`）后精确归零（`0.000`）。
2. **第九轮（`PREREG_A v4`，6 线 `1,419,744` 细胞、`6,148` 基因，LOLO `K_conf = 3`）**：在发展集 `Replogle K562/RPE1` 锁定、排除低靶点 `KOLF2.1J`（`N = 18`）后，盲态评测 `nadig_hepg2`（`N = 799`）、`nadig_jurkat`（`N = 934`）与 `vcc2025_h1`（`N = 80`）。按 Bonferroni `97.5%` 置信区间与 `3/3` 一致性，**Group A（`Hybrid_ProjShared`）与 Group B（`TMVS`）均未确认（`0/4` 标准）**；共波动结构在边际方差外无额外可检测增量（`R14`，`2/3`）。
3. **第十轮（`PREREG_C`，4 条全新单 sgRNA 数据集 `4,264,773` 细胞，**`K_conf = 2 confirmatory cell lines`**）**：`H9` 因高 MOI 导致单 sgRNA 稀疏与超出 `60 GB` 预算在下载前排除（`P2`），`calu3`（`N = 8`）与 `wtc11_i3neuron`（`N = 2`）作次要线，确认集由 `dld1`（`N = 86`）与双供体 `cd4t_rest`（`N = 302`，`284` 跨供体两半有效）组成。在主掩码 `resp_top50` 下，基线倍数映射模型 `Fold-Change Transfer`（`FCT`，同属 VCC 2026 #1/#82 缩放类别）通过 **`H-C1`（`2/2`）**（`+0.1401`、`+0.0560`，经 `S7` 控制 `ln n_cells` 与敲低深度后截距保持 `+0.1400`、`+0.0607`）和并行 **`H-C3` 非劣效（`12/12`）**，但**`H-C2` 未确认（`0/2`）**（两条线均未同时超过 6 个错配对照上界与 100 次基因置换阈值），因此 **Group C 主论断在 `PREREG_C` 下未确认（NOT CONFIRMED）**。

### 2.2 第二段：事后基线消融与掩码选择偏倚发现（Post-Hoc R10-E）
- **明确标注 *Post-Hoc Exploratory***：双团队独立消融（`510` 项核对最大差 `1.11e-16`）推翻了最初“纯变换解释 `FCT`”的错觉，揭示了关键掩码分叉：在 `resp_top50`（按目标观测 $|\Delta^{\mathrm{obs}}_t|$ 选前 50 基因）下，不含目标线基因差异信息的常数基线 `E1` 复现 `FCT` 的 `86.4%–115.5%` 并抬高 `perm100`；在无偏全非靶基因掩码 `all_nontarget` (`S6`) 下，`E1`/`E2`/`E6`/`E7` 与置换对照塌缩至近零或负值，真实 `FCT` 在 `8/8` 有效线保持正值且配对差 `FCT - E1` 的 `97.5% CI` 在 `8/8` 线上均 `> 0`（`dld1 +0.1969`、`cd4t_rest +0.0426`）。

### 2.3 第三段：第十一轮全基因口径前瞻验证与 `Q7` 强制披露（Round 11 Descriptive Results）
- **严格包含 `PREREG_D Q7` + `P1` 全部必报要点（依 `R11_VERDICT.md` §1, §3–5, §7–8）**：
  1. **`K_conf = 1` 与 `P1` 强制声明**：在预锁 `80 GB` 预算与 `3+1` 候选池中，`hcasmc` 因公开 `.h5`/`.h5ad` 矩阵不含细胞级别 guide 分配而在评估前按 `Q3` 可行性排除；`telohaec` 经严格单 sgRNA 过滤（`180,574` 细胞）因该 GWAS 文库敲低较弱仅余 `N = 9 < 50` 个共有可重复靶点未过确认门槛；`feng_ipsc`（`17` 条非 KOLF iPSC 细胞系，`N = 93`）预指定为次要线，致使合格确认性细胞系仅余 `cd4t_stim48` 一条（**`K_conf = 1`**）。遵循预注册 `P1` 规则：**Insufficient confirmatory cell lines (`K_conf < 2`); confirmatory hypothesis testing was not performed（确认性线不足，未做确认性检验），全部第 11 轮结果均为描述性结果。**
  2. **`cd4t_stim48` 描述性结果（`N = 503`，`422` 两半有效，`|G_eval| = 6,131`，`rho_AB_orth = +0.0784`）**：在 `all_nontarget` 口径下，**`H-D1` 通过**（`FCT +0.0532 [+0.0337, +0.0731]`）、**`H-D2` 不通过**（`FCT` 下界 `+0.0337` 低于错配对照 `vcc2025_h1` 上界 `+0.0572`，虽胜过其余 5 个错配线 `<= 0.035` 及置换 `Q0.975 = -0.0129`；在 `8 源线` 敏感性 `N = 1,406` 中仍以 `+0.0403 <= +0.0443` 输给 `vcc2025_h1`）、**`H-D3` 通过**（`FCT - E1 = +0.0735 [+0.0561, +0.0917]`）、**`H-D4` 通过**（`H-D4a E1(resp_top50) = +0.0522 [+0.0257, +0.0792]`；`H-D4b E1(top50) - E1(all) = +0.0725 [+0.0477, +0.0982]`）。
  3. **`cd4t_stim48` 独立性与技术限制（`Q2/Q7`, `Addendum A2`）**：必须明确披露 **`cd4t_stim48` 与第十轮 `cd4t_rest` 共享同一实验室（Marson lab）、同一原代人 `CD4+ αβ T` 细胞谱系、以及同一 sgRNA 文库**（仅供体与 48 小时 TCR 刺激状态为新），且两供体共用 24 条测序通道使 `r_orth_cal` 偏低（偏保守）。
  4. **次要线与跨轮综合**：在次要线 `feng_ipsc`（`N = 93`）上，`FCT +0.1212 [+0.1043, +0.1371]` 描述性通过 `H-D2`、`H-D3`（`+0.0682`）、`H-D4` 及细胞系特异对照 `H-D5`（`FCT - Mis(KOLF2.1J) = +0.0650 [+0.0528, +0.0784]`）。跨轮汇总（R9 `K_conf = 3`、R10 `K_conf = 2`、R11 `K_conf = 1` 描述性）并经过三方复算一致（R11 `48` 点估计差 `<= 1.1e-16`）：`FCT` 在全基因口径下稳定胜过常数基线 `E1`、`resp_top50` 稳定奖励常数基线 `E1`（`H-D4`），但 **细胞系特异性（胜过全部错配基线）仍未确立**（两条 T 细胞线均败于 `vcc2025_h1`），且全基因去幅度 `r_orth_cal` 效应量很小（仅 `0.05–0.12`），绝大部分去幅度线特异响应仍未被解释。

---

## 3. 论文总体结构与章节骨架（Manuscript Structure）

1. **Title & AI Governance / Role Disclosure**（人类角色统一写作 `Human project lead` / 人类项目负责人）
2. **Abstract**（严格三段式：R9+R10 Confirmatory / R10-E Post-Hoc / R11 Descriptive w/ Mandatory Q7 Disclosures）
3. **Section 1: Introduction**
   - 1.1 跨细胞系零样本 CRISPRi 预测困局（守恒响应 vs 线特异斜率；`VCC 2026 Wrap-Up`、`CASTER`、`Molina & Zhang 2026`、`arXiv:2608.00152v3`）。
   - 1.2 从 `CIPHER` 与 `TMVS` 到 `FCT`（`CellularAI`/`AtlasShift` 基线倍数映射类）。
   - 1.3 四重方法学陷阱：`R18` 单标量幅度混杂、`R4/R14` 边际方差 vs 共波动、`R3/R12` 缺位错配/置换零模型、以及 `resp_top50` 结果选基因偏倚（对比 *Nicol, Shivakumar & Irizarry 2026* 共享噪声偏倚）。
   - 1.4 三轮预注册盲态设计概览（R9 `PREREG_A v4` $\to$ R10 `PREREG_C` $\to$ R10-E 事后消融 $\to$ R11 `PREREG_D v2`）。
4. **Section 2: Methods & Pre-Registration Governance**
   - 2.1 三轮共 `14` 项人类 CRISPRi 数据集筹备、单 sgRNA（`P4`）、严格 NTC-only 规则（排除 `safe-targeting` 与 `unassigned`）、事前排除（R10 `H9` P2；R11 `hcasmc` `Q3`）及 HGNC 跨轮符号桥接
   - 2.2 标准化表达、Addendum-B 可重复性、`P5` 靶点准入、Split-Half 物理轴（含 `Addendum A2` 共用测序通道、`A3` 共享供体及 `P3` 双供体）及机械 `P1` 门槛
   - 2.3 预测器数学定义：Group A/B/`R14`/`R18`、Group C/D `FCT`、常数/错配/置换对照（`E1`–`E7`）、与 `H-C3`/`H-D5`
   - 2.4 去幅度指标 `r_orth_cal`、IUT + 2-claim Bonferroni `97.5%` 置信区间、`R19` 三层分解、`S7`（`P6`）混杂回归与多方独立复算（`C3`）
   - **2.5 预注册偏差集中披露（Deviations Across Rounds）**：`R9 D1`（`span3` vs `rank1_q1`）、`R10 D-R10-1`/`D-R10-2`、R10 `13:15` 主控事后更正、`R11 Addendum A1–A4` 与 `DEVIATIONS_R11.md`（`D1` `hcasmc`、`D2` `telohaec` 弱敲低核查、`D3` `K_conf=1` 中止确认、`D4` 冻结前报告层修改）
5. **Section 3: Round 9 Results — Covariance and Marginal-Variance Transfer (`PREREG_A v4`, `K_conf = 3`)**
6. **Section 4: Round 10 Results — Prospective Audit of `FCT` Under `resp_top50` (`PREREG_C`, `K_conf = 2`)**
7. **Section 5: Round 10-E Post-Hoc Discovery — Selection Bias in `resp_top50` and Genome-Wide Divergence (`[事后探索性 / Post Hoc]`)**
   - `BASELINE_ABLATION.md` & `AUDIT_R10E_REPORT.md` 完整数值（`resp_top50` vs `all_nontarget`, `8/8` `FCT - E1 > 0`），解释 `13:15` 更正与回看 R10 `S6`。
8. **Section 6: Round 11 Results — Whole-Transcriptome Evaluation Under `PREREG_D v2` (`K_conf = 1`, Strictly Descriptive Only)**
   - 6.1 `P1` 中止原因（`hcasmc` `Q3`、`telohaec` `N=9` 与 `D2` 经验反思）及 `cd4t_stim48` 独立性/共用通道披露（`Q2/Q7`, `Addendum A2`）
   - 6.2 `cd4t_stim48`（`r9_6` 主源池及 `holdout_r11_src8` 8 源非独立敏感性）：`H-D1` 通过、`H-D2` 因 `vcc2025_h1` 不通过、`H-D3` 通过
   - 6.3 `feng_ipsc`（`N=93` 次要线）：`H-D2`/`H-D3`/`H-D5`（vs `KOLF2.1J` `+0.0650`）与 `Addendum A3` 克隆/通道天花板披露；`telohaec`（`N=9` 附表/区间不解读）
   - 6.4 `H-D4` 选择陷阱全线复现（`3/3` 描述性通过；含 `feng_ipsc` `E1(top50) > FCT`）及三方交叉核对（`max diff = 1.1e-16`）
9. **Section 7: Cross-Round Synthesis, Prior-Art Boundaries, and Prospective Infrastructure**
   - 7.1 依据 `R11_VERDICT.md §7` 的跨轮四点综合（`H-D3` 稳定复现、`H-D4` 稳定复现、对全部错配线特异性未确立尤其两 T 细胞线均败于 `vcc2025_h1`、极小效应量 `0.05–0.12`）
   - 7.2 与前作文献对比（`PRIOR_ART_A–D`：`CGC`、`Molina & Zhang`、`CASTER`、`Nicol et al. 2026`、`cell-eval2 #351`、`Booeshaghi & Pachter 2021`、`AMM-SimWMag`、`arXiv:2608.00152v3`）
   - 7.3 `A-T4` VCC 2026 前瞻推理工具链与 `cell-eval2 v0.22.0` 本地复现（`HUMAN_STEPS.md`, `CELL_EVAL2_REPRO.md`）
10. **Section 8: Limitations, Data/Code Availability, & References**
11. **Section 9（大纲末尾独立节）: `事后探索`**
    - *依 `EXPLORATORY_REVIEW_R11.md §D` 与 `round11/audit/exploratory_check/REPORT.md` 仅列 3 条经独立复核的事后结果（`SA-8` `9/12`、`N_5` `FCT` 胜过任意常数基线组合 `12/12`、未预注册配对 `FCT − Mis_s` `11/12`）与 1 条机制局限（无噪声模型在 `all` 掩码错 `8/12` 线符号）。*

---

## 4. 完整图表清单（Figure & Table Inventory）

> 注意：主文图表中**不包含** `exploratory/null_cone_theorem`（含 `N_5` / 一阶泰勒启发式推导）的图表；`SA-8` 现象列于大纲末尾 §9。

### 4.1 主图清单（Main Figures 1–6）

| 编号 | 图题与核心展示内容 | 数据来源 | 性质与状态 |
| :--- | :--- | :--- | :--- |
| **Figure 1** | **Three-Round Pre-Registered Design, Data Attrition, and `r_orth_cal` Geometry** | R9 (`6` 线) + R10 (`4` 线 + `H9`) + R11 (`3` 线 + `hcasmc`) | **Design Overview**：(A) 三轮盲态门控流程（R9 $\to$ R10 $\to$ R10-E $\to$ R11）；(B) R11 候选池漏斗（`hcasmc` Q3 排除、`telohaec` `N=9`、`cd4t_stim48` `K_conf=1` 触发 P1）；(C) `r_orth_cal` 正交消幅示意图。 |
| **Figure 2** | **Amplitude Confounding (`R18`), Three-Tier Decomposition (`R19`), and `resp_top50` Ceiling Inflation (`S6`)** | R9 (`6` 线) + R10 (`4` 线) | **Confirmatory & Diagnostic**：(A) `GlobalScale` 原始 `r_cal` vs `r_orth_cal = 0.000`；(B) `R19` 堆叠分解；(C) `S6` $\rho^{\perp}_{AB}$ 从 `resp_top50` 到 `all_nontarget` 的下降。 |
| **Figure 3** | **Round 9 Confirmatory Results: `Hybrid_ProjShared` (Group A) and `TMVS` (Group B) Not Confirmed** | R9 (`K_conf = 3` + 3 排除线) | **Confirmatory**：`r_orth_cal`/`slope` vs 错配与置换零分布、`R14` 增量、`H-A4/B4` 劣效森林图、以及 `S1–S4b` & `S2b`（秩 1 `q_1`）。 |
| **Figure 4** | **Round 10 Prospective Audit (`resp_top50`): `FCT` Passes `H-C1` and `H-C3` Yet Fails `H-C2`** | R10 (`K_conf = 2`: `dld1`, `cd4t_rest` + 2 次要) | **Confirmatory**：`FCT` vs 6 错配源均值与 `perm100`、`H-C3` 6 口径非劣效、以及 `S7` `ln(n_cells)` / `|k_p|` 控制回归与 `T1–T3` 分层。 |
| **Figure 5** | **Post-Hoc Discovery (`R10-E`): Selection Bias in `resp_top50` vs. Unbiased `all_nontarget` Divergence** | `round10/exploratory/` & `audit/exploratory/`（10 线） | **Post-Hoc Exploratory（须红字标注 *post hoc*）**：(A) `resp_top50` vs `all_nontarget` 下 `FCT`, `E1`, `E2`, `E3`, `E5`, `E6`, `E7`, `perm100` 双掩码对照；(B) 配对差 `FCT - E1 [97.5% CI]` 在两掩码下对比（`all_nontarget` `8/8 > 0`）。 |
| **Figure 6** | **Round 11 Descriptive Whole-Transcriptome Outcomes (`K_conf = 1`, P1 No Adjudication)** | R11 (`cd4t_stim48`, `feng_ipsc`, `telohaec`; `r9_6` & `src8`) | **Descriptive Only (`P1`: 未做确认性检验)**：(A) `cd4t_stim48` (`N=503` / `src8 N=1,406`) 在 `all_nontarget` 下 `FCT` vs 错配基线（突出显示仅输给 `vcc2025_h1` `+0.0572`）与 `perm100`（`-0.0129`）；(B) `feng_ipsc` `H-D5`（`FCT` vs `Mis(KOLF2.1J)` `+0.0650`）；(C) `H-D4` 验证 `E1(resp_top50)` 相对 `E1(all_nontarget)` 在 `3/3` 线上全面抬升。 |

### 4.2 主表与附表清单（Main Tables 1–7 & Supplementary Tables）

| 编号 | 表题与核心内容 | 性质与来源状态 |
| :--- | :--- | :--- |
| **Table 1** | **Unified Dataset Inventory Across Rounds 9, 10, & 11**（R9 6 线 + R10 4 线 + `H9` P2 排除 + R11 `telohaec`, `hcasmc` Q3 排除, `cd4t_stim48`, `feng_ipsc`：保留单 sgRNA 细胞、NTC 数、Split-Half 轴、`\|G_eval\|`、可重复与评估靶点数 `N`） | 数值全部完备（来自 R9 Table 1 + `R10_QC_COUNTS.md` + `R11_QC_COUNTS.md`） |
| **Table 2** | **Round 9 Main Confirmatory Matrix (`PREREG_A v4`, `K_conf = 3`)**（Group A/B、`GlobalScale`、`R14`、`H-A1..A4` / `H-B1..B4`） | 数值完备 |
| **Table 3** | **Round 10 Main Confirmatory Matrix (`PREREG_C`, `K_conf = 2`, `resp_top50`) & `H-C3` / `S7`** | 数值完备 |
| **Table 4** | **Cross-Round `R19` Decomposition & `S6` Noise-Ceiling Drop (`R9` + `R10`)** | 数值基本完备（定稿可补齐 R10 2 条次要线行） |
| **Table 5** | **`[Post-Hoc Exploratory]` `R10-E` Baseline Ablation (`E1–E7`) & Paired `FCT − E1` Across 10 Lines** | 数值完备（`BASELINE_ABLATION.md` + `AUDIT_R10E_REPORT.md`） |
| **Table 6** | **Round 11 Descriptive Results (`PREREG_D v2`, `K_conf = 1`, Strictly Descriptive Reporting Under `P1`)**（`cd4t_stim48` [`N=503` 及 `src8 N=1,406`]、`feng_ipsc` [`N=93`]、`telohaec` [`N=9`]：`FCT`, `E1`, `FCT - E1`, `max Mis97.5_hi`, `Perm100 Q0.975`, `H-D1..H-D4` 描述性状态, `H-D5 vs KOLF2.1J`, 含 `Q2/Q7`/`Addendum A2/A3` 注解） | 数值全部完备（`R11_VERDICT.md` + `REPORT_holdout_r11_r9_6_log1p.md` + `REPORT_holdout_r11_src8/REPORT_holdout_r11_r9_r10_8_log1p.md` + `AUDIT_R11_REPORT.md`） |
| **Table 7** | **Cross-Round Descriptive Synthesis Matrix (`R10` `S6` + `R11` `all_nontarget`)**: 对比 `dld1`, `cd4t_rest`, `cd4t_stim48`, `feng_ipsc`, `telohaec` 上 `FCT`、`E1`、`FCT-E1`、与 `vcc2025_h1` 错配对照差值、以及 `H-D4` 一致性（直观展示两条 T 细胞线均败于 `vcc2025_h1`） | 数值全部完备（`R11_VERDICT.md §7` + `R10_VERDICT.md`） |
| **Supp. Tables S1–S6** | `S1`: R9 `S1–S5` + `S2b`；`S2`: R10 `H-C3` 6 口径与 `S1/S5`；`S3`: R9/R10/R10-E/R11 多方交叉核对汇总（含 R11 `48` 点估计差 `<= 1.1e-16`）；`S4`: `PRIOR_ART_A–D`（Q1–Q5 完整文献对比表）；`S5`: R11 全模型详表（`r9_6` 与 `src8` 含 `TMVS`/`CIPHER`/`Hybrid`/`E2..E7`）；`S6`: `A-T4` `cell-eval2 v0.22.0` VCC 2025 复现分与格式合规规范 | 完备（注：正表/附表暂不纳入 `SA-8` 机制外推及 `N_5` 表，见 §9） |

---

## 5. 方法、第十一轮设计与预注册偏差细纲（Methods & Deviations）

### 5.1 R9–R11 数据工程、单 sgRNA 过滤与 Split-Half 设计
- **R9 (`6` 线 `1,419,744` 细胞)** 与 **R10 (`4` 线 `4,264,773` 质控细胞，`H9` P2 排除)**：沿用冻结统一预处理（library-size $10^4$, $\operatorname{log1p}$, batch/unit NTC 中心化, Addendum B `n_cells >= 25` + `AD BH < 0.05 >= 50` + 20-PC Energy `p < 0.001`）。
- **R11 候选池与锁定规程（`PREREG_D v2`、`76.338 GB <= 80.000 GB`）**：
  1. `cd4t_stim48`（`GSE314342 Stim48hr`、`CE0008678` D3 vs `CE0010866` D2、`24` 共用测序通道、`MT <= 5%`、`n_genes >= 500`）：`3,911,606` 保留细胞（`3,744,735` targeting + `166,871` NTC，`78,168 / 88,703`），`12,458` 有限基因，`|G_eval| = 6,131`，`11,730` 去重靶点 $\to$ `7,301` 可重复 $\to$ `655` 共享 $\to$ **`503`** 准入（`422` 两半有效）。**按 `Q2`/`Q7` 明确披露与 `cd4t_rest` 共享实验室、细胞类型与文库，按 `Addendum A2` 披露两供体共用 `CD4i_R2L25..R2L48` 测序通道压低 `r_orth_cal`（`rho_AB_orth = +0.0784`）。**
  2. `feng_ipsc`（Feng et al. 2026 `TargetedScreen`，Figshare `27989294 v2`，pre-QC 剔除 `kolf_2` & `kolf_3` `133,631` 细胞与 `unassigned` `470,561` 细胞）：`17` 条非 KOLF HipSci iPSC 线（`9` 供体），`557,672` 保留细胞（`550,406` targeting + `7,266` NTC），`7,595` 基因，`|G_eval| = 5,326`，`444` 去重 $\to$ `135` 可重复 $\to$ `111` 共享 $\to$ **`93`** 准入；按 `Addendum A3` 披露同供体克隆分入两半，仅作次要评估及 **`H-D5`** 检验。
  3. `telohaec`（Schnitzler et al. 2024 `IGVFDS7943IIWZ`/`GSE210681`）：按 `max_umi_ct >= 4` & `X2nd < 4` 提取单 sgRNA，严格仅用 `non-targeting`（排除 `3,273` `safe-targeting` 与 `-TSS*`/增强子），`180,574` 保留细胞（`178,359` targeting + `2,215` NTC，`10` 奇对 `10` 偶通道），`7,324` 有限基因，`|G_eval| = 4,869`，`2,198` 去重但仅 `17` 可重复 $\to$ **`9`** 准入（`< 50` 未过确认门槛）。
  4. `hcasmc`（`IGVFDS6248NDUW` `IGVFFI5999NXAF.h5` `42.885 GB`）：依 **`Q3` 评估前技术可行性闸门**排除（实测 CSR 矩阵 `955,600 x 62,710` 仅含 RNA ENSG 特征、`0` 个 guide 分配，见 `DEVIATIONS_R11.md D1`）。

### 5.2 第十一轮检验定义与多重校正理由（`PREREG_D v2 §4`）
- **主掩码与对比掩码**：`all_nontarget` 为确认主掩码；`resp_top50` 专为测试 `H-D4` 与 `SA-7`；`SA-8` (`resp_top50_split`) 为预指定非判定性描述分析（现象与后续讨论归入 §9）。
- **多重点检校正与 `97.5%` CI（`Q1`, `Q5`, `Q8`）**：R11 共设 **`2` 项独立确认论断**：Claim 1（主生物论断 `H-D1 ∩ H-D2 ∩ H-D3`，三者须同时满足，属 Intersection–Union Test，内部无需再平分 $\alpha$）与 Claim 2（方法学论断 `H-D4 = H-D4a ∩ H-D4b`，亦为 IUT）。对两项论断做 Bonferroni 校正得 $\alpha_{\mathrm{claim}} = 0.025$，对应双侧 **`97.5%` 置信区间**（`n_boot = 4,000`, `seed = 20261123`）。因 `K_conf = 1 < 2`，按 `P1` 机械中止确认裁决，全部输出作描述性报告。

### 5.3 全部预注册偏差与过程记录（Deviations Log: R9 + R10 + R11）
1. **`R9 D1`**：`PREREG_A v4 §4.2` 文本秩 1 $q_1$ vs 冻结代码 `1ca9f60f` 三维 $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$（`S2b` 敏感性证明判定不变；`A-T3` 在该点实现非独立）。
2. **`R10 D-R10-1` / `D-R10-2` & `13:15` 更正**：`H-C3` 代码命名对照；`H-C1` 额外保护条件；`R10_VERDICT.md` `13:15` 公开更正主控初期关于 `E1` 的错觉。
3. **`R11 Addendum A1–A4` 与 `DEVIATIONS_R11.md`（均在代码冻结与查看指标之前记录）**：
   - **`A1`**：沿用 `PREREG_C` 对 `H-D1` 增加 `valid_frac >= 0.95` 与 `rho_AB_orth > 0`（实测 `1.000` 与 `+0.0784 > 0`）。
   - **`A2`/`A3`**：`cd4t_stim48` 共用 24 条通道偏保守压低 `r_orth_cal`；`feng_ipsc` 同供体克隆跨半。
   - **`D1` (`hcasmc`)**：原 `Q3` 字面写作“内存或格式无法按块读取”，实际 `.h5` 可按块读出（峰值 `0.306 GB`）但含 **0** 个 guide 特征且无单细胞 guide 分配，主控查实 `87.29 GB` `.h5ad` 亦无 guide 分配，依文件结构排除。
   - **`D2` (`telohaec`)**：非流程 bug；实测靶基因敲低中位数仅 `46%`（vs `Feng 78%`、`DLD-1 86%`）、前 `5,000` 基因显著 DE 中位数为 `0`，致 `17/2,198` 可重复、`N = 9`；说明基于必需基因富集文库估算的 `20%–40%` 经验比例不适用于弱效应 GWAS 非必需基因文库。
   - **`D3` (`K_conf = 1`)**：机械执行 `P1`，未事后修改门槛或将 `feng_ipsc` 升为确认线。
   - **`D4`**：冻结前小改报告与校验脚本以支持 `K_conf < 2` 时照常输出全部线描述值及 `H-D5` 表，未改动任何指标或假设。

---

## 6. 局限性（Combined Limitations: R9 + R10 + R11 Descriptive）

1. **第十一轮因 `K_conf = 1` 未能完成前瞻确认检验，且唯一达标线独立性有限（`Q2`, `Q7`, `Addendum A2`）**：`hcasmc` 因缺少细胞级 guide 分配排除（`D1`）、`telohaec` 因文库效应微弱仅 `N = 9`（`D2`），使第十一轮按 `P1` 中止确认裁决。唯一达标的 `cd4t_stim48` 虽使用全新供体与刺激态，但与 R10 `cd4t_rest` **同属 Marson 实验室、同为原代 `CD4+ αβ T` 细胞系、并使用同一 sgRNA 文库**；且供体共用 24 条通道拉高噪声天花板、压低 `r_orth_cal`（`rho_AB_orth = +0.0784`）。
2. **`FCT` 跨线特异性的细胞系局限与 `vcc2025_h1` 错配反超之谜（`R11_VERDICT §7`）**：在全基因 `all_nontarget` 下，`dld1`（R10 `S6`）与次要线 `feng_ipsc`（`N = 93`；在 `H-D5` 上超越同为 iPSC 的 `KOLF2.1J` `+0.0650`）虽胜过全部错配基线，但**两条 T 细胞线（R10 `cd4t_rest` 差 `0.002`；R11 `cd4t_stim48` 差 `0.024`，8 源差 `0.004`）均输给了 `vcc2025_h1`（H1 胚胎干细胞）基线**。为什么 H1 基线对 T 细胞具有预测力目前尚无解释，无法确立普适的“胜过全部错配细胞系”特异性。
3. **极小的可解释方差（`Q10`, `R11_VERDICT §7`）**：即使通过 `H-D1` 和 `H-D3` 的全基因 `FCT`，在 R10/R11 达标线上去幅度相关也仅 `r_orth_cal = 0.05–0.14`（R11 `cd4t_stim48` `+0.0532`, `feng_ipsc` `+0.1212`），绝大部分去幅度后的细胞系特异扰动响应仍未得到解释。
4. **必需基因文库与非必需 GWAS 文库的可重复率不可外推（`D2`）**：基于 `Replogle`/`DLD-1` 估算的 `20%–40%` 共享靶点保留率在 `TeloHAEC_2285`（`46%` 敲低中位数）上降至 `< 1%`（`9 / 391`），未来跨谱系实验规划必须先按元数据与先导敲低深度核验文库效应强度，而非套用比例。
5. **模型同源性与 R9 3D 投影非独立性**：三轮主代码（`A-T2`）与独立审计（`A-T3`）均基于同一 AI 模型家族（`Antigravity`/Gemini），虽有主控 `Claude`（Anthropic）另写代码核对了 R9 `TMVS`、R9 `FCT`、R10 确认集及 R11 全部 `48` 点估计（差 `<= 1.1e-16`），但 R9 Group A（`Hybrid_ProjShared`）未由跨模型主控复算，且 R9 中 `A-T3` 复刻了 3D 投影。
6. **外部 VCC 2026 斜率公式未公开、测试集真值未释放及文献查新边界（`R6`, `A-T4`, `PRIOR_ART_D–E`）**：`cell-eval2 v0.22.0` 未公布 `0.02–0.13` 跨线斜率计算式；`gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` 截至 2026-11-22 为空，`AC16`/`MOLM13`/`U-87MG` 受条款限制暂无法开放评测；所有先行文献结论受公开索引范围限制（`PRIOR_ART_A–E`）。
7. **事后跨线迁移表现的子集依赖性（`EXPLORATORY_REVIEW_R11.md §D`, `A-T3` 按线分块复核）**：`FCT` 优于直接 `log1p` 迁移，只在 `expr05` 子集上稳健（`10/10` 线、按线分块置换并 `Bonferroni×3` 校正 `p = 0.0059`；而在 `all` 子集上按线分块非参数检验 `p = 1.0000`）；`Q5a` 跨线守恒论断与无噪声一阶泰勒机制推导均不写入正文。

---

## 7. 待办与合稿推进清单（Manuscript Assembly Checklist）

- [x] **Round 9–11 全部锁定评测、三方复核与事后独立核查**：`R10_VERDICT.md`（含 `13:15` 更正）、`BASELINE_ABLATION.md`、`AUDIT_R10E_REPORT.md`、`PRIOR_ART_D.md`（Q1–Q5）、`R11_VERDICT.md`、`EXPLORATORY_REVIEW_R11.md`（含 `§D`）、`PRIOR_ART_E.md` 与 `round11/audit/exploratory_check/REPORT.md` 均已落盘。
- [x] **大纲 §9 事后探索定稿收口**：已按 `EXPLORATORY_REVIEW_R11.md §D` 缩减为 3 条可报告事后发现与 1 条局限（数字全部采用 `A-T3` 锁定口径复算）；若拟将配对特异性或 `SA-8` 升级为确认性结论，需先由 `A-T1` 筛查全新未用数据并另行开启预注册前瞻检验。
- [ ] **全稿合规终检**：人类角色统一使用 `Human project lead`（绝不使用其他首席研究员英文缩写或头衔）；全部数值核对回源 CSV/JSON 与裁决文档；在主控确认范围后再输出完整草稿或 PDF。

---

## 8. 核心参考文献目录（已验证引文）

- **Round 9 文献**：Adduri et al. (2026) *Cell* (`10.1016/j.cell.2026.08.004`); Arc Institute VCC 2026 Wrap-Up (`2026-11-19`) & `cell-eval2` (`v0.22.0`); Dhaka et al. (2026) `HyperMap` (`10.64898/2026.04.23.720505`); Huang et al. (2026) `CGC` (`10.64898/2026.09.30.755756`); Kawsar et al. (2026) `CASTER` (`10.64898/2026.11.15.766815`); Kuznets-Speck et al. (2025/2026) `CIPHER v2` (`10.1101/2025.06.27.661814`); Luan et al. (2026) (`10.64898/2026.10.20.760514`); Miller et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03307-w`); Molina & Zhang (2026) (`10.64898/2026.07.24.740459`); Nadig et al. (2025) *Nat. Genet.* (`10.1038/s41588-025-02169-3`, `GSE264667`); Nourreddine et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03199-w`); Replogle et al. (2022) *Cell* (`10.1016/j.cell.2022.05.013`); Shoeibi & Yousefi (2026) *arXiv*:`2608.00152v3`; Sun et al. (2026) `AMM-SimWMag` *Genes* (`10.3390/genes17070816`); Viñas Torné et al. (2025) `Systema` *Nat. Biotechnol.* (`10.1038/s41587-025-02777-8`).
- **Round 10 & R11 数据集文献**：Yeung et al. (2026) `DLD-1` (`10.64898/2026.07.10.737863`, `GSE337988`); Zhu et al. (2026) `CD4+ T` *Cell* (`10.1016/j.cell.2026.08.002`, `PMID:42664972`, `GSE314342`); Tian et al. (2021) *Nat. Neurosci.* (`10.1038/s41593-021-00862-0`, `GSE152988`/Zenodo `13350497`); Sunshine et al. (2023) *Nat. Commun.* (`10.1038/s41467-023-41788-4`, `GSE208240`/Zenodo `13350497`); Gogate et al. (2026) `H9` (`IGVFDS9825FOBP`, `PMC13629240`, `PMID:42827575`); Schnitzler et al. (2024) `TeloHAEC_2285` *Nature* (`IGVFDS7943IIWZ`, `GSE210681`, `PMC10921916`, `PMID:38326615`); `HCASMC_971` (`IGVFDS6248NDUW`); Feng et al. (2026) `TargetedScreen` (`PMC12903452`, `PMID:41330380`, Figshare `27989294 v2`).
- **`PRIOR_ART_D.md` (Q1–Q5) 关键方法与基准对比引文**：Nicol, Shivakumar & Irizarry (2026) (`10.64898/2026.05.07.723486`); Illumina `CellularAI`; `AtlasShift`; Xaira `X-Cell`; Adduri et al. (2025) `Arc State` (`10.1101/2025.06.26.661135`); Yan et al. (2026) `PRiMeFlow` (`arXiv:2604.13986v2`); Booeshaghi & Pachter (2021) *Bioinformatics* (`10.1093/bioinformatics/btab085`); Kviman et al. (2025/2026) (`10.1101/2025.03.12.642799`); Liu et al. (2025/2026) (`arXiv:2511.16954`); Panagopoulos et al. (2026) (`arXiv:2610.26267`); Ahlmann-Eltze & Huber (2023) (`10.1038/s41592-023-01814-1`); Ahlmann-Eltze et al. (2025) (`10.1038/s41592-025-02772-6`); Vollenweider & Bühlmann (2026) `SBB` (`10.64898/2026.04.20.719650`); Schäfer et al. (2026) `scPertEval` (`10.64898/2026.07.23.740433`); Wei et al. (2026) `scPerturBench` (`10.1038/s41592-025-02980-0`); Li et al. (2026) *Sci. Adv.* (`10.1126/sciadv.aed3414`); Dibaeinia et al. (2026) (`10.64898/2026.02.04.703804`); Ghajari et al. (2026) (`10.64898/2026.11.04.764088`).

---

## 9. 事后探索

> **口径与隔离说明（依据 `round11/status/EXPLORATORY_REVIEW_R11.md` 第 D 节）**：本节以下全部数字均取自 `A-T3` 按锁定 `P5` 口径独立复算的 `round11/audit/exploratory_check/REPORT.md`（原 `RESULTS_13LINE_NULL_CONE_AUDIT.md` 因未执行 `P5` 靶基因自检过滤、在 `11/12` 条线上用大了 `N` 而作废）。机制推导、解药式夸大表述及 `Q5a` 表达形式守恒论断均**不写入论文正文**。

1. **SA-8 现象（`9/12`）**：在预指定 `SA-8`（`resp_top50_split`，一半细胞选基因、另一半评测并双向平均）及 12 线事后复现中，常数基线 `E1` 与置换得分仍保持正值；`E1(split) − E1(all)` 点估计在 `12/12` 线为正，且联合交叉掩码 `97.5% CI > 0` 在 **`9/12`** 条线成立（如 `cd4t_stim48` `+0.0730 [+0.0384, +0.1098]`、`N = 503`/`n_c = 422`，`feng_ipsc` `+0.1938 [+0.1480, +0.2405]`、`N = 93`，`telohaec` `+0.2663 [+0.0196, +0.5715]`、`N = 9`；而 `cd4t_rest` `+0.0383 [-0.0061, +0.0864]`、`N = 302`/`n_c = 284` 与 `nourreddine_kolf21j` `+0.1211 [-0.0158, +0.2564]`、`calu3` `+0.1048 [-0.0658, +0.2672]` 不显著）。`PRIOR_ART_E.md` 与主控抽查未发现既有文献报告「分半选基因后残差指标仍奖励不含目标信息的变换」，目前作为事后探索发现与待解释现象。
2. **`FCT` 胜过任意常数基线组合（`N_5` 下 `12/12`，并写明按构造为零项）**：在 `N_5`（半偏相关 `r_cone_sp`）下，`FCT` 在 `12/12` 条线保持正非锥信号（`+0.0553` 至 `+0.1899`，非锥比例 $\eta \in [0.741, 0.934]$）并 **`12/12` 胜过单次随机置换 `Perm100 Q97.5`**（`Perm0..99` 非锥能量中位数 $\eta \in [0.6308, 0.7133]$、中位得分 `-0.0010` 至 `+0.0007`）。须写明：**`GlobalScale`（$\eta \in [8.92\times 10^{-16}, 7.30\times 10^{-15}]$）、`b_S`（`Source_Grand_Mean`，$\eta \in [8.36\times 10^{-16}, 5.54\times 10^{-15}]$）以及 `E1`/`E2`/`E6`/`E7`（$|r_{\mathrm{cone\_sp}}| \le 0.0001$，$R^2_{\mathcal{N}_5} \ge 0.992784$）在 `N_5` 下为零或近零是按数学构造如此**，不能用作该指标的经验验证证据；`N_5` 仅作为事后分析工具。
3. **配对比较 `FCT − Mis_s` 的 `CI > 0`（`11/12`，未预注册；不配对时是 `7/12`，两条 T 细胞线都没过）**：在 `N_5` 下按跨全部 6 个错配基线的**配对联合自助法**（`FCT − Mis_s` `97.5% CI > 0`）通过 **`11/12`** 条线（仅 `nourreddine_kolf21j − Mis(k562)` 未过 `+0.0525 [-0.0004, +0.0955]`；`cd4t_rest` 为 `+0.0345 [+0.0261, +0.0432]`，`cd4t_stim48` 为 `+0.0141 [+0.0064, +0.0218]`），但**该配对检验未预注册**，仅属事后结果；而在**不配对**的严格口径（`FCT ci975_lo > max_s Mis_s ci975_hi`）下仅 **`7/12`** 通过，**两条 T 细胞线都没过**（`cd4t_rest` `+0.0351 <= +0.0360` 差 `-0.0009`；`cd4t_stim48` `+0.0559 <= +0.0747` 差 `-0.0188`）。

- **局限**：机制模型在 `all` 掩码上预测错了 `8/12` 条线的符号（用真实 $\mu$ 的无噪声模拟在 `resp_top50_split` 上 `12/12` 预测为正，但在 `all` 掩码上同样 `12/12` 预测为正，错过真实观测中 `8/12` 条线的负号，如 `cd4t_stim48` 模拟 `+0.6084` vs 实测 `-0.0203`、`replogle_rpe1` 模拟 `+0.1738` vs 实测 `-0.0934` 且幅度偏大 `3–15` 倍），所以机制仍未解释。
