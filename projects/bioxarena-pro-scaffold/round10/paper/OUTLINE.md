# 第九轮 + 第十轮（含 R10-E 事后修正）+ 第十一轮合稿论文骨架（`round10/paper/OUTLINE.md`）

- **定位与核心叙事四部曲（Round 9 $\to$ Round 10 $\to$ R10-E Post-Hoc $\to$ Round 11 Confirmatory）**：
  1. **第九轮（`PREREG_A v4`，6 线 LOLO `K_conf = 3`）**：预注册的协方差共波动转移（`Hybrid_ProjShared`）与边际方差缩放（`TMVS`）方法**均未确认**（`0/4` 标准）；共波动在边际方差之外未提供可检测增量（`R14`）。
  2. **第十轮（`PREREG_C`，4 条全新单 sgRNA 人类线 `K_conf = 2`）**：按 VCC 2026 冠军同类的基线倍数缩放模型 `Fold-Change Transfer`（`FCT`）在预注册主掩码 **`resp_top50`** 下通过去幅度正信号（`H-C1 2/2`）与并行非劣效（`H-C3 12/12`），但**未通过细胞系特异性检验（`H-C2 0/2 -> 未确认`）**。
  3. **第十轮事后消融与 `13:15` 主控更正（`R10-E`，严格标注事后探索性 *post hoc*）**：`A-T2` 与 `A-T3` 用正式代码在 **`resp_top50`** 与 **`all_nontarget` (`S6`)** 双掩码下完成 `E1–E7` 消融并交叉复核（`510` 项最大差 `1.11e-16`），推翻此前“仅由变换本身解释”的初步错觉——**在 `resp_top50` 下按目标观测 $|\Delta^{\mathrm{obs}}|$ 选基因制造了严重的选择偏倚**，使无目标基线信息的常数模型 `E1` 也能复现 `86.4%–115.5%` 的 `FCT` 分数并抬高置换阈值（正是 `H-C2` 失败之根源）；一旦切换到无偏的全非靶基因口径 **`all_nontarget`**，`E1/E2/E6/E7` 与置换对照全部塌缩至近零或负值，而真实 `FCT` 在全部 `8/8` 样本量充足（`N >= 18`）的细胞系上保持正值且严格胜过 `E1`（`FCT - E1` `97.5% CI > 0` 占 `8/8`）。
  4. **第十一轮（`PREREG_D`，`[STRICT PLACEHOLDER — 仅预留骨架，不写结论]`）**：以 **`all_nontarget` 为预注册主掩码**，在完全未动用过的新人类细胞系上独立前瞻检验三条共同主假设（`H-D1..D3`，Bonferroni `98.33%` CI）及方法学假设 `H-D4`（证明 `resp_top50` 掩码系统性奖励无基线信息变换）。
- **锁定与版本锚点**：
  - **Round 9**：`PREREG_A v4`（commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`），`CODE_FREEZE`（`1ca9f60f49adf1ac7cce8266e40533273411db30`）。
  - **Round 10 & R10-E**：`PREREG_C`（commit `dc4b4a5`，sha256 `0fb48aeb08aa6befe9583030f33ec668155e9e937616d291ea0c211e0d1bf4c6`），`QC_LOCK`（`86872ed741c3aa681a7a9091e60d03bd9922038e`），`CODE_FREEZE_R10`（`f8d2ebfc238d83a0c3a6c86696950ed292bf90e9`，22/22 冻结文件校验 `OK`），`R10_VERDICT.md`（含 `2026-11-22 13:15` 更正），`round10/exploratory/BASELINE_ABLATION.md`，`round10/audit/exploratory/AUDIT_R10E_REPORT.md`，`PRIOR_ART_D.md`（Q1–Q3 已完成；Q4–Q5 待补充）。
  - **Round 11**：依据 [`R11_第十一轮行动指令.md`](../../controller/directives/R11_第十一轮行动指令.md) 与 [`R11_DATA_SCREEN.md`](../shared/R11_DATA_SCREEN.md)（`R11_SCREEN_DONE`，仅元数据；零表达数据下载）。`PREREG_D` 尚待起草与锁定。
- **AI 治理与作者披露声明（必保留）**：人类负责人 Zijian Zeng 无正规生物学背景；主控 `Claude`（Anthropic）负责预注册审核、结果核验及独立复算（并在 R10-E 主动记录及更正了初期事后直觉判断）；执行智能体 `Antigravity`（Google Gemini 模型）负责数据整合与元数据筛查（`A-T1`）、主分析与基线消融（`A-T2`）、独立复核与 `C3` 交叉核对（`A-T3`）、VCC 2026 前瞻管线/合规核查与合稿骨架（`A-T4`）及 API 文献/源码审计（`A-T5`）。

---

## 1. 标题候选（Title Candidates）

> **措辞原则**（沿用 Round 9 `PAPER_REVIEW_R1.md` **B3**）：禁止使用 *"Rejects"*、*"Refutes"*、*"Proves"* 等过强表述；区分已锁定确认判定与事后探索机制，并为 Round 11 结果预留分支后缀。

1. **四部曲机制+预注册总括型（主推，适应 R11 任意判定分支）：**  
   *Outcome-Selected Gene Masks Obscure Baseline-Driven Specificity in Zero-Shot CRISPRi Perturbation Prediction: Two Confirmatory Audits, Post-Hoc Ablation, and a Prospective Genome-Wide Test*
2. **偏方法学与评测陷阱警示型：**  
   *Separating Amplitude, Selection Bias, and Baseline Specificity in Zero-Shot CRISPRi Generalization Across Human Perturb-seq Atlases*
3. **贴合 VCC 2026 基准争议型：**  
   *Why Zero-Shot Line-Specific CRISPRi Benchmarks Fail or Misattribute Baseline Transfer: Multi-Round Pre-Registered Audits From Co-Fluctuation to Fold-Change Scaling*
4. **预留 Round 11 判定收束后缀（仅在 `PREREG_D` 完成后按实选择）：**  
   - *... Confirms Genome-Wide Fold-Change Baseline Specificity on Unseen Human Cell Lines*（若 `H-D1..H-D3` 全票通过且 `H-D4` 通过）  
   - *... Fails Unbiased Pre-Registered Genome-Wide Specificity Despite `H-D4` Mask-Bias Confirmation*（若 `H-D1..H-D3` 未通过）

---

## 2. 摘要骨架（Abstract Skeleton）

> **排版与结构约束**：成稿控制总长紧凑，**严格分成三大块**——（1）R9 + R10 已锁定确认结果；（2）R10-E 事后探索与 `13:15` 更正机制（明确标注 *post hoc*）；（3）R11 前瞻检验槽位（明确声明结论空白，待 `PREREG_D` 留出验证完成后填入）。

### 2.1 第一段：第九轮与第十轮已锁定确认性结果（Confirmatory Results: Round 9 & Round 10）
1. **研究背景与幅度混淆（`R18`）**：利用未扰动对照零样本预测细胞系特异 CRISPRi 响应在 Arc VCC 2026 中表现停滞（跨系斜率 `0.02–0.13` 对生物重复 `1.02`）。我们证明单凭每细胞系单一标量乘子（`Transfer_Mean_GlobalScale`）即可制造高原始残差相关（R9 确认集 `-0.126` 至 `+0.250`，全 6 线 `-0.316` 至 `+0.352`；R10 `dld1 +0.240`，`cd4t_rest -0.262`），而在单扰动掩码内正交去除共享跨线响应 $\hat{\Delta}^{\mathrm{TMean}}_p$ 后（`r_orth_cal`）精确坍缩至 `0.000`。
2. **第九轮（`PREREG_A v4`，6 线 `1,419,744` 细胞、`6,148` 基因，LOLO `K_conf = 3`）**：以 `Replogle K562/RPE1` 为开发集、排除 `KOLF2.1J`（`N = 18`），在盲态留出线 `nadig_hepg2`（`N = 799`）、`nadig_jurkat`（`N = 934`）与 `vcc2025_h1`（`N = 80`）上按 Bonferroni `97.5%` 置信区间与 `3/3` 一致性检验，**Group A（`Hybrid_ProjShared`，协方差共波动）与 Group B（`TMVS`，边际方差缩放）均未确认（`0/4` 标准）**：均未达 `slope > 0.13`（`0/3`）、均在 `H1` 上负向失利（`H-A2`/`H-B2: 2/3`）、均未击败错配对照与 100 次基因置换（`0/3`）且对 `Transfer_Mean` 全部劣效（`3/3`）。共波动结构在边际方差外无额外可检测增量（`R14`，`2/3` 未一致性通过）。
3. **第十轮（`PREREG_C`，4 条全新单 sgRNA 人类数据集 `4,264,773` 质控细胞，**`K_conf = 2 confirmatory cell lines`**）**：针对与 VCC 2026 #1（`CellularAI`）及 #82（`AtlasShift`）同属基线倍数映射类的 `Fold-Change Transfer`（`FCT`），预锁 6 步盲态流程；`H9` 神经元因高 MOI 导致单 sgRNA 靶点稀疏且超出 `60 GB` 预算在下载前排除（`P2`），`calu3`（`N = 8`）与 `wtc11_i3neuron`（`N = 2`）作次要报告，`dld1`（`N = 86`）与双供体 `cd4t_rest`（`N = 302`，`284` 跨供体两半有效）进入确认裁决。在预注册主掩码 `resp_top50` 下：**`H-C1` 确认（`2/2`）**（`dld1 +0.1401 [+0.0870, +0.1917]`、`cd4t_rest +0.0560 [+0.0171, +0.0938]`，且经 `S7` 控制 `ln n_cells` 与敲低深度 `|k_p|` 后截距仍为正 `+0.1400`、`+0.0607`）；然而 **`H-C2` 未确认（`0/2`）**（两条线 `CI97.5_lo` 均未同时超出 `6` 错配均值上界 `+0.1092`/`+0.0761` 与置换 `p97.5 +0.1179`/`+0.0554`），故**主科学论断「`FCT` 捕捉到目标线特异信号」在 `PREREG_C` 下未确认**。并行 **`H-C3` 非劣效通过（`12/12`）**（`cd4t_rest` `cent_de50_r` `+0.0305 [+0.0177, +0.0435]` 严格更优），三方复算精确吻合（`max <= 6.66e-16`）。

### 2.2 第二段：事后探索性发现与掩码选择偏倚机制（Post-Hoc Exploratory `R10-E` — Explicitly Demarcated）
- **[事后探索性（post hoc），不改变 R10 确认判定]**：针对“无目标基线信息模型是否即可产生同样增益”这一疑问，双团队独立消融（`E1–E7` 对比 `FCT`，`510` 项核对最大差 `1.11e-16`）揭示了**关键的掩码依赖分叉**：
  1. 在 **`resp_top50`**（按真实 $|\Delta^{\mathrm{obs}}_t|$ 选前 50 基因）下，不含目标线基因差异信息的常数中位数基线 `E1` 确能复现 `FCT` 去幅度得分的 **`86.4%–115.5%`**（`dld1` `+0.1567` 对 `+0.1401`，配对差跨越零；`cd4t_rest` `+0.0495` 对 `+0.0560`），同时 100 次置换中位数虚高（`+0.0892`、`+0.0347`）——表明最初由表观 `E1 ≈ FCT` 引出的“纯变换无基线增益”解读并不成立，实为按目标结果选基因导致的选择偏倚。
  2. 在无选基因偏倚的 **`all_nontarget` (`S6`)** 全非靶基因掩码下，所有无目标基线变体（`E1`、`E2`、`E6` $\ln\overline{\mathrm{FC}}$、`E7` $\overline{\mathrm{FC}}-1$）及 `perm100` 均塌缩至零附近或负值（`dld1` `E1 = -0.0614`，`perm100 p97.5 = -0.0352`；`cd4t_rest` `E1 = +0.0130`，`perm100 p97.5 = +0.0131`），而真实 `FCT` 在全部 **`8/8` 样本量充足细胞系（`N >= 18`）中保持稳健正值**（`dld1 +0.1355`、`cd4t_rest +0.0556`、R9 六线 `+0.1081..+0.1848`），且配对差 `FCT - E1` 的 `97.5% CI` 在 **`8/8` 线上全部严格大于 `0`**（`dld1 +0.1969 [+0.1673, +0.2253]`、`cd4t_rest +0.0426 [+0.0248, +0.0605]`）。回查 R10 预注册次要分析 `S6`，`dld1` 在全基因口径下已满足特异性（`+0.1082 > +0.0793` / `-0.0352`），而 `cd4t_rest` 仅差 `0.002`（`+0.0376 <= +0.0395`）。
  3. **文献核查定位（`PRIOR_ART_D.md`）**：虽然 `CGC`、`Molina & Zhang`、`CASTER` 与 `VCC 2026 Wrap-Up` 已证明许多深度跨系模型的总体增益可被无情境均值基线解释（Q3），但此前**未见研究**（Q1/Q2 **没有找到**）系统测试 `FCT` 的“常数基线替换消融”，亦无人证明**仅凭 `log1p(c·FC)−log1p(c)` + 几何平均即能在结果筛选（`resp_top50`）掩码下凭空制造 `86%–116%` 虚假 `r_orth_cal` 得分**。

### 2.3 第三段：第十一轮前瞻验证（`[ROUND 11 PROSPECTIVE PLACEHOLDER — ZERO CONCLUSIONS]`）
- *To test this post-hoc hypothesis prospectively without circularity, Round 11 (`PREREG_D`) locks **`all_nontarget`** as the pre-registered primary mask across brand-new, unseen human CRISPRi cell lines selected strictly from metadata screening (`R11_DATA_SCREEN.md`; zero prior evaluation) before downloading expression matrices. Using Bonferroni-adjusted two-sided **`98.33%` CIs** (`[0.833%, 99.167%]`) over three co-primary hypotheses (`H-D1`: `FCT r_orth_cal > 0`; `H-D2`: `FCT > max(6 mismatched control means, perm100)` in `all_nontarget`; `H-D3`: paired `FCT - E1 > 0`) plus methodological secondary hypothesis `H-D4` (`E1(resp_top50) > 0` and `E1(resp_top50) - E1(all_nontarget) > 0`), `[INSERT K_conf, LINE NAMES, N_eval_targets, PASS/FAIL COUNTS FOR H-D1..H-D4, AND FINAL VERDICT STRICTLY AFTER PREREG_D HOLDOUT & C3 AUDIT COMPLETE]`.*

---

## 3. 论文总体结构与章节骨架（Manuscript Structure）

1. **Title & AI Governance / Role Disclosure**
2. **Abstract**（严格三段式：R9+R10 Confirmatory / R10-E Post-Hoc Discovery / R11 Prospective Placeholder）
3. **Section 1: Introduction**
   - 1.1 跨细胞系零样本 CRISPRi 预测困局：守恒主成分高相关 vs. 细胞系特异斜率 `0.02–0.13`（`VCC 2026 Wrap-Up`、`CASTER`、`Molina & Zhang`、`arXiv:2608.00152v3`）。
   - 1.2 从线性响应理论（`CIPHER`）与边际方差（`TMVS`）到基于基线表达的倍数缩放（`FCT` / `CellularAI` / `AtlasShift`）。
   - 1.3 零样本评测与基线验证中的四重方法学陷阱：`R18` 全局标量幅度、`R4/R14` 边际方差 vs. 共波动、`R3/R12` 缺位错配/置换零模型、以及 **`resp_top50` 结果选基因偏倚与细胞数纠缠（`S6`, `S7`）**。
   - 1.4 三轮递进盲态研究架构（R9 `PREREG_A v4` $\to$ R10 `PREREG_C` $\to$ R10-E 事后分解 $\to$ R11 `PREREG_D`）。
4. **Section 2: Methods & Pre-Registration Governance**
   - 2.1 R9（6 线）+ R10（4 线 + `H9` P2 下载前排除）数据整编、单 sgRNA（`P4`）、双供体 `P3` 与跨轮 HGNC 符号桥接（附 R11 候选数据准入规范与零表达筛查机制）
   - 2.2 预处理、Addendum-B 可重复性、Split-Half 物理单元定义与运行时盲态闸门（`PREREG_COMMIT`, `CODE_FREEZE`）
   - 2.3 全部预测模型与负对照数学定义：Group A/B（`Hybrid_ProjShared`, `TMVS`, `GlobalScale`, `R14`, `MismatchedCov/Var`, `SrcCorrTgtVar`）、Group C（`FCT`, `MismatchedMean`, `Perm`）、**`R10-E` 消融家族（`E1, E2 p10/p90, E3, E5, E6, E7`）**
   - 2.4 评估终点与校正统计：`r_orth_cal`、联合重抽样 Bonferroni 区间（R9/R10 `97.5%`; R11 三共同主假设 `98.33%`）、`H-C3` 非劣效、`R19` 三层分解与 `S7`（`P6`）混杂回归
   - 2.5 独立多组复算协议（`A-T2`, `A-T3`, 主控 `Claude`）
   - **2.6 预注册偏差与过程修正记录（Deviations from Pre-Registration）**：`R9 D1`（`span3` vs `rank1_q1`）、`R10 D-R10-1`、`R10 D-R10-2`、$\mu_t$ `0` vs `1e-8`、`H9` P2 原因校订及 R10 `13:15` 主控事后解读更正说明
5. **Section 3: Round 9 Results — Covariance and Marginal-Variance Transfer (`PREREG_A v4`, `K_conf = 3`)**
   - `R18` 标量混淆证实、Group A 与 Group B 确认性未通过（`0/4`）、`R14` 增量阴性、`R19`/`S6` 诊断、全套敏感性（含 `S2b`）与二次探索 `FCT` 线索
6. **Section 4: Round 10 Results — Prospective Verification of `FCT` Under `resp_top50` (`PREREG_C`, `K_conf = 2`)**
   - `H9` 预排除与机械准入、HGNC 挽回（`83->86`, `301->302`）、`H-C1 2/2 PASS` vs **`H-C2 0/2 FAIL`（主论断未确认）**、`H-C3 12/12 PASS`（及 `S1` 降落）、`S7` 细胞数对照、`SA-3` R9 锚点表现
7. **Section 5: Round 10-E Post-Hoc Discovery — Unmasking Selection Bias in `resp_top50` and Genome-Wide `FCT` Specificity (`[事后探索性，已正式复核]`)**
   - *全节显式标注 `Post-Hoc Exploratory`（不改变 R10 `PREREG_C` 未确认裁决）：*
   - 5.1 为什么 `resp_top50` 会诱发错觉：按 $|\Delta^{\mathrm{obs}}_t|$ 选基因令 `E1` 复现 `86.4%–115.5%` 的 `FCT` 并抬高 `perm100`，使配对差 `FCT - E1` 在 5 条大样本/R9 线上跨越零（解释 `13:15` 更正的来龙去脉）。
   - 5.2 切换至无选基因偏倚的 `all_nontarget` (`S6`) 掩码：`E1/E2/E6/E7` 与 `perm100` 集体塌缩，真实 `FCT` 在 **`8/8` 有效线保持正值且 `FCT - E1` 97.5% CI 全胜（`8/8`）**；重新解读 R10 冻结 `S6` 中 `dld1`（已过特异性）与 `cd4t_rest`（差 `0.002`）的表现，及低基线权重的非线性压低机理假说。
   - 5.3 双团队 `R10-E` 交叉复核（`510` 项 `<= 1.11e-16`）与 `PRIOR_ART_D.md` 边界定位（Q1/Q2 填补，Q3 引用前作边界）。
8. **Section 6: Round 11 Prospective Genome-Wide Validation (`PREREG_D`) (`[STRICT PLACEHOLDER]`)**
   - *仅写设计占位框架与候选屏筛元数据（[R11_DATA_SCREEN.md](../shared/R11_DATA_SCREEN.md)），**不写入任何 Round 11 结果或结论**：*
   - 6.1 候选新细胞系元数据与功效预估（如 `#1 TeloHAEC_2285` `1.47 GB` 动脉内皮静息态预计 `75–150` R9 靶点、`#2 HCASMC_971` `42.89 GB` `.hdf5` 冠脉平滑肌预计 `40–80` 靶点、可选 `Feng 2026` 43 株非 KOLF iPSC；明确放弃高 MOI `WTC-11 EC`、`TeloHAEC_GW` 与超预算 `WTC-11 CM`；结合 R10 `cd4t_rest` `0.002` 窄边界说明样本量要求）。
   - 6.2 `PREREG_D` 确认性判定表占位（`H-D1`, `H-D2`, `H-D3` 在 `all_nontarget` 下，`98.33%` Bonferroni CI；支持 `8 源线` 敏感性与 `S7`、`H-C3` 实用性）。
   - 6.3 `H-D4` 方法学验证占位（前瞻实证 `E1(resp_top50) > 0` 与 `E1(resp_top50) - E1(all_nontarget) > 0`）。
9. **Section 7: Discussion, Synthesis, and Prospective Infrastructure**
   - 7.1 四阶段实证对 VCC 2026 与基准设计的重塑意义（为何 `all_nontarget` + 正交消幅 + 常数基线对照 + `n_cells` 回归缺一不可）
   - 7.2 与既有文献的关系矩阵（`PRIOR_ART_A–D`：`CGC`, `Molina & Zhang`, `CASTER`, `CIPHER v2`, `Systema`, `Miller 2026`, `AMM-SimWMag`, `arXiv:2608.00152v3`, `cell-eval2` Jensen 亏空, `Booeshaghi 2021`, `arXiv:2511.16954` + `[PENDING: PRIOR_ART_D Q4 & Q5]`）
   - 7.3 `A-T4` VCC 2026 前瞻推理与 `cell-eval2 v0.22.0` 本地复现基础设施（[HUMAN_STEPS.md](../../round9/p2/vcc_linespecific/src/prospective/HUMAN_STEPS.md)、[CELL_EVAL2_REPRO.md](../../round9/p2/vcc_linespecific/src/prospective/CELL_EVAL2_REPRO.md)，含合规边界与双阶段 `non-targeting` 拼接桥）
10. **Section 8: Limitations**（R9 + R10 + R10-E + R11 Placeholder + 文献边界）
11. **Section 9: Data, Code, AI-Governance Disclosure, and References**

---

## 4. 完整图表清单（Figure & Table Inventory）

### 4.1 主图清单（Main Figures 1–7）

| 编号 | 图题与核心展示内容 | 数据来源 | 性质标注与状态 |
| :--- | :--- | :--- | :--- |
| **Figure 1** | **Multi-Round Pre-Registered Blinding Architecture & Geometric Orthogonalization** | R9 (`6` 线) + R10 (`4` 线 + `H9` P2) + R11 (`PREREG_D` 框架) | **已定稿设计**：(A) 四阶段锁链（R9 LOLO $\to$ R10 `PREREG_C` $\to$ R10-E 事后诊断 $\to$ R11 `PREREG_D` 未见过细胞系）；(B) 通道/双供体 `P3` Split-Half；(C) `r_orth_cal` 几何去幅图。 |
| **Figure 2** | **Whole-Line Amplitude Confounding (`R18`), Three-Tier Decomposition (`R19`), and `resp_top50` Ceiling Inflation (`S6`)** | R9 (`6` 线) + R10 (`4` 线) | **已定稿（Confirmation & Diagnostic）**：(A) `GlobalScale` 原始 `r_cal` vs `r_orth_cal = 0.000`；(B) `R19` 堆叠柱状图（`all_nontarget` 下 `L1 = 1.6%–51.3%`、确认线 `11.3%–34.8%`）；(C) `S6` $\rho^{\perp}_{AB}$ 在 `resp_top50` 到 `all_nontarget` 的断崖下降。 |
| **Figure 3** | **Round 9 Confirmatory Results: `Hybrid_ProjShared` (Group A) and `TMVS` (Group B) Not Confirmed** | R9 (`K_conf = 3` + 3 排除线) | **已定稿（Confirmatory）**：(A) 对比 `5` 错配对照组与 `perm100` 分布（`H-A1..A3`, `H-B1..B3`）；(B) `R14` 增量与 `SrcCorrTgtVar`；(C) `H-A4/B4` 劣效图；(D) `S1–S4b` & `S2b`（秩 1 `q_1`）敏感性。 |
| **Figure 4** | **Round 10 Prospective Audit (`resp_top50`): `FCT` Passes `H-C1` and `H-C3` Yet Fails `H-C2` Specificity** | R10 (`K_conf = 2`: `dld1`, `cd4t_rest` + `calu3`, `i3Neuron`) | **已定稿（Confirmatory）**：(A) `FCT` vs `6` 错配源均值（`+0.1092`, `+0.0761`）与 `perm100`（`+0.1179`, `+0.0554`）；(B) `H-C3` `12/12` 非劣与 `S1` `4/12` 劣效对比；(C) `S7` `r_orth_p` 对 `ln(n_cells)` / `|k_p|` 控制回归与 `T1/T2/T3` 三分位。 |
| **Figure 5** | **Post-Hoc Discovery (`R10-E`): Selection Bias in `resp_top50` Masks Genome-Wide `FCT` Baseline Specificity** | `round10/exploratory/` & `audit/exploratory/`（R9 6 线 + R10 4 线） | **已核实数值，图面必须显式加注「事后探索性（Post Hoc）」红框**：(A) `resp_top50` vs `all_nontarget` 双面板哑铃/分组柱图对比 `FCT`、`E1(Med)`、`E2(P10/P90)`、`E3`、`E5`、`E6`、`E7` 与 `perm100`（直观展示无基线模型仅在 `resp_top50` 下高企而在 `all_nontarget` 下集体崩塌）；(B) 配对差 `FCT - E1 [97.5% CI]` 在两个掩码下的森林图对比（`resp_top50` 在 5 大线跨零 vs `all_nontarget` `8/8` 有效线全部严格 `> 0`）；(C) R10 冻结 `S6` (`all_nontarget`) 错配对照与置换回看（`dld1` 通过 vs `cd4t_rest` 差 `0.002`）。 |
| **Figure 6** | **`[ROUND 11 PROSPECTIVE PLACEHOLDER]` Genome-Wide Confirmatory Validation (`PREREG_D`, `all_nontarget`)** | **[待 R11 新线完成后绘制]** | **预留空图槽（禁止预填结论）**：(A) 全新人体细胞系在 `all_nontarget` 主掩码下 `FCT` vs `6` 错配、`perm100` 与 `E1` 配对差检验（`H-D1..H-D3`, Bonferroni `98.33%` CI）；(B) `H-D4` 检验 `E1(resp_top50)` vs `E1(all_nontarget)` 配对差确认选择偏倚；(C) `8 源线` 敏感性与 `S7`/`H-C3`。 |
| **Figure 7** | **Multi-Round Methodological Comparison & VCC 2026 Prospective Readiness** | `PRIOR_ART_A–D` + `A-T4` `cell-eval2` | **概念/验证图**：(A) 跨方法机制对比象限图（加性均值 `PertMean`/`β_p` vs 非加性 `log1p·FC` 缩放 vs `resp_top50` 选择偏倚）；(B) 公开 VCC 2025 数据集 `cell-eval2 v0.22.0` Signal vs Null 验证及时间戳锁定流程。 |

### 4.2 主表与附表清单（Main Tables 1–7 & Supplementary Tables）

| 编号 | 表题与核心内容 | 状态 |
| :--- | :--- | :--- |
| **Table 1** | **Unified Dataset Inventory Across Rounds 9, 10, & Planned R11 Candidate Metadata (`R11-0`)**（前 10 套已用细胞系完整参数 + `H9` P2 排除 + R11 元数据候选项如 `TeloHAEC_2285`、`HCASMC_971`、`Feng 2026`，严格标明 R11 仅元数据筛查） | R9 & R10 数值已备；R11 行暂仅列元数据预估，待 `PREREG_D` 锁定后更新实际入组计数 |
| **Table 2** | **Round 9 Main Confirmatory Matrix (`PREREG_A v4`, `K_conf = 3`)**（`Hybrid_ProjShared`, `TMVS`, `GlobalScale`, `R14`, `H-A1..A4` / `H-B1..B4`） | 数值完备 |
| **Table 3** | **Round 10 Main Confirmatory Matrix (`PREREG_C`, `K_conf = 2`, `resp_top50`) & `H-C3` / `S7` Summary**（`FCT`, `max mismatch`, `perm100`, `H-C1/H-C2/H-C3`, `S7` OLS intercept & Spearman, `SA-3` R9 anchors） | 数值完备 |
| **Table 4** | **Cross-Round `R19` Three-Tier Decomposition & `S6` Noise-Ceiling Comparison (`R9` + `R10`)** | R9 全 6 线 + R10 `dld1`/`cd4t_rest` 备齐（定稿时可由 `summary_holdout_r10_log1p_primary.json` 补齐 `calu3`/`i3Neuron` `R19` 行） |
| **Table 5** | **`[事后探索性 / Post Hoc]` `R10-E` Baseline Ablation (`E1–E7`) Across `resp_top50` vs. `all_nontarget` (10 Cell Lines)**（必标「事后探索性（post hoc）」） | **数值已由 `A-T2` 产出并经 `A-T3` 独立复核（`1.11e-16` 精度一致），直接摘自 `BASELINE_ABLATION.md` 表 1–2** |
| **Table 6** | **`[事后探索性 / Post Hoc]` Paired Difference `FCT − E1` (`97.5% CI`) & Frozen `S6` Specificity Audit**（必标「事后探索性（post hoc）」） | **数值已备齐（来自 `BASELINE_ABLATION.md` 表 3 与 `R10_VERDICT.md` 13:15 更正：`dld1 +0.1969 [+0.1673, +0.2253]`，`cd4t_rest +0.0426 [+0.0248, +0.0605]`；R9 全 6 线 `+0.0989..+0.2444`）** |
| **Table 7** | **`[ROUND 11 PROSPECTIVE PLACEHOLDER]` `PREREG_D` Confirmatory Outcomes (`H-D1..H-D4`, `all_nontarget`, `98.33%` CIs)** | **严格留空占位**（仅保留表头：`Target Line | Role | N | FCT | Max Mis CI98.33_hi | Perm100 p97.5 | E1 | FCT - E1 [98.33% CI] | H-D1/D2/D3 | H-D4 E1_top50 - E1_all | C3 Diff`） |
| **Supp. Tables S1–S6** | `S1`: R9 `S1–S5` + `S2b`；`S2`: R10 `H-C3` 6 口径全文 + `S1`/`S5`；`S3`: 两轮+`R10-E` 双/三方 `C3` 交叉核对总表（R9 `456+90`，R10 `349+364`，R10-E `510`）；`S4`: `PRIOR_ART_A–D` 矩阵（含 `[PENDING: Q4 & Q5]`）；`S5`: `A-T4` `cell-eval2 v0.22.0` VCC 2025 复现分；`S6`: `[PLACEHOLDER: R11 8-source sensitivity & S7/H-C3]` | `S1–S3, S5` 完备；`S4` 待补 Q4/Q5；`S6` 待 R11 完成后填入 |

---

## 5. 方法、机制解读与预注册偏差细纲（Detailed Methods, R10-E Diagnostics, & Deviations）

### 5.1 R9 & R10 数据与模型底座（复用并保持与冻结代码一致）
- **Round 9（6 线 `1,419,744` 细胞，`6,148` 基因）**与 **Round 10（4 线 `4,264,773` 单 sgRNA+NTC 细胞，`H9` 因 P2 下载前排除；`dld1` 通道 Split-Half、`cd4t_rest` 双供体 `CE0006864`/`CE0008162` Split-Half `P3`；HGNC 跨轮桥接挽回 `83->86` 与 `301->302` 靶点）** 完整细节沿用两轮冻结规范。
- **预测器族与对比基线**：
  - R9 Group A (`Hybrid_ProjShared`), Group B (`TMVS`), `R18` (`Transfer_Mean_GlobalScale`), `R14` (`Hybrid_on_TMVS`), R10 Group C (`FCT`, 归属 `CellularAI`/`AtlasShift` 映射类) 与对应错配/置换对照。
  - **R10-E 事后探索基线消融族（完全固定 `round10/src/` `f8d2ebfc`，不改动任何冻结文件；`4,000` 次联合 bootstrap，`seed = 20261122`，`resp_top50` 与 `all_nontarget` 双掩码）**：
    - `E1_ConstMed`: $c = \operatorname{median}(\mu_t)$，$\hat{\Delta}^{\mathrm{E1}}_{p,g} = \operatorname{log1p}(c\,\overline{\mathrm{FC}}_{p,g}) - \operatorname{log1p}(c)$（所有基因共享同一目标总体尺度常数，零逐基因信息）；
    - `E2_ConstP10` / `E2_ConstP90`: 常数分别取 $\mu_t$ 的 10% 与 90% 分位数；
    - `E3_SrcMeanCtrl`: 用源线 libnorm 对照均值平均 $\bar{\mu}_{\mathcal{S}}$ 替代 $\mu_t$；
    - `E5_RankMatchedMu`: 保留目标 $\mu_t$ 经验分布，但按 $\bar{\mu}_{\mathcal{S}}$ 排序重新分配基因位置；
    - `E6_LogFCBar`: $c \to +\infty$ 渐近极限 $\ln(\overline{\mathrm{FC}}_{p,g})$；
    - `E7_FCBarMinus1`: $c \to 0^+$ 渐近极限 $\overline{\mathrm{FC}}_{p,g} - 1$。

### 5.2 第十轮 `13:15` 探索性更正与选择偏倚机理（写入方法/结果衔接小节）
1. **初期错觉与正式更正流程**：主控最初（12:50）在自写验证脚本中仅观察了主掩码 `resp_top50` 下的 `E1/E3/E5`，曾假定“由于 `E1` 复现 `FCT` 的 `86%–112%`，`FCT` 优势主要源于变换本身而非目标基线”。随后 `A-T2` 与 `A-T3` 使用正式冻结模块在 `resp_top50` **和** `all_nontarget` (`S6`) 两种掩码下完整复算 `E1–E7`，并在 13:15 正式推翻该初期解读。
2. **已验证的掩码分叉核心数据（`BASELINE_ABLATION.md` & `AUDIT_R10E_REPORT.md`，*post hoc*）**：
   - **`resp_top50` 下（证实偏倚）**：`E1` 在大样本确证系与 R9 线复现 `FCT` 的 `86.4%–115.5%`（`dld1`: `FCT +0.1401` vs `E1 +0.1567` [`111.8%`], `FCT - E1 = -0.0166 [-0.0708, +0.0396]`; `cd4t_rest`: `+0.0560` vs `+0.0495` [`88.4%`], `diff +0.0065 [-0.0136, +0.0270]`; `nadig_jurkat`: `+0.1047` vs `+0.1209` [`115.5%`], `diff -0.0162 [-0.0326, +0.0006]`; `K562`: `+0.0039`; `RPE1`: `-0.0120`），且 `E6`（`dld1 +0.1492`, `cd4t +0.0498`, `HepG2 +0.2349`, `H1 +0.3117`）与 `Perm100 Med`（`+0.0892`, `+0.0347`, `+0.1264`, `+0.2101`）同步虚高。
   - **`all_nontarget` 下（剥离偏倚）**：`E1`、`E2`、`E6`、`E7` 与 `perm100` 全部塌缩至零附近或负值（`dld1`: `E1 = -0.0614`, `E6 = -0.0713`, `E7 = -0.0541`, `perm100 p97.5 = -0.0352`; `cd4t_rest`: `E1 = +0.0130`, `E6 = +0.0090`, `E7 = +0.0169`, `perm100 p97.5 = +0.0131`; `HepG2`: `E1 = -0.0470`; `Jurkat`: `E1 = -0.0782`; `K562`: `-0.0377`; `RPE1`: `-0.0934`）。相反，真实 **`FCT` 在全部 `8/8` 有效线上保持显著正值**（`dld1 +0.1355 [+0.1082, +0.1649]`，`cd4t_rest +0.0556 [+0.0376, +0.0741]`，R9 六线 `+0.1081` 到 `+0.1848`），且 **配对差 `FCT - E1` 的 `97.5% CI` 在全部 `8/8` 线上均严格 `> 0`**（`dld1 +0.1969 [+0.1673, +0.2253]`，`cd4t_rest +0.0426 [+0.0248, +0.0605]`，R9 `+0.0989` 到 `+0.2444`）。`FCT` 点估计亦在 `8/8` 线上全面超越 `E3` 与 `E5`。
   - **回看 R10 预注册次要分析 `S6`（`all_nontarget`）**：`dld1` 达到特异性（`FCT CI_lo +0.1082 > max mismatch CI_hi +0.0793` [`replogle_rpe1`] 且 `> 事后 perm100 p97.5 -0.0352`），而 `cd4t_rest`（`CI_lo +0.0376`）距错配均值上界 `+0.0395`（`vcc2025_h1`）仅差 `0.002`。

### 5.3 第十一轮 `PREREG_D` 方法框架占位（严格未运行、不写结果）
- **盲态锁定顺序**：`PREREG_D` 提交锁定 $\to$ `R11_GO.txt` $\to$ `A-T1` 下载与仅计数 QC $\to$ 确认线名单锁定 $\to$ 代码冻结 $\to$ 留出评测与独立复算。
- **候选数据集与功效预估占位（来自 `R11_DATA_SCREEN.md`）**：
  - 优先未用新谱系：`TeloHAEC_2285`（`IGVFDS7943IIWZ`/`GSE210681`，静息内皮细胞，`1.47 GB`，`214,449` 单 sgRNA 细胞，`>2,200` NTC，预计 R9 `75–150` 靶点）；`HCASMC_971`（`IGVFDS6248NDUW`，血管平滑肌 `ZIM3-KRAB`，`.hdf5 42.89 GB`，预计 R9 `40–80` 靶点）；必要时评估 `Feng 2026 TargetedScreen` 非 KOLF iPSC（`4.48 GB`，预计 `60–115` 靶点，但属已知 iPSC 状态）。鉴于 `cd4t_rest` 在 `all_nontarget` 仅差 `0.002`，`PREREG_D` 须仅凭 R10 `dld1`/`cd4t_rest` 已知 CI 宽度（不碰新数据）给出所需最小靶点数 `N` 的功效预估。
- **主掩码与四条预注册假设**：
  - 主掩码固定为 `all_nontarget`（`resp_top50` 降为次要）；主源池为 R9 `6` 线（附加 `R9+R10 dld1/cd4t_rest` `8` 源线敏感性）；`4,000` bootstrap，`seed = 20261123`，`K_conf >= 2` 且全票一致。
  - **三条共同主假设（Bonferroni 双侧 `98.33%` CI `[0.833%, 99.167%]`）**：`H-D1`（`FCT r_orth_cal CI98.33_lo > 0`）、`H-D2`（`FCT CI98.33_lo > max(错配对照均值 CI_hi, all_nontarget perm100 p97.5)`）、`H-D3`（同一组 bootstrap 下配对差 `FCT - E1 CI98.33_lo > 0`）。
  - **方法学次要假设 `H-D4`（单独报告）**：`E1(resp_top50) CI_lo > 0` **且** 配对差 `E1(resp_top50) - E1(all_nontarget) CI_lo > 0`。
  - **其余同步保留**：置换零分布同时在 `all_nontarget` 与 `resp_top50` 下计算；保留 `S7`、`H-C3` 非劣效、`K_conf < 2` 自动退回描述性报道。

### 5.4 与预注册的偏差汇总（Deviations Across Rounds）
1. **`R9 D1`**：`PREREG_A v4 §4.2` 文本秩 1 $q_1$ vs. 冻结代码 `1ca9f60f` 三维 $\operatorname{span}\{\mathbf{1}, \bar{b}_S, a_S\}$（`S2b` 秩 1 敏感性证实判定不变；`A-T3` 在该点复刻 3D 投影的非独立性已披露；R10 `P7` 明确锁定 3D 锚点）。
2. **`R10 D-R10-1`**：`H-C3` `ctrl_disc_ownmasked`/`panel` 论文命名 vs. 冻结代码 `ctrl_disc_raw`/`raw_panel`（数学完全一致）。
3. **`R10 D-R10-2`**：`H-C1` 代码额外保护条件（`valid_frac >= 0.95`, `rho_ab_orth > 0`），两条确认线均严格通过（`1.0`, `> 0.64`），无判定分歧。
4. **`R10` 微小代码差异与文书修正**：`fct_from_fc_bar` 目标 $\mu_t$ 下限 `0` vs `1e-8`（影响 `~1e-8`）；`H9` P2 排除原因校订为双重理由（单 sgRNA 靶点稀疏 + 60 GB 预算）；**`R10_VERDICT.md` `13:15` 更正记录（公开记录主控初期事后解读被正式双掩码消融推翻的过程）**。

---

## 6. 局限性（Combined Limitations: R9 + R10 + R10-E + `[R11 Placeholder]`）

1. **执行智能体同源性与 R9 3D 投影非独立性**：两轮执行 `A-T2` 与审计 `A-T3` 均属 `Antigravity`（Gemini）；虽主控 `Claude`（Anthropic）另写代码独立核对了 R9 `TMVS`、R9 `HepG2/H1 FCT`、R10 全确认集与初步消融原型，但 R9 Group A（`Hybrid_ProjShared`）未由跨模型主控单独复算，且 `A-T3` 复刻了 R9 的 3D 投影。
2. **确认样本量与谱系受限（R9 `K_conf = 3`, R10 `K_conf = 2`）及 `cd4t_rest` 跨供体天花板（`P3`）**：R10 最终仅含 `DLD-1` 与 `CD4+ T` 两条确认线，`cd4t_rest` 因双供体+双批次分割引入真实生物学供体变异（$\rho^{\perp}_{AB} = 0.6479$，`284/302` 半有效），恰解释了其在 `S6` 上以 `0.002` 窄幅失利；该局限亦将通过 `[PLACEHOLDER: R11 新谱系如 TeloHAEC/HCASMC]` 缓解。
3. **可重复靶点过滤范围与 `S1` 全部靶点表现落差**：`PREREG_A`/`C` 聚焦于 Addendum-B 可重复强效应单 sgRNA 扰动；当 R10 `S1` 纳入非可重复靶点时，`FCT` 相关减半且 `H-C3` 在质心口径出现 `4` 项劣效。
4. **`R10-E` 的事后探索属性及 `all_nontarget` 低信噪比代价**：尽管 `all_nontarget` 完全消除了 `resp_top50` 选择偏倚，全非靶基因矩阵中混杂大量未变化/低信号基因，令噪声天花板大幅走低（`dld1 0.8277 -> 0.2533`、`cd4t_rest 0.6479 -> 0.1468`），对去衰减检验的渐近稳定性提出更高样本要求——此正需要第 11 轮预注册 `PREREG_D` 独立验证。
5. **外部 VCC 2026 斜率公式闭源与测试集真值未公开（`R6`, `A-T4`）**：`cell-eval2 v0.22.0` 未公布 `0.02–0.13` 斜率计算函数且 `10.64898/2026.11.17.767214` 暂未见索引；`gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` 截至 2026-11-22 为空，`AC16`/`MOLM13`/`U-87MG` 对照与真值受条款与释放进度限制暂不可用于开放研究验证。
6. **文献查新边界与待补 Q4/Q5（`PRIOR_ART_D.md`）**：已完成 Q1–Q3（确认 FCT 常数基线消融与 `log1p` 几何平均制造虚假 `r_orth_cal` 增益未见前作，而上下文均值解释零样本模型已有 `CGC`、`Molina & Zhang`、`CASTER` 证实）；**`[PLACEHOLDER: 待 A-T5 完成 R11_LIT_DONE（Q4 关于按观测响应选基因造成选基因偏倚的文献；Q5 关于跨细胞类型扰动在 FC 空间 vs log1p 差值空间守恒性的文献）后补全定位]`**。
7. **`[ROUND 11 LIMITATIONS PLACEHOLDER]`**：留空待 `PREREG_D` 执行完成后按实补充（含任何 R11 预处理偏差、实际 `K_conf`、以及如启用高 MOI 或新数据集产生的特定限制）。

---

## 7. 待完成与合稿交接核对清单（Checklist Before Drafting Final Combined Paper）

- [x] **A-T2 (`R10E_ABLATION_DONE`)** & **A-T3 (`R10E_AUDIT_DONE`)**：`BASELINE_ABLATION.md` 与 `AUDIT_R10E_REPORT.md` 已完成并交叉复核通过（`510` 项最大差 `1.11e-16`）；`R10_VERDICT.md` `13:15` 更正已落地。
- [x] **A-T5 (`R10E_LIT_DONE`)**：`PRIOR_ART_D.md` Q1–Q3 已完成。
- [x] **A-T1 (`R11_SCREEN_DONE`)**：`round10/shared/R11_DATA_SCREEN.md` 元数据筛查已完成。
- [ ] **A-T5 (`R11_LIT_DONE`)**：完成 `PRIOR_ART_D.md` Q4（观测选基因造成的选择偏倚）与 Q5（FC 空间 vs `log1p` 差值空间跨系守恒性）$\to$ 补全本文 §1、§6 与 Supp. Table S4。
- [ ] **A-T3 + A-T1 + A-T2 (`PREREG_D` 全流程执行)**：`round11/prereg/PREREG_D.md` 锁定 $\to$ `R11_GO.txt` $\to$ 下载与计数 QC $\to$ `K_conf` 确认线锁定 $\to$ 代码冻结 $\to$ `H-D1..H-D4` 留出验证与 `C3` 复算 $\to$ **再填充本文摘要 §2.3、Section 6、Figure 6、Table 7 与 Supp. Table S6**。

---

## 8. 核心参考文献目录（已核验引文 + 占位符）

- **Round 9 已有验证文献**：Adduri et al. (2026) *Cell* (`10.1016/j.cell.2026.08.004`); Arc Institute VCC 2026 Wrap-Up (`2026-11-19`) & `cell-eval2` (`v0.22.0`); Dhaka et al. (2026) `HyperMap` (`10.64898/2026.04.23.720505`); Huang et al. (2026) `CGC` (`10.64898/2026.09.30.755756`); Kawsar et al. (2026) `CASTER` (`10.64898/2026.11.15.766815`); Kuznets-Speck et al. (2025/2026) `CIPHER v2` (`10.1101/2025.06.27.661814`); Luan et al. (2026) (`10.64898/2026.10.20.760514`); Miller et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03307-w`, `PMID:42823510`); Molina & Zhang (2026) (`10.64898/2026.07.24.740459`); Nadig et al. (2025) *Nat. Genet.* (`10.1038/s41588-025-02169-3`, `GSE264667`); Nourreddine et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03199-w`); Replogle et al. (2022) *Cell* (`10.1016/j.cell.2022.05.013`); Shoeibi & Yousefi (2026) *arXiv*:`2608.00152v3`; Sun et al. (2026) `AMM-SimWMag` *Genes* (`10.3390/genes17070816`); Viñas Torné et al. (2025) `Systema` *Nat. Biotechnol.* (`10.1038/s41587-025-02777-8`).
- **Round 10 数据集与元数据出处**：Yeung J et al. (2026) `DLD-1` (*bioRxiv*, `10.64898/2026.07.10.737863`, `GSE337988`); Zhu X et al. (2026) `CD4+ T` (*Cell*, `10.1016/j.cell.2026.08.002`, `PMID:42664972`, `GSE314342`); Tian R, Kampmann M et al. (2021) `WTC11 i3Neurons` (*Nat. Neurosci.*, `10.1038/s41593-021-00862-0`, `PMID:34031600`, `GSE152988`/Zenodo `13350497`); Sunshine S, Hein MY, Weissman JS et al. (2023) `Calu-3` (*Nat. Commun.*, `10.1038/s41467-023-41788-4`, `PMID:37803001`, `GSE208240`/Zenodo `13350497`); Gogate et al. (2026) `H9` (`IGVFDS9825FOBP`, `PMC13629240`, `PMID:42827575`).
- **`PRIOR_ART_D.md` 新增核查文献（Q1–Q3 已验证，Q4–Q5 待补充）**：Illumina `CellularAI` (`GitHub Illumina/CellularAI` & `2026` Blog); `AtlasShift` (`GitHub kaipengm2/Virtual-Cell-Challenge-2026`); Xaira `X-Cell`; Adduri et al. (2025) `Arc State` (`10.1101/2025.06.26.661135`); Yan et al. (2026) `PRiMeFlow` (`arXiv:2604.13986v2`); Booeshaghi & Pachter (2021) *Bioinformatics* (`10.1093/bioinformatics/btab085`); Kviman et al. (2025/2026) (`10.1101/2025.03.12.642799`); Liu et al. (2025/2026) (`arXiv:2511.16954`); Panagopoulos et al. (2026) (`arXiv:2610.26267`); Ahlmann-Eltze & Huber (2023) (`10.1038/s41592-023-01814-1`); Ahlmann-Eltze et al. (2025) (`10.1038/s41592-025-02772-6`); Csendes et al. (2025) (`PMC12016270`); Li et al. (2025) (`PMC12202205`); `PertEval-scFM` (`10.1101/2024.10.02.616248`); Kernfeld et al. (2025) (`10.1186/s13059-025-03840-y`); Kedzierska et al. (2025) (`10.1186/s13059-025-03574-x`); Bendidi et al. (2024) (`arXiv:2410.13956`); Wu et al. (2024/2025) `PerturBench` (`arXiv:2408.10609`); Li et al. (2026) *Sci. Adv.* (`10.1126/sciadv.aed3414`, `PMC13557097`); Dibaeinia et al. (2026) (`10.64898/2026.02.04.703804`, `PMC12919078`); Ghajari et al. (2026) (`10.64898/2026.11.04.764088`).
- **Round 11 潜在候选数据文献（来自 `R11_DATA_SCREEN.md` 元数据）**：Schnitzler et al. (2024) *Nature* `TeloHAEC_2285` (`IGVFDS7943IIWZ` / `GSE210681`, `PMC10921916`, `PMID:38326615`); `HCASMC_971` (`IGVFDS6248NDUW`); Feng et al. (2026) *Cell Genom.* (`PMC12903452`, `PMID:41330380`, Figshare `27989294`); Sivakumar et al. (2025) `WTC-11 CM` (`IGVFDS6332VCTO`, `PMC11785042`); `[PENDING: PRIOR_ART_D Q4/Q5 & R11 最终入组线文献]`.
