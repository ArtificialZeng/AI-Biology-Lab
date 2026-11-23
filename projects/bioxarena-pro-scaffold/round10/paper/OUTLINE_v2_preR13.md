# 第九至十二轮合稿论文骨架（`round10/paper/OUTLINE.md`）

> **先前版本存档**：第十一轮末版已完整另存为 [`round10/paper/OUTLINE_R11.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/paper/OUTLINE_R11.md)。

- **定位与核心叙事（Round 9 $\to$ Round 10 $\to$ R10-E / R11 Post-Hoc Provenance $\to$ Round 12 `K_conf = 1` Confirmatory Core $\to$ Round 13 Placeholder）**：
  1. **第九轮（`PREREG_A v4`，6 线 LOLO `K_conf = 3`）**：预注册的协方差共波动转移（`Hybrid_ProjShared`）与边际方差缩放（`TMVS`）方法**均未确认**（`0/4` 标准）；共波动在边际方差之外未提供可检测增量（`R14`）；单一标量缩放即可制造原始残差假高相关而在正交去幅度 `r_orth_cal` 下精确归零（`R18`）。
  2. **第十轮（`PREREG_C`，4 条单 sgRNA 人类线 `K_conf = 2`）**：在预注册主掩码 **`resp_top50`** 下，基线倍数缩放模型 `Fold-Change Transfer`（`FCT`）通过去幅度正信号（`H-C1 2/2`）与并行非劣效（`H-C3 12/12`），但**未通过细胞系特异性（`H-C2 0/2 -> NOT CONFIRMED`）**。
  3. **第十轮事后消融与第十一轮（`R10-E` & `PREREG_D v2`，`K_conf = 1` 仅作描述性报告）**：双团队复算发现 `resp_top50` 按目标 $|\Delta^{\mathrm{obs}}_t|$ 选基因制造了选择偏倚；第十一轮在全非靶基因口径 **`all_nontarget`** 下因 `hcasmc` 缺 guide 分配与 `telohaec`（`N = 9`）未过门槛而使 `K_conf = 1 < 2`，按 `P1` 中止确认性检验、全部结果仅作描述性报告（唯一达标线 `cd4t_stim48` 在不配对口径仍败于 `vcc2025_h1`）。随后对 R9–R11 共 `12` 条线开展的事后探索检查（[`REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/REPORT.md)）发现：**两项新现象形成第十二轮两个论断的假设来源**——(i) 在全非靶基因 `all_nontarget` 下**配对对比** `FCT − Mis_s` 可解决两条 T 细胞线不配对检验功效不足的问题；(ii) 即使在分半独立选基因（`SA-8`，`resp_top50_split`）下，不用目标线基线的预测器 `E1_ConstMed` 仍持续高企。
  4. **第十二轮（`PREREG_E`，全文确认性核心；`K_conf = 1`，仅 `A549` CRISPRko，推广性有限）**：在全新跨细胞系与跨扰动方式测试集 **`a549_ko`（`GSE345058` CRISPRko，`N = 445`）** 上对上述两项来自 R9–R11 事后观察的假设进行盲态前瞻检验，**论断 1（`H-E1`，靶线特异性）与论断 2（`H-E2`，`SA-8` 偏倚）均获三方一致确认**：
     - **确认性核心（论断 1）**：按 `PREREG_E §8` 严格表述为**一项预注册的前瞻检验：一个只用基线的变换（`FCT`）在从未用过的细胞系和扰动方式上、在去幅度的残差指标（`r_orth_cal`）上胜过 `6` 个错配细胞系基线**（`a549_ko` `6/6` 配对差 `97.5% CI > 0`，并胜置换对照）。同时必须写明：**`FCT` 不是最好的预测器（`Source_Grand_Mean` `+0.405` > `CIPHER_Linear` `+0.282` > `Hybrid_on_TMVS` `+0.237` > `TMVS` `+0.222` > `FCT` `+0.190`）——论断 1 讲的是特异性（目标线自身对照谱带来特异信息），不是预测精度。** 必须引用 **Kawsar et al. 2026（`CASTER`）**，不声称首次做正确与错配对照的对比。
     - **次要方法学观察（论断 2，不作核心发现）**：即使去掉共享反应、用独立细胞分半选 top-50 基因，不用目标线基线的 `E1_ConstMed` 仍从全基因近零（`-0.020`）升至 `+0.257`（`H-E2b +0.277`），提示检验靶线特异性应使用全非靶基因。必须引用 **Nicol et al. 2026 v2** 与 **Heidari et al. 2026** 并写明各自已做工作；该现象仅限 `r_orth_cal` 空间、**机制未知**（R11 真实 $\mu$ 无噪声模型在全基因上预测错了 `8/12` 线符号）。
     - **单线与偏离披露**：因偏离 `D3`（按 `PREREG_E §3.2–§3.3` 将 `HEK293` NTC 扣除改为按半扣除），`HEK293` 评测靶点从先前 `SUPERSEDED_v1` 的 `62` 个（`K_conf = 2`）降至 `16 < 50`，机械退为次要线（`SECONDARY_ONLY`），最终仅 `A549` CRISPRko 一条确认性细胞系（**`K_conf = 1`，推广性有限**）。
  5. **第十三轮（`GSE327644` 原代人巨噬细胞，占位小节）**：正在预注册阶段；大纲仅留占位小节，不写任何预期结果。

- **锁定与版本锚点**：
  - **Round 9**：`PREREG_A v4`（commit `4cc89ebae67f0ed0748b80f994d26dcc9c2b95ef`），`CODE_FREEZE`（`1ca9f60f49adf1ac7cce8266e40533273411db30`）。
  - **Round 10 & R10-E**：`PREREG_C`（`dc4b4a5`），`QC_LOCK`（`86872ed`），`CODE_FREEZE_R10`（`f8d2ebfc`，`FREEZE_MANIFEST_R10.txt` `22/22 OK`），`R10_VERDICT.md`（含 `13:15` 更正），`BASELINE_ABLATION.md`，`AUDIT_R10E_REPORT.md`，`PRIOR_ART_D.md`。
  - **Round 11**：`PREREG_D v2`（`72fc0c937358b089e78664b464f9db5aafb790b1`，`sha256 1a054498...`）+ `PREREG_D_ADDENDUM_A.md`，`CODE_FREEZE_R11`（`cc834d863105c66874acb4c27cc214fd79704955`，`FREEZE_MANIFEST_R11.txt` `41/41 OK`），`DEVIATIONS_R11.md`（`D1–D4`），`CONFIRMATORY_LINES_R11.txt`（`K_conf = 1`），`R11_VERDICT.md`，`EXPLORATORY_REVIEW_R11.md`（含 `§D`），`PRIOR_ART_E.md`，`round11/audit/exploratory_check/REPORT.md`。
  - **Round 12（确认性核心）**：`PREREG_E`（commit `4c606f4732ea080768e27757c689d40d7892cd3d`，`sha256 d10f92ced16a4519c43ea2b92b8b7ee881ab3d1764e39ffd10580473e4108924`），`CODE_FREEZE_R12`（commit `11d91c72025264b838af091abcf0eac71588eb91`，`FREEZE_MANIFEST_R12.txt` `46/46 OK`），`R12_RESULTS_COMMIT`（`5232aa5466d1fa3dca57a8bbe160a20dc04ab46f`），`DEVIATIONS_R12.md`（`D1–D4`），`RULING_R12_PREFREEZE.md`，`RULING_R12_HEK293_UNITS.md`，`CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt`（`K_conf = 2` 旧版记录），`CONFIRMATORY_LINES_R12.txt`（`K_conf = 1`），[`R12_VERDICT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/R12_VERDICT.md)（`2026-11-23 07:47`），`PRIOR_ART_F.md`。
  - **Round 13（在研/预注册）**：[`PRIOR_ART_CONCLUSION_R13.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/status/PRIOR_ART_CONCLUSION_R13.md)（`2026-11-23 08:25`），`PRIOR_ART_G.md v2`，`GSE327644` 原代巨噬细胞（正处预注册阶段）。
- **AI 治理与作者披露声明（必保留；人类角色仅称 `Human project lead` / `人类项目负责人`）**：人类项目负责人 Zijian Zeng（Human project lead Zijian Zeng，无正规生物学背景）；主控 `Claude`（Anthropic）负责预注册审核、结果核验、裁决与偏差审查及独立复算；执行智能体 `Antigravity`（Google Gemini 模型）负责数据整合（`A-T1`）、主分析与消融（`A-T2`）、独立复核与 `C3` 交叉核对（`A-T3`）、合稿骨架与基准校验（`A-T4`）及 API 文献/源码审计（`A-T5`）。

---

## 1. 标题候选（Title Candidates）

> **措辞原则**：去掉已被第十二轮推翻的“细胞系特异性未确认”等表述；不得夸大预测精度或新颖性；**每个候选标题均须明示仅有一条确认性细胞系（`K_conf = 1`，推广性有限）**，且不使用过度绝对化的英文动词。

1. **全周期写实主推型（突出跨模态单线确认 + 推广性有限 + 方法学审计）：**  
   *Prospective Confirmation of Target-Baseline Specificity Over Mismatched Controls in One Unseen CRISPRko Cell Line (`A549`, $K_{\mathrm{conf}}=1$, Limited Generalizability) Following Multi-Round Pre-Registered Audits*
2. **特异性 vs 精度辨析 + 单线界定型：**  
   *Separating Specificity From Accuracy in Zero-Shot Perturbation Prediction: Single-Line Prospective Confirmation ($K_{\mathrm{conf}}=1$, Limited Generalizability) of Fold-Change Transfer Over Mismatched Baselines and Persistent Split-Half Selection Inflation*
3. **四轮预注册因果脉络型：**  
   *Multi-Round Pre-Registered Audits Confirm Cross-Modality Baseline Specificity Over Mismatched Controls in `A549` CRISPRko ($K_{\mathrm{conf}}=1$, Limited Generalizability) While Exposing Split-Half Selection Bias Under Orthogonal Residuals*

---

## 2. 摘要骨架（Abstract Skeleton）

> **摘要硬约束（依据 `PREREG_E §2.2, §7.2, §8`、`R12_VERDICT.md`、`PRIOR_ART_CONCLUSION_R13.md`）**：  
> (1) 必须写明**只有一条确认性细胞系（`A549` CRISPRko，`K_conf = 1`，推广性有限）**；  
> (2) 必须写明**第十二轮两个论断的假设均源于第 9–11 轮的事后观察**；  
> (3) 论断 1 严格按 `PREREG_E §8` 表述，**引用 `Kawsar et al. (2026)`（`CASTER`），不声称首次**，并明确写明 **`FCT` 不是最好的预测器（`Source_Grand_Mean` `+0.405 > FCT +0.190`），论断 1 讲的是特异性而非精度**；  
> (4) 论断 2 **仅作为方法学观察，不作核心发现**，引用 **`Nicol et al. (2026) v2`** 与 **`Heidari et al. (2026)`**，写清机制未知且仅限 `r_orth_cal`。

### 2.1 第一段：背景与第 9–11 轮预注册排查（Post-Hoc Hypothesis Provenance）
1. **评估混淆与早期否定（Rounds 9–10）**：单靠目标未扰动对照零样本预测 CRISPR 扰动响应在 VCC 2026 中斜率偏低。首先证实每细胞系单一标量缩放（`Transfer_Mean_GlobalScale`）即可制造高原始相关性，而在单扰动内对共享跨线均值 $\hat{\Delta}^{\mathrm{TMean}}_p$ 正交消幅（`r_orth_cal`）后精确归零（`0.000`；`R18`）。在第九轮（`PREREG_A v4`，`K_conf = 3`：`nadig_hepg2`/`nadig_jurkat`/`vcc2025_h1`）中，`Hybrid_ProjShared` 与 `TMVS` 均未确认；第十轮（`PREREG_C`，`K_conf = 2`：`dld1`/`cd4t_rest`）在 `resp_top50` 下检验基线倍数缩放变换 `FCT`，因未能同时越过 6 个错配基线上界而未确认。
2. **第 10–11 轮掩码偏倚、描述性审计与假设形成**：双团队消融（`R10-E`）与第十一轮（`PREREG_D v2`，因 `hcasmc` 缺 guide 与 `telohaec` `N = 9` 而致 `K_conf = 1` 仅作描述性报告）揭示 `resp_top50` 严重抬高只用标量基线中位数、不使用目标线基线谱的预测器 `E1_ConstMed`。对 R9–R11 共 12 条线开展的事后探索（`REPORT.md`）提出两条待验证假设：(i) 全非靶基因口径 **`all_nontarget`** 下**配对差 `FCT − Mis_s`** 在 11/12 线上优于 `Mis_s`；(ii) 即使分半独立选基因（`SA-8`，`9/12` 线 `E1(split) − E1(all) > 0`），`E1_ConstMed` 仍被抬高。

### 2.2 第二段：第十二轮前瞻确认（确认性核心：`K_conf = 1`，`A549` CRISPRko，推广性有限）
1. **单线门控与偏差披露**：在完全未见过的第十二轮候选池（`PREREG_E`）上，`HEK293`（`XuCao2023.h5ad`）因偏差 **`D3`** 将 NTC 扣除单位由错误单批次（旧 `SUPERSEDED_v1` `62` 靶点）更正为按半扣除后仅余 `16 < 50` 个评测靶点而退出确认性线、归为次要描述线，故**仅有 `a549_ko`（`Liu et al. 2026`, `GSE345058`, CRISPRko, `N = 445`, `|G_eval| = 6,121`, $\rho^\perp_{AB} = 0.789$）一条确认性细胞系（`K_conf = 1`），因此所有确认性结论仅依赖单一细胞系，推广性有限。**
2. **确认性核心发现——论断 1（`H-E1`，靶线特异性确认；三方复核最大差 $\le 1.9\times 10^{-16}$）**：承接 `Kawsar et al. (2026)`（`CASTER`）在神经网络中对单错配对照做逐靶点替换的先例（本文不声称首次），本项预注册的前瞻检验显示，在 **只用基线的变换（`FCT`）** 下，该变换在 **从未用过的细胞系（`A549`）和扰动方式（`CRISPRi` 源 $\to$ `CRISPRko` 靶）** 上，在去幅度残差指标 **`r_orth_cal`**（`all_nontarget`）上得分为 **`+0.1898 [+0.1742, +0.2055]`**（`97.5%` Bonferroni CI），**在配对对比中全部 `6/6` 胜过 6 个错配细胞系基线**（差值 `+0.0796 [+0.0693, +0.0900]` vs `nadig_jurkat` 至 `+0.1869 [+0.1693, +0.2049]` vs `nadig_hepg2`），并胜过 100 次基因置换（`Q0.975 = +0.0018`）。**然而，`FCT` 并非最好的预测器：同一指标下，不用任何目标线信息的 `Source_Grand_Mean` 达 `+0.405`（`+0.4051`），`CIPHER_Linear` `+0.282`、`Hybrid_on_TMVS` `+0.237`、`TMVS` `+0.222` 均高于 `FCT` 的 `+0.190`；论断 1 证明的是目标线自身对照谱带来了特异信息（特异性），而非预测精度领先。**

### 2.3 第三段：方法学副观察（论断 2，非核心发现）与先占边界、第十三轮进展
- **方法学观察——论断 2（`H-E2` / `SA-8` 确认，非核心发现）**：在 `A549` 上，保留源线平均倍数变化但**不用目标线自身基线谱的预测器 `E1_ConstMed`** 在全基因上近零（`-0.0201`），但在独立细胞分半选 top-50 基因下仍跃升至 **`+0.2565 [+0.2195, +0.2981]`**（跨掩码差 `+0.2766 [+0.2412, +0.3164]`，保留约 `90%` 同半选基因虚高）。参照已报道拆半 DEG 评估及有信息基线正分现象的 **Nicol, Shivakumar & Irizarry (2026) v2**（`Split DEG + control` 下 `ctl-mean` 归零而 `nonctl-mean` 仍呈正相关）与多轴基准工作 **Heidari et al. (2026)**（`CrossSplit`/`*_vs_topdeg`/`define_trivial`/`rigorous_pearson_delta_pert`/`baseline_logfc`），本文仅将其作为一条**机制未知**、严格局限于正交残差指标 `r_orth_cal` 的**方法学观察**（不外推至 DEG 指标）：检验靶线特异性须在全部非靶基因上执行。后续第十三轮（`GSE327644` 原代人巨噬细胞）正在预注册中。

---

## 3. 论文总体结构与章节骨架（Manuscript Structure）

1. **Title & AI Governance / Role Disclosure**（`K_conf = 1` 与推广性有限进标题；`Human project lead Zijian Zeng`）
2. **Abstract**（三段式：R9–R11 假设来源 $\to$ R12 确认性核心 `H-E1` + `K_conf = 1` + `FCT` 非最佳预测器 + `CASTER` $\to$ 论断 2 非核心方法学观察 + `Nicol v2` / `Heidari 2026` + 机制未知 + R13 预注册占位）
3. **Section 1: Introduction**
   - 1.1 跨细胞系零样本 CRISPR 扰动预测的核心追问：靶线自身对照表达是否携带超越共享均值与错配对照的特异信息？
   - 1.2 方法学陷阱梳理：`R18` 单标量幅度混杂、`R14` 共波动 vs 边际方差、以及从同半选基因到独立分半选基因的掩码抬高（衔接 *Nicol et al. 2026 v1/v2*、*Heidari et al. 2026*、*Miller et al. 2026*）。
   - 1.3 领域已有对照基准：*Kawsar et al. (2026, CASTER)* 在部分筛选适配设定下的正确 vs 单邻系错配对照、*Heidari et al. (2026)* 的 `baseline_logfc`、与本文盲态前瞻 `all_nontarget` 配对 IUT 检验（`CRISPRi -> CRISPRko`）之区别（**明确不声称首次**）。
   - 1.4 四轮预注册递进与假设来源透明声明（`E6`）：R9–R10 原设确认检验 $\to$ R10-E/R11 事后观察催生 `H-E1` 与 `H-E2` $\to$ R12 盲态前瞻检验 `PREREG_E` $\to$ R13 正在预注册。
4. **Section 2: Methods, Pre-Registration Governance, & Deviations**
   - 2.1 R9–R12 人类 CRISPRi / CRISPRko 数据全集（16 项人类数据筹备，严格单 sgRNA `P4` 与严格 NTC 扣除规则，`Liu et al. 2026 GSE345058`, `Xu et al. 2023 GSE218566` 等）。
   - 2.2 标准化、HGNC 映射（含 `A549` `QARS` `591` 细胞缺失映射不入 `G_eval` 说明）、`P5` 靶点准入与独立分半物理/确定性规则（`A549` 32 条 `Lane` 严格嵌套 `Half A/B`；`HEK293` 2-stage 规则触发 Stage 2 `seed=42` 分层对半）。
   - 2.3 预测器精确定义与机制边界修正：`FCT`、`Source_Grand_Mean`、`Transfer_Mean_VarScaled`、`CIPHER_Linear`、`Hybrid_ProjShared`、`Hybrid_on_TMVS`、6 个 `NegCtrl_FCT_MismatchedMean(s)`、`Perm100`、以及**严格限定 `E1_ConstMed` 定义为「包含源线平均倍数变化 $\overline{\mathrm{FC}}_{p,g}$、但省略目标细胞系自身对照基线谱 $\mu_{t,g}$ 的迁移预测器」**（明确区分于完全零信息对照）。
   - 2.4 评估体系与检验规程：`r_orth_cal`（`all_nontarget` 与 `SA-8` `resp_top50_split`）、IUT + 2-Claim Bonferroni `97.5%` 置信区间（`seed = 20261201`, `n_boot = 4,000`）、R12 非对称门控规则（`A549` 必需 + `HEK293` 条件准入 + `K_conf = 1` 强制推广性声明）及仅报告非确认 `§5` 分析（含 `N_5`、未配对 `H-D2`、`SA-Cas9`、`SA-Src8`）。
   - **2.5 全部预注册偏差与裁决集中记录（重点含 R12 `D1–D4` 与 `SUPERSEDED_v1`）**：
     - R9 `D1`（3D vs rank-1）；R10 `D-R10-1/2` & `13:15` 更正；R11 `D1–D4`（`hcasmc` `Q3`, `telohaec` 弱敲低 `N=9`, `P1` 中止）。
     - **R12 `D1–D4`（均于任何指标计算前依据计数与元数据决定）**：`D1`（`HEK293` 禁用继承自 R11 代码的 `20%` MT 默认上限，依 `PREREG_E §3.2` 保留 `98,274` 细胞；`.obs` `percent_mito` 中位数 `25.81%`）；`D2`（AnnData `attrs['_index']` 与 `A549` `int32` indptr 溢出修正）；**`D3`（`HEK293` Stage 2 NTC 扣除单位由误设单一 `batch_0` 更正为按半 `half_1/half_2`，致其评测靶点从旧版 [`CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt) 的 `62` 个 [`K_conf = 2`] 降至 `16 < 50`，机械退为 `SECONDARY_ONLY`，最终 `K_conf = 1`）**；`D4`（Stage 2 排列划分细节对齐下载前官方代码）。
5. **Section 3: Rounds 9–10 Confirmatory Audits (`Hybrid`, `TMVS`, `FCT` Under `resp_top50`)**
   - 3.1 R9（`K_conf = 3`）：`Hybrid_ProjShared` 与 `TMVS` 未确认（`0/4`）；`R18` 与 `R14` 诊断结果。
   - 3.2 R10（`K_conf = 2`：`dld1`, `cd4t_rest`）：`resp_top50` 下 `H-C1`、`H-C3` 通过而 `H-C2` 未确认（`0/2`）。
6. **Section 4: Rounds 10-E & 11 — Mask Selection Bias, `K_conf = 1` Descriptive Outcomes, & Post-Hoc Provenance for Round 12**
   - 4.1 `R10-E` 双掩码消融（`resp_top50` vs `all_nontarget`，`8/8` `FCT − E1 > 0`）与 R11（`PREREG_D v2`, `K_conf = 1` `cd4t_stim48` 按 `P1` 仅作描述性报告；`Q7` 必报项：`hcasmc`、`telohaec N=9`、`cd4t_stim48` 共享 Marson lab/谱系/文库及共用通道、`feng_ipsc H-D5`）。
   - 4.2 **R9–R11 事后汇总（明确标注为 Round 12 两个论断的事后假设来源，非确认性证据）**：采自 `A-T3` [`REPORT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/audit/exploratory_check/REPORT.md) 的 12 线锁定口径复算——`r_orth_cal` 下 `FCT − Mis_s` 配对 `11/12` 过线（不配对 `7/12`，两 T 细胞线不配对失败）、`SA-8` `E1(split) − E1(all) > 0` 在 `9/12` 线成立、`N_5` `12/12`（并明示按构造为零对照）、及无噪声真实 $\mu$ 模型在 `all` 掩码预测错 **`8/12`** 线符号（因此机制未解）。
7. **Section 5: Round 12 Confirmatory Results — Cross-Modality Validation in `A549` CRISPRko (`PREREG_E`, `K_conf = 1`, 推广性有限)**
   - 5.1 门控筛选（`a549_ko` `N = 445`, `|G_eval| = 6,121`; `hek293_xucao` `D3` 后 `N = 16` 退为次要线，`K_conf = 1`）与**三方独立复算一致性**（`A-T2`, `A-T3`, 主控 `r12_check.py`：`90` 项最大差 `1.9e-16`；`177` 点估计 `1.1e-16`；`264` CI `8.3e-17`；`SA-Src8` `195` 点估计 `1.7e-16`；全部判定 `100%` 一致）。
   - 5.2 **确认性核心——论断 1（`H-E1` 靶线特异性确认）与「特异性 $\neq$ 精度」辨析**：
     - `a549_ko` `H-E1a` `+0.1898 [+0.1742, +0.2055]`；`H-E1b` 对全部 `6` 个错配源线配对差均 `97.5% CI > 0`（最小估计 `+0.0796 [+0.0693, +0.0900]` vs `nadig_jurkat`、最小下界 `+0.0682` vs `replogle_rpe1`，最大 `+0.1869` vs `nadig_hepg2`）；`H-E1c` 胜置换（`0.1742 > +0.0018`）。
     - **核心边界对比**：在完全同一指标（`a549_ko` `all_nontarget` `r_orth_cal`）下，`Source_Grand_Mean` **`+0.4051`**、`CIPHER_Linear` **`+0.2817`**、`Hybrid_on_TMVS` **`+0.2373`**、`TMVS` **`+0.2217`** 均显著高于 `FCT` **`+0.1898`**；证实 **`FCT` 不是最好的预测器，论断 1 只证明目标线自身未扰动基线带来真实跨系/跨模态特异信息，而非预测精度最高。**
   - 5.3 **非核心方法学观察——论断 2（`H-E2` / `SA-8` 确认，机制未知，局限于 `r_orth_cal`）**：
     - `a549_ko` 上 `E1(all_common) = -0.0201`，而在独立分半 `resp_top50_split` 下 `E1(split) = +0.2565 [+0.2195, +0.2981]`（`H-E2a`）、跨掩码差 `+0.2766 [+0.2412, +0.3164]`（`H-E2b`），相对同半 `resp_top50`（`+0.2866`）保留约 `90%` 抬高；全基因 `FCT − E1 = +0.2099` 清晰分离而在 `resp_top50` 下塌缩至 `+0.0207 [-0.0078, +0.0516]`。明确界定机制未知且不外推至 DEG 指标。
   - 5.4 **仅报告敏感情境与次要线（`PREREG_E §5`）**：`SA-Cas9`（`N = 445`，`FCT +0.1950`, `H-E2b +0.2735` 全部通过）、`SA-Src8`（`N = 528`, `8/8` 配对差均过，最小下界 `+0.0372` vs `dld1`, `FCT +0.1563`）、`hek293_xucao`（`N = 16`, `rho_AB_orth = 0.3005`, 论断 1 描述性 `6/6` 通过 `FCT +0.1561`，论断 2 因 `N=16` 低功效未过 `+0.0339 [-0.0868, +0.1797]`，不能据之否定 HEK293 上存在 `H-E2`）、不配对 `H-D2`（A549 通过 `+0.0504`；HEK293 不通过 `-0.0487`）、以及 `N_5` 报告表（明示 `Source_Grand_Mean`/`GlobalScale` 按构造精确为零 `0.0000`, `E1..E7` $R^2_{\mathcal{N}_5} \ge 0.9998$）。
8. **Section 6: 第十三轮占位（Round 13 Pre-Registration Placeholder）**
   - **`GSE327644` 原代人巨噬细胞（正处预注册阶段）**：本小节仅为占位，待预注册锁定、代码冻结、评测及三方核验完成后再填入结果；**当前不写任何预期结果或预先推断**。
9. **Section 7: Discussion & Prior-Art Boundaries (`PRIOR_ART_CONCLUSION_R13.md`)**
   - 7.1 为何特异信息仅占极小方差（$R^2 \approx 3.6\%$ on A549）而 `Source_Grand_Mean`（`+0.405`）显著高于 `FCT`（`+0.190`）？区分“学习全局源线跨扰动相对响应结构”与“靶线自身基线条件化”。
   - 7.2 与先占文献的精确划分：
     - **论断 1 vs `Kawsar et al. 2026 (CASTER)` & `Heidari et al. 2026`**：不声称首次；界定 `CASTER`（部分筛选神经网络单错配控制替换，零样本逊于 `PerturbMean`）及 `baseline_logfc`（单模态 OOD）与本文全源配对 IUT + 跨模态 `CRISPRi -> CRISPRko` 前瞻检验的重合与差异。
     - **论断 2 vs `Nicol et al. 2026 v2` & `Heidari et al. 2026`**：如实缩小新颖性，说明 `Nicol v2` 已提出 `Split DEG + control` 并报告 `nonctl-mean` 维持正分（`Adamson 0.732`, `ReplogleRPE1 0.492`；且在标准 `PearsonΔ` 下 `top-20` 均值 `0.259` 与全基因 `0.255` 相近，而 `ctl-mean` 归零）、`Heidari 2026` 已涵盖 `CrossSplit`/`*_vs_topdeg`/`define_trivial`（`trivial`/`non_trivial`/`non_significant`）/`rigorous_pearson_delta_pert`/`baseline_logfc`；本文补充的是在正交去均值残差 `r_orth_cal` 下不用目标基线的 `E1_ConstMed` 从近零（`-0.020`）到 `+0.257` 的跃升现象，定位为机制未知的方法学观察。
   - 7.3 跨谱系/跨模态 `K_conf = 1` 推广边界及 R10–R11 两 T 细胞线与 R12 `A549` 的对比反思。
10. **Section 8: Combined Limitations, Checklist & References**
11. **Appendix / Archive (§9): 第 9–11 轮事后探索基准存档（Hypothesis Provenance, Non-Confirmatory）**

---

## 4. 完整图表清单（Figure & Table Inventory）

### 4.1 主图清单（Main Figures 1–7）

| 编号 | 图题与核心展示内容 | 数据来源 | 性质与状态 |
| :--- | :--- | :--- | :--- |
| **Figure 1** | **Four-Round Pre-Registered Design, Attrition Funnel, & `r_orth_cal` Orthogonal Geometry** | R9 (`6`) + R10 (`4+1`) + R11 (`3+1`) + R12 (`A549`, `HEK293 D3`) | **Design Overview**：(A) R9 $\to$ R10 $\to$ R10-E/R11（假设生成）$\to$ R12 前瞻确认 $\to$ R13 占位流程；(B) 各轮筛选与 `HEK293` `D3` 退组路径；(C) `r_orth_cal` 消幅与配对对照几何。 |
| **Figure 2** | **Confounding Diagnostics Across Early Rounds: `R18` Scalar Scaling & `resp_top50` Inflation** | R9 + R10 (`BASELINE_ABLATION.md`) | **Confirmatory & Post-Hoc Diagnostic**：(A) `GlobalScale` 正交归零；(B) `resp_top50` vs `all_nontarget` 下 `E1`/`E2`/`E6`/`E7` 与 `FCT - E1` 分化（`8/8 > 0`）。 |
| **Figure 3** | **Rounds 9–10 Unconfirmed Baselines: `Hybrid_ProjShared`, `TMVS`, & R10 `FCT` (`resp_top50`)** | R9 (`K_conf=3`) + R10 (`K_conf=2`) | **Confirmatory**：R9 `0/4` 与 `R14`；R10 `H-C1`/`H-C3` 通过与 `H-C2` `0/2` 失败。 |
| **Figure 4** | **Round 11 Descriptive Whole-Transcriptome Outcomes (`K_conf = 1`, P1 Aborted) & 12-Line Provenance** | R11 (`cd4t_stim48`, `feng_ipsc`, `telohaec`) + `REPORT.md` | **Descriptive & Post-Hoc Provenance**：(A) R11 `P1` 仅作描述性报告（两 T 线不配对输给 `vcc2025_h1`）；(B) R9–R11 `12` 线配对 `FCT − Mis_s`（`11/12`）与 `SA-8`（`9/12`）事后模式（生成 R12 假设）；(C) 无噪声机制模型 `all` 掩码错 `8/12` 符号。 |
| **Figure 5** | **Round 12 Confirmatory Core (`A549` CRISPRko, `K_conf = 1`): Claim 1 (`H-E1`) Paired Specificity Across 6 Mismatched Sources** | `round12/results/` & `round12/audit/` (`a549_ko` 主图，含 `SA-Src8`/`SA-Cas9`/`HEK293`) | **Confirmatory Core (Single Line, `K_conf = 1`)**：(A) `a549_ko` `FCT` vs 6 错配基线与 `Perm100` 点估与 `97.5% CI`；(B) `H-E1b` `6/6` 配对差森林图（`min +0.0796`, 标明跨模态 `CRISPRi -> CRISPRko`）；(C) `HEK293`（`N=16` 次要）、`SA-Cas9`、`SA-Src8`（`8/8`）平行展示。 |
| **Figure 6** | **Specificity Is Not Best Accuracy: `FCT` (`+0.190`) vs Stronger Target-Agnostic Source Predictors (`Source_Grand_Mean +0.405`)** | `a549_ko` (`REPORT_holdout_r12_r9_6_log1p.md`) | **Mandatory Safeguard Comparison**：按同一指标 `r_orth_cal`（`all_nontarget`）横向排布 `Source_Grand_Mean (+0.405)`、`CIPHER_Linear (+0.282)`、`Hybrid_on_TMVS (+0.237)`、`TMVS (+0.222)`、`FCT (+0.190)`、`Mis_s`（`-0.003..+0.110`）与 `E1 (-0.020)`，直观证明论断 1 讲特异性而非最高精度。 |
| **Figure 7** | **Methodological Observation (Non-Core): Claim 2 (`H-E2` / `SA-8`) Persistent Constant-Baseline Elevation Under Split-Half Selection** | `a549_ko` (`+0.2565`, `H-E2b +0.2766`) + `hek293` (`N=16` 低功效) + `R11 12` 线背景 | **Confirmed Methodological Caveat (`r_orth_cal` only, Mechanism Unknown)**：(A) `E1` 在 `all_nontarget (-0.020)`、`resp_top50 (+0.287)`、`resp_top50_split (+0.257)` 阶梯跳升（`~90%` 保留）；(B) `FCT − E1` 在全基因（`+0.210`）与 `top50`（`+0.021` NS）对比；标注与 *Nicol v2* / *Heidari 2026* 边界。 |

### 4.2 主表与附表清单（Main Tables 1–9 & Supplementary Tables）

| 编号 | 表题与核心内容 | 性质与来源状态 |
| :--- | :--- | :--- |
| **Table 1** | **Unified Dataset Inventory Across Rounds 9, 10, 11, & 12**（`16` 项数据集：包含 R12 `a549_ko` `605,484` 细胞/`45,320` NTC/`N=445`/`\|G_eval\|=6,121` 与 `hek293_xucao` `98,274` 细胞/`2,757` NTC/`D3` 前 `62` vs 后 `16` 靶点、以及 R13 `GSE327644` 占位标注） | 数值完备（R9–R11 表 + `R12_QC_COUNTS.md` + `DEVIATIONS_R12.md`） |
| **Table 2** | **Round 9 & Round 10 Confirmatory Matrices (`PREREG_A v4`, `PREREG_C`)** | 数值完备 |
| **Table 3** | **Cross-Round `R19` Decomposition, `S6`, & `R10-E` Baseline Ablation (`10` Lines)** | 数值完备 |
| **Table 4** | **Round 11 Descriptive Outcomes (`PREREG_D v2`, `K_conf = 1`, P1 Aborted) & 12-Line Post-Hoc Provenance (`SA-8 9/12`, `Paired 11/12`, `Noiseless 8/12 Sign Errors`)** | 数值完备（`R11_VERDICT.md` + `round11/audit/exploratory_check/REPORT.md`） |
| **Table 5** | **Round 12 Confirmatory Matrix (`A549` CRISPRko, `K_conf = 1`, `r_orth_cal`)**: 逐项列出 `H-E1a`（`+0.1898`）、全部 6 项 `H-E1b` 配对差（`+0.0796..+0.1869`）、`H-E1c`（`Perm Q975 +0.0018`）、`H-E2a`（`+0.2565`）、`H-E2b`（`+0.2766`）及 `valid_frac = 1.000`, `rho_AB_orth = 0.789` | 数值完备（`R12_VERDICT.md` + `REPORT_holdout_r12_r9_6_log1p.md`） |
| **Table 6** | **Specificity vs Accuracy Benchmark on `A549` (`all_nontarget` & `resp_top50`)**: 并列展示 `Source_Grand_Mean` (`+0.4051`)、`CIPHER_Linear` (`+0.2817`)、`Hybrid_on_TMVS` (`+0.2373`)、`TMVS` (`+0.2217`)、`FCT` (`+0.1898`)、`Hybrid_ProjShared` (`+0.1554`)、`E5` (`+0.1022`)、`E3` (`+0.0794`) 与 `E1`/`E2`/`E6`/`E7` | 数值完备（`REPORT_holdout_r12_r9_6_log1p.md`） |
| **Table 7** | **Round 12 Non-Confirmatory Sensitivity & Secondary Line (`SA-Cas9`, `SA-Src8`, `HEK293 N=16`, Legacy Unpaired `H-D2`)** | 数值完备（`R12_VERDICT.md §4` + `src8` 报告） |
| **Table 8** | **Prior-Art Scope & Boundary Comparison Matrix (`PRIOR_ART_CONCLUSION_R13.md`)**: 对照 `Kawsar et al. 2026 (CASTER)`、`Nicol et al. 2026 v1/v2`、`Heidari et al. 2026 (scPerturb-eval)`、`Miller et al. 2026`、`Liang et al. 2026 (SafeDelta)` 各自已做内容、未能读取部分（bioRxiv 403）及与我方差异 | 完备（`PRIOR_ART_CONCLUSION_R13.md` + `PRIOR_ART_G.md v2`） |
| **Supp. Tables S1–S7** | `S1–S4`: R9–R11 敏感性与文献对比；`S5`: 四轮多方独立复算误差汇总（含 R12 主控/A-T2/A-T3 最大差 `1.9e-16`）；`S6`: R12 `N_5` `r_cone_sp` 完整表（逐行标明 `ZERO_BY_EXACT_CONSTRUCTION` 与 `SUPPRESSED_BY_CONE_CONSTRUCTION`）；`S7`: `a549_ko` `QARS` 与 `HEK293` MT 分布及 `SUPERSEDED_v1` vs `v2` 计数对照 | 完备（`REPORT_holdout_r12_r9_6_log1p.md` + `DEVIATIONS_R12.md`） |

---

## 5. 方法、第十二轮设计与预注册偏差细纲（Methods & Deviations）

### 5.1 四轮预注册脉络与第十二轮核心设计（`PREREG_E`）
- **核心假设来源（`PREREG_E §0.1 E6`）**：在 R9–R10 原主假设先后未确认、R11 因 `K_conf = 1` 中止确认后，研究人员于 R11 数据集（12 条线）上开展事后探索并提出两项待检验假设：**`H-E1`**（用配对 `FCT − Mis_s` 取代功效不足的不配对 `H-D2`）与 **`H-E2`**（`SA-8` 拆半选基因下 `E1_ConstMed` 仍明显为正）。第十二轮锁定 [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md)（`4c606f4`）在**完全未见**的 R12 候选数据上作盲态前瞻检验。
- **跨扰动模态与不对称准入规则（`PREREG_E §0.2, §2.2`）**：
  - 源池固定为 R9 `r9_6` `6` 条 `CRISPRi` 细胞系（`6,148` 共有基因）；主目标候选 **`a549_ko`**（Liu et al. 2026 `GSE345058`, A549 肺腺癌 `CRISPRko` CROP-seq dox-Cas9，`10.997 GB`，`32` 通道按 `batch1+3` vs `batch2+4` 嵌套划分两半与 `Lane` 级别 NTC 扣除）既是全新谱系又是全新扰动模态（`CRISPRko`）。设 `a549_ko` 为必需线（`>= 50` 靶点 & `>= 3,000` 基因），允许 **`K_conf >= 1`** 执行判定，但强制任何仅 `K_conf = 1` 的结果必须在题目与摘要中声明“**单线确认、推广性有限**”。
  - 条件候选 **`hek293_xucao`**（Xu et al. 2023 `XuCao2023.h5ad`, `Zenodo 13350497 / GSE218566`, `dCas9-KRAB-MeCP2 CRISPRi`, `0.337 GB`）：先过 `E3` 原始整数 UMI 与引导标签可行性检查（实测 **`PASS`**：`98,315` 单 sgRNA 细胞、`2,758` NTC across `15` `NO-TARGET_*`），再经 `P5`（须 `n_eval_targets >= 50` 且 `>= 3,000` 基因才入 `K_conf`，否则 `SECONDARY_ONLY`）。**严禁按敲低 `> 40%` 事后筛细胞/sgRNA。**
- **两项确认论断与多重校正（`PREREG_E §4`）**：`IUT + 2-Claim Bonferroni`（$\alpha_{\mathrm{claim}} = 0.025$，双侧 `97.5%` 联合扰动自助法区间 `[1.25%, 98.75%]`，`seed = 20261201`, `n_boot = 4,000`, `n_perm = 100`）。
  - **论断 1（`H-E1` 确认性核心，`all_nontarget`）**：`H-E1a`（`FCT` `97.5% CI_lo > 0`、`valid_frac >= 0.95`、`rho_ab_orth > 0`）$\cap$ **`H-E1b`（对全部 `6` 条 `r9_6` 错配基线逐一计算配对差 $\Delta^{\mathrm{pair}}(s') = r^{\mathrm{orth\_cal}}(\mathrm{FCT}) - r^{\mathrm{orth\_cal}}(\mathrm{Mis}_{s'})$，每条均要求 `97.5% CI_lo > 0`）** $\cap$ `H-E1c`（`FCT CI_lo > Perm100 Q_0.975`）。
  - **论断 2（`H-E2` / `SA-8` 非核心方法学观察）**：`H-E2a`（`resp_top50_split` 下 `E1_ConstMed` `97.5% CI_lo > 0`）$\cap$ `H-E2b`（公共有效扰动上联合交叉掩码差 $E1_{\mathrm{split}} - E1_{\mathrm{all}}$ `97.5% CI_lo > 0`）。明确 `E1_ConstMed` 包含源线平均倍数变化 $\overline{\mathrm{FC}}_{p,g}$、仅不用目标线对照基线谱 $\mu_{t,g}$。

### 5.2 第九至十二轮全部预注册偏差与 `HEK293` 退组过程（必须进方法与局限）
- **R9–R11 偏差回顾**：`R9 D1`（3D vs rank-1 投影）、`R10 D-R10-1/2` & `13:15` 公开更正、`R11 D1–D4`（`hcasmc` 0 guide 排除、`telohaec` 弱敲低 `N=9`、`P1` 中止 `K_conf = 1`）。
- **R12 偏差 `D1–D4`（[`DEVIATIONS_R12.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/DEVIATIONS_R12.md)，全部在任何数据上计算任何评测指标之前，仅依据计数与元数据决定）**：
  1. **`D1`（`11-23 01:20`，裁定 B）**：`HEK293` 主分析不设线粒体比例上限。`PREREG_E §3.2` 仅规定 `>= 300 UMI` 与 `>= 200` 基因、保留全部 QC 合格单 sgRNA 细胞；沿用 R11 默认 `20%` 上限会删掉 `76.7%` 细胞（`.obs` `percent_mito` 5/25/50/75/95% 为 `12.89% / 19.91% / 25.81% / 32.32% / 42.42%`），文本优先故关闭 MT 过滤，亦不新增 `MT <= 20%` 敏感性分支。
  2. **`D2`（`11-23 01:20`，裁定 A）**：读取代码修正 `AnnData` `attrs['_index']`（避免 `HEK293` 索引列与条形码排序丢失）及 `A549` `indptr` `int32` 溢出（`3.60e9 > 2^31 - 1` 转 `np.int64`），不改算法。另记 `A549` `QARS`（`591` 细胞）因冻结 HGNC 表缺少映射而排除，不影响 `r9_6`、`G_eval` 与 `445` 评测靶点。
  3. **`D3`（`11-23 04:00`，裁定 2 —— 影响 `K_conf` 归属的关键修正）**：`PREREG_E §3.2–§3.3` 原文规定 `HEK293` 无 `.obs` 重复列时进入 Stage 2 并按**半（`replicate/half for HEK293`）**作为 NTC 扣除单位。最初共享代码误用单一 `batch_0`（全库 `2,757` NTC 统一扣除），使 `HEK293` 得到 `96/200` 可重复、`62` 个评测靶点（`>= 50`）并曾一度将 `HEK293` 写入 **`CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt`**（`K_conf = 2`）。经 `A-T3` 独立读取核对发现不符后，在计算任何预测或相关指标之前，将 Stage 2 扣除单位改回与原文一致的 `half_1` / `half_2` 按半扣除；官方修正版与 `A-T3` 独立重建均证实 `HEK293` 降为 **`19/200` 可重复、`16` 个评测靶点（`|G_eval| = 6,126`，`< 50`）**，遂将 `v1` 归档为 `SUPERSEDED_v1`，按 `§2.2` 机械改列 `HEK293` 为 **`SECONDARY_ONLY`**，使 `K_conf` 从 `2` 降至 **`1`**（[`CONFIRMATORY_LINES_R12.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/CONFIRMATORY_LINES_R12.txt)）。
  4. **`D4`（`11-23 04:00`，裁定 2 §2）**：针对 Stage 2 `seed=42` 分层 50/50 划分未写明的细节（单一生成器 vs 每层重建、键排序 `Python` 字符串顺序、`(cell_barcode, row_idx)` 排序、`perm[:n//2]` 归 A 半），锁定采用下载前已写好的官方代码规则，不看任何指标。

---

## 6. 局限性汇总（Combined Limitations: Rounds 9–12 & Pre-Submission Audit）

1. **单线确认性与跨谱系推广边界（`K_conf = 1` + `D1–D4` + `SUPERSEDED_v1` 退出）**：第十二轮虽然在 `A549` CRISPRko（`N = 445`）上同时确认了论断 1 和论断 2，但由于偏离 **`D3`** 修正 `HEK293` Stage 2 NTC 扣除单位（从错误单一 `batch_0` `62` 靶点 [`CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/CONFIRMATORY_LINES_R12_SUPERSEDED_v1.txt) 依预注册原文改为按半扣除后降至 `16 < 50` 退出确认线），全轮确认仅依赖 **唯一一条确认性细胞系（`a549_ko`，`K_conf = 1`）**，跨更多原代与稳态细胞系的推广性有限。此外，此前 R10 `cd4t_rest` 与 R11 `cd4t_stim48` 两条 T 细胞线在全基因不配对口径均败于 `vcc2025_h1`（且与第十一轮 `K_conf = 1` 中止并存），跨正常/原代谱系的普及性仍待后续补充。同时 `HEK293` 受 `~2,155` 浅测序与无 MT 上限（`D1`，中位 `25.81%`）限制，$\rho^\perp_{AB} = 0.30$、仅有 `16` 个靶点，`H-E2` 极低功效（不能据之否定该现象）。
2. **`FCT` 不是最好的预测器且总体可解释方差很小（特异性 $\neq$ 精度）**：即使在第十二轮取得确认的 `A549` 上，`FCT` 的 `r_orth_cal` 也仅为 **`+0.190`**（$R^2 \approx 3.6\%$），明显低于无需目标线自身基因基线信息的 **`Source_Grand_Mean`（`+0.405`）**、`CIPHER_Linear`（`+0.282`）、`Hybrid_on_TMVS`（`+0.237`）与 `TMVS`（`+0.222`）。`H-E1` 仅证明目标线自身对照表达带有相对于错配对照的特异信号，**不代表 `FCT` 是高精度或实用最优模型**。
3. **假设事后 provenance 与单数据集外部有效性限制（`E6`）**：第十二轮的两个论断（配对 `FCT − Mis_s` 与 `SA-8` 偏倚）均产生于看过第 9–11 轮数据之后的事后归纳；因此第 9–11 轮内部的 `11/12` 配对表现或 `9/12` `SA-8` 均属事后发现，真正意义上的盲态前瞻确认迄今只有 `A549` CRISPRko 这 `1` 个独立数据集。
4. **论断 2 机制未知、作用域严格受限及机理模型在 `all` 掩码上 `8/12` 线符号错误**：论断 2 仅作为**方法学观察**（不作核心发现），其适用范围严格限于本研究的正交去共享均值残差空间 **`r_orth_cal`**，不能外推至 Nicol 等人分析的 DEG 类指标。底层机制至今未知：第十一轮用真实 $\mu$ 分布进行的无噪声模拟虽然在分半 `top-50` 预测全为正，但在全基因 **`all` 掩码上同样全预测为正（`12/12`），从而预测错了实测中 `8/12` 条细胞系的负号**（如 `cd4t_stim48` 模拟 `+0.6084` vs 实测 `-0.0203`、`replogle_rpe1` 模拟 `+0.1738` vs 实测 `-0.0934`，幅度偏大 `3–15` 倍）。此外，即便在 `A549` 上，`SA-8` 分半也仅消除了约 `10%` 的虚高，而约 `90%`（`0.277 / 0.307`）仍存在且原因未明。
5. **先占文献检索的平台壁垒（`HTTP 403` 投稿前必补人工 PDF 核读）**：按 [`PRIOR_ART_CONCLUSION_R13.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/status/PRIOR_ART_CONCLUSION_R13.md)，**Nicol et al. (2026) v2** 与 **Heidari et al. (2026)** 的 `bioRxiv` 全文均受网站机器人验证拦截（`HTTP 403`，主控与 `A-T5` 均不绕过验证），当前仅核实了二者的公开代码仓、notebook 输出、v1 `PMC13192888.1` 全文及摘要元数据。**投稿前必须有人手动下载并通读这两篇 PDF 正文**，确认正文中未另行讨论「`E1` 式基线在 `top-k` 与全基因之间的跳升」。
6. **跨线表达形式对比的事后局限与非独立工程边界**：在 R11 探索核查中，经 `A-T3` 按线分块置换并 `Bonferroni×3` 校正后，**`FCT` 优于直接 `log1p` 迁移，只在 `expr05` 子集上稳健**（`10/10` 线，`p = 0.0059`；`all` 子集 `p = 1.0000`），而 `Q5a` 表达形式守恒论断不成立也不写入正文。此外仍存在模型同源性（主代码 `A-T2` 与审计 `A-T3` 同基于 Gemini，虽 R12 已由 Anthropic `Claude` 主控独立复算全部关键指标最大差 $\le 1.9\times 10^{-16}$）、外部 VCC 斜率公式未公开与测试集未开放等限制。

---

## 7. 待办与合稿推进清单（Manuscript Assembly Checklist）

- [x] **历史档保留**：现版 R11 大纲已另存为 [`OUTLINE_R11.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/paper/OUTLINE_R11.md)。
- [x] **Round 9–12 全部锁定评测、三方复核与 R13 文献定论落盘**：`R10_VERDICT.md`、`R11_VERDICT.md`、`round11/audit/exploratory_check/REPORT.md`、`R12_VERDICT.md`、`DEVIATIONS_R12.md`（`D1–D4` & `SUPERSEDED_v1`）、`PRIOR_ART_CONCLUSION_R13.md`、`PRIOR_ART_G.md v2` 均已纳入本大纲。
- [ ] **投稿前必做人工核读（`PRIOR_ART_CONCLUSION_R13.md`）**：由人类手动下载并通读 **Nicol et al. (2026) v2** 与 **Heidari et al. (2026)** 的 `bioRxiv` PDF 全文，确认无遗漏对 `top-k` vs 全基因跳升的讨论。
- [ ] **Round 13（`GSE327644` 原代人巨噬细胞）**：待预注册锁定、计算与三方复核全部完成后再填入 §6 占位小节。
- [ ] **全稿合规与发布审批**：保持人类角色称谓仅用 `Human project lead` / `人类项目负责人`；**不生成未经用户逐次同意的新闻稿、推文或任何对外文字**；不改动 `README`；定稿前全面核对回源 CSV/JSON 与冻结哈希（R10 `22/22`, R11 `41/41`, R12 `46/46`）。

---

## 8. 核心参考文献目录（已验证引文）

- **Round 12–13 关键方法与先占必引文献**：
  - **Kawsar, Dip, Paul, Noman, Zhang & Moni (2026)**, *Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line* (`CanPert / CASTER`), `doi:10.64898/2026.11.15.766815`（论断 1 必引；`t1` 神经网络单邻系错配对照先例）。
  - **Nicol, Shivakumar & Irizarry (2026)** v1/v2 (`doi:10.64898/2026.05.07.723486`, `PMC13192888.1` + `github.com/phillipnicol/systema`; v2 bioRxiv 正文 `HTTP 403` 待人工 PDF 核读)（论断 2 必引；`Split DEG + control` 与 `nonctl-mean` 持续正分）。
  - **Heidari, Karimpour, Srivatsa & Montazeri (2026)**, *Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward* (`scPerturb-eval`), `doi:10.64898/2026.02.14.705879`, `github.com/LCBB-UT/scPerturb-eval`（v1 bioRxiv 全文 `HTTP 403` 待人工 PDF 核读）（论断 1/2 必引；`CrossSplit`, `*_vs_topdeg`, `define_trivial`, `rigorous_pearson_delta_pert`, `baseline_logfc`）。
- **Round 12 数据集文献**：
  - **Liu et al. (2026)**, `A549` dox-Cas9 CROP-seq `CRISPRko` (`GSE345058`, `doi:10.64898/2026.06.01.728087`)（仅作数据来源引用）。
  - **Xu, Cao et al. (2023)**, `HEK293` PerturbSci-Kinetics (`GSE218566` / `Zenodo 13350497` `XuCao2023.h5ad`, `doi:10.1038/s41587-023-01948-9`)。
- **Round 9 文献**：Adduri et al. (2026) *Cell* (`10.1016/j.cell.2026.08.004`); Arc Institute VCC 2026 Wrap-Up (`2026-11-19`) & `cell-eval2` (`v0.22.0`); Dhaka et al. (2026) (`10.64898/2026.04.23.720505`); Huang et al. (2026) `CGC` (`10.64898/2026.09.30.755756`); Kuznets-Speck et al. (2025/2026) `CIPHER v2` (`10.1101/2025.06.27.661814`); Miller et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03307-w`); Molina & Zhang (2026) (`10.64898/2026.07.24.740459`); Nadig et al. (2025) *Nat. Genet.* (`10.1038/s41588-025-02169-3`, `GSE264667`); Nourreddine et al. (2026) *Nat. Biotechnol.* (`10.1038/s41587-026-03199-w`); Replogle et al. (2022) *Cell* (`10.1016/j.cell.2022.05.013`); Shoeibi & Yousefi (2026) `arXiv:2608.00152v3`; Viñas Torné et al. (2025) `Systema` *Nat. Biotechnol.* (`10.1038/s41587-025-02777-8`).
- **Round 10–11 数据集与领域对照文献**：Yeung et al. (2026) `DLD-1` (`10.64898/2026.07.10.737863`, `GSE337988`); Zhu et al. (2026) `CD4+ T` *Cell* (`10.1016/j.cell.2026.08.002`, `GSE314342`); Tian et al. (2021) (`10.1038/s41593-021-00862-0`); Sunshine et al. (2023) (`10.1038/s41467-023-41788-4`); Schnitzler et al. (2024) (`GSE210681`); Feng et al. (2026) `TargetedScreen` (`PMID:41330380`); `IGVFDS6248NDUW`; Liang et al. (2026) `SafeDelta` (`10.64898/2026.10.05.756906`); Huang (2026) (`arXiv:2610.32962v2`); Huang et al. (2026) `AmbiModBench` (`arXiv:2609.32527v1`); Gao (2026) (`arXiv:2611.10658v1`); Vollenweider & Bühlmann (2026) `SBB` (`10.64898/2026.04.20.719650`); Schäfer et al. (2026) `scPertEval` (`10.64898/2026.07.23.740433`); Panagopoulos et al. (2026) (`arXiv:2610.26267v1`); Booeshaghi & Pachter (2021) (`10.1093/bioinformatics/btab085`); Ahlmann-Eltze & Huber (2023) (`10.1038/s41592-023-01814-1`); Townes et al. (2019) (`10.1186/s13059-019-1861-6`).

---

## 9. 第 9–11 轮事后探索归档（Non-Confirmatory Provenance for Round 12）

> **口径与隔离说明（依据 [`EXPLORATORY_REVIEW_R11.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/status/EXPLORATORY_REVIEW_R11.md) 第 D 节）**：以下全部 R9–R11 数字均取自 `A-T3` 按锁定 `P5` 口径独立复算的 `round11/audit/exploratory_check/REPORT.md`；它们在 R9–R11 内仅为**事后探索与第十二轮假设来源**，直到第十二轮 `PREREG_E` 才获得针对 `A549` 单线（`K_conf = 1`）的盲态前瞻确认。

1. **SA-8 现象（R9–R11 `9/12` 线）**：`E1(split) − E1(all)` 点估计 `12/12` 为正、`97.5% CI > 0` 在 `9/12` 线成立（`cd4t_stim48` `+0.0730 [+0.0384, +0.1098]`、`cd4t_rest` `+0.0383 [-0.0061, +0.0864]` 不显著），直接催生第十二轮 `H-E2`。
2. **`FCT` 胜过任意常数基线组合（`N_5` 下 `12/12`）**：R11 与 R12 报告表均写明 `GlobalScale`（$\eta < 10^{-14}$）、`b_S`（$\eta < 10^{-14}$）及 `E1`/`E2`/`E6`/`E7`（$R^2_{\mathcal{N}_5} \ge 0.9928$）在 `N_5` 下为零或近零**是按数学构造如此**，仅作事后/非确认分析工具（R12 `a549_ko` `FCT N_5 = +0.1992`）。
3. **配对比较 `FCT − Mis_s`（R9–R11 `11/12` 事后 vs R12 `A549` 前瞻确认）**：R9–R11 未预注册事后配对 `FCT − Mis_s` 在 `N_5`（及 `r_orth_cal`）下 `11/12` 通过（而严苛不配对口径仅 `7/12`，两条 T 细胞线均未过），直接促成第十二轮把 `r_orth_cal` 配对检验 `H-E1b` 写入 `PREREG_E` 并于 `A549` 获得 `6/6` 前瞻确认。

- **机制局限（延续至 R12）**：R11 用真实 $\mu$ 的无噪声模拟在 `all` 掩码上 `12/12` 预测为正，预测错了 `8/12` 条线真实观测到的负号（如 `cd4t_stim48` `+0.6084` vs `-0.0203`），所以底层机制仍未解释。
