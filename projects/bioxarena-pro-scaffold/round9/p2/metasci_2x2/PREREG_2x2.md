# PREREG_2x2：同厂商与跨厂商科研代码审计 $2 \times 2$ 对照实验预注册方案

- **版本与状态**：Pre-Registration Draft v1.0（设计冻结阶段；尚未调用任何付费 LLM API，尚未在任何留出基准上执行受测审计）
- **撰写角色**：AI-Biology-Lab Round 9 Phase 2 智能体 `A1`（元科学·2×2 审计设计）
- **人类 PI / 授权人**：Zijian Zeng ｜ **主控**：Claude
- **对应背景文献**：Dong & Chen (2026, `arXiv:2608.28631v1`, CrossAudit) 与 Dong (2026, `arXiv:2611.23250v1`)；关联参考 Song (2026, `arXiv:2610.01471v2`)、Baumgärtner & Gurevych (2026, `arXiv:2601.12910v3`, SciCoQA)、Vaccaro (2026, `arXiv:2606.11217v1`)
- **配套冻结工件（均位于 `round9/p2/metasci_2x2/`）**：
  - [LITERATURE_AUDIT.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/LITERATURE_AUDIT.md)：`2608.28631` 与 `2611.23250` 逐章已做/未做核查清单
  - [CONSTITUTION_2x2.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/CONSTITUTION_2x2.md)：冻结版双分级科研代码审计规则手册
  - [audit_response_schema.json](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/audit_response_schema.json)：强制连续风险分与规则引用 JSON Schema
  - [matched_fpr_eval.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/matched_fpr_eval.py)：锁定评估器（Neyman-Pearson 匹配误报率、$\text{pAUROC}_{[0,0.20]}$、阈外校准检验、含阈值重估的任务簇自助法）
  - [power_and_cost_sim.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/power_and_cost_sim.py)、[complementarity_power.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/complementarity_power.py)、[post_repair_power.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/post_repair_power.py) 及 [outputs/](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/outputs/)：先验功效、I 类错误、覆盖率与调用成本仿真

---

## 1. 科学问题与从先驱文献继承的精确空白

`Dong & Chen (2026, arXiv:2608.28631 §3.2, §4.4)`、`Dong (2026, arXiv:2611.23250 §D.6–D.7, §F)` 与 `Song (2026, arXiv:2610.01471 §4.1, Limitations)` 已在单生成器环境下发现同/跨厂商路由存在检测差异，但同时明确指出以下四个核心混杂尚未解决（详见 `LITERATURE_AUDIT.md`）：

1. **单生成器混杂（Generator-Auditor Confounding）**：`Dong (2026b)` 与 `Song (2026)` 均仅使用单一生成器（分别为 `claude-haiku-4-5-20251001` 与 `Claude Opus 4.6`），而 `CrossAudit (§4.4)` 使用的是确定性 Python 模板生成器。单一生成器下，**「同厂商 vs 跨厂商」与「审计器自身能力/严格度主效应」完全不可分**。
2. **误报率 / 操作点不匹配（No Matched FPR）**：`CrossAudit (§4.4)` 中 `gpt-5.1` 对干净增量误拦截率为 `10/10`（100%）而 Claude 为 `0/10`；`Dong (2026b §B.2, Fig. 3)` 中 `self-strong`（`claude-sonnet-4-6`）在 4 条离散裁决规则下从 `3.3%` FPR 直接跳至 `34.0%` FPR，在跨厂商审计器所在的 `16.0%` FPR 区间完全空白。
3. **采样参数与调用工具链混杂（Sampling & Harness Confounding）**：`Dong (2026b §A.2, §D.6 Limit 1)` 同厂商传 `temperature=0`、跨厂商不传采样参数；科学实验中同厂商改走本地 CLI；`Song (2026)` 混用 Codex CLI、Gemini API 与 Claude Code。
4. **科学多文件代码与修后外部结局缺失**：`Dong (2026b §4.1–§4.3)` 的科研代码任务局限于单步函数且未测量修复后的科研测试通过率；`SciCoQA (arXiv:2601.12910 §7)` 虽涵盖多文件论文–代码对齐，但其 `543` 例合成差异由 `GPT-5` 生成（作者承认使 GPT-5 家族在合成集上占优），且未测误报率。

本预注册实验旨在通过**对称 $2 \times 2$（生成器厂商 $\times$ 审计器厂商）因子交叉 + 双阶段洁净/植入缺陷构建 + 确定性匹配误报率（Matched-FPR）比较 + 同构无工具采样控制 + 修后确定性验证结局**，在预先锁定所有判据的前提下检验跨厂商审计假说。

---

## 2. 预注册假说与双向可发表裁决规则（Hypotheses & Decision Criteria）

设生成器厂商 $G \in \{V_1, V_2\}$、审计器厂商 $A \in \{V_1, V_2\}$（核心 2×2 为 `Anthropic Claude` 与 `OpenAI GPT`；可选扩展至 3×3 / 4×4 纳入 `Google Gemini` 与 `DeepSeek`）。定义**对称跨厂商交互对比量（Cross-vs-Same Interaction Estimator）**：
$$\Delta_{\text{cross}-\text{same}}(Y) = \frac{1}{2}\Big[\big(Y_{1,2} - Y_{1,1}\big) + \big(Y_{2,1} - Y_{2,2}\big)\Big],$$
同时定义两个方向上的**简单效应（Simple Effects）**：$S_1 = Y_{1,2} - Y_{1,1}$（在 $G=V_1$ 上的跨减同）与 $S_2 = Y_{2,1} - Y_{2,2}$（在 $G=V_2$ 上的跨减同）。

### 2.1 主假说与次级假说

- **H1（匹配误报率下的对称跨厂商检出优势 — Co-Primary Endpoints）**：
  - **H1a**：在 $\text{FPR}=10\%$ 精确匹配下，跨厂商审计比同厂商审计具有更高的真阳性检出率 $\Delta_{\text{cross}-\text{same}}(\text{TPR@10\%FPR}) > 0$。
  - **H1b**：在实用低误报区间 $\text{FPR} \in [0, 0.20]$ 的归一化部分 ROC 曲线下面积上，$\Delta_{\text{cross}-\text{same}}(\text{pAUROC}_{[0,0.20]}) > 0$。
- **H2（定位与机制真实性 — Secondary Endpoint E3）**：
  - 在 $\text{FPR}=10\%$ 匹配门槛下，不仅看整文件是否拦下（"flagging"），更用确定性 AST/文件/符号匹配检查是否准确命中植入缺陷所在文件与函数/行区，检验 $\Delta_{\text{cross}-\text{same}}(\text{TLR@10\%FPR}) > 0$（回应 `Dong 2026b §D.4` *"Flagging is not naming"*）。
- **H3（闭环一轮修复后的外部确定性科研通过率 — Co-Primary/Key Endpoint E4）**：
  - 将各 $(G, A)$ 单元在 $\text{FPR}=10\%$ 门槛（或实际 5 折交叉验证校准阈值 $\hat{\tau}^{\text{cv}}_{g,a}$）下被拦截的样本交回原生成器 $G$ 进行 **1 轮受控修复**（遵守防作弊修复守卫），在**外部确定性测试集**上测量净通过率增量 $\Delta_{\text{cross}-\text{same}}(\text{Pass}_{\text{post-repair}}) > 0$（同时追踪干净样本误拦截修复后的回归率 $\text{Regress}_{\text{clean}}$）。
- **H4（等预算 2-Call 盲区正交性 — Secondary Endpoint E5）**：
  - 在严格匹配单次调用边际 $\text{FPR}=10\%$ 的前提下（消除 `Song 2026 §5.3` 检出率差异混杂），跨厂商配对的缺陷重叠度低于同厂商两次独立阅读：$J_{\text{within}}(A_{\text{same},1}, A_{\text{same},2}) - J_{\text{cross}}(A_{\text{same}}, A_{\text{cross}}) > 0$；且在**总组合 $\text{FPR}=10\%$ 匹配**下，2 次调用的融合（均值评分融合 $\frac{1}{2}(s_1 + s_2)$）对比 `1 同 + 1 跨` 与 `2 同` 的检出率差。
- **H5（模型能力非对称 vs 真正跨厂商对称效应的定性交互检验）**：
  - 分解审计器主效应 $\hat{\beta}_{\text{aud}} = \frac{1}{2}[(Y_{1,2}-Y_{1,1}) - (Y_{2,1}-Y_{2,2})]$ 与交互效应 $\Delta_{\text{cross}-\text{same}}$。

### 2.2 双向裁决与「不利/零结果」发表承诺（遵循 P2_RULES #8）

对两个 Co-Primary 指标（$\text{TPR@10\%FPR}$ 与 $\text{pAUROC}_{[0,0.20]}$）以及 $E_4$（修后通过率）在 **Holm-Bonferroni 校正**后评估双侧 95% 阈值重估任务簇自助法置信区间（4,000 次重采样）及对应双侧自助法 $p$ 值：
1. **支持对称跨厂商盲区假说**：至少一项 Co-Primary 指标在 Holm 校正后 $p < 0.05$（且 95% CI 严格 $> 0$），**并且**两个生成方向上的简单效应 $S_1, S_2$ 不呈显著反向冲突（即非单一审计器在两个生成器上绝对碾压导致的非对称假象）。
2. **推翻对称假设 / 证实能力主效应主导**：若 $\Delta_{\text{cross}-\text{same}}$ 的 95% CI 包含 0（或未能通过 Holm 校正），而 $\hat{\beta}_{\text{aud}}$ 显著非零且 $S_1 \cdot S_2 < 0$（同一家强审计器在两组生成器上同时胜出），**必须如实报告「跨厂商优势未确认」**，将前人观察归因于「模型能力/操作点非对称」而非同源盲区。
3. **检出–修复解耦（Goodhart / 修复回归）**：若审计检出率显著提升但 $E_4$ 修后通过率增量 $\le 0$ 或干净回归率抵消收益，如实报告**「审计报警不自动转化为科学通过率提升」**，禁止临场更换指标或事后过滤任务补救。

---

## 3. 2×2 实验架构与防混杂执行双轨（Execution Tracks）

### 3.1 因果识别矩阵

| 生成器 $G$ \ 审计器 $A$ | 审计器 $V_1$（如 Anthropic Claude） | 审计器 $V_2$（如 OpenAI GPT） | 行内比较 |
|---|---|---|---|
| **生成器 $V_1$（Claude）** | $(V_1, V_1)$ 同厂商 ($Y_{11}$) | $(V_1, V_2)$ 跨厂商 ($Y_{12}$) | $S_1 = Y_{12} - Y_{11}$ |
| **生成器 $V_2$（GPT）** | $(V_2, V_1)$ 跨厂商 ($Y_{21}$) | $(V_2, V_2)$ 同厂商 ($Y_{22}$) | $S_2 = Y_{21} - Y_{22}$ |
| **交互对比 $\Delta_{\text{cross}-\text{same}}$** | — | — | $\frac{1}{2}(S_1 + S_2)$ |

> **关键代数性质**：若数据生成模型为 $Y_{ijm} = \mu + \alpha_i^{\text{gen}} + \beta_j^{\text{aud}} + \gamma \mathbb{I}(i \neq j) + u_m + \varepsilon_{ijm}$，则 $(S_1 + S_2)/2 = \gamma + \bar{\varepsilon}_{\text{int}}$ **精确消去**生成器难易主效应 $(\alpha_2 - \alpha_1)$ 与审计器能力/严格度主效应 $(\beta_2 - \beta_1)$。

### 3.2 执行双轨与混杂控制

我们在审阅本机 `/Users/mac/AgentHub/src/adapters/`（`claude.ts`、`codex.ts`、`deepseek.ts`、`antigravity.ts`）源码后确认：AgentHub 原生支持 `newThread` 上下文隔离，但**其适配器未透传 `temperature` 参数**，且默认注入 `HUB_ROLE_PROMPT` 与 `mcp__agenthub` 工具，同时分别经由四种不同 CLI/App-server 工具链调用。若混用会重现 `Dong (2026b §D.6)` 的工具链与采样混杂。因此本方案设定严格双轨：

- **Track A（主确认轨 — 同构直连 API，零工具链/采样混杂）**：
  - **统一执行器**：单一 Python 轻量调度器直接调用各厂商官方 HTTPS API（或兼容端点），**无外部文件读取工具、无代码执行工具、无 MCP 工具**（遵循 `CrossAudit §4.1` toolless defense）。
  - **严格三重隔离（I1）**：
    1. *参数隔离*：不同厂商独立权重（冻结快照 ID 记录于运行清单）；
    2. *上下文隔离*：每次读取为全新无历史对话，审计器仅接收围栏化的任务规范文本 + 候选代码文件，绝不接收生成器对话轨迹、推理链或生成厂商身份标识；
    3. *标识盲化*：提示词与代码注释剥离任何模型/厂商名称（如 `Generated by ...`）及随机文件名指纹。
  - **匹配采样超参**：所有 $(G, A)$ 单元统一设定 `temperature = 0.3`（或 $0.2$，在 Pilot 前一次性冻结）、`top_p = 1.0`、`max_output_tokens = 4096`、推理深度（如 reasoning effort）固定为各旗舰模型等价中/高档，每件样本独立读取 $K = 4$ 次（主推荐档）。
- **Track B（可选现实对照轨 — 本机 AgentHub CLI 订阅模式）**：
  - 仅在 API 预算受限或专门评估现实 CLI 工具链效应时启用（使用人类已有的 Claude Code / Codex / Antigravity / DeepSeek 订阅）。
  - **强制前置改造**：必须禁用 `mcp__agenthub`（防止审计器跨线程偷看生成器对话）、在只读临时目录运行以隔绝工作区泄漏、每次调用强制 `newThread: true`，并在报告中**显式标注「CLI 工具链与采样参数与厂商身份绑定」**，不与 Track A 混合回归。

---

## 4. 任务集、确定性真值与「匹配孪生」人工植入错误协议

### 4.1 公开科研代码任务底库（均已核实许可证与数据量，严禁 `xatlas/Orion`）

为同时覆盖**可执行数值科研代码**、**生物医学建模代码**与**多文件论文–代码对齐**，且彻底规避 `SciCoQA` 合成集由 `GPT-5` 生成带来的数据污染（`arXiv:2601.12910 §7`），任务池由三部分组成（均不触碰任何禁区数据）：

1. **主库 S1（可执行科研计算任务 — `SciCode-Verified`, `Apache-2.0`）**：
   - 来源：`huggingface.co/datasets/shhu2001/SciCode-Verified`（`arXiv:2608.04975`）。经实测 `data/problems_test.jsonl` 为 `951,750 bytes`（`64` 个独立科学问题，`290` 个子步骤，覆盖物理、化学、生物、材料、数学与统计），`test_data_cleaned.h5` 为 `1.03 GiB`（远低于全队 `60 GB` 总下载上限）。每题含经过专家双审修复的规范描述、函数头、标准解答与确定性 `np.allclose` 数值单元测试。
2. **主库 S2（生物医学建模 CPU 任务 — `BioXArena-Pro-Public` 非 `xatlas` 子集, `CC BY-NC 4.0` + `MIT grade.py`）**：
   - 来源：`huggingface.co/datasets/mbzuai-ai4bio/BioXArena-Pro-Public`（`arXiv:2611.13994`）中 `chemical-biology`（8 任务）与 `phenotype-disease`（8 任务），共 `16` 个独立研究任务（下载总计约 `587 MB`，**严格排除 `perturbation/xatlas-orion-perturb`**）。利用本地生物单元切分验证集与官方 MIT `grade.py` 提供确定性评分（不向 Codabench 提交）。
   - **合并独立问题簇**：S1（`64` 题）+ S2（`16` 题）天然构成 **$M = 80$ 个完全独立科研问题簇**。若启用高功效档（$M = 120$），再从 **S3（`SciCoQA`, `CC BY 4.0`, `3.49 MiB`）** 的真实/非 CS 学术仓库池（Physics 51、Statistics 40、EE 114、Q-Bio/Math 等，**只取干净原始代码 `original_code` 与论文章节，不使用其 GPT-5 生成的 `changed_code`**）中抽取 `40` 个独立开源仓库，使总独立任务簇达 **$M = 120$**。
3. **补充静态对照层（`SciCoQA` Real Split, 92 例真实差异, 72 篇论文, `CC BY 4.0`）**：
   - 作为非 LLM 生成的第三方人类科研真实误差基准（GitHub issues 42 + 复现论文 50），用于验证审计器在多文件长上下文真实科研差异下的文件定位率（匹配 `relevant_code_files`）与未标定基线对比（不计入 $2\times 2$ 交互主推断）。

### 4.2 两阶段构造：「对称洁净基线 + 第三方确定性植入错误」

若仅收集 $G_1$ 与 $G_2$ 各自自然写错的代码，会导致 $G_1$ 的天然 bug 与 $G_2$ 的天然 bug **难度、类型完全不匹配**，无法检验 $2\times 2$ 交互。因此本方案采用**匹配孪生植入（Matched-Twin Planted Defect）**设计：

1. **阶段一：对称洁净基线生成与双门控验证（Clean Baseline Generation & Gate）**
   - $G_1$ 与 $G_2$ 在相同任务规范下各自独立生成完整可运行解答（最多允许 $R_{\text{gen}} = 2$ 次尝试，即带 1 次基于可见运行报错的自我修正以保证较高任务通过率；每题需生成 $C = 2$ 个独立洁净实现或互不重叠的子步骤/变体以降低 $\text{FPR}=10\%$ 阈值估计方差）。
   - **进入匹配池准入条件**：仅当某一任务 $m$ 上 **$G_1$ 与 $G_2$ 双方均成功产出通过全部确定性测试的洁净实现**，方将该任务纳入配对分析池（不对称仅单方通过的任务归入「单边探索性副集」而不混入主 $2\times 2$ 对比，防止难度选择偏倚）。
   - **洁净样本自身干净度门控（Clean-Standard Gate, 继承 `Dong 2026b §5 / §D.5` 教训）**：所有标为 `clean`（$D=0$）的样本必须额外通过 Tier-0 `DCL` 检查（语法、依赖、输入输出契约、无隐藏异常吞没、规范与代码签名完全一致），确保误报率测的是**审计器的误报**而非实验构造缺陷。

2. **阶段二：供应商中立的确定性错误植入（Vendor-Neutral Planted Defect Injection）**
   - **严禁使用 $G_1$ 或 $G_2$（或任何待测同家族模型）生成植入错误**，以免重演 `SciCoQA` 中 GPT-5 偏爱 GPT-5 植入缺陷的偏差！
   - 错误由**确定性 Python AST 变换脚本 + 人类 PI/第三方脚本核定模板库**植入每个洁净基线 $C_{m,g}$，覆盖 `CONSTITUTION_2x2.md` 中的 **5 类科研实质缺陷**（每任务选取 $D = 2$ 个来自不同缺陷家族的植入变体，全部 5 类在 $M$ 个任务间分层均衡轮转）：
     1. `ALGO_MATH_LOGIC` (`CA-SCI-001`)：如 L2 改 L1 范数、梯度/势能符号反转、积分边界差一、归一化分母缺失；
     2. `STAT_LEAK_METRIC` (`CA-STAT-001`)：如全集 fit scaler/PCA 后再切分训练测试集、特征选择泄漏目标标签、指标 `max(0, r)` 截负、跨组泄漏；
     3. `SPEC_CODE_MISMATCH` (`CA-SPEC-001`)：省略规范要求的正则项/收敛条件/激活函数，或加入改变科学语义的未声明裁剪/激活（对应 `SciCoQA Fig. 1`）；
     4. `NUM_PROVENANCE_FABRICATION` (`CA-SCI-002`)：跨计算环节与简要结果声明一致引入 $\times 10^3$ 比例错误或 $5\%$ 偏差，或硬编码常数值覆盖模型输出；
     5. `DATA_UNIT_SCHEMA` (`CA-DATA-001`)：未屏蔽 `-999` 哨兵值、部分行混用单位、特征列错位。
   - **植入缺陷双确定性验证门（Defect Validity Gate）**：每个植入样本 $C_{m,g}^{\text{fault}}$ 必须同时满足：
     - (a) 可正常导入并执行，无语法/运行时崩溃；
     - (b) 通过**最小可见冒烟测试**（如返回形状与类型正确、非全 NaN，排除一眼即破的崩溃错误）；
     - (c) **确定性未能通过隐藏/外部科研真值测试**（如 `SciCode-Verified` 数值断言或 `BioXArena-Pro` `grade.py` 指标显著低于洁净基线或报错），且记录**精确植入文件路径、函数符号与行跨度**作为 $E_3$ 定位真值。

---

## 5. 匹配误报率（Matched-FPR）比较方法与统计推断

### 5.1 单项审计输出聚合与连续 ROC 构建

1. 每件样本 $x$ 在每个 $(G, A)$ 单元下在**全新上下文**中被独立读取 $K$ 次（推荐 $K = 4$，统一温度 $T=0.3$），输出需严格符合 `audit_response_schema.json`（不合规或无规则引用者记为格式失败 `audit_failure` 并单独统计失效率；其 `overall_block_score` 在主敏感度分析中分别按「最保守 0 分」与「剔除重采」双口径报告）。
2. 对每个样本 $x$ 计算跨 $K$ 次有效读取的平均阻塞风险分：
   $$\bar{s}_{g,a}(x) = \frac{1}{K_{\text{eff}}} \sum_{k=1}^{K_{\text{eff}}} s_{g,a,k}(x) \in [0, 100].$$
3. **预设防坍缩回退机制（Pilot-Gated Graded-Strictness Fallback）**：在进入正式批前先运行 $M_{\text{pilot}} = 10$ 题校准检查。若任一审计器在干净或缺陷样本上的 $\bar{s}_{g,a}(x)$ 出现极端二值坍缩（如中位离散唯一取值数 $< 5$，且在 $\text{FPR} \in [0.03, 0.20]$ 区间内不存在阈值点、出现 $> 5$ 个百分点的间隙），则自动在正式批中启用 `Dong (2026b §D.7)` 建议的三档严格度提示词重问（`lenient` / `balanced` / `strict`，总读取分配为 $2 + 2 + 2$ 次），合并分段 ROC 包络曲线；否则仅使用单档连续评分 + $K=4$ 读取（节省 33% 审计成本）。

### 5.2 精确确定性 Neyman-Pearson 匹配阈值（零随机抖动）

对于每个单元 $(g, a)$ 与目标误报率 $\alpha \in \{0.05, 0.10\}$（主基准 $\alpha = 0.10$），在干净样本得分集合 $\mathcal{C}_g$（样本数 $N_0 = M \cdot C$）上按 `matched_fpr_eval.py` 确定性求解临界阈值 $(\tau_{g,a}^*, \theta_{g,a}^*)$（其中 $\theta_{g,a}^* \in [0, 1]$ 为同分平局确定性权重），使得：
$$\frac{1}{N_0}\sum_{x \in \mathcal{C}_g} \Big[\mathbb{I}(\bar{s}_{g,a}(x) > \tau_{g,a}^*) + \theta_{g,a}^* \mathbb{I}(\bar{s}_{g,a}(x) = \tau_{g,a}^*)\Big] = \alpha \quad \text{（精确恒等）}。$$
对每件缺陷样本 $x \in \mathcal{D}_g$ 赋予匹配分数示性值 $z_{g,a}(x; \alpha) = \mathbb{I}(\bar{s}_{g,a}(x) > \tau_{g,a}^*) + \theta_{g,a}^* \mathbb{I}(\bar{s}_{g,a}(x) = \tau_{g,a}^*)$，得到 $\widehat{\text{TPR}}_{g,a}(\alpha) = \frac{1}{N_1}\sum_{x \in \mathcal{D}_g} z_{g,a}(x; \alpha)$。同时计算 **$\text{pAUROC}_{[0, 0.20]}$**（在实用可接受误报率区间 $[0, 0.20]$ 内的归一化 ROC 面积，对单点分位波动更稳健）与全曲线 **$\text{AUROC}$**。

### 5.3 阈外实际校准验证（5-Fold OOF Calibration Check）

鉴于实际部署中并不知道新任务池的精确 $\tau^*$，我们对每个 $(g, a)$ 单元额外执行 **5 折任务簇交叉验证**（`oof_threshold_calibration_check`）：在 4 折干净样本上拟合 $\hat{\tau}^{\text{train}}(\alpha)$ 并应用于留出折的干净与缺陷样本，报告真实留出 $\widehat{\text{FPR}}_{\text{OOF}}$ 是否落在 $[0.05, 0.15]$ 容差内以及对应的 $\widehat{\text{TPR}}_{\text{OOF}}$。

### 5.4 统计推断与我们通过仿真发现的陷阱规避

> **方法学自检发现（已由 `power_and_cost_sim.py` 前期验证）**：若先在全样本干净集上拟合 $\hat{\tau}_{g,a}^*$ 固定不动，再把 $z_{g,a}(x)$ 视作固定二值结果对任务簇差值 $d_m$ 做符号翻转置换检验，**会导致零假设下的 I 类错误严重膨胀至 `0.300`（30%）**！原因是有限干净样本上的分位阈值估计误差 $(\hat{\tau} - \tau^*)$ 作为共享偏移注入所有 $M$ 个任务，被固定阈值置换检验误判为真实信号。

因此，本预注册**禁止**在匹配阈值上使用固定阈值符号翻转检验，**唯一锁定「含阈值重估的任务簇自助法（Threshold-Re-Estimated Task-Cluster Bootstrap）」**作为主推断方法：
- 在全部 $4,000$ 次按任务簇 $m \in \{1, \dots, M\}$ 有放回重采样中，**在每一次重采样 $b$ 内部对全部 4 个 $(G, A)$ 单元重新估计 $\hat{\tau}_{g,a}^{*(b)}(\alpha)$、$\text{TPR}^{*(b)}$、$\text{pAUROC}^{*(b)}$ 与 $\text{AUROC}^{*(b)}$**，从而把干净样本阈值估计方差与缺陷样本抽样方差同时正确传导至 95% CI 与双侧自助法 $p$ 值。
- **跨指标多重检验控制**：对两个 Co-Primary 指标（$\text{TPR@10\%FPR}$ 与 $\text{pAUROC}_{[0,0.20]}$）采用 **Holm-Bonferroni** 校正控制族错误率 $\alpha_{\text{FWER}} = 0.05$；$E_4$（修后外部通过率）采用独立确认或并入 3 项 Holm 家族（以主方案明确为准，见下表）。

---

## 6. 样本量、功效分析、I 类错误与调用成本估计（全数可追溯至脚本输出）

以下数字**全部由本地脚本计算并输出至 `outputs/`**（`power_and_cost_sim.py` 种子 `20261121`，250 次 MC $\times$ 200 次簇内重估自助法；基线 $\text{TPR@10\%FPR}=0.35$，强审计器非对称偏移 logit $= 0.80$ 即未校准前误报率约 `5% vs 25%`，生成器非对称 logit $= 0.30$，任务内 $\text{ICC}=0.20$，跨审计器配对相关 $=0.35$）：

### 6.1 样本量与统计功效对照表（来源：[simulation_summary.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/outputs/simulation_summary.md)）

| 方案档位 | 独立任务数 $M$ | 干净/任务/厂 $C$ | 植入缺陷/任务/厂 $D$ | 零假设下 I 类错误 (`TPR10` / `pAUC20` / `AUROC`) | 95% CI 覆盖率 (`TPR10` / `pAUC20` / `AUROC`) | $\Delta=0.10$ 功效 (`TPR10` / `pAUC20` / `AUROC`) | $\Delta=0.12$ 功效 (`TPR10` / `pAUC20` / `AUROC`) | $\Delta=0.15$ 功效 (`TPR10` / `pAUC20` / `AUROC`) |
|---|---:|---:|---:|---|---|---|---|---|
| 最小可行试探档 | 40 | 1 | 2 | `0.008 / 0.036 / 0.044` | `0.992 / 0.964 / 0.956` | `0.124 / 0.244 / 0.448` | — | `0.272 / 0.480 / 0.680` |
| 标准入门档 | 60 | 1 | 2 | `0.012 / 0.048 / 0.088` | `0.988 / 0.940 / 0.900` | `0.204 / 0.364 / 0.556` | — | `0.376 / 0.632 / 0.844` |
| **推荐主确认档** | **80** | **2** | **2** | `0.020 / 0.040 / 0.056` | `0.980 / 0.968 / 0.948` | `0.360 / 0.636 / 0.852` | `0.524 / 0.736 / 0.924` | `0.772 / 0.912 / 0.984` |
| **高功效顶刊档** | **120** | **2** | **2** | `0.040 / 0.064 / 0.064` | `0.960 / 0.936 / 0.948` | `0.592 / 0.812 / 0.952` | `0.708 / 0.868 / 0.976` | — |
| 超高分辨档 | 160 | 2 | 2 | `0.020 / 0.032 / 0.052` | `0.980 / 0.964 / 0.944` | `0.700 / 0.864 / 0.972` | — | — |

**功效审计结论**：
1. 为什么 $M=40, C=1$ 不足？在 $N_0=40$ 干净样本上估计 $10\%$ 分位阈值本身带来约 $\pm 0.08\text{–}0.10$ 的 SD，即使 $\Delta=0.15$ 时 `TPR@10%` 功效也仅 `27.2%`！
2. 提升至 $C=2$ 洁净变体且 $M \in \{80, 120\}$ 后：
   - **$M=80, C=2, D=2$（推荐档）** 在 $\Delta=0.10$ 时 `AUROC` 功效达 **`85.2%`**，`pAUROC[0,0.20]` 在 $\Delta=0.12\text{–}0.15$ 达 **`73.6%–91.2%`**；在 Endpoint $E_4$（一轮修后通过率差 $+0.10$）上功效达 **`92.6%`**（$+0.08$ 为 `77.2%`，零假设 T1E=`0.054`，见 [post_repair_power.json](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/outputs/post_repair_power.json)）；在等预算 2-call 匹配 FPR 互补性上（当低温同厂商重复相关性 $\rho_w=0.85$、跨厂商 $\rho_c=0.45$ 时）检测 Jaccard 下降（$\Delta J \approx 0.249$）功效达 **`95.1%`**、2-call 联合增益（`+8.0` pp）功效达 **`69.6%`**（$\rho_c=0.35$ 时 **`81.2%`**，见 [complementarity_power.json](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/outputs/complementarity_power.json)）。
   - **$M=120, C=2, D=2$（高功效档）** 将 `pAUROC[0,0.20]` 在 $\Delta=0.10$ 处推过 80% 线（**`81.2%`**），$E_4$ 在 $+0.08$ 处达 **`90.3%`**，2-call 联合增益在 $\rho_c=0.45$ 处达 **`85.1%`**。

### 6.2 调用次数与指示性成本估计（来源：[cost_estimates.json](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/metasci_2x2/outputs/cost_estimates.json)）

按单次调用经验足迹假设：生成约 `2.5k` 输入 / `1.5k` 输出 tokens（含平均 1.5 次生成尝试确保双厂商过 Clean Gate）；审计约 `5.0k` 输入 / `1.0k` 输出 tokens；1 轮修复约 `6.0k` 输入 / `1.5k` 输出 tokens。

| 方案档位 | 厂商数 $V$ | 目标任务 $M$ (生成池) | 每任务 $(C,D)$ | 读取 $K$ | 总制品数 | Pilot 调用 | 阶段一生成调用 | 阶段二审计调用 | 阶段三修复调用 | **总 LLM 调用数** | 估算 Token (In / Out M) | 经济模型 (~$0.5/$2.5) | 中等旗舰 (~$3/$15) | 高阶旗舰 (~$5/$25) |
|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|
| `2x2_Minimal` | 2 | 40 (55) | (1, 2) | 3 | 240 | 152 | 165 | 1,440 | 144 | **1,901** | 9.18 / 2.07 M | $9.76 | $58.54 | $97.57 |
| `2x2_Standard` | 2 | 60 (80) | (1, 2) | 4 | 360 | 190 | 240 | 2,880 | 216 | **3,526** | 17.17 / 3.77 M | $18.01 | $108.05 | $180.08 |
| **`2x2_Recommended`** | **2** | **80 (105)** | **(2, 2)** | **4** | **640** | **190** | **630** | **5,120** | **320** | **6,260** | **29.97 / 6.75 M** | **$31.86** | **$191.16** | **$318.60** |
| **`2x2_HighPower`** | **2** | **120 (155)** | **(2, 2)** | **4** | **960** | **228** | **930** | **7,680** | **480** | **9,318** | **44.66 / 10.04 M** | **$47.43** | **$284.58** | **$474.30** |
| `3x3_Extension` | 3 | 80 (115) | (2, 2) | 3 | 960 | 405 | 1,035 | 8,640 | 720 | **10,800** | 52.02 / 11.70 M | $55.26 | $331.56 | $552.60 |
| `4x4_AllLocalVendors` | 4 | 80 (130) | (2, 2) | 3 | 1,280 | 700 | 1,560 | 15,360 | 1,280 | **18,900** | 91.73 / 20.35 M | $96.74 | $580.44 | $967.40 |

> **注**：若采用 `Track B`（AgentHub 本机订阅模式），边际 API 现金成本接近 `$0`，但受制于各桌面订阅每小时/每日消息速率限制（如 `6,260` 次调用按单条 15–30 秒串行需约 26–52 小时，4 通道并行约 8–15 小时且可能触及订阅限流），且引入 CLI 工具链混杂（详见 `HUMAN_DECISIONS.md`）。

---

## 7. 操作执行流程与硬停机门槛（Execution Protocol & Stopping Gates）

1. **Step 0 — 锁定并提交预注册**：在运行任何生成或审计 LLM 调用前，由主线执行者提交 `PREREG_2x2.md`、`CONSTITUTION_2x2.md`、`audit_response_schema.json` 与 `matched_fpr_eval.py`，并在外部运行日志记录提交 commit 哈希。
2. **Step 1 — 数据准备与确定性门控验证（0 LLM 费用）**：由 `A-T1` 在 `p2/shared/DATA_SOURCES.md` 登记 `SciCode-Verified`、`BioXArena-Pro-Public`（非 `xatlas`）与 `SciCoQA` URL、许可证、大小与 SHA-256。
3. **Step 2 — 10 题 Pilot 门控（停机检验 1）**：
   - 检查 Clean Gate 通过率（要求双方生成器在 $\le 2$ 次尝试内联合通过率 $\ge 55\%$；否则调整子步骤难度池并在 Amendment 记录，**不读取任何审计器结果**）。
   - 检查 Pilot 上 4 个单元的格式合规率（`I3` 格式失败率必须 $< 10\%$）与分数分散度（检验是否需要触发 3 档严格度提示词回退）。
   - 检查是否发生任何预算超支（实际单次 token 超预估 $1.8\times$ 时立即暂停并由人类调整 $K$ 或 $M$）。
4. **Step 3 — 正式批生成、确定性注入与 2×2 审计（全盲态/顺序随机化）**：
   - 所有样本 ID 随机乱序发送，每次审计附带精确 `constitution_commit`；原始请求/响应（包含 `system_fingerprint` / 服务端返回模型 ID / token 用量）完整落盘保存（补齐 `CrossAudit I2` 第二层 *raw-exchange preservation*）。
5. **Step 4 — 锁定分析与受控 1 轮修复**：
   - 直接调用冻结的 `matched_fpr_eval.py` 输出主表格，随后对 matched 10% FPR 拦截样本运行 1 轮带防作弊修复守卫（禁止删断言、跳过测试、`try/except` 吞错误）的修复调用并执行确定性单元测试。

---

## 8. 文献核实记录（2026-11-21 通过 arXiv / HuggingFace 公开 API 核验）

> **arXiv API 使用声明**：本设计使用 arXiv 公开 API 检索并下载论文全文（遵循 https://info.arxiv.org/help/api/index.html 使用条款，请求间隔 $\ge 3$ 秒），读者与复现者应自行核查各检索论文与数据集的原始许可证限制。

| 文献 / 数据集 | 标识符 / URL | 核实方式 | 核实要点 |
|---|---|---|---|
| Dong & Chen (2026). *CrossAudit: A Git-Native, Cross-Vendor Audit Loop for Agentic Science*. | https://arxiv.org/abs/2608.28631 | arXiv API + PDF 全文精读 | 2 位作者，2026-08-05，19 页；确认 §3.2/§4.2/§4.4/§5–6 原文边界 |
| Dong (2026). *What Moves the Ceiling on Automated Audit: from General Code to Scientific Code, Results and Data*. | https://arxiv.org/abs/2611.23250 | arXiv API + PDF 全文精读 | 1 位作者，2026-09-24，27 页；确认 §3–4/§A–G 全文限制与统计方法 |
| Song (2026). *When Does a Second Model Help? Cross-Model Review in LLM Verification*. | https://arxiv.org/abs/2610.01471 | arXiv API + PDF 全文精读 | 1 位作者，2026-10-01 (`v2` 10-06)，16 页；确认 §4–5 与 Limitations |
| Baumgärtner & Gurevych (2026). *SciCoQA: Quality Assurance for Scientific Paper--Code Alignment*. | https://arxiv.org/abs/2601.12910 | arXiv API + PDF 全文精读 + HF API | ACL 2026；`UKPLab/scicoqa`，`CC BY 4.0`，`3,657,481 bytes`（real 92, synthetic 543, pooled 129） |
| Hu, Huang, Deng & Chen (2026). *SciCode-Verified*. | https://arxiv.org/abs/2608.04975 | arXiv API + HF API | 2026-08-05；`shhu2001/SciCode-Verified`，`Apache-2.0`，`problems_test.jsonl` (`951,750 bytes`, 64 题 / 290 子步) |
| Li et al. (2026). *BioXArena-Pro*. | https://arxiv.org/abs/2611.13994 | arXiv API + 本地仓库核对 | 2026-11-10；`CC BY-NC 4.0` 数据 + `MIT` 评分脚本，64 题 / 8 域 |
| Vaccaro (2026). *Preregistration for Experiments with AI Agents*. | https://arxiv.org/abs/2606.11217 | arXiv API | ICML 2026 Spotlight；AI 实验自由度与预注册模板 |
| Panickssery, Bowman & Feng (2024). *LLM Evaluators Recognize and Favor Their Own Generations*. | https://arxiv.org/abs/2404.13076 | arXiv API | NeurIPS 2024；自识别因果影响自偏好基线 |
| Chen et al. (2025). *Do LLM Evaluators Prefer Themselves for a Reason?* | https://arxiv.org/abs/2504.03846 | arXiv API | 区分合理能力偏好与有害错误自偏好 |
| Roytburg et al. (2026). *Are LLM Evaluators Really Narcissists?* | https://arxiv.org/abs/2601.22548 | arXiv API | ICML 2026；评估者质量基线修正自偏好 |
| Tian et al. (2024). *SciCode: A Research Coding Benchmark Curated by Scientists*. | https://arxiv.org/abs/2407.13168 | arXiv API | NeurIPS 2024；80 题 / 338 子步源基准 |
| Bansal & Balaji (2026). *OpenDiscoveryTrace*; Son et al. (2025). *SPOT*; Thomas, Gligoric & Shah (2026). | `2609.09203`, `2505.11855`, `2606.27687` | arXiv API | 过程轨迹审计、低召回科研验证与 LLM 预注册防 p-hacking |
