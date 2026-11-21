# 文献全文精读审计：arXiv:2608.28631 与 arXiv:2611.23250 的已做与未做清单

本文档逐节对照两篇核心先驱文献的**全文**（2026-11-21 通过 arXiv API 下载并全文精读），并交叉核对 `arXiv:2610.01471v2`（Song, 2026）与 `arXiv:2601.12910v3`（SciCoQA），精确列明其**已完成工作**与**尚未解决的开放问题**，作为 `PREREG_2x2.md` 实验设计的直接依据。

---

## 1. 文献一：Dong & Chen (2026), *CrossAudit: A Git-Native, Cross-Vendor Audit Loop for Agentic Science* (`arXiv:2608.28631v1`, 2026-08-05)

### 1.1 已做工作（附章节、图表与实测数字）

| 维度 | 具体已完成内容 | 引用章节 / 图表 |
|---|---|---|
| **问题定义与定位** | 提出代理科学中的「同源监督问题（same-source supervision problem）」：生成器与审查器同家族/同厂商时面临自偏好偏差与潜在相关盲区，且监督记录常留存在不可回放的平台临时日志中。与 The AI Scientist [1]、AI co-scientist [4]、PoLL [18]、非正式跨厂商建议及 `statcheck` [19]/PPS [20] 做结构对比。 | §1–§2, Table 1 (pp. 1–4) |
| **协议形式化（八项不变量）** | 提出基于双仓库分权（Science Repo vs Audit Repo）与人类撰写版本化规则手册（Constitution）的 CrossAudit 协议，定义 8 项不变量：`I1` 异构性（参数、上下文、凭证三重隔离）、`I2` 账本完备性（区分 3 层，当前达第 1 层 parsed-record inspectability）、`I3` 引用有效性（未引用规则或无法校验宪法版本直接 ESCALATE）、`I4` 确定性优先（DCL 脚本先于 LLM 且不可被模型豁免）、`I5` 有界修订（默认 `max_rounds = 3` 防震荡）、`I6` 分级中断（仅 `BLOCKER` 与 `ADVISORY` 两级，避免严重度通胀谈判）、`I7` 回执绑定（SHA-256 清单与版本绑定）、`I8` 失效关闭准入（无有效回执拒绝准入）。定义 4 态终止图与人类发起的单次争议通道（`CA-META-001`–`004`）。 | §3.1–§3.4, Figs. 1–3, Table 2 (pp. 4–8) |
| **开源参考实现与自审计** | 基于 GitHub Actions 与 Python 构建参考实现（`github.com/dongzhaohe321418-lab/crossaudit`, MIT）。对本仓执行 3 轮跨厂商审计（记录于 `audits/`），修复作用域仅由 SHA 推导、回执校验、`I3 > I4` 裁决优先级及围栏提示词防注入等接线漏洞。 | §4.1, §4.5 (pp. 8–9, 14–15) |
| **计算化学真实部署案例（无对照观测）** | 在钙钛矿添加剂筛选管道中实际运行变体（Claude 生成器 + Codex-CLI 审计器）：第 1 阶段 7 轮（每轮发现数从 `7→4→3→2→1→1→0` 收敛，14 个不同问题中 12 个行为确认关闭）；截至 2026-08-03 第 2 次冻结累计 21 轮、39 次发现事件（23 个不同 ID，20 个确认关闭），其中 `11/39` 标记为 `DETERMINISTIC`、`28/39` 为 `JUDGMENT`，单轮审计中位耗时 `6` 分钟（范围 `3.7–18.2` 分钟）。作者明确将其定性为无对照组的「存在性证明（existence proof）」，且指出 3 处与参考协议的偏离（工具型 CLI 审计器、外置控制器状态、4 级严重度）。 | §4.2, Fig. 4 (pp. 9–11) |
| **第二实现教训（产品线）** | 基于 `crossaudit_v3`（v3.2.0）提炼 3 条协议细化：准入保障是分层级而非开关；对话式路由器属于可信基内部；领域规则包缺失时 `I4` 机械保障变薄弱，界面需诚实显示「未测量」而非「0」。 | §4.3, Table 3 (pp. 11–12) |
| **探索性植入缺陷小试（30 增量 / 43 缺陷）** | 用**确定性 Python 脚本**生成 30 个合成增量（10 个干净 + 20 个含缺陷，跨 10 类共 43 个植入缺陷），对比 DCL（`17/43`）、Claude 子智能体（`38/43` 冻结映射 / `43/43` 事后裁定，`0/10` 干净增量拦截）与 `gpt-5.1`（`T=0`，`41/43` 宽松 / `37/43` 严格，但对 **10/10 个干净增量全部拦截**，发出 26 条误报 BLOCKER）。使用问题簇自助法与符号翻转检验（宽松差异 `+3/43`, 95% CI `[-2, +8]`, `p=0.46`；严格差异 `-1/43`, `[-6, +4]`, `p=1.0`）及 2000 次置换零基线（`4.8 / 22.4 / 31.8`）。披露跨厂商自审推翻了原先的「盲法」宣称，并记录 1 项 5 增量单位换算算术失败案例（`SMOKE-FINDINGS.md`）。 | §4.4, Table 3 (pp. 12–14) |
| **威胁模型与例外管理假说** | 系统梳理共享语料相关错误、提示词注入、宪法 Goodhart 博弈、共谋、操作者篡改、自声明厂商身份（未做加密证明）与保密性边界；提出从 $O(\text{increments})$ 降至 $O(\text{escalations})$ 的人类注意力经济学假说。 | §5–§7 (pp. 15–18) |

### 1.2 未做部分（含原文明确声明的开放问题与章节定位）

1. **完全未实例化「同厂商 vs 跨厂商」对称比较（§4.4, p. 13）**  
   - §4.4 语料由**确定性脚本生成器**而非 LLM 生成，原文明确声明：*"The corpus's author is the deterministic generator... no arm audits work that its own model family produced, and the own-generation setting in which self-preference is demonstrated [5] is not instantiated here."*  
   - 原文仅测了两家审计器读同一非 LLM 语料，**未形成「生成器厂商 × 审计器厂商」2×2 交叉**。
2. **未测量 `I1` 的「盲区去相关（blind-spot decorrelation）」核心猜想（§3.2, p. 4; §4.4, p. 14）**  
   - §3.2 (`I1`) 明言：*"It is likewise designed to decorrelate vendor-idiosyncratic failure modes, an effect this paper does not measure (§4)."*  
   - §4.4 明确总结：*"a single model pair, run once, cannot separate a vendor effect from a prompt, temperature, or sampling effect."*
3. **未在匹配误报率（matched FPR）下对比审计器（§4.4, pp. 13–14; §6, p. 18）**  
   - `gpt-5.1` 对干净增量误拦截率高达 `10/10`（100%），而 Claude 为 `0/10`。作者指出这种情况下 `gpt-5.1` 等同于无条件拦截桩（"reliability-without-validity failure [9]"），原始召回率无法裁决审计器优劣，且 10/10 误拦截会摧毁 $O(\text{escalations})$ 经济学前提。
4. **工具链（harness）与采样参数混杂（§4.4, pp. 12–14）**  
   - Claude 组通过带会话采样的 subagent 运行，GPT-5.1 组通过 CI runner 以 `temperature=0` 运行，且每组仅单次运行。
5. **真实部署无对照组、未测跨厂商审计对修复后外部科研结果的影响（§4.2, pp. 10–11; §4.4, p. 13）**  
   - 真实部署缺乏无审计或同厂商对照组；人工审阅注意力耗时未记录；组件级试验未测量经审计反馈修改后的真实可执行科研指标变化；`I2` 未保存原始模型请求/响应；`I1` 厂商身份仅自声明。

---

## 2. 文献二：Dong (2026), *What Moves the Ceiling on Automated Audit: from General Code to Scientific Code, Results and Data* (`arXiv:2611.23250v1`, 2026-09-24)

### 2.1 已做工作（附章节、图表与实测数字）

| 维度 | 具体已完成内容 | 引用章节 / 图表 |
|---|---|---|
| **审计分层架构刻画** | 刻画已发布的 CrossAudit 系统 3 档配置（`general`、`science` 含 `number→source` 溯源与单位/收敛检查、`research` 含文献工具引用核验）、确定性裁决优先级、防篡改修复守卫（禁止删除断言/跳过测试/新增宽泛异常捕获）与加密摘要回执。 | §2 (pp. 2–3) |
| **通用代码天花板测量（HumanEval+ / MBPP+）** | 冻结单一生成器 `claude-haiku-4-5-20251001` 生成 540 题 × 2 遍，筛选「通过可见测试但未通过隐藏测试」的 **110 个缺陷样本**（56 题）与 **150 个正确样本**（143 题）。测量：  - **多次独立阅读收益（$K=1\to 8$）**：默认跨厂商审计器（`gpt-5.6-terra`）联合召回率从 `10.7%` 升至 `30.0% [20.0, 40.7]`（代价为 `16.0% [10.1, 22.3]` 正确代码误报），第 8 次仍增 `+1.93` pp（未达 `<1.0` 平台线）；3 家族 20 次联合达 `48.2%`（误报 `36.0%`）。  - **模型路线与裁决规则扫频**：更强同厂商模型（`claude-sonnet-4-6`, `T=0`）在 BLOCKER 规则下召回率 `-26.4` pp（`[-37.3, -15.6]`，误报 `3.3%` vs `16.0%`），同规模同厂商（`claude-haiku-4-5`, `T=0`）为 `-12.7` pp（`[-25.0, -0.9]`）；但在事后 `any-finding` 规则下方向反转（`self-strong` `59.1%` vs `cross` `34.5%`，误报 `34.0%` vs `19.3%`）。高推理跨厂商模型 `gpt-6-astra`（$K=4$）达 `32.7%`（误报 `10.7%`）。  - **规则手册干预**：加一句「寻找可见测试未覆盖之处」，在 56 个单轮修订缺陷上使 $K=1$ 标记率增 `+26.8` pp（正确代码 `+12.5` pp `[-3.6, +26.8]`），但**一轮修订后隐藏测试通过率仅 `+5.4` pp（`[-0.9, +12.1]`, 未显著）**；在 110 样本 $K=8$ 上召回率从 `30.0%` 升至 `60.9%`，代价是误报飙升 `+20.7` pp 至 `36.7%`。  - **漏检残差归因与澄清干预**：3 家族 20 读漏检的 57 个缺陷中，事后双模型评审认为 `44/57`（`77.2%`）属「规范未唯一确定（specification-undetermined）」；补充规范后诊断从 `0/32` 升至 `9/32`（簇符号翻转 `p=0.0625`，全集中于未与隐藏 Oracle 矛盾的 9 例）。 | §3, §A–§B, Figs. 1, 3, 4, Table 1 (pp. 3–5, 10–20) |
| **三项预注册科研审计研究** | **1) 科学代码（SciCode, §4.1/§C.1）**：由 `claude-haiku-4-5` 生成 81 个缺陷实例（30 题）与 150 个正确实例（53 题），审计器不可见测试。整体问题表述下 $K=8$ 召回率 `87.7%`、正确代码误报 `68.7%`；改写任务文本聚焦单步后，误报降至 `36.0%`（`-32.7` pp），召回 `72.8%`（`-14.8` pp, `p=0.0625`）。生成器同族 CLI 路线召回比跨厂商低 `46.9` pp。  **2) 科学结果（74 SciCode 实例, §4.2/§C.2）**：49 例报告不一致（vs 日志/results.json）+ 25 例跨文件一致伪造（$\times 10^3$ 或 $5\%$）+ 74 例干净对照。确定性 `science` 配置（`number→source`）检出 `49/49` 不一致与 `0/25` 一致伪造（`0` 误报）；LLM 审计器（4 读）检出 `23/25` 一致伪造（事后核查 21 例合理，McNemar `p=2.4e-7`），但对干净样本误报 `12/74`（事后发现 `8` 例实为基准构造缺陷、`4` 例模型算术/语法误判）。  **3) 科学数据（4 份 UCI CC BY 4.0 物理表, 280 样本, §4.3/§C.3）**：140 干净 + 140 含 7 类注入缺陷。文档派生确定性校验器检出 `92/140`（`0` 误报）；跨厂商审计器检出 `65/140`（干净误报 `4/140 = 2.9%`）；二者联合 `112/140`（`80.0%`）。模型补齐了混用单位与换列（校验器漏检中抓回 `18` 例），但漏掉全部 `20/20` 打乱目标列与 `20/20` 重复行。同族 CLI 模型在此任务彻底失效（误报 `67/140` 干净数据，至少 48 例存在事实性行列计数幻觉）。 | §4.1–§4.3, §C.1–§C.3, Fig. 2 (pp. 4–7, 20–22) |
| **方法学审计与透明化** | 提出「洁净样本必须先经审计标准本身清洁核验」原则；通过模拟揭示样本率 `<0.10` 时问题簇百分位自助法覆盖率降至 `~0.80`（§A.3）；记录所有研究经 `gpt-6-astra` 3–21 轮独立审查与拒收历史。 | §5, §A.3–§A.4, §D–§G (pp. 7–8, 11–27) |

### 2.2 未做部分（含原文明确列举的限制与后续计划）

1. **单一生生成器（Single Generator），无法区分「同厂商特性」与「审计模型强弱」（§3, §D.6 Limit 5, §D.7, §F p. 27）**  
   - 全部程序和结果仅由单一模型 `claude-haiku-4-5-20251001` 生成。  
   - 原文 §F (p. 27) 明确划清边界：*"We do not measure self-preference. A frozen generator's code audited by routes that differ in model, operating point and sampling cannot separate preference from performance, so we report the same-vendor contrast as a contrast between routes, on one substrate."*  
   - §B.2 (p. 15) 明确留白：*"whether it is a property of same-vendor auditing is open."*  
   - §D.7 (p. 25) 明确建议：*"A second generator would show whether the residual's character belongs to the benchmark or to one model's failures."*
2. **未在匹配误报率（Matched Operating Point / Matched FPR）下比较审计器（§D.6 Limit 2, §D.7, Fig. 3）**  
   - §D.6 Limit 2 原话：*"No auditor is compared with another at a matched false-positive rate. Our sweep re-grades archived findings at four deterministic decision rules, under which one family has no point between 3.3% and 34.0% false positives. We built no randomised policy between them, so the paper cannot say whether the families differ at the same operating point."*  
   - §D.7 点出最需要的下一步实验：*"The experiment this paper most needs is to re-ask, instructing each family at graded severity and reading afresh, which can place the families at a common false-positive rate and decide whether the auditor contrast survives there."*
3. **采样温度与调用渠道不对称（§A.2, §B.2, §C.1, §D.2, §D.6 Limit 1）**  
   - 通用代码实验中，跨厂商组不传采样参数（8 次独立阅读输出波动，联合召回增 `+19.3` pp），同厂商组传 `temperature=0`（输出几乎确定，8 次仅增 `+1.8` 至 `+3.0` pp）。科学代码/数据实验中，同厂商改走带额外上下文的 CLI 工具链而非同构 API。§D.6 Limit 1 明言：*"Every cross-family union-of-K contrast on general-purpose code inherits this, and it is the confound we would fix first."*
4. **无先验样本量与最小可检测效应（MDE）功效计算（§A.3 p. 11, §D.6 Limit 3, §E p. 26）**  
   - 原文 §E (p. 26) 明言：*"Reading depth K, the number of arms, and the size of every intervention population were set by cost rather than by a power calculation, and no sample size was justified in advance by a minimum-detectable-effect argument."*
5. **科学场景未测多文件/跨模态论文–代码对齐，也未测科学任务修后外部通过率（§D.6 Limit 4, §E p. 25）**  
   - SciCode 为单步独立函数；结果与数据研究为小规模合成包装；仅通用代码 1 项规则研究测过 1 轮修后隐藏测试通过率（`+5.4` pp, 不显著），科学代码未测闭环修复后的确定性通过率增益。
6. **无外部人工裁决，且评审模型 `gpt-6-astra` 存在「既做选手又做裁判」的自评冲突（§A.4, §D.6 Limit 6, §E p. 26）**  
   - 原文要求对涉及 `astra` 排名的结论保留悬置，直至非参评模型或外部确定性裁判核验。

---

## 3. 关联基准与互补文献对照（为什么 `2×2` 设计必须规避其已知陷阱）

1. **Song (2026, `arXiv:2610.01471v2`)**：测试 30 份韩语制品 × 150 个人工错误 × 10 种审查条件（900 次会话），同样**仅用单一生成器**（`Claude Opus 4.6`），且混淆了审查模型能力、CLI/API 工具链（Codex CLI vs Gemini API vs Claude Code）、提示语言（韩/英混杂）与信息限制，作者在 Limitations 中显式指出：*"reviewer capability and model identity are confounded... a factorial design is left to future work"*，且明确提醒不同审查器的 `Jaccard` 差异受检出率高低影响（§5.3）。
2. **Baumgärtner & Gurevych (2026, `arXiv:2601.12910v3`, SciCoQA)**：提供 `CC BY 4.0` 公开多文件论文–代码对齐数据集（92 例真实 + 543 例合成，129 例 pooled；总下载 `3.49 MiB`）。但注意两点不可直接照搬的缺陷：  
   - 其 `543` 例合成差异由 `GPT-5` 生成，论文 §7（Table 10）实证 `GPT-5` 在合成子集召回率虚高（`70.0%` vs 真实子集 `41.3%`），存在**数据生成模型亲和污染**；  
   - 其原基准未控制误报率（不同模型每篇预测条数 `# Preds.` 从 `1.6` 到 `6.7` 不等，Table 10），且主要依赖 `GPT-OSS 20B` 作为 LLM 裁判（本身存在自偏好提示，§5）。
   - **因此本设计吸收 `SciCoQA` 的分类体系与干净原始仓库/真实无厂商污染子集，同时以「双厂商对称生成干净基线 + 第三方确定性脚本/人工植入科学缺陷 + 确定性可执行真值与精确文件/符号定位」消除 LLM 裁判偏差与 GPT-5 生成污染。**

---

## 4. 总结映射：本 `PREREG_2x2.md` 精确填补的 6 项空白

| 开放空白（来源） | 本 `2×2` 预注册方案的直接解法 |
|---|---|
| 单生成器混淆「同/跨厂商」与「模型强弱」(`2608.28631 §4.4`, `2611.23250 §D.6/§F`) | **对称 $2\times 2$ 交互估计量** $\Delta_{\text{cross-vs-same}} = \frac{1}{2}[(Y_{12}-Y_{11}) + (Y_{21}-Y_{22})]$ 精确消去生成器与审计器主效应，并辅以双向一致性检验。 |
| 误报率/操作点不匹配 (`2608.28631 §4.4`, `2611.23250 §D.6-D.7`) | 输出连续阻塞风险分 `overall_block_score` $\in [0, 100]$（$K=4$ 读平均）+ **确定性 Neyman-Pearson 精确匹配 FPR**（$\alpha \in \{0.05, 0.10\}$）+ $\text{pAUROC}_{[0,0.20]}$ + 必要时分层严格度重问回退，外加 5 折交叉验证阈外样本实测 FPR 校准检验。 |
| 采样参数与工具链混杂 (`2608.28631 §4.4`, `2611.23250 §D.6 Limit 1`) | 主确认轨迹（Track A）强制同构无工具围栏 API、同温度（$T=0.3$）、同上下文隔离；同时单独对比 AgentHub 本机 CLI 轨迹（Track B）并显式量化工具链效应。 |
| 仅「标记」未核实「定位准确度」与「修后真实科研通过率」(`2611.23250 §D.4, §B.3`) | 增加确定性文件/符号定位真阳率（TLR@10%FPR）及**受控 1 轮修复后的外部确定性测试通过率增量**（含防作弊修复守卫）。 |
| 无先验功效与阈值估计方差分析 (`2611.23250 §A.3, §E`) | 专门针对「含干净样本阈值重估」的 2×2 任务簇自助法做 Monte Carlo 功效与 I 类错误模拟（见 `power_and_cost_sim.py`），揭示并规避固定阈值置换检验导致的 ~30% I 类错误膨胀陷阱。 |
| 基准干净样本未自洽洁净与 LLM 裁判偏差 (`2611.23250 §D.5`, `2601.12910 §5/§7`) | 双阶段**洁净门控 + 确定性缺陷门控**（干净通过所有测试与 DCL 审查，植入缺陷通过冒烟测试但必败于隐藏/外部科研真值），零 LLM 裁判。 |
