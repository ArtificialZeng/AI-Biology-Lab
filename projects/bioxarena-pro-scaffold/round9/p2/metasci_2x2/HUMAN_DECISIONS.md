# HUMAN_DECISIONS：2×2 同厂商 vs 跨厂商审计实验人类决策清单

**决策对象**：Zijian Zeng（人类 PI）｜ **协同审核**：Claude（主控）  
**状态**：等待人类审批。A1 严格遵守任务卡红线，**未调用任何付费 API、未向 AgentHub 其他智能体发送消息、未下载大型测试集、未运行 git**。

---

## 决策 1：选择实验执行通道与账号授权（关键混杂 vs 现金预算权衡）

根据对本机 `/Users/mac/AgentHub/src/adapters/`（`claude.ts`、`codex.ts`、`deepseek.ts`、`antigravity.ts`）的源码只读核查，四种桌面智能体在技术控制能力上与直连 API 存在本质差异：

| 执行轨道 | 使用本机账号/凭证 | 上下文隔离 (`I1`) | 温度/采样控制 | 工具链 (`harness`) 一致性 | 原始 API 报文留痕 (`I2`) | 现金成本 | 方法学定位 |
|---|---|---|---|---|---|---|---|
| **Track A：直连统一 HTTPS API**（**推荐主确认轨**） | Anthropic / OpenAI / Google / DeepSeek API Key | 每次新建无状态请求，100% 隔离 | **可精确锁定**（如全线 `T=0.3`） | **完全同构**（单一无工具围栏 Python 脚本） | **完整落盘**（含 `model` 指纹与精确 token） | 按量计费（见决策 3） | 彻底剥离 `Dong 2026b §D.6` 与 `Song 2026` 的采样及工具链混杂，**适合顶刊主表** |
| **Track B：AgentHub 本机订阅 CLI**（低现金备选） | 已有 Claude Code / Codex / DeepSeek (`codewhale`) / Antigravity 桌面订阅 | 支持 `newThread`，但需先关闭 `mcp__agenthub` 与本地文件工具 | **不支持传 `temperature`**（四个适配器均无温度参数） | **混杂**（分别经 `claude -p`、`codex app-server`、`codewhale`、`agentapi` 运行，且默认注入 `HUB_ROLE_PROMPT`） | 仅能录得适配器事件流，缺统一 API 指纹 | 边际 API 现金接近 **$0**（受订阅速率限流约束） | 只能测「端到端厂商 CLI 产品组合」，**无法将模型族与工具链/默认采样分离** |
| **Track A+B：Pilot 探路 + 正式直连** | 先用少额 API 或 Track B 跑 Pilot | 同上 | 正式批全锁定 | 正式批全同构 | 正式批完整 | 先低成本验证后再放行正式预算 | 推荐折中方案：用 10 题 Pilot 先核验 JSON 合规率与分数分布，再放行全量费用 |

### 请人类决定（勾选一项）：
- [ ] **批准 Track A（直连 API 主确认轨）** — 请明确授权可调用的 API 环境变量（例如 `ANTHROPIC_API_KEY`、`OPENAI_API_KEY`；若扩展 3×3/4×4 则加 `GEMINI_API_KEY`、`DEEPSEEK_API_KEY`）。
- [ ] **仅批准 Track B（AgentHub 订阅轨，零新增 API 预算）** — 若选此项，需人类允许执行前由工程人员对 AgentHub 启动临时沙箱配置：禁用 `mcp__agenthub` 工具与文件写入权限、每轮 `newThread: true`、隔离空工作目录、并按订阅限流节流（预计 `6,260` 次串行需 26–52 小时，4 路并行约 8–15 小时）。
- [ ] **批准 Track A+B 混合**（例如先跑 Pilot `M=10` 检查格式与方差，再由人类二次确认正式批预算）。

---

## 决策 2：锁定参测厂商组合与具体模型快照（Model Selection）

### 2.1 矩阵规模选择
- [ ] **核心 $2 \times 2$ 矩阵（推荐首选）**：集中全部统计功效于两家主导基线厂商（如 **Anthropic Claude** $\times$ **OpenAI GPT/Codex**），与 `CrossAudit`（`arXiv:2608.28631 §4.2`）和 `Dong (2026b §3)` 直接正面对话。
- [ ] **扩展 $3 \times 3$ 矩阵**：加入 **Google Gemini**（`Claude × GPT × Gemini`）。
- [ ] **全本机四厂商 $4 \times 4$ 矩阵**：`Claude × GPT × Gemini × DeepSeek`（注意网格规模从 4 格增至 16 格，总调用量增至约 3 倍）。

### 2.2 能力档位与精确模型快照（必须在启动前锁定，防止能力跨代混杂）
请人类从以下两档中选择每家厂商的**具体冻结模型 ID**（建议生成器与审计器同代对齐，例如全部使用旗舰级或全部使用高性价比级；**禁止**一家用旗舰而另一家用轻量版）：

| 厂商 | 经济/高效档（~$32–$47 跑完推荐 $2\times 2$） | 旗舰/主确认档（~$191–$474 跑完推荐/高功效 $2\times 2$） | 人类最终锁定快照 ID（待填） |
|---|---|---|---|
| **Anthropic (Claude)** | e.g. `claude-haiku-4-5-20251001`（或当前最新 Haiku） | e.g. `claude-sonnet-4-6` / `claude-opus-4-x` | `____________________` |
| **OpenAI (GPT/Codex)** | e.g. `gpt-5-mini` / `gpt-5.4-mini` | e.g. `gpt-5.4` / `gpt-5.6-terra` / `gpt-6-astra`（按账号实际可用最新快照） | `____________________` |
| **Google (Gemini)** *(若选 3×3/4×4)* | e.g. `gemini-2.5-flash` / `gemini-3.8-flash` | e.g. `gemini-2.5-pro` / `gemini-3.1-pro` | `____________________` |
| **DeepSeek** *(若选 4×4)* | `deepseek-chat` / `deepseek-reasoner` | `deepseek-reasoner` (最新可用快照) | `____________________` |

> **关键说明**：`SciCoQA` 原论文发现 `GPT-5` 家族在其自建合成集上有偏高召回（`arXiv:2601.12910 §7`），本设计已**剔除** `SciCoQA` GPT-5 合成 diff，仅采用第三方确定性脚本植入，确保任何选用模型均可公平对比。

---

## 决策 3：选择样本量功效档位与硬预算上限（Budget Cap）

以下调用次数与成本直接追溯自 `outputs/power_simulation.json`、`outputs/post_repair_power.json` 与 `outputs/cost_estimates.json`：

| 档位 | 配置 $(V, M, C, D, K)$ | 功效摘要（$\Delta=0.10$ 时 `pAUROC20` / `AUROC` / 修后通过率 $E_4$） | 总 LLM 调用次数 | 估算 Token (In / Out) | 经济档估算 | 中等旗舰估算 | 高阶旗舰估算 |
|---|---|---|---:|---|---:|---:|---:|
| 档位 1：`2x2_Standard` | $V=2, M=60, C=1, D=2, K=4$ | `36.4% / 55.6% / 83.3%`（`TPR10` 仅 `20.4%`） | `3,526` | 17.17 / 3.77 M | $18.01 | $108.05 | $180.08 |
| **档位 2：`2x2_Recommended`（推荐平衡档）** | **$V=2, M=80, C=2, D=2, K=4$** | **`63.6% / 85.2% / 92.6%`**（$\Delta=0.12$ 时 `pAUC20` `73.6%`，$\Delta=0.15$ 时 `91.2%`） | **`6,260`** | **29.97 / 6.75 M** | **$31.86** | **$191.16** | **$318.60** |
| **档位 3：`2x2_HighPower`（顶刊高功效档）** | **$V=2, M=120, C=2, D=2, K=4$** | **`81.2% / 95.2% / 97.3%`**（即使 $\Delta=0.10$ 仍破 80% 功效线） | **`9,318`** | **44.66 / 10.04 M** | **$47.43** | **$284.58** | **$474.30** |
| 档位 4：`3x3_Extension` | $V=3, M=80, C=2, D=2, K=3$ | 三厂商交叉对称验证 | `10,800` | 52.02 / 11.70 M | $55.26 | $331.56 | $552.60 |
| 档位 5：`4x4_AllLocalVendors` | $V=4, M=80, C=2, D=2, K=3$ | 全机四厂商对称验证 | `18,900` | 91.73 / 20.35 M | $96.74 | $580.44 | $967.40 |

### 请人类决定：
1. **选定样本量档位**：`____________________`（推荐 `2x2_Recommended` 或 `2x2_HighPower`）。
2. **硬预算上限（Hard Cap）**：最多允许消耗 **USD `$________`**（或订阅消息上限 `________` 条）；执行脚本须在达到 `50%` 与 `90%` 预算时输出账单汇总，超限立即自动停机。
3. **Pilot 放行授权**：是否批准先拨付 Pilot（`M=10` 题，约 `190` 次调用，中阶旗舰约 `$5.80` / 经济档约 `$0.97`）用于校准检验？ [ ] 是 / [ ] 否

---

## 决策 4：公开数据下载授权（需指派 `A-T1` 登记 `DATA_SOURCES.md`）

按 `P2_RULES.md` 第 1、5 条，所有公开数据集在下载前必须由 `A-T1` 在 `round9/p2/shared/DATA_SOURCES.md` 登记 URL、许可证、大小与 SHA-256，全队总量不超过 60 GB，**严禁接触 X-Atlas / Orion 及其衍生数据**：

| 数据集 | 许可证 | 大小 | 用途 | 人类是否批准下载 |
|---|---|---|---|---|
| `UKPLab/scicoqa` | `CC BY 4.0` | `3.49 MiB` (`3,657,481 bytes`) | S3 干净多文件科研代码库 + `Real` 92 例真实差异对比层 | [ ] 批准 |
| `shhu2001/SciCode-Verified` (`problems_test.jsonl` + `test_data_cleaned.h5`) | `Apache-2.0` | `951.75 KB` + `1.03 GiB` | S1 `64` 题 / `290` 子步可执行科研代码真值（若只做轻量方案，可讨论仅用能自包含执行的子题或由人类批准 1.03 GiB 下载） | [ ] 批准全量 (`1.03 GiB`) / [ ] 仅 JSONL |
| `mbzuai-ai4bio/BioXArena-Pro-Public` 非 `xatlas` 域（`chemical-biology` 8 + `phenotype-disease` 8） | `CC BY-NC 4.0`（数据）+ `MIT`（本地 `grade.py`） | `587.24 MB` (`273.55 + 313.70 MB`) | S2 `16` 项非 `xatlas` 生物建模 CPU 任务（本地生物单元切分，**绝不提交 Codabench，绝不碰 `xatlas-orion-perturb`**） | [ ] 批准 |

> 注：以上全部批准合计约 **1.62 GB**，远低于 `60 GB` 全队红线。

---

## 决策 5：预注册提交锁定与启动指令

1. **在任何付费模型调用执行前**，请主控指派唯一 git 写入者（`对话 19f7514f`）按 `p2/status/COMMIT_REQUEST_A1.txt` 提交冻结 `round9/p2/metasci_2x2/` 下的 `PREREG_2x2.md`、`CONSTITUTION_2x2.md`、`audit_response_schema.json`、`matched_fpr_eval.py` 与仿真脚本/输出。
2. 人类填写完本清单中「执行轨（决策 1）」「模型快照（决策 2）」「档位与硬预算上限（决策 3）」「数据下载授权（决策 4）」后，方可执行 Pilot 与正式批脚本。
