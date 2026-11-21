# CONSTITUTION_2x2（冻结版科研代码审计宪法 v1.0）

遵循 `CrossAudit`（`arXiv:2608.28631 §3.3`）规则手册设计规范：每条规则具备**稳定编号**、**两级严重度**（`BLOCKER` 或 `ADVISORY`，遵守 `I6`）、**作用范围**、**可判定接受标准（decidable criterion）**与**依据来源**。审计器不得脱离规则自创 `BLOCKER`（`CA-META-003`），亦不得跳过可验证缺陷；所有规则在预注册提交后即冻结，实验周期内不得修改。

---

## 1. 元规则（Audit Meta-Rules，遵循 CrossAudit I3 / I6）

| Rule ID | Severity | Scope | Decidable Criterion | Rationale |
|---|---|---|---|---|
| `CA-META-001` | `BLOCKER` | All audited increments | 缺失必需的证据/声明路径即视为发现，不得默认通过（`missing evidence is a finding, never an implicit pass`）。 | 防止沉默放行（`arXiv:2608.28631 §3.3`）。 |
| `CA-META-002` | `BLOCKER` | Audit output JSON | 输出必须包含有效 JSON、引用的 `constitution_commit` 与对应规则 ID；空引用/畸形输出由确定性控制器转交 `ESCALATE`（`I3`）。 | 杜绝无据空洞通过报告。 |
| `CA-META-003` | `BLOCKER` | Audit findings | 仅可针对明确标注 `BLOCKER` 的规则发出 `BLOCKER`；风格偏好或非判定性改进建议一律降级为 `ADVISORY`。 | 防止审计器自由裁量误报膨胀。 |
| `CA-META-004` | `ADVISORY` | Beyond-evidence queries | 超出当前代码与规范文本证据边界的问题应标为不可判定并升级/存疑，不得凭空猜测通过。 | 对应 `Dong 2026b §B.4` 未定规范边界。 |

---

## 2. 确定性检查层优先规则（DCL Tier-0，遵循 CrossAudit I4）

> `DCL` 脚本先于任何 LLM 运行；任一 `DCL` 硬失败直接触发 `BLOCKED`，任何模型不可豁免。

| Rule ID | Severity | Scope | Decidable Criterion | Source |
|---|---|---|---|---|
| `CA-DCL-001` | `BLOCKER` | Python/script syntax & imports | 代码可无语法错误解析（`ast.parse`）且所有声明依赖在冻结环境内可导入；禁止超出白名单的隐藏断言跳过或 `try/except: pass` 吞错误。 | `arXiv:2611.23250 §2` repair guard |
| `CA-DCL-002` | `BLOCKER` | Visible smoke checks & provenance | 若附带运行报告/日志，报告中的每个数值须匹配其声明日志/JSON 字节跨度（`number→source`），且不得违背基本物理/数学边界（如概率 $\notin [0,1]$、NaN/Inf 输出）。 | `arXiv:2611.23250 §2, §4.2`; `AUDIT_LEDGER.csv` #4, #5, #10, #13 |

---

## 3. 五类科学实质缺陷规则（模型审计核心域，对齐 SciCoQA 与 AUDIT_LEDGER）

| Rule ID | Severity | Defect Family | Decidable Criterion | Grounding Source |
|---|---|---|---|---|
| `CA-SCI-001` | `BLOCKER` | `ALGO_MATH_LOGIC` | 实现的数学公式、范数/符号、积分/迭代边界、矩阵维数或算子顺序与规范/论文描述存在实质性逻辑矛盾，导致在合法输入上输出与规范预期不一致。 | `SciCoQA Table 1` (`Difference` / `Algorithm`); `SciCode-Verified` |
| `CA-STAT-001` | `BLOCKER` | `STAT_LEAK_METRIC` | 存在训练–验证/测试划分泄漏（如全集标准化/特征选择后再切分、跨分组样本混合）、评估口径篡改（如擅自 `max(0, r)` 截断或使用非对称基线）、或统计检验独立性假设违背。 | `AUDIT_LEDGER.csv` #3, #11, #14, #17; `BioXArena-Pro` |
| `CA-SPEC-001` | `BLOCKER` | `SPEC_CODE_MISMATCH` | 规范/论文中明确要求的关键步骤、激活/正则项、条件分支或收敛判据在代码中缺失（Code Omission），或代码加入了改变科学行为的未声明组件（Paper Omission）。 | `SciCoQA Table 1` (`Code Omission` / `Paper Omission`) |
| `CA-SCI-002` | `BLOCKER` | `NUM_PROVENANCE_FABRICATION` | 代码或计算链路中存在硬编码伪造数值、跨环节一致性的量纲/比例偏差（如 $\times 10^3$ 或常数偏移未与物理/数学推导对齐）、或报告声明步骤在代码中并未实际执行。 | `arXiv:2611.23250 §4.2`; `AUDIT_LEDGER.csv` #5, #10 |
| `CA-DATA-001` | `BLOCKER` | `DATA_UNIT_SCHEMA` | 输入/输出变量单位混用、列映射错位、未处理哨兵缺失值（如 `-999`）或未按文档声明的量纲/容差约束处理边界数据。 | `arXiv:2611.23250 §4.3`; `CrossAudit §4.4` |
| `CA-SCI-003` | `ADVISORY` | Style / readability / minor optimization | 代码命名风格、非语义性重构、或未改变规范定义科学输出的性能建议。**不得用于 BLOCK。** | `CrossAudit §3.2` (`I6`) |

---

## 4. 连续评分校准说明（Matched-FPR Calibration Contract）

除遵循上述二值 `BLOCKER`/`ADVISORY` 规则外，审计器必须在每次输出中给出连续阻塞风险分 `overall_block_score` $\in [0, 100]$：
- `0–19`：未发现违背上述 `BLOCKER` 规则的证据（仅含零条或多条 `ADVISORY`）；
- `20–49`：存在可疑边界线索或局部证据不足，但尚不足以确认违反特定规则；
- `50–79`：高度怀疑违反至少一条 `BLOCKER` 规则并能指明文件/符号位置；
- `80–100`：确证违反至少一条 `BLOCKER` 规则，具备可复现代码路径与数学/统计机制说明。

**关键纪律**：评测对比时，**不使用厂商自带默认 `verdict == BLOCKED` 阈值**（防止重演 `Dong 2026b` 中 3.3% vs 16.0% 误报率不可比问题），而是使用经 $K$ 次独立重读平均分 $\bar{s}(x)$ 上的 **Neyman-Pearson 确定性匹配误报率阈值** $\hat{\tau}_{g,a}^*(\alpha)$ 与 $\text{pAUROC}_{[0, 0.20]}$ 进行跨厂商比较。
