# A-T3 独立预注册终审报告：`PREREG_I.md` (`v5-FINAL-PRE-LOCK`)

- **审计智能体**：`A-T3`（A 队 · 预注册与统计审计）
- **审计时间**：`2026-11-26 00:20 +08:00`
- **审计对象（仅正文 + 洁净室 `round16/audit/`；严格未读 `round16/src/`）**：
  - [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md)（`284` 行，`v5-FINAL-PRE-LOCK`，SHA-256 `43d242e460fd64db2e4af1e628e886cf1d94fb4ccd1e95925928a99cb96ee4a8`）
  - [`round16/audit/`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/)（[`run_audit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py), [`audit_r16_state_adapter.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py), [`audit_r16_raw_controls.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py), [`audit_r16_sampling.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_sampling.py), [`audit_r16_metrics.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_metrics.py), [`audit_r16_adjudicate.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_adjudicate.py)）
- **洁净室与合规声明**：严格遵守仅读 `PREREG_I.md` 正文与自有 `round16/audit/` 源码；未读取或导入 `round16/src/`；未下载 `GSE295214` 或 `adata_pred.h5ad`；未读取或推断任何真实单细胞预测/指标；未修改 `round16/prereg/`、`round16/audit/`、`round16/src/`、`paper/` 或 Rounds 9–15 冻结文件。

---

## 0. 总体审计结论摘要

| 检查维度 | 总结论分类 | 关键计数 | 说明 |
|---|---|---|---|
| **1. `v5` 内部自洽性（§2 常量 / §3 裁决 / Group 3 `GSE295214`）** | **无阻塞性冲突；Gate 3 可按现稿锁库** | `0` 项「锁定前必须改」；`2` 项「可以锁定后作为补记」；`3` 项「无问题」 | `§2.1–§2.3` 与 `§3.1–§3.3` 核心数学/常量/分层/`ci975 == 0.0` 并列优先级完全闭环自洽；`§2.4` Group 3 无数据窥探且 `build_gse295214_parts` 在 Gate 4 下载前冻结。 |
| **2. `PREREG_I v5` vs. 现有 `round16/audit/` 一致性** | **须在 Gate 3 锁库后按 `v5` 对齐代码（本次未改）** | `7` 项对齐清单（含 `2` 项语义级、`5` 项常量/接口/新增功能） | 现有 `round16/audit/` 仍含 Gate-2 前默认与 `v4` 命名（如 `ST-HVG` 默认 `libnorm1e4_log1p`、`median_ntc_umi` 取全量 NTC vs 抽样 512、`SE-600M` 启用 `HGNCBridge` 回退 vs `v5 L81` 精确匹配、Group 3 尚未实现等）；均须在锁库后、Gate 4 冻结前统一更新。 |
| **3. 残留“看了数据才能决定 / 视情况”漏洞检查** | **无确证性 `H-I1`/`H-I2` 后验漏洞；存在次要/校验条款宽松度** | `0` 项「锁定前必须改」；`5` 项「可以锁定后作为补记」（均限于 `§1.3`/`§2.4` 敏感性/`§4`）；`2` 项「无问题」 | 主裁决（`H-I1`, `H-I2`, Group 1–3 核心池与门控）无主观自由度；仅 `§1.3 adata_pred.h5ad` 检验细胞系与阈值、`§2.4.4` 次要 `full r9_6` 靶点集、`§4.3 Perm100`/`CellularAI`、`§4.4 G4` 配对方式待在 Gate 4 补记/代码冻结中固定。 |

---

## 1. 检查一：`PREREG_I v5` 内部自洽性（§2 常量、§3 裁决优先级与平局规则、Group 3 `GSE295214`）

### 1.1 `§2` 锁定常量与 `§3` 裁决规则交叉核对 — **【无问题】**
1. **核心评估规模与合格线数一致**：
   - `§2.1.1`（[`PREREG_I.md#L65-L67`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L65-L67)）：`|G_model| = 2,000`（`0` HGNC 冲突），`|G_eval^R16|` 锁定为 `nadig_hepg2=1,843`、`nadig_jurkat=1,843`、`a549_ko=1,837`、`dld1=1,832`、`cd4t_rest=1,830`、`feng_ipsc=1,575`。
   - `§2.1.2`（[`PREREG_I.md#L70-L74`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L70-L74)）、`§2.3`（[`PREREG_I.md#L105-L114`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L105-L114)）与 `§3.3`（[`PREREG_I.md#L218-L226`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L218-L226)）三处核对完全一致：
     - Group 1（$K_1 = 2$）：`nadig_hepg2 = 740`（`740 hv`; `118 masked`）、`nadig_jurkat = 873`（`873 hv`; `137 masked`）；
     - Group 2（$K_2 = 3 \ge 2$：$K_{2,\mathrm{ko}} = 1$，$K_{2,\mathrm{i}} = 2$）：`a549_ko = 403`（`403 hv`; `75 masked`）、`dld1 = 75`（`75 hv`; `7 masked`）、`cd4t_rest = 245`（`230 hv`; `46 masked`）；退化保护 `|G_2^qual| < 2` 按预定不触发；
     - `Q4` 降级与保留控制：`feng_ipsc = 88`（`88 hv`; `13 masked`，`|G| = 1,575`）严格降级为 `SENSITIVITY_ONLY` 靶线，同时作为 `unseen_by_checkpoint` 异源对照保留于 $\mathcal{U}_{\mathrm{local}} = \{\text{a549\_ko}, \text{dld1}, \text{cd4t\_rest}, \text{feng\_ipsc}\}$；`§0.2`、[`#L74`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L74)、[`#L105-L114`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L105-L114)、[`#L178`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L178) 与 [`#L222`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L222) 无任何矛盾。
     - 检验力表 `§3.3` 中的 $\delta / \mathrm{FCT}$ 与剩余安全裕度经逐项算术复核均自洽（如 `a549_ko`: `0.04/0.1898 = 21.1%`、`0.04/0.0796 = 50.3%`；`dld1`: `28.6%–29.5%`；`cd4t_rest`: `71.9%`；`feng_ipsc`: `33.0%`；`vcc2025_h1`: `0.04 - 0.0299 = 0.0101`，约 `0.66 SE`；`hek293_xucao`: `h_{0.95} > 0.04`，TOST 检验力 `0%`）。
2. **两分层匹配池与 4 检查点重叠披露自洽**：
   - Group 1 每线 `seen = 1`（对侧 Nadig 系）、`unseen = 4`（$\mathcal{U}_{\mathrm{local}}$）；Group 2 每线 `seen = 2`（Nadig 双系）、`unseen = 3`（$\mathcal{U}_{\mathrm{local}} \setminus \{t\}$）。`§3.1` [`#L176-L178`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L176-L178) 与 `§2.3` 表严格一致；`§0.2` [`#L28`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L28) 与 `§4.3.2` [`#L245`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L245) 对 4-ckpt 平均中的 `3/4` 训练重叠披露及单折 LOO 交叉核对也相互呼应。

### 1.2 并列规则（`ci975 == 0.0`）与 4-Way TOST / 线级 8 档优先级 — **【无问题】**
经形式化检查 [`PREREG_I.md#L180-L209`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L180-L209)：
1. **单对 4 分类优先级（6 步级联）互斥且完备**：
   - `Step 1`: `valid_frac < 0.95` 或 `97.5%` CI 非有限 $\implies$ `INCONCLUSIVE`（且强制 `within_delta = False`）。
   - `Step 2`: `ci975_lo == 0.0` 或 `ci975_hi == 0.0` $\implies$ `INCONCLUSIVE`（先于 `MATCHED_BETTER`/`MISMATCHED_BETTER`/`EQUIVALENT`，哪怕 `within_delta == True` 如 `[0.0, 0.03]` 亦严格归入 `INCONCLUSIVE`，完全消除了闭区间/半开区间边界冲突）。
   - `Step 3/4`: `ci975_lo > 0.0`（`MATCHED_BETTER`）与 `ci975_hi < 0.0`（`MISMATCHED_BETTER`）在百分位自举单调性 $q_{1.25\%} \le q_{98.75\%}$ 下不可能重叠。
   - `Step 5/6`: 当 `97.5%` CI 严格跨越 `0.0` 时，`within_delta == True`（要求 `valid_frac >= 0.95` 且 `(ci95_lo > -0.04 and ci95_hi < 0.04)` 开区间严格成立）$\implies$ `EQUIVALENT`，其余（包括精确贴着 `±0.04` 边界的 `ci95_lo == -0.04` 或 `ci95_hi == +0.04`）$\implies$ `INCONCLUSIVE`。
2. **线级 `H-I1` 8 档裁决与 `H-I2` 3 档裁决无歧义**：
   - `matched_arm_valid`（含 `rho_ab_orth > 0.0`）在第 2 优先级先于 `LINE_EQUIVALENT_NULL_SIGNAL`（第 6）与 `LINE_MISMATCHED_BETTER`（第 7），防止 `rho_ab_orth <= 0` 或低有效自举比例误产出 `EQUIVALENT` 或 `MISMATCHED_BETTER`。
   - `matched_positive_ci975 == False` 时即便所有对比均 `MATCHED_BETTER`（例如 matched=0, mismatched=-0.08），因第 3/4 档要求 `matched_positive_ci975 == True`，也会安全落入第 8 档 `MIXED_OR_INCONCLUSIVE`。
   - `H-I2`（[`#L207-L209`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L207-L209)）三分类对称，`0.0` 边界（`ci975_lo <= 0.0 <= ci975_hi`）自动归为 `INDISTINGUISHABLE`。

### 1.3 Group 3（`GSE295214`）读取、分半、准入、`FCT`/`STATE` 对照池与单独报告机制检查
`§2.4`（[`PREREG_I.md#L125-L160`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L125-L160)）已锁定下载审批门（#L129）、`38,606` RNA 整数校验与 `6,638` Guide 的 `fit_poisson_gaussian_thresholds_from_hist(..., min_thr=5, max_thr=20, default_thr=7)` 单导过滤（#L132-133）、NTC=`non-targeting_*`（372 条，#L134）、QC `n_umis >= 300 & n_genes >= 200 & max_pct_mt=None`（#L135）、85 个 GSM 通道按字典序奇/偶分入 Half A / Half B 及组内 `< 15` NTC 中心化回退（`min_batch_nt=15`, 禁止跨半分享，#L136）、`16×32` 分层抽样（#L138）、先计数后算指标的 `R16_GSE295214_QC_COUNTS.md` 门控（`N_eval >= 50`, `|G_eval^R9| >= 3,000`, `per-half NTC >= 15`，#L140-144）、主 `FCT` 对照池排除 `vcc2025_h1`（`r9_6 \ {vcc2025_h1}` 5 系，#L147）、`zeroshot/ensemble4` + `seen (n=2)`/`unseen (n=4)` 分层（#L150-152）、独立数据冻结 `FREEZE_MANIFEST_R16_DATA.txt`（#L159）与强制 `H1`/Replogle-Weissman 披露（#L156）。

尽管主体逻辑非常完整，审计发现 **2 个次要实现细则未在 `§2.4` 文字中显式写死（均不构成锁前阻塞，原因：下载仅在 Gate 3 锁库 + Gate 4 `build_gse295214_parts` 代码冻结且获人工下载批准之后才发生）**：
1. **【可以锁定后作为补记】`§2.4.4 L148` 次要敏感性 `full r9_6`（含 `vcc2025_h1`）是否固定沿用主 `N_eval^R16` 靶点集**：
   - `§2.4.3`（[`#L141`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L141)）定义 `N_eval^R16(gse295214_h1)` 须在 `r9_6 \ {vcc2025_h1}` 至少 1 个源中可重复；而 `§2.4.4`（[`#L148`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L148)）增加“完整 `r9_6`（含 `vcc2025_h1`）作为次要谱系重叠敏感性”时，未写明该敏感性是**严格固定在主池锁定的同一个 `N_eval^R16` 靶点子集上**（仅让 `vcc2025_h1` 参与 `FCT` 权重与 `Transfer_Mean` 拟合，保证 `H-I2` 严格配对可比），还是重新与 `vcc2025_h1` 扩张求交。建议在 Gate 4 补记中固定为**严格冻结在主分析同靶点集 `N_eval^R16` 上**。
2. **【可以锁定后作为补记】`§2.4.2 L133-L137` Guide-to-Target 解析规则、重名基因去重及 `P5` 数值阈值**：
   - `§2.4.2`（[`#L133-L137`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L133-L137)）明确写了 `non-targeting_*`（372 条），但未把剩余靶向 guide（6,266 条）如何解析为基因名（分隔符 `_`/`-`、`-TSS*`/安全靶点/双靶点排除规则）、`38,606` 个 10x RNA feature 若存在同名 gene symbol 时如何择一保留（如按最大 NTC 均值或 `select_canonical_genes`）、以及与 `PREREG_H` 一致的 `P5` 阈值（最低单扰动细胞数/Energy `p`/AD DE 基因数）直接内联写出，而是依赖 Gate 4 前冻结的 `build_gse295214_parts`（[`#L137`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L137)）。因 Gate 4 冻结严格早于下载，不会引入数据窥探；建议在 Gate 4 代码冻结与测试文档中把上述解析与去重确定性规则写明。
3. **【无问题】其余 `§2.4` 核心规则（下载前置门、`r9_6 \ {vcc2025_h1}` 主源池、`85` GSM 奇偶分半、`min_batch_nt = 15` 禁止跨半分享、只计数的 `R16_GSE295214_QC_COUNTS.md` 门控、不与 Group 1/2 合并及 `H1`/Replogle-Weissman 双重披露）**：均可不经修改直接按 `v5` 执行。

---

## 2. 检查二：`PREREG_I v5` 与现有 `round16/audit/` 实现的逐项差异清单

> **严格遵守“以 `v5` 为准，锁定后再改代码”**：本次审计未修改 `round16/audit/*.py`。下表逐一列出 Gate 3 `PREREG_I_COMMIT.txt` 锁库后、`CODE_FREEZE_R16.txt` 前须在 `round16/audit/` 中修复对齐的所有差异（均属于 Gate 4 代码对齐必改项）。

| # | `round16/audit/` 文件与行号 | 当前实现状态 | `PREREG_I v5` 锁定要求 | 锁后修改建议与风险 |
|---|---|---|---|---|
| **A1** | [`audit_r16_raw_controls.py#L128`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L128), [`#L228`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L228), [`#L411`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L411), [`#L525`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L525), [`#L593`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L593), [`#L725`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L725); [`run_audit_r16.py#L448`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py#L448) | **按全库所有通过 QC 的 NTC 细胞**计算 `median_ntc_umi`（如 `np.median(all_umis[all_umis > 0.0])`），并作为 `target_median_umi` 覆盖 `preprocess_hvg_expression` | [`PREREG_I.md §2.1.4#L83`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L83)：`sc.pp.normalize_total(adata)`（**每个被抽样对照池**全转录组总 UMI 的中位数，`target_sum=None`）$\to$ `sc.pp.log1p` | **高一致性影响（`STATE_HVG`）**：锁库后须改用抽样 `512` 细胞自身全转录组 `cell_total_umis` 的非零中位数（不传入全库 NTC 中位数覆盖），并与 `round16/src/` 三方对照核实 `scanpy.pp.normalize_total` 一致性。 |
| **A2** | [`audit_r16_state_adapter.py#L394-L451`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L394-L451), [`#L489-L495`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L489-L495), [`#L663-L676`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L663-L676), [`#L986-L992`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L986-L992) | `map_genes_to_se600m_pe` 先做精确匹配，再用 `HGNCBridge` 补救旧/别名基因（且 [`align_hvg_expression_to_model#L303-L352`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L303-L352) 亦含 `HGNCBridge` 补救） | [`PREREG_I.md §2.1.4#L81`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L81) 与 [`§2.1.5#L89`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L89) 明确规定 `SE-600M` 使用 **`exact gene-symbol string lookup` 映射到 `protein_embeddings.pt`**，对应精确映射下界（`9,022 / 8,391 / 13,766 / 11,714 / 11,535 / 7,360`） | **高一致性影响（主模型 `STATE_SE`）**：锁库后须按 `v5 L81` 禁用 `SE-600M` 的 `HGNCBridge` 补救（仅保留精确字符串查表及确定性重名去重）；同时在 Gate 4 与主线核对 `ST-HVG` [`§2.1.4#L83`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L83) 的 `2,000` 列提取是否也应采用精确匹配而非 `HGNCBridge`。 |
| **A3** | [`audit_r16_state_adapter.py#L183`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L183), [`#L694`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L694) | `preprocess_hvg_expression` 与 `CleanRoomArcStateAdapter.__init__` 默认参数为 `hvg_preproc_mode="libnorm1e4_log1p"`（仅 `run_audit_r16.py` 传入 `"official_median_log1p"`） | [`§1.2#L47`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L47), [`§2.1.4#L83`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L83), [`§4.3.1#L244`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L244)：`official_median_log1p` 为主，`libnorm(1e4)` 为敏感性 | 锁库后将两处函数/类默认值直接改为 `"official_median_log1p"`。 |
| **A4** | [`run_audit_r16.py#L78-L112`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py#L78-L112), [`#L319-L332`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py#L319-L332)；及各模块 docstring | `prereg_version = "v4-GATE2-PENDING-V5"`，`locked_n_eval_exact = None`，模块头注释仍写 `PREREG_I.md v4`；`expected_n_eval_r16_min` 含 `feng_ipsc: 50` 却未写入 `740, 873, 403, 75, 245, 88` 精确锁死与 `half_valid` 计数 | [`§2.3#L109-L118`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L109-L118)：`v5-FINAL-PRE-LOCK`，精确锁定 6 线 `N_eval^R16`、`half_valid` 及 4 条局部敏感性系（`cd4t_stim48=401`, `a549_ko_cas9pos=403`, `a549_myllia=20`, `hek293_xucao=16`） | 锁库后填入 `locked_n_eval_exact` 与 `half_valid` 校验，更新版本标识与注释为 `v5`。 |
| **A5** | [`audit_r16_adjudicate.py#L307-L308`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_adjudicate.py#L307-L308), [`#L270-L321`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_adjudicate.py#L270-L321) | 退化分支字面量为 `"DOWNGRADED_SINGLE_STRATUM_DESCRIPTIVE_ONLY"`；`adjudicate_group_r16` 尚无独立的 `group3_prospective` 输出段 | [`§3.1#L203`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L203)：字面量为 `SINGLE_STRATUM_DESCRIPTIVE_ONLY`；[`§2.4.5#L154-L156`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L154-L156) & [`§3.1#L202`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L202)：Group 3 严格独立汇报并附带 `H1` 与 Replogle/Weissman 披露 | 锁库后将退化状态字面量改回 `SINGLE_STRATUM_DESCRIPTIVE_ONLY` 并增加 `group3_prospective` 独立汇总分支。 |
| **A6** | [`audit_r16_raw_controls.py#L104-L115`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L104-L115), [`#L850-L854`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_raw_controls.py#L850-L854); [`audit_r16_sampling.py#L37-L44`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_sampling.py#L37-L44) | 尚无 `gse295214_h1`（`85` GSM `.h5`）读取/NTC 抽样/`QC_COUNTS` 门控实现；且 `load_line_eval_reference_r16` 对所有非 Nadig 线固定使用全部 6 系 `r9_6` | [`§2.4#L131-L159`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L131-L159)：须在下载前实现并通过合成 `.h5` 单元测试的 `GSE295214` 读取器；主 `FCT` 池强制为 `r9_6 \ {vcc2025_h1}`（5 系），`full r9_6` 仅作次要敏感性 | Gate 3 锁库后至下载批准前在 `round16/audit/` 补齐合成单测验证的 Group 3 模块及 `r9_6 \ {vcc2025_h1}` 主参考构筑。 |
| **A7** | [`run_audit_r16.py#L74`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py#L74), [`#L505`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/run_audit_r16.py#L505); [`audit_r16_state_adapter.py#L774-L788`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L774-L788), [`#L869`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L869), [`#L1086-L1126`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/audit_r16_state_adapter.py#L1086-L1126) | (1) 仅计算 `2,058` 维**中心化**质心余弦+L2，缺 **`2,048` 维**与非中心化余弦；<br>(2) `custom_batch_idx` 仅为标量 `int`，尚不支持逐细胞 `gem_group` 向量映射（Group 1 `G4` 2×2 与 `§1.3` `(M, B_M)` 必需）；<br>(3) 输出为 `audit/results/{model}_{ckpt}__{line}.npz`，尚缺对主线 `round16/results/arrays/{observed,preds_*}.npz` 导出/重算、Group 2 分折分解、`ST-HVG libnorm+Perm100`、`G4` 面板、实用精度指标（Pearson/DE-50/cosine）及 `adata_pred.h5ad` 检验的调度接口 | [`§1.3#L55`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L55), [`§2.1.5#L89`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L89), [`§4.3–§4.4#L243-L266`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L243-L266), [`§5.5#L281-L282`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L281-L282) | 锁库后逐项补全 `audit/` 适配器逐细胞 `batch` 向量、`2,048`+`2,058` 质心诊断、主线 `results/arrays/*.npz` 独立重算器及全套次要诊断脚本；同时核对 `bypass_zero`（设 `batch_encoder=None` vs 显式输入 `0 ∈ R^128` 向量）与主线 `round16/src/` 数值对齐。 |

---

## 3. 检查三：残留“看了数据才能决定 / Gate 5 时再定 / 视情况”漏洞检查

对 [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) 全文逐条筛查可能引入后验选择自由度的条件条款：

1. **【可以锁定后作为补记】`§1.3 L55` & `§5.5 L282`：Arc 归档 `adata_pred.h5ad` 验证对象与通过阈值未锁死**
   - 原文写 *"at least one Group 1 primary line (`nadig_hepg2` or `nadig_jurkat`) archived `adata_pred.h5ad`"*，既允许 Gate 5 后在 `hepg2` 与 `jurkat` 中任选其一，也未规定如何处理 Arc 原生 eval 对照细胞集合与我方 `512` 细胞抽样集合之间的差异及何种偏离算校验失败。
   - **评价与建议**：此项属于适配器工程保真核查，不影响 `H-I1`/`H-I2` 统计裁决，可不阻断 Gate 3 锁库；建议锁定后在 Gate 5 前补记明确：**两条 Group 1 主线均检或至少锁定先检 `nadig_hepg2`**，并用同一对照细胞子集（或规定相关性下界）作为通过判据。

2. **【可以锁定后作为补记】`§2.4.4 L148` & `§2.4.2 L133-L137`：Group 3 次要 `full r9_6` 靶点边界与 guide/重名解析**
   - 见第 1.3 节：因严格位于下载前（Gate 4 冻结 `build_gse295214_parts`），无看数据再调的可能；建议在锁库后、下载批准前补记固定“次要 `full r9_6` 严格锁定在与主池相同的 `N_eval^R16` 靶点集上评估”。

3. **【可以锁定后作为补记】`§4.3.1 L244` & `§2.1.4 L83`：`STATE_HVG` 的 `Perm100` 置换单元与 `G_model` 列提取规则**
   - `Perm100`（[`#L47`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L47), [`#L123`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L123), [`#L244`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L244)）未写明置换对象（输入 `2,000` HVG 基因标签 vs 输出基因列）、随机种子及适用细胞系；`ST-HVG` `2,000` 列提取也未明示精确字符串 vs `HGNCBridge`（会影响 `AARS`/`DARS`/`H2AFZ`/`TAZ` 等 4 个符号）。因 `STATE_HVG` 属次要对照且 Gate 4 冻结前可统一三方代码，不构成主假设漏洞。

4. **【可以锁定后作为补记】`§4.3.3 L246`：`CellularAI_OpenCode_PublicRefit` 标注为 `"optional and non-blocking"`**
   - 存在“视后续进度决定做或不做”的操作自由度。由于被明确限制为仅描述性（[`#L49`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L49), [`#L246`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L246)），不进入 `H-I1`/`H-I2`、不触发多重性问题；但若最终未跑，须在状态档案中客观记录工程原因而不得选择性略去负面结果。

5. **【可以锁定后作为补记】`§4.4.3 L260-L265`：`G4` 2×2 分解的 `(B_M, B_S)` 配对对齐方式**
   - Group 1 未写明 `(M, B_S)` 中对方细胞系 512 个原生 `gem_group` 标签按抽样顺序一对一配对还是重抽；Group 2/3 未明文指定 `(b_1 = 0, b_2 = 17)` 中 `B_M := 0, B_S := 17` 方向以及跨多 `S` 汇总口径。属 `§4` 机理分解诊断，建议在 Gate 4 代码冻结中固定为“Group 1 按抽样后 0..511 位置对齐；Group 2/3 固定 `B_M=0, B_S=17` 并附 `(17,0)` 对称检查”。

6. **【无问题】`§2.4.1 L129-L130` Group 3 人类下载批准条件**：
   - 是否批准下载发生在**任何 `.h5` 下载和计数之前**；一旦获批，准入由纯计数门槛（`N_eval >= 50`, `|G_eval^R9| >= 3,000`, `NTC >= 15`）机械判定，无数据后验选择空间。

7. **【无问题】`§0.3–§0.4` & `§2.3 L120` 外部扩展边界**：
   - `VCC 2026` 伴随分析（目标 `< 2026-11-27` 锁库前）、`Occam`、以及 `replogle_k562_essential` / `replogle_rpe1` / `vcc2025_h1`（~45 GB）均被隔离在独立附录或单独人类批准之后，严禁混入本协议 Group 1/2/3 主判决。

---

## 4. `A-T3` 给主控的 Gate 3 锁库与 Gate 4 对齐操作建议

1. **Gate 3 锁库判决**：建议**直接锁定 `PREREG_I.md v5`（签发 `PREREG_I_COMMIT.txt`，满足 `< 2026-11-27` 时间要求）**；`v5` 核心确证性设计无「锁定前必须改」的统计漏洞或内部数字矛盾。
2. **锁定后立即行动（Gate 4 `CODE_FREEZE_R16.txt` 前）**：
   - 补记固定第 1.3 节与第 3 节列出的 5 项次要/接口规范（重点：`SE-600M` `v5 L81` 精确查表 vs `ST-HVG` `HGNCBridge`、`Perm100` 定义、`G4` `B_M/B_S` 配对、`full r9_6` 固定同靶点集、`adata_pred.h5ad` 检验线与阈值）；
   - `A-T3` 同步将 [`round16/audit/`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/audit/) 中上述 **A1–A7** 共 7 项对齐到 `v5` 与 Gate-4 冻结细则，完成后与主线 `round16/src/` 及主控三方完成纯合成数值一致性检验（`<= 1e-10`、`100%` verdict 一致）。
