# Round 16 独立统计学与代码一致性审核报告 (`A-T3`)

- **审核智能体:** `A-T3`（A 队 · 预注册与统计审计）
- **审核时间:** `2026-11-24 +08:00`
- **主控参照:** [`PREREG_I_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_AUDIT_R1.md)（`I1–I14`，含 `11:15` v3 复核）与 [`MAIN.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/MAIN.txt)（`c4370cfc`）
- **审核对象及当前 SHA-256（只读核验）:**
  - [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md)（`v3-DRAFT`, `168` 行, `29,598 B`, `sha256: 3f4b0edfd88034ee6ab857043837a19373157d1b8e48aa97f8cf1a1ef5b1c755`）
  - [`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py)（`228` 行, `9,291 B`, `sha256: d6cd9daa3dc2667d814e0c1954e171ba53b9d12aca3e23ed77b477e0b841e2f7`）
  - [`test_unit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/test_unit_r16.py)（`342` 行, `17,373 B`, `sha256: c6464aaaa97dc2ee7885c7e13ea61f4dd9818971a4d22d1cb3bf6c29c32f8c45`；独立执行 `python3 -m pytest -v round16/src/test_unit_r16.py`：`8 passed in 4.30s`）
- **盲法与范围合规声明:** 未打开任何其余 `round16/src/*.py` 文件，未下载任何模型权重，未安装 `arc-state`，未克隆 `CellularAI`，未打开存档 `adata_pred.h5ad`，未对真实细胞执行推断，未改动 R9–R15 任何冻结文件，未运行 `git`。
- **总体审核结论:** **有条件通过（需先完成下述 6 项「必须改」并在 Gate 2 结构检查中补齐 `TBD_GATE2` 后，方可进入 Gate 3 锁定 `PREREG_I_COMMIT.txt` 与 Gate 4 代码冻结）**。`v3` 对 `I14` 的分层本地错配池改造及时且方向正确，Bonferroni 双论断与 TOST 配对逻辑主体成立，但代码与文本之间仍存在分半相关键名歧义、`within_delta` 未受 `valid_frac` 门控且与线级确认边界模糊、Gate 2 退化与组级汇总规则缺失、`512` 对照采样/4 检查点重叠/协变量处理规范不完整、以及低 $N$ 队列 TOST 检验效能未披露等问题。

---

## 1. 必须改（Gate 3 锁定与 Gate 4 冻结前硬性阻塞项）

### M1. 代码与预注册键名不一致：`adjudicate_r16.py` 检查 `"rho_ab"` 而非正交分半相关 $\rho^\perp_{AB}$ (`rho_ab_orth`)
- **证据定位:**
  - `PREREG_I.md` §3 L102、L112 明确规定匹配臂前提为 ``valid_frac >= 0.95`, $\rho^\perp_{AB} > 0$ ``；
  - [`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py#L160) 实际读取 `rho = float(orth_suite.get("rho_ab", float("nan")))`；
  - [`test_unit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/test_unit_r16.py#L165) `test_5` 构造夹具 `orth_pos = {"ci975_lo": 0.05, "valid_frac": 1.0, "rho_ab": 0.6}`。
- **统计风险:** 本项目 R10–R14 所有结果表均区分非正交分半相关 `rho_half_ab`（对应 raw 相关分母）与正交分半相关 `rho_half_ab_orth`（$\rho^\perp_{AB}$，对应 `r_orth_cal` 分母，例如 R12 开发复核中 `cd4t_rest` 在 `all_nontarget` 下 $\rho^\perp_{AB}=+0.147$ 远低于非正交分半相关）。键名写成 `"rho_ab"` 极易与非正交 $\rho_{AB}$ 混淆；若上游 `compute_orth_suite_r16` 实际输出 `"rho_ab_orth"`，`orth_suite.get("rho_ab")` 会静默返回 `NaN` 导致 `matched_positive=False`，而若输出非正交 $\rho_{AB}$ 则会放行正交分半不可复现的结果。
- **修改要求:** 在 `adjudicate_r16.py` L160、`test_unit_r16.py` L165 及 `PREREG_I.md` 中统一显式采用正交分半相关键名 **`rho_ab_orth`**（或在冻结前由 `test_7` 断言上游输出字典键与 $\rho^\perp_{AB}$ 一一对应），并在 `test_5` 新增负例验证 `rho_ab_orth <= 0` 或 `valid_frac < 0.95` 时必定阻断 `LINE_MATCHED_SPECIFIC`。

### M2. `within_delta` 的单双侧文本错位、未受 `valid_frac >= 0.95` 门控，且 `within_delta = True` 与 `LINE_MATCHED_SPECIFIC` 存在语义冲突
- **证据定位:**
  - `PREREG_I.md` §3 L105–106 将 `MATCHED_BETTER` / `MISMATCHED_BETTER` 的 `within_delta` 分别写成单侧条件（`ci95_hi < 0.04` / `ci95_lo > -0.04`），而 [`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py#L55-L70) 始终按双侧 TOST 区间计算 `tost_within = bool(finite_95 and lo95 > -float(delta) and hi95 < float(delta))` 并直接赋给 `"within_delta": tost_within`；
  - 更重要的是，`classify_pair_tost_r16` L57–70 在 `vf < MIN_BOOT_VALID_FRAC` (`0.95`) 或 `not finite_975` 时虽把 `category` 设为 `"INCONCLUSIVE"`，却**没有重置 `within_delta`**（即使有效 bootstrap 仅 10% 且 `ci975` 为 `NaN`，只要剩余样本偶然落在 `(-0.04, 0.04)` 内，函数仍返回 `within_delta=True`）；
  - 实质冲突：`PREREG_I.md` §3 L95 定义 $\delta = 0.04$ 为“**最小有实际意义的目标对照特异性效应量**”；但在现行 [`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py#L103-L104) 与 [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md#L112) 下，若大样本线（如 $N \approx 800\text{--}1000$）出现估计值 `+0.015`、`95% CI [+0.006, +0.024]`、`97.5% CI [+0.004, +0.026]`（即 `test_5` 的 `mb_small`），它既显著大于 0 又显著低于实际意义边界 $\delta=0.04$，代码却仍授予最高确认标签 **`LINE_MATCHED_SPECIFIC`**（v2 曾有 `(with ci95_hi >= 0.04 or reported alongside within_delta)` 占位注释，v3 直接删掉了且未定义线级区分）。
- **修改要求:**
  1. 在 `adjudicate_r16.py` 中强制 `within_delta = bool(vf >= MIN_BOOT_VALID_FRAC and finite_975 and tost_within)`；`PREREG_I.md` §3 L105–106 统一改写为精确双侧 `ci95_lo > -0.04 and ci95_hi < +0.04`（并说明有限有效区间嵌套性）。
  2. 在 `PREREG_I.md §3` 锁定前明确 `within_delta = True` 的线级处理：至少应在分层与线级摘要中显式区分 **`LINE_MATCHED_SPECIFIC`**（不仅统计显著 > 0，且各配对不满足亚阈值等价，或至少不存在全配对 `within_delta=True`）与 **`LINE_MATCHED_SPECIFIC_SUB_DELTA`**（统计显著 > 0 但效应明确小于 $\delta=0.04$），并在返回字典导出 `all_within_delta` / `any_within_delta`。

### M3. `LINE_EQUIVALENT_NO_CONTROL_USE` 未核验匹配臂有效性/零信号平庸等价，且缺组级汇总规则与 Gate 2 `G_2_qual` 退化保护
- **证据定位:** [`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py#L105-L110)（L185–188）与 `PREREG_I.md` §3 L113–114。
- **统计与逻辑缺陷:**
  1. **`EQUIVALENT` 的两种截然不同机制被混为一谈：** `STRATUM_EQUIVALENT` 和 `LINE_EQUIVALENT_NO_CONTROL_USE` 仅检查所有配对满足 `EQUIVALENT`，完全不检查匹配臂 `orth_suite` 的 `valid_frac >= 0.95`、`rho_ab_orth > 0` 以及 `matched_positive_ci975`！如果 `STATE_SE` 在该细胞系完全崩塌（`r_orth_cal ≈ 0` 或 `< 0`，或 $\rho^\perp_{AB} \le 0$），匹配与错配同样给出零信号，也会落入 `LINE_EQUIVALENT_NO_CONTROL_USE`；而真正的科学命题（呼应 Arc 2026 结项文）是“**模型具有非零超均值信号（`matched_positive_ci975 = True`），却不随目标对照改变（`EQUIVALENT`）**”。二者统计事实与机制含义完全不同，不可混用同一标签。同理 `LINE_MISMATCHED_BETTER` 与 `EQUIVALENT` 都应要求 `orth_suite` 自身 `valid_frac >= 0.95` 且分半/CI 有效。
  2. **Gate 2 交集缩水导致 `MISSING_STRATUM` 退化陷阱:** Group 2 候选池目前仅 4 条线（`a549_ko`, `dld1`, `cd4t_rest`, `feng_ipsc`），且尚未在 Gate 2 与 `ST-SE-Replogle` 的 `pert_onehot_map.pkl`（多为 Replogle 必需基因）求交集，也未扣除 `SE-600M` 预训练重叠线。一旦 Gate 2 后 $|\mathcal{G}_2^{\mathrm{qual}}| = 1$（例如仅 `a549_ko` 留下），则该唯一幸存线的 `unseen_by_checkpoint = G_2_qual \ {t} = ∅`，触发 `MISSING_STRATUM`（L115 / L181），直接导致**整个 Group 2 无法做任何线级判定（0/1）**；若 $|\mathcal{G}_2^{\mathrm{qual}}| = 0$，则 Group 1（`nadig_hepg2`, `nadig_jurkat`）的 `unseen` 也变空，**全研究 0 条线可判定**。另一方面，Group 1 本身 `seen_by_checkpoint` 仅有 `n_seen = 1`（`hepg2` 仅对比 `jurkat` 一条异系线，反之亦然），单对通过即判 `STRATUM_MATCHED_SPECIFIC`，脆弱度远高于多源池。
  3. **缺少组级/研究级判定规则:** §3 定义了“配对级”和“线级”结果并禁止跨线池化，但**未预注册 Group 1（`K_1=2`）、Group 2（`K_2<=4`）与 `a549_ko`（CRISPRko）vs CRISPRi 线的最终综合结论准则**（例如按各线计数 `k/K` 报告、是否要求全线一致、或 `SE-600M` 重叠后 $K_2$ 不足 2 时是否降级为探索性/分层描述性报告）。
- **修改要求:**
  - 在 `LINE_EQUIVALENT_NO_CONTROL_USE` 中至少拆分为（或并列输出布尔标志 distinguishing）**`LINE_EQUIVALENT_NONZERO_SIGNAL`**（`matched_positive_ci975=True` 且双层 `EQUIVALENT`）vs **`LINE_EQUIVALENT_NULL_SIGNAL`**（匹配臂 `r_orth_cal <= 0` 或不可区分于零，且双层 `EQUIVALENT`），并对 `orth_suite` 自身 `valid_frac < 0.95`（及 `EQUIVALENT` 分母不稳定如 `rho_ab_orth <= 0`）设为 `INCONCLUSIVE` 门控。
  - 在 §2.2 & §3 预注册 Gate 2 后 $|\mathcal{G}_2^{\mathrm{qual}}| \le 1$ 时的事前确定性兜底方案（例如：当 $|\mathcal{G}_2^{\mathrm{qual}}| \le 1$ 时，`unseen` 退化已证明主要线不足，明确写明是否改以本地既有完整原始计数且无泄漏的敏感性线——如仅用于错配控制对照源而不作主确认目标线的 `a549_myllia` 等——或事先声明直接降级为“单层描述性报告，不授予双层线级结论”），并明确写死 Group 1 (`K_1=2`, $n_{\mathrm{seen}}=1$) 与 Group 2 (`K_2<=4`, CRISPRko 与 CRISPRi 分别汇总) 的研究级报告准则。

### M4. TOST $\delta = 0.04$ 统计检验效能（Power）在中等/小 $N$ 队列存在严重分化——敏感性线甚至数学上不可能通过
- **基于本仓库 R11–R14 实测配对自助法（`n_boot=4,000`, `all_nontarget`, `r_orth_cal`）的经验效能审计:**
  令配对差 `95% CI` 半宽为 $h_{0.95} = (\mathrm{ci95\_hi} - \mathrm{ci95\_lo})/2$，近似标准误 $\widehat{\mathrm{SE}} \approx h_{0.95}/1.960$（当仅有 `97.5% CI` 时按 $h_{0.95} \approx h_{0.975}\times 1.960/2.2414$ 换算）。在真零假设 $\Delta=0$ 下，单对获得 `EQUIVALENT` 的接受域为 $|\hat{\Delta}| < 0.04 - h_{0.95}$；一旦 $h_{0.95} \ge 0.04$，单对效能严格为 **0%**。多源联合（如 Group 2 全 $5$ 对同时通过）进一步压低全线一致 `LINE_EQUIVALENT_NO_CONTROL_USE` 的联合概率。各线基于历史 `FCT - Mis_s`（或 `FCT - E1`）数据的实测效能评估如下：

| 组别 / 角色 | 细胞系 | 已观测对齐 $N_{\mathrm{eval}}$ ($\rho^\perp_{AB}$) | 历史实测配对 `95% CI` 半宽 $h_{0.95}$ | $\Delta=0$ 剩余边界 $0.04 - h_{0.95}$ ($\mathrm{SE}$ 倍数) | TOST $\delta=0.04$ 效能评估与 Gate 2 风险 |
|---|---|---|---|---|---|
| **Group 2 主线** | `a549_ko` | `445` (`+0.789`, R12) | **`0.0073 – 0.0155`**（R12 6 源实测 `FCT-Mis`） | `0.0245 – 0.0327` (`3.1–8.8 SE`) | **高效能（单对及联合均 >95%）**；前提是 Gate 2 与 `P_model` 交集后保留 $N \gtrsim 250$。 |
| **Group 2 主线** | `cd4t_rest` | `302` (`+0.147`, R12 dev) | **`0.0075 – 0.0202`**（R12 dev 6 源 `h975` 换算；`FCT-E1` 实测 `0.0162`） | `0.0198 – 0.0325` (`1.9–8.5 SE`) | **中高但偏脆**：因 $\rho^\perp_{AB}$ 偏低（`+0.147`），`nadig_hepg2` 对半宽已达 `~0.0202`；若 Gate 2 `P_model` 将 $N$ 压至 $\le 150$，效能会明显下滑。 |
| **Group 2 主线** | `feng_ipsc` | `93` (`+0.557`, R11/R12 dev) | **`0.0110 – 0.0170`**（R12 dev 6 源换算；R11 `Mis_kolf21j` 实测 `0.0111`，`FCT-E1` `0.0214`） | `0.0230 – 0.0290` (`2.6–5.1 SE`) | **仅在 $N\approx 93$ 维持时较高**；$N=93$ 距 $N \ge 50$ 门槛很近，若 Gate 2 交集跌到 $N=50\text{--}60$，$h_{0.95}$ 扩至 `0.023–0.030`，联合效能降至中等以下。 |
| **Group 2 主线** | `dld1` | `86` (`+0.253`, R11 dev) | **`0.0174`**（R11 dev `Mis_kolf21j`）；`FCT-E1` 实测 **`0.0255`** | `0.0145 – 0.0226` (`1.1–2.5 SE`) | **效能偏低（边缘线）**：$N=86$ 极接近 $50$ 门槛，若其他错配对或 `STATE_SE` 噪声接近 `0.0255`，单对在 $\Delta=0$ 仅剩 `~1.1 SE`（~50–70%），5 对联合通过概率明显偏低，容易输出 `MIXED_OR_INCONCLUSIVE`。 |
| **Group 1 主线** | `nadig_hepg2` / `nadig_jurkat` | 本地对齐 `797` / `930` | （R10 `dev_r9` CSV 未报 `r_orth_cal` 配对 `FCT-Mis`） | 按 $N\propto 800\text{--}930$ 尺度缩放预期高于 `a549_ko` | **预期高效能**，前提是 Gate 2 后 $N_{\mathrm{eval}} \gtrsim 400$ 且 $\rho^\perp_{AB} \ge 0.2$。 |
| **可选扩展** | `vcc2025_h1` | `80` (R12 dev) | 换算约 `0.0154 – 0.0299`（`replogle_rpe1` $h_{975}=0.0341$） | 最低仅 `0.0101` (`~0.66 SE`) | 若批准扩展，`N=80` 对最噪错配源效能极低。 |
| **敏感性线** | `cd4t_stim48` | `503` (R11) | `0.0137 – 0.0154`（R11 `Mis_kolf21j` / `E1`） | `0.0246 – 0.0263` (`3.1–3.8 SE`) | **高 TOST 效能**（仅作非独立敏感性参考）。 |
| **敏感性线** | `a549_myllia` | `32` (`+0.388`, R14) | **`0.0233 – 0.0337`**（R14 6 源实测；`FCT-E1` `0.0369`） | 最差仅 `0.0063` (`0.37 SE`) | **严重欠效能（Single <30%，多对联合 <10%）**。 |
| **敏感性线** | `hek293_xucao` | `16` (`+0.301`, R12) | **`0.0403 – 0.0847`（全部 > 0.04！）** | **$< 0$（接受域为空！）** | **数学上 0% 效能**：在 $N=16$ 且与 `FCT` 类似方差下，即使 $\hat{\Delta}=0$，半宽已超过 $\delta=0.04$，永远不可能判 `EQUIVALENT`。 |

- **修改要求:** 在 `PREREG_I.md §3` 强制加入**事前 TOST 检验效能分层披露**，明确声明：（1）$N_{\mathrm{eval}}^{\mathrm{R16}} \ge 50$ 门槛仅保证 bootstrap 数值稳定与极大差异检测，**不能**保证 $\delta=0.04$ 下 TOST `EQUIVALENT` 拥有 $80\%$ 效能；（2）一旦 Gate 2 交集把 `dld1`、`feng_ipsc` 等压近 $50\text{--}90$，或对 `a549_myllia`/`hek293_xucao` 跑敏感性检验，落入 `INCONCLUSIVE`/`MIXED_OR_INCONCLUSIVE` 应解释为“**小样本效能不足**”，不得误读为“排除了等价”。

### M5. `STATE` 推断协议仍存在关键欠规范与混杂漏洞（`512` 对照抽样、4 检查点不对称、与 `batch_encoder=true` 冲突）
- **证据定位:** `PREREG_I.md` §1.2 L41、§2.1 L57–61、§2.2 L82–83、`PREREG_I_AUDIT_R1.md` §4 L88。
- **具体缺口与统计学风险:**
  1. **`16 × 32 = 512` 对照抽样细则未写死（必须改）：** `PREREG_I.md` L82–83 只写了 `512` 细胞/系与 `16 × 32` 集合，却未锁定：（a）每个控制源的 `512` 个 NTC 是全扰动共享同一份 `16×32` 划分还是按扰动重排（共享同一划分会在扰动间引入共有控制采样方差相关性）；（b）若某系（或过滤/QC 后如 `a549_ko`）NTC 细胞不足 `512`，是有放回补抽还是按 `32` 整除截断；（c）是否按批次/供体分层抽样（例如 `cd4t_rest` 跨 50 个 `.h5` 文件 / 多通道，`a549_ko` 跨 4 文件）；（d）**匹配对照基线自重叠偏差**：若 `STATE_SE`（`predict_residual=true`）将那 512 个采样目标 NTC 的均值直接用于残差预测或偏移校正，而观测值 `delta_obs` 使用同一批总体 NTC 计算残差，则匹配臂可能与 `delta_obs` 的控制基线噪声产生弱正/负关联，而错配臂则完全独立；需要在 Gate 2 检查 `predict_residual=true` 的精确数学定义（预测值是“相对输入 32 细胞的残差”还是“预测后验表达再减去对照均值”），并记录各线总 NTC 数量及抽样占比 $512 / N_{\mathrm{NTC}}$。
  2. **Group 2 4 检查点均值 vs Stratum A 训练重叠披露缺失（必须改）：** 主控在 [`PREREG_I_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_AUDIT_R1.md#L88) 已要求“锁定前补一句披露：第 2 组用 4 个检查点的平均，`hepg2` 和 `jurkat` 分别各被其中 3 个见过”，但目前 `PREREG_I.md` v3（`3f4b0edf`）正文尚**未写入这一句准确机制说明**！更重要的是，Group 1 用单检查点 LOO，Group 2 用 4 检查点集成平均（会人为压缩预测随机方差、提高 Group 2 相对 Group 1 的信噪比），且在 Stratum A 中，匹配臂 $t$ 对 4/4 检查点均为未见 OOD 分布，而错配对照 `nadig_hepg2`/`nadig_jurkat` 对 3/4 检查点是训练内分布（甚至若 `zeroshot/hepg2` 对 `nadig_hepg2` 是唯一 LOO 检查点，却被其余 3 个见过的检查点掩盖）。因此：必须在 `PREREG_I.md` §0.3 / §2.2 / §3 明确补充该 `3/4` 披露，并在 `§4.3` 分检查点拆解中预先指明“对照严格 LOO 单检查点（如用 `zeroshot/hepg2` 评估 `mis=nadig_hepg2`）以彻底隔离训练集污染”。
  3. **`batch_encoder=true` 相同协变量的潜在人为错配/压制风险（必须改）：** `PREREG_I.md` §2.1 L60 要求匹配与错配运行使用完全相同的 `batch/cell_type` 协变量映射。这一约束在“未见系均传入中性/零向量协变量”时是合理的；但如果 Gate 2 发现 `ST-SE-Replogle` 实际读取 `cell_type_onehot_map.pkl` 且按以下任一方式工作，就会产生偏倚：（a）若对 Group 1 `nadig_hepg2`/`jurkat` 误传了可识别的协变量前缀、或未知批次统一被 fallback 映射为某一个训练批次 ID（如索引 0），会让模型在所有运行中强依赖该固定批次编码，从而**压制模型对输入对照表达式的敏感度（伪造 `EQUIVALENT`）**；（b）反之，如果匹配臂使用了目标细胞系真实协变量 token，则在错配臂仍沿用同一 token 相当于给模型输入“**目标系标签 + 异系细胞表达**”的矛盾输入。**必须要求 Gate 2 源码/权重结构检查不仅读取 `pkl` 映射表，还要在合成张量上运行两项零真实细胞冒烟核验：① 在固定相同中性协变量下变换输入对照 `X_ctrl` 时，`ST-SE` 输出确实发生非零变化；② 确认 Group 1 与 Group 2 在零样本推理时传给 `batch_encoder` 的协变量路径完全对等（均为零向量或 `None` 等中性处理）。**

### M6. 对 R9–R14 结果已知前提下的事后偏倚（Hindsight Bias）防控仍需补强披露
- **审计核查的既知事实:** `PREREG_I.md` v3（L19, L69–79）及 [`PREREG_I_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_AUDIT_R1.md)（`I11c`）已要求标记为“预注册的回顾性模型审计”。但基于 R10–R12 实际报告核查，仍存在 3 处高敏感事后偏倚点未在 `PREREG_I.md` 中充分透明化：
  1. **$\delta = 0.04$ 锚定强线导致低基线组“相对宽松度失衡”:** $\delta = 0.04$ 明确取自 R12 唯一确认线 `a549_ko` 的最小 `FCT - Mis` 差值（`+0.0796` vs `nadig_jurkat`）的一半。这仅占 `a549_ko` `FCT` 得分（`+0.1898`）的 **`21%`**，占 `dld1`（`+0.1401` R10 / `+0.1355` R11 dev）的 **`~29–30%`**、占 `feng_ipsc`（`+0.1212`）的 **`33%`**，却占 **`cd4t_rest`（`FCT=+0.0556` R12 dev）的 `72%`！**（R12 dev 中 `cd4t_rest` 有 3/6 条 `FCT - Mis` 本身就在 `+0.0347..+0.0398 <= 0.04`）。若事先不写明这种相对比例差异，后续在 `cd4t_rest` 上即使 `STATE_SE` 表现极差或仅剩微弱信号，也极易因高绝对门槛而被笼统解释为 `EQUIVALENT`。
  2. **`H-I2` 中 `FCT` 强弱事先已知带来的解读偏倚:** R10 中 `FCT` 在 `dld1`、`cd4t_rest` 未通过非配对 `H-C2` 门（`0/2 未确认`），`cd4t_stim48` 在 R11 也未通过 `H-D2`（且因共享 Marson CD4 供体/文库依赖性降为敏感性），而 `a549_ko` 与 `feng_ipsc` 通过。因此 `H-I2`（`STATE_SE - FCT`）并非在盲态下比拼未知基线，而是已知混合了高分基线（`a549_ko +0.190`）与弱基线（`cd4t_rest +0.056`）。必须在 §3 预注册：**禁止通过 `cd4t_rest` 胜过/打平 `FCT` 来掩盖在强 `FCT` 线（`a549_ko`）或整体特异性上的失败**。
  3. **入池普查证明（防“挑选候选线”嫌疑）:** 建议（并作为合规必要披露）在 `PREREG_I.md §2.2` 明文确认 `{nadig_hepg2, nadig_jurkat} + {a549_ko, dld1, cd4t_rest, feng_ipsc}` 构成本地同步盘/SSD 上**全部**已有原始单细胞计数、独立、且 R9–R14 既往目标数 $\ge 50$ 的队列全集（仅 `cd4t_stim48` 因 R11 已锁定的非独立性规则排除，`a549_myllia`/`hek293_xucao` 因 $< 50$ 靶点排除，`nourreddine_kolf21j` 与 `k562`/`rpe1`/`vcc2025_h1` 因本地无完整原始单细胞计数排除）。

---

## 2. 建议改（提升稳健性与工程严谨度）

1. **补充 `test_unit_r16.py` 边界与无效路径覆盖度（强烈建议）:**
   - 目前 [`test_unit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/test_unit_r16.py#L149-L247) `test_5` 与 `test_7` 尚未测试以下分支：
     - `valid_frac < 0.95` 与 `NaN` 在 `classify_pair_tost_r16` 中强制返回 `INCONCLUSIVE`（且修复 M2 后验证 `within_delta == False`）；
     - 严格开区间边界相等（`ci975_lo == 0.0`、`ci95_lo == -0.04`、`ci95_hi == +0.04` 均必须落入 `INCONCLUSIVE`）；
     - `orth_suite` 自身 `ci975_lo <= 0`、`valid_frac < 0.95` 或 `rho_ab_orth <= 0` 时阻断 `LINE_MATCHED_SPECIFIC`；
     - `adjudicate_claim_i1_line` 在 `seen_sources=()`（空 Stratum A）时返回 `MISSING_STRATUM`；
     - `adjudicate_claim_i2_line` 的 `FCT_SUPERIOR_TO_MODEL` 分支、以及 `valid_frac < 0.95` / `NaN` 强制归入 `INDISTINGUISHABLE`（`test_7` L320–327 目前仅测了 `MODEL_SUPERIOR_TO_FCT` 和 `INDISTINGUISHABLE`）；
     - 将 `test_7` 中的真实调用产物 `orth_good` 直接检查是否包含 `rho_ab_orth` / 正确传递分半相关字段。
2. **明确核实并修正 `PREREG_I.md` §2.2 本地路径前缀与 `nourreddine_kolf21j` 事实记录:**
   - 经只读文件系统核验：`nourreddine_kolf21j` 在 [`round9/p2/shared/lines/nourreddine_kolf21j/`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/shared/lines/nourreddine_kolf21j/) 确仅有已筛 `6,148` 基因的 `ctrl_cells.npz` (`1.1 GiB`) 与 `pseudobulk_delta.npz` (`46 MiB`)，无全转录组原始单细胞 `.h5ad`，排除决定准确无误。
   - `round12/shared/raw/a549_ko/`、`round10/shared/raw/{dld1, cd4t_rest}/`、`round11/shared/raw/feng_ipsc/`、`round14/shared/raw/gse337804` 均确认在本地 SSD；`nadig2025`（`GSE264667_hepg2_raw_singlecell_01.h5ad` 5.2 GiB 与 `GSE264667_jurkat_raw_singlecell_01.h5ad` 8.7 GiB）及 `replogle2022`（仅归一化 bulk `.h5ad`：`K562_essential_normalized_bulk_01.h5ad` 76 MiB、`rpe1_normalized_bulk_01.h5ad` 91 MiB、`K562_gwps_normalized_bulk_01.h5ad` 357 MiB）实际挂载在 `/Volumes/雷雳传输/A顶会论文同步盘/A_BIG_MAC顶会论文写稿/AI自主研究/gemini/external_data/`，而非项目根目录相对路径 `./external_data/`。建议在 `PREREG_I.md §2.2` 表中注明同步盘绝对/配置前缀，避免后续脚本找不到文件。
3. **增加对原始 Pearson / 余弦和 `SE-600M` 表征非退化的敏感性理智检查:**
   - 若因 `SE-600M` 基因名映射缺失率过高（如 Ensembl ID vs HGNC symbol 在 `19,790` 基因表中大量未命中）导致匹配与错配输入的 `SE-600M` 嵌入都接近常数向量，也会人为压低差异。建议 Gate 2 除检查 `G_model`、`P_model` 外，锁定各输入文件基因 ID 映射到 `SE-600M` `19,790` 表的覆盖率（例如记录映射基因数），并在报告中并行核对匹配与错配 `SE-600M` 对照嵌入的中心余弦相似度及欧氏距离是否显著非零。

---

## 3. 无问题（经逐行核对确认数学与程序合规项）

1. **双共主要论断 Bonferroni 与 TOST 等价边界配对严格精确:**
   - 仅保留 `H-I1` 与 `H-I2` 两项共主要论断、仅设 `STATE_SE` 为主要模型（`I5`），按两检验 Bonferroni 校正 $\alpha_{\mathrm{claim}} = 0.05/2 = 0.025$：
     - 优效/劣效双侧检验采用 `97.5%` 百分位自助法置信区间（`[1.25%, 98.75%]`）；
     - TOST 等价检验原假设 $H_0: |\Delta^{\mathrm{pair}}| \ge \delta$ 拆为两个单侧 0.025 检验，由 Berger 交并检验（IUT）定理精确等价于双侧 **`95.0%` CI（`[2.5%, 97.5%]`）** 严格落入 $(-\delta, +\delta)$ 内部（即 `ci95_lo > -delta and ci95_hi < +delta`）。置信水平配对无误。
2. **4 类配对互斥与完备性（在有效自助法样本下）:**
   - 在 `classify_pair_tost_r16`（`adjudicate_r16.py` L57–66）中，`vf < 0.95` 或非有限 CI 率先判 `INCONCLUSIVE`；余下按 `lo975 > 0` $\to$ `MATCHED_BETTER`，`hi975 < 0` $\to$ `MISMATCHED_BETTER`，`tost_within` $\to$ `EQUIVALENT`，`else` $\to$ `INCONCLUSIVE`。由于同一有序自助法样本满足 `lo975 <= lo95 <= hi95 <= hi975`，前两分支不可能同时成立，后两分支补集覆盖全部有限区间，严格**互斥且完备**。
3. **`I14` 两层错配池隔离护栏严谨有效:**
   - `adjudicate_claim_i1_line`（L140–157）强行校验证 `seen_sources` 与 `unseen_sources` 互不相交，且其并集与 `paired_mismatched['per_source']` 键集合完全相等（防漏报、防混入）；`test_6` 验证非本地线（`nourreddine_kolf21j`、`replogle_k562_essential`、`replogle_rpe1`、`vcc2025_h1`）被自动拒绝。
4. **`I1–I4`, `I7`, `I9–I13` 认知诚实性与盲法隔离到位:**
   - `r_orth_cal`（`all_nontarget`）已恢复为唯一主要终点（`I2`）；
   - `r_dorth_cal` 引用 `Systema` (`doi:10.1038/s41587-025-02777-8`) 并降为 `§4.1` 描述性分析，几何零空间“性质 P1”迁入单元测试 `test_2`（`<= 1e-14`，`I3`）；
   - “机制假设 M2”（`I4`）与 `CellularAI_OpenCode_PublicRefit`（`I7`）均降为 `§4` 描述性分析且不做确认性判定；
   - `r_orth_cal` 与冻结 `metrics_r12.py` 数值一致性在 `test_1` 验证至 `<= 1e-12`；
   - 鲜明新随机种子 `20261601` 未与 R9–R15 共享；
   - Gate 1（查新 + 人类许可）$\to$ Gate 2（仅权重词表结构检查，零真实细胞）$\to$ Gate 3（锁定 `PREREG_I_COMMIT.txt`）$\to$ Gate 4（`CODE_FREEZE_R16.txt`）$\to$ Gate 5（实细胞推断 + 锁后 Group 1 `adata_pred.h5ad` 核验）时序逻辑清晰，且 `test_4` 验证了未锁定/未冻结前禁止推断与零训练泄漏拦截；
   - [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md)、[`adjudicate_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/adjudicate_r16.py)、[`test_unit_r16.py`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/src/test_unit_r16.py) 中不存在违背 `I1` 的优先权宣称，也未改动 R12–R15 冻结文件。
