# Round 15: Arc Virtual Cell Challenge (VCC) 2026 预发布可行性评估报告

- **任务智能体**：`A-T1`（A 队 · 数据工程）
- **核查时间窗口**：`2026-11-23T08:30:21Z` – `2026-11-23T08:42:35Z` UTC（北京时间 **`2026-11-23 16:30 – 16:43 +08:00`**；最终平台复检时间戳 **`2026-11-23T08:41:41Z`** / `16:41:41 +08:00`）
- **合规声明**：严格遵守只读免登录边界。**未下载任何表达矩阵**，未注册 `virtualcellchallenge.org` / `vcc-cli` 账号，未接受任何 Terms of Use / EULA，未向外部发送任何邮件/issue/评论，未修改 `PREREG` 或 `A-T5` 产物。所有 2025 `.h5ad` 结构检查均通过 HTTP `Range` 请求仅读取头文件元数据（单文件拉取 `338 KB – 829 KB`）。

---

## 一、发布状态全平台核查（截至 `2026-11-23T08:42:35Z`）

| 平台 / 入口 | 查询地址 / 接口 | 查询时间 (UTC / +08:00) | 核查结果 |
|---|---|---|---|
| **GCS Bucket（`2026/` 前缀）** | [`storage.googleapis.com/.../prefix=virtual-cell-challenge/2026/`](https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?delimiter=/&prefix=virtual-cell-challenge%2F2026%2F) | `08:38:50Z` / `16:38:50` | **为空**：`prefixes = None`, `items = []`（`0` 个文件/子目录） |
| **GCS Bucket（`virtual-cell-challenge/`）** | [`storage.googleapis.com/.../prefix=virtual-cell-challenge/`](https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?delimiter=/&prefix=virtual-cell-challenge%2F) | `08:38:50Z` / `16:38:50` | **仅含 `2025/`**：`prefixes = ['virtual-cell-challenge/2025/']`（`adata_Training/Validation/Test.h5ad`、`gene_names.csv`、`pert_counts_*.csv` 上传 `timeCreated` 均为 `2025-12-16T17:08–17:36Z`，外加 `FASTQ/`）；根目录仅 `scbasecount/`、`tahoe100M/`、`virtual-cell-challenge/` |
| **HuggingFace (`arcinstitute`)** | `huggingface.co/api/{datasets,models,spaces}?author=arcinstitute` | `08:41:41Z` / `16:41:41` | **`0` 份 VCC 2026 数据/模型**：共 `22` 个数据集（最新修订 `2026-10-28T21:52:33Z`，VCC 相关仅 `arcinstitute/PIE_arc_vcc_25`）、`28` 个模型、`0` 个 Space；全站 `search=virtual-cell-challenge` 仅返回第三方 `Arcticbun/2025_Virtual_Cell_Challenge_Test_Data`（`2026-05-18`，2025 测试集） |
| **NCBI GEO (`gds`) & BioProject** | `eutils.ncbi.nlm.nih.gov` (`gds`, `bioproject`, `pubmed`) | `08:41:41Z` / `16:41:41` | **`0` 份 VCC 2026 记录**：`gds` 检 `"Cal-33" AND "MOLM13"` = `0` 条；`"Virtual Cell"[Title] AND 2026[PDAT]` = `0` 条；`bioproject` `"Arc Institute" AND 2026[PDAT]` 共 `6` 条（最新 `PRJNA1540732` / `GSE349934` 为 HEK293T ecCLIP，均非 VCC 2026）；`PubMed 42648290` 为 2026-09-17 *Cell* 评述 (`10.1016/j.cell.2026.08.004`) |
| **Zenodo & Figshare** | `zenodo.org/api/records` & `api.figshare.com/v2/articles` | `08:41:41Z` / `16:41:41` | **`0` 份 Arc VCC 2026 数据**：Zenodo `"Virtual Cell Challenge"` 仅 `4` 条命中（2025 `17342530`、社区脚本 `18004721` 等）；Figshare 仅返回 10 条无关条目 |
| **`virtualcellchallenge.org`** | `/datasets`, `/rules`, `/terms`, `/faq`, `/evaluation`（Next.js JS chunks） | `08:30–08:42Z` / `16:30–16:42` | **仅限注册参赛者登录后下载**：对照包（`controls-val` 2026-08-20 释出；`vcc_2026_controls_final.zip` 约 `615 MiB` 2026-10-22 释出）须登录 + **Persona 人脸自拍验证** + 同意 `/terms` 后通过网页或 `vcc datasets download` 获取；投稿已于 `2026-11-05 23:59 UTC` 截止 |
| **`arcinstitute.org/news`** | Newsroom 与 `2026-11-19` Wrap-Up / `2026-10-01` Behind the Data | `08:30–08:40Z` / `16:30–16:40` | **最新博文停留在 `2026-11-19`**（*Virtual Cell Challenge 2026: Winners and Wrap-Up*），其后**无任何数据开放新公告**；官方**未公布 2026 测试集真值或完整筛选数据的具体开放日历日期**，仅文末引导 *"Explore our open datasets in the Virtual Cell Atlas"* |
| **关联代码与预印本核查** | GitHub `ArcInstitute/cell-eval2`、`PyPI vcc-cli`、`api.biorxiv.org` | `08:38–08:42Z` / `16:38–16:42` | `cell-eval2==0.22.0` 与 `vcc-cli==0.3.0` 公开，已释出 `configs/vcc2026.yaml`、`vcc2026-metrics-brief.md/.pdf` 及 `vcc2026-metrics.pdf`；Metrics 预印本 `10.64898/2026.11.17.767214` 尚在 bioRxiv API 索引中（`status: no posts found`） |

- **历史释放节奏参考（非官方承诺）**：2025 年结项博文发布于 **`2025-12-06`**，完整 H1 训练集/验证集/测试集 `.h5ad` 于 **`2025-12-16 17:08–17:36 UTC`**（**结项后第 10 天**）上传至 `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/`。2026 年 Wrap-Up 发布于 **`2026-11-19`**，截至 **`2026-11-23`** 刚过 4 天，公共桶仍为空。

---

## 二、实验设计与数据结构审查

### 2.1 官方披露的 VCC 2026 实验设计

| 维度 | VCC 2026 官方设计详情（源自 `2026-10-01` Behind the Data、`2026-11-19` Wrap-Up、`/datasets`、`vcc2026-metrics-brief.md`） | 与 BioXArena `PREREG` 契约兼容性评估 |
|---|---|---|
| **细胞系构成（6 条）** | - **验证组 (`A/B/C`)**：`Cal-33`（头颈鳞癌）、`HeLa`（宫颈腺癌）、`Jurkat`（T 细胞白血病）<br>- **测试组 (`D/E/F`)**：`AC16`（心肌细胞）、`MOLM13`（急性髓系白血病）、`U-87MG`（胶质母细胞瘤） | **⚠️ `Jurkat` 必须排除**：R9 源细胞系集合已包含 `nadig_jurkat`（同源 Jurkat，不能作独立目标线）。其余 **`Cal-33, HeLa, AC16, MOLM13, U-87MG`（至多 5 条）** 均未在 R9–R14 使用 |
| **扰动方式与机制** | **单基因 CRISPRi 敲低**（全 6 线使用相同 sgRNA 序列，先经高 CRISPRi 活性阈值筛选，窄范围低 MOI 递送、统一采样时间点；完整筛选每基因 2 条 guide，竞赛集挑选其中 **"the better-performing of its two guides"** 即每基因 1 个构体）。*注：2025 文档明文写明使用 `dCas9-KRAB` (`Gilbert et al. 2014`) 双 guide 盒；2026 公开博文与预印本摘要未明文标注具体阻遏结构域名称（如 `KRAB` vs `Zim3`），仅说明 CRISPRi + 定制 sgRNA 捕获探针* | **合规（CRISPRi 低 MOI 单扰动）**，且竞赛筛选每细胞对应单一选定构体；但需注意筛选偏倚（见第五节） |
| **测序与预处理管线** | **10x Genomics Flex**（固定细胞探针杂交，**18,533** 基因公共轴） + **Ultima Genomics UG100**（每线 `>800B` reads，中位 **20,000 UMIs/细胞**），由 Arc 开发的 **`cyto`**（`10.64898/2026.01.21.700936`）处理 | 深度远高于普通 scRNA-seq（中位 20k UMI/细胞），无极端 FACS 表型门控；但与 `vcc2025_h1` 共享同实验室 Flex 工艺（见第五节） |
| **双规模数据体制** | - **完整未筛选屏幕**（尚未开源）：每线 **`>10,000`** 个 CRISPRi 敲低 + **`100`** 条 NTC（每线 `~50,000` 个对照细胞，6 线累计 `>3,000万` 细胞，中位 **`500` 细胞/扰动**）<br>- **竞赛精简子集**：验证/测试各 **`300` 靶点**（两面板重合 **`129`** 个，测试新增 **`171`** 个），筛除未达 **`>= 80%` 中位靶基因敲低率** 的扰动并统一下采样至 **`400` 细胞/扰动** + **`46` 条 NTC × `400` = `18,400` 对照细胞**（平均每条线与对应控制） | 两套规模差异极大：完整屏幕若公开将大幅拓宽基因覆盖；`300` 靶点竞赛集则受严苛敲低阈值与配额下采样制约 |
| **NTC 与原始整数 UMI** | - 对照包 `context_A..F.h5ad` 已确认提供 **`target_gene == "non-targeting"`** 及引导列 **`ntc_id`**（`46` 条非靶向 sgRNA × `400` 细胞 = `18,400` 细胞）<br>- `/datasets` 与 `configs/vcc2026.yaml` 明确规定竞赛与对照均为 **raw counts（非负全整数原始 UMI，未 log 归一化，单细胞上限 `10^6`）** | 对照充足（`18,400` > `>=15`/half）；整数原始 UMI 符合 `X` 契约 |

### 2.2 通过 2025 `.h5ad` HTTP `Range` 头文件对照检查（推测并核对缺口）

本轮使用 `h5py` + HTTP `Range` 对 GCS `2025/` 全部三个 `.h5ad` 只读拉取头部结构（未读取完整 `/X` 矩阵）：
1. **`adata_Training.h5ad`**（总大小 `15,482,497,461 B`，仅拉取 `338,357 B` / `6` 次 `Range` 请求）：`/X` shape `[221273, 18080]`；`obs = ['_index', 'batch', 'guide_id', 'target_gene']`；`var = ['_index', 'gene_id']`；`n_batches = 48`（`Flex_1_01 .. Flex_3_16`）。
2. **`adata_Validation.h5ad`**（总大小 `6,928,967,541 B`，仅拉取 `829,301 B` / `13` 次 `Range` 请求）：`/X` CSR `float32` shape `[98927, 18080]`，前 20 个存储值全部为精确整数（`[1.0, 4.0, ...]`）；`has_raw = False`、`layers = []`；`obs = ['target_gene', 'guide_id', 'batch']`；`guide_id` 记录完整构体（`31` 条 `non-targeting_*|non-targeting_*` 等）；`n_batches = 48`（`Flex_1_01 .. Flex_3_16`）。
3. **`adata_Test.h5ad`**（总大小 `11,950,739,168 B`，仅拉取 `577,248 B` / `9` 次 `Range` 请求）：结构完全一致（shape `[170846, 18080]`、`48` 个 Flex 批次、`obs = ['target_gene', 'guide_id', 'batch']`）。

> **⚠️ 对 2026 Split-Half 通道/批次列的关键核查缺口**：在 2026 竞赛方公开展示的对照文件 schema（`context_A.h5ad`）中，`.obs` 仅列出了 `target_gene`、`context`、`ntc_id`，**未出现类似 2025 Atlas 文件的 `batch`（`Flex_x_xx` 捕获通道）列**，且 `vcc2026-metrics-brief.md` 第 8 节显示官方参考上线使用 **`400` 细胞的 5 次随机半数切分**（random half-splits）。因此，未来 VCC 2026 正式开源到 Atlas 时，**必须先用 HTTP `Range` 检查公开 `.h5ad` 是否保留了 `Flex` 通道 `batch` 列**——若保留则可行按通道/批次不相交切分 + per-half NTC 中心化；若被匿名剔除，则只能回退到 Stage 2 分层随机对半切分（同 R12 `hek293_xucao` 分支）。

---

## 三、靶点名单公开状态与 `u6_rep` 通过率推算

### 3.1 是否可在不登录情况下获取靶点名单？

**不能。** 截至 `2026-11-23T08:42Z`：
- 验证/测试面板的 **`300` 靶点名单 (`pert_counts.csv`)** 仅打包于需登录 + Persona 自拍验证下载的 **`vcc_2026_controls_*.zip`** 内，完整 **`>10,000` 靶点列表**亦未公布。
- 免登录公开页面 `/datasets` 的代码框仅展示 4 个首字母示例基因（`ABCD1`, `ACLY`, `ADNP`, `AGO1`），其余全部处于未公开状态。
- **直接结论**：当前**无法**针对 VCC 2026 真实靶点直接计算 `u6_rep`、`u6_rep ∩ R9_COMMON` 及预测过门禁数量。

### 3.2 可靠参考基线：`vcc2025_h1` 实测（本轮脚本直接重算）与 2026 公开线索

1. **R9 冻结参考集合**：`u6_rep = 2,117`、`R9_COMMON = 6,148`、`u6_rep ∩ R9_COMMON = 1,698`（其中已纳入 `vcc2025_h1` 自身的 `204` 个可复现靶点）。为避免自包含偏差，排除 `vcc2025_h1` 后的其余 5 条 R9 源系参考集为 **`u5_rep = 2,004`**、**`u5_rep ∩ R9_COMMON = 1,610`**。
2. **2025 H1 实测基准（基于 `round9/p2/shared/lines/vcc2025_h1/pseudobulk_delta.npz` 直接实算）**：
   - 2025 同样为 **`300` 靶点**（Train `150` + Val `50` + Test `100`），总通过 `half_valid = 300/300`，通过 R9 可复现门禁 **`reproducible = 204 / 300 = 68.0%`**（Train `104/150 = 69.3%`、Val `34/50 = 68.0%`、Test `66/100 = 66.0%`）。
   - 若仅统计与另 5 条 R9 源系（`u5_rep`）重合的靶点：**总重合 `95 / 300`（其中 `u5_rep ∩ R9_COMMON` 为 `84`）**，通过可复现门禁 **`91 / 95 = 95.8%`**（`u5_rep ∩ R9_COMMON` 内为 **`80 / 84 = 95.2%`**！包含自身 `u6_rep` 则为 `208` 重合/`172` 在 `R9_COMMON`、`204`/`168` 可复现）。可见只要 Flex 深测序（`>400–1,000` 细胞/扰动）叠加官方敲低质控筛选，在先验重合靶点上的存活率（`~95%`）可高于 R12 `A549`（`445 / 530 = 83.96%`），但注意 `2025` 测试集平均细胞数远超 `400`。
3. **2026 Wrap-Up 公开重合线索与场景估算（非实算，待公布核验）**：
   - Wrap-Up（`2026-11-19`）明确披露：测试集 `300` 靶点中有 **`277 / 300`（`92.3%`）** 出现在 Replogle 等 (2022) 全基因组 K562 CRISPRi 数据中，`298 / 300` 出现在 X-Atlas/Orion 中。需注意 `277` 是对 Replogle 全基因组（约 `9,867` 靶点）而非 `u6_rep`（`2,117`）的重合，真实落入 `u6_rep ∩ R9_COMMON` 需待名单公布后实算。
   - **情景 A（若仅公开 `300` 靶点竞赛精简集）**：假设重合率接近 2025 基准或略高（例 `u6_rep ≈ 95–150`；假设落入 `u6_rep ∩ R9_COMMON ≈ 80–125`）：
     - 若按 **`20–40%`** 常规经验：`19–60`（`R9_COMMON` 内 **`16–50`**）——处于 **`< 50` 边缘或刚及格**；
     - 若按 **`HEK293 9.76%` (`16/164`)** 保守校准：**`9.3–14.6`**——明显 **`< 50` 淘汰**；
     - 若按 **`A549 83.96%` (`445/530`)** 乐观深测序校准（考虑到 `400` 细胞 + `>=80% KD` 强筛选）：**`79.8–125.9`**（`R9_COMMON` 内 `67–105`）——**有望跨过 `>= 50` 门禁**。
   - **情景 B（若同时公开 `>10,000` 靶点完整屏幕）**：几乎必覆盖绝大多数 `u6_rep ∩ R9_COMMON`（`1,698` 中的绝大部分，例 `1,200–1,600+`），即便套用最严苛 **`9.76%`** 保守比率（`117–156+`）或 **`20–40%`**（`240–640+`）亦可**稳健超过 `50`**。

---

## 四、许可证与合规使用状况

1. **2025 明确先例（已核实）**：Arc Virtual Cell Atlas 页面（`https://arcinstitute.org/tools/virtualcellatlas`）明文标注 `vcc 2025` 全套训练/验证/测试数据：**`© 2025. This work is openly licensed via CC0 1.0`**，允许无限制科研下载、衍生处理与公开发表，已合入 R9 `vcc2025_h1`。
2. **2026 当前受限状态（硬阻断）**：
   - `https://virtualcellchallenge.org/terms`（Effective Date `August 4, 2026`）与 `/rules`（Revised `September 16, 2026`）把下载数据定义为受版权保护的 *"Content / Services"*，明文禁止 **"Use your account for any purpose other than participating in the Challenge"**，并禁止未经明确许可 **"Use, copy, reproduce, modify, ... publish, ... distribute, ... create derivative works based on, or otherwise exploit for any purpose any of our Services"**，且规则优先于 FAQ。
   - 2026 赛事已于 **`2026-11-05 23:59 UTC`** 结束投稿、**`2026-11-19`** 结项。**截至目前，`2026` 年 6 条细胞系尚无 `CC0 1.0` / `CC BY` 等开源许可证落地**，在官方未公开发布到 Atlas 并标明开放许可证前，**严禁注册/借用参赛账号将竞赛包用作赛外独立评估集**。

---

## 五、重大技术与方法学风险

1. **与源系 `vcc2025_h1` 的同实验室 / 同探针工艺相似性风险**：
   - `vcc2025_h1`（R9 已有源系）与 VCC 2026（目标评估候选）均产自 **Arc Institute Technology Center + 10x Genomics Flex（探针杂交捕获）+ Ultima UG100 + `cyto` + CRISPRi**，而其余 5 条 R9 源系多为标准 poly-A 3'/5' scRNA-seq。即使细胞系完全独立（`Cal-33, HeLa, AC16, MOLM13, U-87MG`），shared Flex 探针检测特性（及同机构建库操作/筛选管线）可能引入 **技术平台相关性（pseudo-replication / 平台校准偏差）**，从而使 `vcc2025_h1` 对 2026 目标线的转移斜率或相关性呈现非生物性虚高。未来若立项预注册（`PREREG`），**必须预设 `leave-vcc2025_h1-out` 敏感性对照**。
2. **竞赛集 `>= 80%` 中位靶基因敲低预筛选与 `400` 细胞均衡抽样的选择偏倚（Selection Bias）**：
   - 2026 竞赛集主动剔除了所有中位靶基因敲低 `< 80%` 的扰动，并且对每基因仅选取双 guide 中效果更好的单条构体，再统一下采样至 `400` 细胞。这会截断真实生物学扰动效应分布的弱/中等效应尾部（survivorship bias），可能系统性抬高 `half_valid` 与 `reproducible` 存活率及响应幅值（slope/效应量），**使其性能估计不能直接外推至全基因组无偏敲低表现**。因此若未来 Arc 同时释放 **`>10,000` 靶点的完整未筛选屏幕**，应在验证中**优先使用完整屏幕（或至少做对照分层）**。
3. **验证与测试面板跨组重合 + 基因轴交集收缩风险**：
   - 若同时取 2026 验证线（`Cal-33, HeLa`）与测试线（`AC16, MOLM13, U-87MG`），注意两套 `300` 靶点面板间有 **`129` 个共有靶点**（而非各自完全独立）且 `U-87MG` 中有 `21` 个靶点已见于 McFaline-Figueroa 2024。
   - 2026 基因轴为 `18,533` 个 Flex 探针基因（2025 为 `18,080`，`vcc2025_h1` 经多批次过滤后留 **`14,656`** 基因），虽通常足够满足 `|G_eval| >= 3,000`，仍需正式 `.h5ad` 释出后逐线核验全批次有限性与表达下限门控。
4. **通道/批次字段缺失风险**：
   - 如 `2.2` 节所述，竞赛对照包不含 `batch`（Flex 通道），若公开后 `.h5ad` 亦缺失 `batch` 列，将无法做硬件通道/批次级不相交 `split-half`，只能退化为分层随机对半切分。

---

## 六、最终结论与重启评估触发条件

- **本轮综合结论**：**`不确定（需等待官方免登录开放发布并完成靶点/批次列与许可证核查）`**。
- **核心依据**：截至 **`2026-11-23T08:42Z` (`16:42 +08:00`)**，VCC 2026 对我们原则上具有**极高潜在价值**（提供除已用 `Jurkat` 外的 `Cal-33, HeLa, AC16, MOLM13, U-87MG` **至多 5 条独立人类细胞系**，低 MOI 单基因 CRISPRi、20k 深度原始整数 UMI、充足真实 NTC），但当前同时受制于 **三项硬缺口**：
  1. `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` 与所有公共仓均为 **`0` 文件**；
  2. 现有唯一下载渠道受 `virtualcellchallenge.org/terms` 禁止赛外使用条款与 Persona 人脸自拍门禁阻断，且尚无 `CC0 1.0` 授权；
  3. **`pert_counts.csv` / `>10,000` 全靶单未公开**，无法实算 `u6_rep ∩ R9_COMMON`，且未知公开版是否保留 Flex 通道 `batch` 列以及是否含完整未筛选屏幕。

### 转为 `可行` 的明确先决触发条件（需全部满足）

1. **免登录开源落地**：Arc Institute 在 `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/` 或 Virtual Cell Atlas 正式发布 VCC 2026 数据，并明文附带允许科研与衍生评估的开放许可证（如 **`CC0 1.0`** 或 `CC BY 4.0`）。
2. **元数据与结构 Range 验证通过**：在不下载全矩阵前提下，通过头文件审查确认包含 `Cal-33, HeLa, AC16, MOLM13, U-87MG`（**硬排除 `Jurkat`**），`/X`（或 `/raw/X`）为非负整数原始 UMI、`target_gene == "non-targeting"` NTC 充足（每 half `>= 15`）、具备构体/guide 标注，并核实是否存在可用于通道不相交切分的 `batch` 列。
3. **靶点重叠与偏倚治理达标**：实算新线对 `u6_rep ∩ R9_COMMON` 与过门禁预估稳健 `>= 50` 且 `|G_eval| >= 3,000`；若官方同时提供 `>10,000` 完整屏幕，则优先采用完整屏幕规避 `>= 80% KD` 预筛选偏倚，并在未来 `PREREG` 中将 **`vcc2025_h1` 留一敏感性检验（leave-one-source-out）** 写入预注册方案。
