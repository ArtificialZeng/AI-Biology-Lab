# Round 14 `GSE337804`（`A549`，独立实验室）可行性核查报告

- **执行日期**：`2026-11-23`
- **执行智能体**：`A-T1`（严格遵守 [R14 指令](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R14_第十四轮_GSE337804可行性指令.md)：仅拉取 GEO/SRA/BioProject 元数据与 `<= 100 MB` 的小文件 `GSE337804_2cl_description.pdf`（`236,820 B`）及 `GSE337804_sgRNA_lib.csv.gz`（`12,261 B`）；**未下载任何表达矩阵**、未注册账号、未同意条款、未执行 `git`）
- **明确结论**：**`可行`**（技术与契约层面具备完整可执行条件；同时必须预先锁定 `PREREG_G` 前置无指标 HDF5/单 sgRNA 整数门控，且明确披露 **`n_eval_targets >= 50` 容错极窄（理论上限 `61`，仅允许最多 `11` 个脱落；乐观按 R12 `A549` 为 `~60.1`，保守体制均不达 `50`）**，仅回答「**同一细胞系 `A549`、独立实验室、不同化学方法/Cas9 与不同 guide 库上的独立数据集复现**」，不主张跨细胞系推广）。

---

## 1. GEO 全部补充文件与样本清单

### 1.1 系列补充文件（`GSE337804`，`Public on Jul 10 2026`，BioProject `PRJNA1491103`）

| 文件名 | 大小（`HTTP HEAD` / 实测） | 本轮是否下载 | 用途与备注 |
|---|---:|---|---|
| `GSE337804_2cl_description.pdf` | `236,820 B` (`236.8 KB`) | **已下载至 scratch 审阅** (`sha256: 6374660318364c2e3d407cecad5904c5fb42ad6b5ee3a36ebbb5cd8cebecd87e`) | 4 页官方实验与数据结构说明（来源文件标注 `2cl_dataset_description.html`） |
| `GSE337804_sgRNA_lib.csv.gz` | `12,261 B` (`12.3 KB`) | **已下载至 scratch 审阅** (`sha256: 96f6eb62449ba276d36591af2f7e0fca47ffedb9f19d711bcc27671d50882df6`) | 全库 `694` 条 sgRNA 清单及克隆质粒测序读数（7 列：`id,symbol,rank,class,entrez,guide,count`） |
| `GSE337804_A549_ch_1.RDS.h5ad` | `563,808,587 B` (`563.8 MB`) | **未下载** (`>100 MB` 表达矩阵) | 通道 1（`batch 1`）Seurat v3 转出 `.h5ad` |
| `GSE337804_A549_ch_2.RDS.h5ad` | `608,829,330 B` (`608.8 MB`) | **未下载** (`>100 MB` 表达矩阵) | 通道 2（`batch 2`）Seurat v3 转出 `.h5ad` |
| **合计** | **`1,172,886,998 B` (`1.173 GB`)** | — | 全部公开补充文件仅此 4 个（无 `.MTX` 文件夹上传至 GEO；PDF 提及的 `HCT116` 未在本 GEO 或 `PRJNA1491103` 释放） |

### 1.2 样本（GSM）与 SRA 测序深度核实（`PRJNA1491103`）

| GSM | 标题与库类型 | 细胞系 / 通道 / 条件 | 平台与化学方法 | SRA（`8` NovaSeq X lanes/GSM） |
|---|---|---|---|---|
| `GSM9860781` | `A549_ch1_WTS`（全转录组 mRNA） | `A549`, `batch: 1`, `Basal`, `Single-cell perturb seq` | 10x Chromium Next GEM `5' HT v2` + Jumpcode `DepleteX` | `SRX34196100`: **`1,178,809,048` reads** (`356.00 Gb`) |
| `GSM9860782` | `A549_ch1_sgRNA`（guide 捕获） | `A549`, `batch: 1`, `Basal`（标题误写 `CRISPIRi`） | 同上 `5' HT v2` | `SRX34196101`: **`77,102,681` reads** (`23.29 Gb`) |
| `GSM9860783` | `A549_ch2_WTS`（全转录组 mRNA） | `A549`, `batch: 2`, `Basal`, `Single-cell perturb seq` | 10x Chromium Next GEM `5' HT v2` + Jumpcode `DepleteX` | `SRX34196102`: **`1,448,829,095` reads** (`437.55 Gb`) |
| `GSM9860784` | `A549_ch2_sgRNA`（guide 捕获） | `A549`, `batch: 2`, `Basal`（标题误写 `CRISPIRi`） | 同上 `5' HT v2` | `SRX34196103`: **`68,538,637` reads** (`20.70 Gb`) |

- 两通道 GEX 总测序量 **`2,627,638,143` reads（`793.55 Gb`）**；按 PDF 预期目标 `36,300` singlets/通道（合计 `~72,600` singlets）折合约 **`~36,200` 原始 GEX reads/singlet**（另叠加 `DepleteX` 降噪增益）。

---

## 2. 计数、单 Guide、NTC 区分、通道与独立性核查

1. **整数原始 UMI 计数**：
   - **元数据证实**：GEO `!Sample_data_processing` 明文记录 *`".H5AD files contain raw single cell mRNA counts, with cell filtering and guide calling"`*；PDF 指出其源自 Seurat v3 对象（`RNA` assay 存 mRNA counts，映射参考为 GRCh38 + GENCODE v39 protein-coding & IG/TR）。
   - **待下载后必查边界**：文件后缀为 `*.RDS.h5ad`（Seurat 转出 `.h5ad`）。部分转出工具会将 log-normalized 值写入 `/X` 而把整数原始 counts 写入 `/raw/X` 或 `/layers/counts`。因本轮严格禁下矩阵，**真实 HDF5 矩阵槽位及非负整数校验须在 `PREREG_G` 锁定后执行前置 `E3` 式无指标检验确认**。
2. **逐细胞单 sgRNA 标签（低 MOI CRISPRko）**：
   - PDF 明确核实为 **低 MOI（`~6% EGFP+`）+ 嘌呤霉素筛选的 `CRISPRn (= CRISPR knockout)` 筛选**（修正 GEO 整体设计中“dual guide”措辞及 Sample 标题“CRISPIRi”误拼：`218` 基因 × `3` sgRNA，单细胞内按单 sgRNA 处理，多 sgRNA 视作 multiplet 剔除）。
   - PDF 说明 `seurat_objects`（即 `*.RDS.h5ad` 来源）**`contain only filtered cells with exact 1 sgRNA called`**，分配结果记在细胞元数据列 **`sgrna_id`** 与 **`sgrna_symbol`** 中（具体保留细胞数与是否全为严格单 guide 需下载后通过 `.obs` 门控核实）。
3. **NTC 标记（`20` `ctrl_NT` vs `20` `ctrl_T` 能否区分；`ctrl_T` 性质）**：
   - **完全可区分，`ctrl_T` 确认为 `safe-targeting`**：`GSE337804_sgRNA_lib.csv.gz` 中两类对照字段严格分离：
     - **`ctrl_NT`（`20` 条，真非靶向对照）**：`class == "ctrl_NT"`、`id == "non-targeting_1".."non-targeting_20"`、`symbol == "non-targeting"`（质粒库总 reads `38,830`，占比 **`1.69%`**，中位 `2,211.5` / min `370` / max `3,153`）。
     - **`ctrl_T`（`20` 条，safe-targeting 对照）**：`class == "ctrl_T"`、`id == "gene_desert_1".."gene_desert_20"`、`symbol == "gene_desert"`（质粒库总 reads `66,105`，占比 **`2.88%`**，中位 `3,493.0` / min `1,083` / max `4,162`）；PDF 明确说明为 ***"20 sequences targeting gene desert regions"***（基因组荒漠切割对照，含 Cas9 DSB DNA 损伤效应）。
   - **确定性规则**：因二者通过 `class`、`id`（`non-targeting_*` vs `gene_desert_*`）及 `symbol`（`non-targeting` vs `gene_desert`）**完全可分**，建线时应 **严格只取 `ctrl_NT` 作为 `NTC`，并将 `ctrl_T`（`gene_desert_*`）彻底排除出 NTC 与靶点**（与 [dataset_readers_r13.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/src/dataset_readers_r13.py) 及 [PREREG_E](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) 规则完全一致）。按质粒比例粗估，`ctrl_NT` 预期总计约 `~1,230` 细胞（每通道约 **`~615` NTC**），远超 `>= 15`/half 要求。
4. **通道 / 重复 / 处理与潜在干扰项**：
   - **天然互斥双通道**：共 **`2` 个独立 10x 通道**（`ch_1` vs `ch_2`），天然构成 `Half A (ch_1)` vs `Half B (ch_2)` 与各自通道内 NTC 中心化扣除（`>=2` batches，无需 `seed=42` fallback）。
   - **单一条件**：两通道均为 `Basal` 未刺激单克隆 Cas9 `A549`（Day 16 收样），**无药物处理或时间点混合**。
   - **特殊工艺（`DepleteX`）与出处边界**：GEX 库使用了 Jumpcode `DepleteX` 消除约 `255` 个高丰度核糖体/线粒体/非变异管家转录本。经核对，**`R9_COMMON_6148` 内仅含 `11` 个 `RPL/RPS/MT-*` 基因**，且 `GSE337804` 交集靶点 `61` 个中 **`0` 个属 `RPL/RPS/MT-*`**，不会阻塞 `|G_eval| >= 3,000` 门槛（实际 `|G_eval|` 与无 `%MT` 上限处理需在 `PREREG_G` 前明文披露，尤其若 `MT-*` 被耗竭需同 R12 HEK293 / R13 Flex 处理不设 `%MT` 上限）。
   - **文献溯源与独立性**：`GSE337804` 贡献作者（Stefanie Fellinger, Adam Krejci, Tilmann Burckstummer, Christopher Carpenter, Hani Goodarzi）对应 *Cell* 2026 论文 `STATE`（Adduri et al., `10.1016/j.cell.2026.07.052`, PMID `42673963`, Myllia Biotechnology × UCSF/Arc），与 R12 `GSE345058`（Liu et al. `10.64898/2026.06.01.728087`, `SA-Cas9`）实验团队、化学平台（`SpCas9` + `10x 5' HT` + `DepleteX`）及载体均独立；且不属 `X-Atlas/Orion`（`HCT116`/`HEK293T`）或 `Arc Virtual Cell Challenge`（`H1`/`AC16`/`MOLM13`/`U-87MG`）。

---

## 3. `sgRNA_lib.csv.gz` 靶点核查与多种比例估算

基于 `GSE337804_sgRNA_lib.csv.gz` 全量 `694` 行 + 冻结 `hgnc_complete_set.txt`（[verify_r10_lines.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/shared/scripts/verify_r10_lines.py) 与 [prepare_r10_lines.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/shared/scripts/prepare_r10_lines.py) 同源函数实算）：

- **靶基因总数**：**`218`** 个独特基因（精确 `3` guides/基因 = `654` 条；质粒总 reads 最小为 `4,573`，无整基因 dropout），**`218 / 218` 全部为 HGNC `approved` 官方符号**（`0` unmapped，原始与规范符号一致）。
- **与 R9 冻结参考集合交集**：
  - `u6_rep` (`2,117`) 交集：**`65`**
  - **`u6_rep ∩ R9_COMMON` (`1,698`) 交集（理论 `n_eval_targets` 最大上限）**：**`61`**（不在 `R9_COMMON` 的 4 个为 `LDB1`, `SPCS2`, `UBE3D`, `WTAP`）
- **对 R12 冻结 `a549_ko` 的无数据泄露对照核验（只查 R12 缓存 `targets/reproducible`，未碰 `GSE337804` `.h5ad`）**：
  - 上述 `61` 个交集靶点中，**`15` 个**曾出现在 R12 `a549_ko` (`GSE345058`) 中，且 **`15 / 15`（`100.0%`）在 R12 中均通过 Addendum-B 可复现门槛**；其余 **`46 / 61`（`75.4%`）** 为 R12 从未测试过的独立新增基因！

### 评测靶点数多口径校准估算（`>= 50` 门槛要求 `>= 81.97%` 存活）

| 校准情景 | 来源依据 | 换算公式 | `n_eval_targets` 估算 | 是否过 `>= 50` 门槛 |
|---|---|---|---:|---|
| **乐观（R12 `A549` 对 `R9C` 比例）** | R12 `a549_ko` 建库清单分母（[R12_DOWNLOAD_PLAN.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_DOWNLOAD_PLAN.md)）：`445 / 452 = 98.45%`（build 后 `445 / 451 = 98.67%`） | `61 × (445 / 452)` | **`60.1`**（build 后口径 **`60.2`**） | **通过**（仅余 `~10–11` 容错边际） |
| **次乐观（R12 `A549` 对 `u6_rep` 比例）** | R12 `a549_ko` 建库清单分母：`445 / 530 = 83.96%`（build 后 `445 / 529 = 84.12%`） | `65 × (445 / 530)` | **`54.6`**（build 后口径 **`54.7`**） | **勉强通过** |
| **通用经验范围（`20%–40%` 规则）** | 跨数据集常规过滤率 | `61 × [0.20, 0.40]` (`65 × [0.20, 0.40]`) | **`12.2–24.4`** (`13.0–26.0`) | **不通过** |
| **保守（R12 `HEK293` 比例）** | R12 `hek293_xucao` 建库清单分母：`16 / 143 = 11.19%` / `16 / 164 = 9.76%`（build 后 `16 / 125 = 12.80%` / `16 / 139 = 11.51%`） | `61 × (16 / 143)` | **`6.8`** (`6.3`；build 后 **`7.8`** / `7.5`) | **不通过** |
| **极端保守（R13 `mac_ctrl` 比例）** | R13 `mac_ctrl`：`3 / 91 = 3.30%`（`3 / 107 = 2.80%`） | `61 × (3 / 91)` | **`2.0`** (`1.8`) | **不通过** |

---

## 4. 可行性总结论与 `PREREG_G` 前置建议

> **总结论：`可行`（推荐推进 `PREREG_G` 同系独立实验室复现检验，附带严格的 `.RDS.h5ad` 前置计数/槽位门控与高淘汰边界披露）**

**理由如下：**
1. **技术与契约全面满足**：开放免登录直接下载（总计 `1.173 GB`）；全 `218` 基因 100% 匹配 HGNC；每基因 `3` guide 质粒均匀分布；元数据确认具备单 sgRNA 标签（`sgrna_id`/`sgrna_symbol`）、原始单细胞 RNA counts（需锁定后核实 `/X` vs `/raw/X` vs `/layers` 精确位置）、以及严格分离的 `20` `ctrl_NT`（non-targeting）与 `20` `ctrl_T`（safe-targeting，排除）；天然具备 `ch_1` 与 `ch_2` 独立双通道分半及单一 `Basal` 无药条件；且 GEX 测序极深（`2.63B` reads + `DepleteX`），`R9` 交集内 `15/15` 共同靶点在 `a549_ko` 中全部可复现、`46/61` 为全新独立基因。
2. **必须事先写明的风险与局限**：
   - `n_eval_targets` 上限仅 **`61`**，虽然 R12 同系 `A549` 曾实现 `98.45%`（`~60` 个），但 `GSE337804` 预计单细胞规模约 `~318` targeting 细胞/靶点（仅为 R12 `GSE345058` 的约 `57%`），一旦有 `>= 12` 个靶点未达 Addendum-B（`>=25` 细胞、`AD >= 50`、`Energy p < 0.001`）即落入 `< 50`（按保守比例仅为 `2.0–26.0`）；
   - `DepleteX` 耗竭高丰度核糖体/线粒体/部分 NVG，故 `PREREG_G` 须明确 **`max_pct_mt = None`（不设线粒体上限）** 及 GENCODE v39 基因宇宙对齐 `G_eval >= 3,000` 门控；
   - **科学解释边界锁定**：即便 `GSE337804` 通过 `§5` 门槛与 Claim 1/2，**仅能声明在「同一细胞系 `A549`、独立实验室、不同 Cas9/5' HT + `DepleteX` 化学平台与 75% 新靶点/新 sgRNA 序列」上的独立数据集复现，不可外推为多细胞系泛化**。
