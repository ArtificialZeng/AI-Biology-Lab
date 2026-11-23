# Round 14 Clean-Room Raw Reader Cross-Check Report (`READER_XCHECK_R14`)

- **执行智能体**：`A-T3`（Team A · 预注册与统计审计）
- **日期**：`2026-11-23 +08:00`
- **核对对象**：`a549_myllia`（`GSE337804`, [shared/lines/a549_myllia](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/shared/lines/a549_myllia) vs [independent_r14_loader.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/audit/independent_r14_loader.py) 干净室独立重建）
- **复现脚本**：[xcheck_r14_reader.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/audit/xcheck_r14_reader.py)（严格使用 `python3 -I`，未修改 `round14/shared/lines/`，**未计算任何预测值、相关系数或统计指标**）
- **最终结论**：**`R14_READER_XCHECK_OK`**（全部槽位检验、细胞/通道/NTC/`gene_desert` 计数、基因集、可重复与评测靶点集合及矩阵数组逐项完全一致；在使用 `prepare_r10_lines` 默认 `chunk_size=8192` 时 `ctrl` 均值与 `delta` 系列数组 `max_abs_diff = 0.0 <= 1e-6`）

---

## 1. 门控与干净室合规确认
- `enforce_r14_gate(phase="reader_xcheck_r14")` 通过：`PREREG_G_COMMIT = d34560d7ae568b299dca15c766c2e356b483c39f`，`R14_LINES_AUTHORIZED = R14`（来自 [R14_LINES_AUTHORIZED.txt](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/status/R14_LINES_AUTHORIZED.txt)）。
- 未读取、引用或导入 `round14/src/`、`round14/results/` 或 `round14/status/A-T2.txt`；`round14/shared/lines/a549_myllia/` 文件时间戳与内容保持原样未动。
- 未执行任何模型预测、残差相关、置信区间或效应量计算。

---

## 2. 原始 `.h5ad` 槽位核查与一致性检验 (`PREREG_G §4.1`)

对 [GSE337804_A549_ch_1.RDS.h5ad](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/shared/raw/gse337804/GSE337804_A549_ch_1.RDS.h5ad) 与 [GSE337804_A549_ch_2.RDS.h5ad](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round14/shared/raw/gse337804/GSE337804_A549_ch_2.RDS.h5ad) 的全量非零值（100% 扫描，零抽样）结果如下：

| 通道 | 候选槽位 | `excluded_layers` | 形状 `(cells × genes)` | 非零元素数 (`nnz`) | 全量 UMI 和 | 最大值 | 负值 / 非有限 / 非整数 | 多槽位全量对比 | 槽位记录对比（独立重建 vs `shared/lines`） |
|---|---|---|---|---:|---:|---:|---:|---|---|
| `ch_1` | `/raw/X`, `/X` | `[]` | `25,702 × 18,711` | `121,101,403` | `363,210,965` | `1,018` | `0 / 0 / 0` | `MULTIPLE_IDENTICAL_INTEGER_SLOTS_VERIFIED` (`/raw/X == /X`) | 独立记录 `chosen_slot="/raw/X"` vs shared `x_slot_used="X"` (`int_slots=["X","raw/X"]`)；两槽位逐元素完全相同 |
| `ch_2` | `/raw/X`, `/X` | `[]` | `27,809 × 18,808` | `130,604,474` | `402,730,416` | `1,059` | `0 / 0 / 0` | `MULTIPLE_IDENTICAL_INTEGER_SLOTS_VERIFIED` (`/raw/X == /X`) | 独立记录 `chosen_slot="/raw/X"` vs shared `x_slot_used="X"` (`int_slots=["X","raw/X"]`)；两槽位逐元素完全相同 |

---

## 3. 逐细胞、分通道、NTC 与 `gene_desert` 排除计数核对 (`§4.2–§4.5, §5`)

| 项目 | 独立重建（`independent_r14_loader`） | `round14/shared/lines/a549_myllia` (`A-T1`) | 是否完全一致 |
|---|---:|---:|---|
| 原始总细胞数 (`ch_1 + ch_2`) | `53,511` (`25,702 + 27,809`) | `53,511` (`25,702 + 27,809`) | 一致 |
| RNA QC 排除数 (`<300 UMI` or `<200 genes`, `max_pct_mt=None`) | `0` | `0` | 一致 |
| 多 guide / 未分配 / 未映射排除数 | `0` | `0` | 一致 |
| **`ctrl_T` (`gene_desert_1..20`) 排除数合计** | **`2,063`** (`ch_1=988`, `ch_2=1,075`) | **`2,063`** (`ch_1=988`, `ch_2=1,075`) | **一致** |
| **保留总细胞数** | **`51,448`** (`targeting=49,355`, `NTC=2,093`) | **`51,448`** (`targeting=49,355`, `NTC=2,093`) | **一致** |
| **`ch_1` (`Half A`, batch=`ch_1`) 保留细胞 / NTC / 靶细胞** | **`24,714`** / **`1,033`** / `23,681` | **`24,714`** / **`1,033`** / `23,681` | **一致** |
| **`ch_2` (`Half B`, batch=`ch_2`) 保留细胞 / NTC / 靶细胞** | **`26,734`** / **`1,060`** / `25,674` | **`26,734`** / **`1,060`** / `25,674` | **一致** |
| 低 NTC fallback 细胞数 (`< 15` NTC/batch) | `0` | `0` | 一致 |

---

## 4. 基因宇宙、可重复靶点集合与评测靶点集合核对 (`§5–§7`)

| 项目 | 独立重建 | `round14/shared/lines/a549_myllia` | 集合与顺序比较 |
|---|---:|---:|---|
| 跨通道交集 `raw_shape`（`n_genes_multifile_dropped`） | `[51,448, 18,480]` (`231`) | `[51,448, 18,480]` (`231`) | 完全一致 |
| `G_final`（HGNC 规范化且两通道有限，`genes.txt`） | **`15,456`**（排除未映射 `379`） | **`15,456`**（排除未映射 `379`） | `genes_exact_order_match = true` |
| 与 `r9_6` 共同基因组对齐 `|G_eval|` | **`6,126`** | **`6,126`** | 完全一致 |
| 去重后靶点组数（`raw 218 -> >=10 cells`；`half_valid`） | `212`（`212 / 212`） | `212`（`212 / 212`） | `targets_exact_order_match = true`; `half_valid` 0 差异 |
| **可重复靶点集合（`reproducible == True`）** | **`81 / 212`** | **`81 / 212`** | **0 差异（完全相同集合）** |
| **评测靶点集合（`eval_targets`, `r9_6` + `self ∈ G_eval`）** | **`32`** (`< 50`, `K_CONF = 0`) | **`32`** (`< 50`, `K_CONF = 0`) | **0 差异（完全相同集合）** |
| `§7` 重叠亚组 `overlap_r12_a549`（数量与集合） | **`4 / 15`** | **`4 / 15`** | `[CREBBP, HDAC3, KDM6A, PTEN]` 完全一致 |
| `§7` 新增亚组 `novel_not_in_r12_a549`（数量与集合） | **`28 / 46`** | **`28 / 46`** | 完整 28 基因集合完全一致 |

- **81 个可重复靶点完整集合（字母序，两边完全相同）**：`AMOTL2, ATRX, ATXN7L3, BRD1, BRD4, BRPF1, CBX3, CDK13, CMTR2, CNOT4, CREBBP, DDX28, DGKD, DHX29, DHX36, EED, EGFR, EHMT1, EP300, FIG4, FOSL1, HDAC1, HDAC3, HOXC13, HSD17B10, INTS14, IRF2, IRF2BP2, KAT6A, KDM1A, KDM3B, KDM6A, KDM8, KLF5, KMT5B, LDB1, MAPK14, MEAF6, MIB1, MRM2, MRM3, MSL1, MTERF4, MTPAP, NCOR2, NF2, NFIB, NME6, NSD1, NSUN4, PCGF1, PIK3R4, POLG, PRORP, PTEN, PTPN11, RAF1, RBBP7, RING1, RNF146, RPUSD3, RPUSD4, RRAGC, SHOC2, SIRT1, SMAD4, SMARCA5, SMARCD1, SRFBP1, SUPT20H, SUPV3L1, SUZ12, SYVN1, TAOK1, TFPT, TRAF3, TRIM28, TRMT10C, TRMT61A, UBA5, UBE3D`。
- **32 个评测靶点完整集合（字母序，两边完全相同）**：`BRPF1, CBX3, CREBBP, DHX29, DHX36, EHMT1, EP300, HDAC1, HDAC3, HSD17B10, INTS14, IRF2BP2, KAT6A, KDM1A, KDM3B, KDM6A, MSL1, MTPAP, PIK3R4, PRORP, PTEN, RAF1, RPUSD3, RPUSD4, SHOC2, SMARCA5, SRFBP1, SUPT20H, SUPV3L1, TFPT, TRMT10C, UBA5`。

---

## 5. `ctrl` 均值与 `delta` 数组逐元素数值对比

| 数组（`.npz` 键） | 形状 | `chunk_size=8192` 最大绝对差 | `<= 1e-6` | 备注（含初始 `chunk_size=4096` 对照发现） |
|---|---|---:|:---:|---|
| `ctrl_cells.npz::ctrl_mean_libnorm` | `(15456,)` | `0.0` | `PASS` | `chunk_size=4096` 与 `8192` 均为 `0.0` |
| `ctrl_cells.npz::ctrl_mean_libnorm_full` | `(15456,)` | `0.0` | `PASS` | `chunk_size=4096` 与 `8192` 均为 `0.0` |
| `ctrl_cells.npz::ctrl_mean_log1p` | `(15456,)` | `0.0` | `PASS` | `chunk_size=4096` 与 `8192` 均为 `0.0` |
| `ctrl_cells.npz::X_log1p` | `(2093, 15456)` | `0.0` | `PASS` | `chunk_size=4096` 与 `8192` 均为 `0.0` |
| `pseudobulk_delta.npz::delta_log1p` | `(212, 15456)` | `0.0` | `PASS` | `chunk_size=4096` 时曾为 `8.252e-07`（因 `process_r10_line` 内部 `float32` 分块累加顺序差异） |
| `pseudobulk_delta.npz::half_a_delta_log1p` | `(212, 15456)` | `0.0` | `PASS` | `chunk_size=4096` 时曾为 `1.427e-06`；改用 `prepare_r10_lines` 默认 `chunk_size=8192` 后完全为 `0.0` |
| `pseudobulk_delta.npz::half_b_delta_log1p` | `(212, 15456)` | `0.0` | `PASS` | `chunk_size=4096` 时曾为 `1.798e-06`；改用 `prepare_r10_lines` 默认 `chunk_size=8192` 后完全为 `0.0` |
| `pseudobulk_delta.npz::knockdown_depth_log1p` | `(212,)` | `0.0` | `PASS` | `chunk_size=4096` 时曾为 `7.451e-08` |
| `pseudobulk_delta.npz::n_cells` | `(212,)` | `0.0` | `PASS` | 全整数完全一致 |
| `pseudobulk_delta.npz::ad_counts` | `(212,)` | `0.0` | `PASS` | 全部 212 靶点 AD 计数完全一致 |
