# Round 13 第二独立细胞系增补筛选报告（`R13_DATA_SCREEN.md`）

- **执行日期**: `2026-11-23`
- **执行智能体**: `A-T1`（纯元数据排查；全程未下载任何 `>100 MB` 文件、未下载任何表达矩阵、未注册账号、未同意条款/EULA、未执行 `git`、未修改 `round12/` 或 `PREREG_E`）
- **基线参照**: 接续 [round11/shared/R12_DATA_SCREEN.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R12_DATA_SCREEN.md) 与 [round12/shared/R12_QC_COUNTS.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_QC_COUNTS.md)
- **已用 / 明令排除集合**: R9–R12 已用系（`K562`, `RPE1`, `HepG2`, `Jurkat`, `H1`, `KOLF2.1J`/`KOLF2.2J`, `DLD1`, `Calu3`, `CD4+ T [rest & stim48]`, `HipSci iPSC`, `TeloHAEC`, `HCASMC`, `WTC11 iPSC/i3Neuron`, `A549`, `HEK293`）、全部 `X-Atlas/Orion`（含 Huang et al. `HCT116`/`HEK293T`）、全部 `VCC`、以及全部受控访问（`dbGaP`/NDA 等）。

---

## 一、总结论

> **结论：没有。**

截至 **`2026-11-23`**，在 `round11/shared/R12_DATA_SCREEN.md` 已排雷基础上，对 **`2026-11-22` 之后**（并全面回扫 2026 下半年及近期预印本）公开释放的 **`GEO`、`IGVF`、`scPerturb`、`ArrayExpress/BioStudies`、`CELLxGENE`、`Zenodo`、`Figshare`、`bioRxiv/medRxiv/EuropePMC`、`HuggingFace`** 元数据核查，以及对两项 R12 边际悬案（**`GSE327644`** 原代人巨噬细胞、`IGVF HUES8`）和 2026 新近高频误报候选（**`GSE301119`** 原代人巨噬细胞、`Perturb-ME` `A375`、`VIPerturb-seq`、`GSE343369` `KOLF2.2J`、`GSE320250` 小鼠 HSC、`Figshare 34002315` HKU 汇集库、`IGVF H9/HFF-1/Hs27` 等）的深度核实表明：

**目前不存在任何可在无受控权限前提下稳定通过 `PREREG` 确认性门槛（`n_eval_targets >= 50`、`|G_eval| >= 3,000`、per-half `NTC >= 15`、低 MOI 单 sgRNA CRISPRi/CRISPRko 全转录组、非方法学严重偏倚）的新增独立未用人类细胞系。**

- 唯一在乐观校准下可触及 `>= 50` 的未用人类原代细胞候选为 **`GSE327644`**（`SLICeVLP` 原代人巨噬细胞 CRISPRko，`569` 靶点），但其 **`u6_rep = 108`（`u6_rep ∩ R9_COMMON = 92`）**：按通用 **`20%–40%`** 经验仅 **`21.6–43.2`（`R9_COMMON` 内 `18.4–36.8`）**、按 R12 保守 **`HEK293` 校准（`9.76%`）** 仅 **`10.5`（`R9_COMMON` 内 `10.3`）**，均不达 `50`；仅按 R12 最乐观 **`A549` 深度校准（`83.96%`）** 才得 **`90.7`**，且元数据明确存在 **单供体（`RG5581`）、仅 4 通道 10x Flex（`Ctrl` 与 `LPS` 混池、靶向探针面板）、未披露单 sgRNA 细胞数与 MOI、且仅提供 `6.58 GB` 原始 `.tar`（无轻量 QC `.h5ad`）** 等重大风险，不足以作为可锁定的第二确认线。

---

## 二、R12 双校准比率定义与基准参数

沿用 R9 冻结参考集合（经 [verify_r10_lines.py](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round10/shared/scripts/verify_r10_lines.py) 与 HGNC 冻结表在本轮复核）：`union6_rep (u6_rep) = 2,117`、`R9_COMMON = 6,148`、`u6_rep ∩ R9_COMMON = 1,698`。

每个候选均并列报告以下三套 `n_eval_targets` 估算（不单列乐观值）：
1. **R10–R11 经验 `20%–40%` 规则**：`0.20 × u6_rep .. 0.40 × u6_rep`（及 `u6_rep ∩ R9_COMMON` 的 `20%–40%`）。
2. **R12 乐观校准（`A549` 深测序/多通道体制）**：R12 实测 `a549_ko` 为 `445 / 530 = 83.96%`（占 `u6_rep`），或 `445 / 452 = 98.45%`（占 `u6_rep ∩ R9_COMMON`；见 [R12_QC_COUNTS.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/shared/R12_QC_COUNTS.md)，保留单细胞总数 `605,484` 细胞 / `999` 靶点、NTC `45,320`）。
3. **R12 保守校准（`HEK293` 浅覆盖 + 半数批减均值体制）**：R12 实测 `hek293_xucao` 为 `16 / 164 = 9.76%`（占 `u6_rep`），或 `16 / 143 = 11.19%`（占 `u6_rep ∩ R9_COMMON`；保留 `98,274` 细胞 / `200` 靶点、NTC `2,757`、`half_1/half_2` per-half NTC 中心化）。

---

## 三、候选数据集与悬案深度核查表

> 下表覆盖两项指定核查标的（`GSE327644`、`IGVF HUES8`）及本轮全部新检出或易混淆的高关注候选；已于 [R12_DATA_SCREEN.md](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round11/shared/R12_DATA_SCREEN.md) 详述排除的数据集不重复列为主体。

| 候选 / 数据集 | 访问路径 | 细胞系 | 扰动方式 | 靶点数 | 单 sgRNA 细胞数 | NTC 数 | `u6_rep` / `u6_rep ∩ R9_COMMON` | 评测靶点估算（`20–40%` \| `A549 84%` \| `HEK293 10%`） | 文件大小 | 结论 |
|---|---|---|---|---:|---|---|---|---|---|---|
| **`GSE327644`** (`SLICeVLP`, Jung et al. *Nat Biotechnol* `2026-08-17`) | [GEO GSE327644](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327644) / PMID `42608566` / `PRJNA1450051` | **原代人巨噬细胞** (`CD14+` 分化；单供体 `RG5581`；未用系) | **CRISPRko** (`Cas9-eVLP + VPX-lenti sgRNA`，非 GEO 特性误标的 `CRISPRi`；10x GEM-X Flex；`Ctrl` 与 `50 ng/ml LPS` 混池) | **`569`** 原始唯一基因 (`569` HGNC；`2,332` targeting guides；论文摘要称 `583`) | **未在元数据/摘要披露**（未解压 `6.58 GB` `.tar`；仅 4 GEO 样本 `Mac_L01..L04` 共 8 个 `.h5`，每样本含 `Mac_Ctrl` 与 `Mac_LPS`） | guide 模板含 **`593`** 条非靶向对照 sgRNA；细胞数未披露 | **`108`** / **`92`** | 按 `u6_rep=108`：<br>`20–40%`: **`21.6–43.2`**<br>`A549 (83.96%)`: **`90.7`**<br>`HEK293 (9.76%)`: **`10.5`**<br>*(按 `R9C=92`：`18.4–36.8` \| `90.6` \| `10.3`)* | `guide_template.txt.gz`: `38,166 B`；`GSE327644_RAW.tar`: `6,579,169,280 B` (`6.58 GB`) | **边际高风险淘汰 / 不可定为可靠二线**（`20–40%` 与 `10%` 体制均 `<50`；仅在最乐观 `84%` 下可达 `~91`，但受单供体、`Ctrl+LPS` 混样、Flex 探针、未核实 MOI 与单 sgRNA 深度多重制约） |
| **`IGVF HUES8` Perturb-seq** (`IGVFDS0956UUQN` 等 97 组) | [IGVF API](https://api.data.igvf.org/search/?type=MeasurementSet&samples.sample_terms.term_name=HUES8)（最新释放 `2026-11-05..06`） | **`HUES8`** hESC（及定型内胚层诱导 `48h`） | `Sp-dCas9-KOX1-KRAB` CRISPRi **高 MOI Enhancer/TSS 调控元件** | 针对开放染色质区/TSS/TF 调控位点，非低 MOI 单基因库 | `0`（**无任何已处理矩阵或 per-cell guide 注释**；仅 `FASTQ + seqspec YAML`） | 无已处理对象（`derived_analysis_sets = None`） | N/A（非编码顺式元件为主；`HUES8` 的 28 个 `AnalysisSet` 仅 ATAC/ChIP/RNA-seq） | N/A（`MOI = 10–12` 违反低 MOI 单 sgRNA `P4` 门槛） | 仅原始 FASTQ（`AnalysisSet = 0`） | **淘汰**（无已处理矩阵/单细胞 guide 分配，且 `MOI = 10` 或 `12` 为高 MOI 非编码增强子筛选） |
| **`GSE301119`** (Wang et al. *iScience* **`2026-11-06`**, `PMC13760980`) | [GEO GSE301119](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301119) / [PMC13760980](https://pmc.ncbi.nlm.nih.gov/articles/PMC13760980/) | **原代人巨噬细胞 `hMΦ`**（`Donor1`, `Donor2`；未用系） | **CRISPRi** (`AdV dCas9-KRAB + lenti sgRNA`) + 平行 **CRISPRa**（排除） | **`203`** 基因（`5` sgRNA/基因；`glmGamPoi` 仅 **`86`** 个有显著靶点 KD） | CRISPRi 合计 **`28,446`** 细胞（含多 guide 按最高 UMI 强制归因），**单 guide 比例仅 `33.7%`（约 `9,586` 个严格单 sgRNA 细胞**，均摊仅约 `~41` 细胞/基因、`~20`/half） | 库内含 **`100`** 条 NT sgRNA（未单独披露单 guide NTC 细胞数） | 靶单 `mmc2.xlsx` 受 PMC 网页 challenge 阻断未直取；`203` 基因全库即使按高估计 `30–40%` 重合仅 `u6_rep ≈ 60–81`（显著 KD 仅 `86` 个更低） | 即使取理想上界 `u6_rep=80`：<br>`20–40%`: **`16–32`**<br>`A549 (84%)`: **`67`** *(受 `~20`/half 浅度阻断不可能达标)*<br>`HEK293 (10%)`: **`7.8`** | 仅 `GSE301119_CRISPRi_seurat5.rds` (`406.8 MB`) 与 `CRISPRa` `.rds` (`345.8 MB`)；无 `.h5/.h5ad` | **淘汰**（单 sgRNA 细胞覆盖极浅 `~9,586` 分配给 `203` 基因+`100` NT，远低于 `n_cells >= 25`/half + AD/Energy 门槛；且仅有 `86` 基因显著 KD、无 `.h5/.h5ad`） |
| **`Figshare 33264225` `Perturb-ME`** (`10.64898/2026.08.17.745330`) | [Figshare 33264225](https://api.figshare.com/v2/articles/33264225) + [Genentech/PerturbME](https://github.com/Genentech/PerturbME) | **`A375`** 黑色素瘤（未用系；IFNγ 处理后测表面 HLA） | **CRISPRko** (`Perturb-CITE-seq`, `MOI < 0.2`)，**FACS 极端表型富集（bottom/top 5% HLA-I）** | 模型保留 **`221`** 基因（`target_modules.csv`；全基因组文库但多数未进模型） | 全 QC `353,091` 细胞；`MOI=1` 单 guide **`277,152`** 细胞（**但 20 个 10x 通道中 18 个为 `HLA-low/high` FACS 分选，仅 `2` 个未分选 `Control`**） | 包含 `HIGH_CTL` / `LOW_CTL` 与未分选对照通道（按通道比例无偏 `Control` 全基因组平均仅 `~1–2` 细胞/基因） | 对 `221` 模型靶点精确计算：**`77`** / **`66`** | 按 `u6_rep=77`：<br>`20–40%`: **`15.4–30.8`**<br>`A549 (83.96%)`: **`64.6`**<br>`HEK293 (9.76%)`: **`7.5`**<br>*(按 `R9C=66`：`13.2–26.4` \| `65.0` \| `7.4`)* | `PerturbME_A375_HLA_MOI1.h5ad`: `1,430.6 MB` | **淘汰（方法学 + 量级双不合规）**（`90%` 细胞经 HLA-I 极端 5% FACS 门控引入强烈碰撞/选样偏倚；且 `20–40%` 与 `10%` 体制均 `<50`） |
| **`Zenodo 18460279` `VIPerturb-seq`** (Bradu et al. `10.64898/2026.02.12.705613`) | [Zenodo 18460279](https://zenodo.org/api/records/18460279) + [Figshare 33719593](https://api.figshare.com/v2/articles/33719593) | 图 3 全基因组与 `vimentin` 筛选为主系 **`K562`**（已用系！）；仅图 2 mixing pilot 混 `K562+HEK293+HAP1` | **CRISPRi** (`GuEST-List` 分裂探针 + 10x Flex) | `manifest`: **`18,866`** (`c>=20`: `17,475`; `c>=50`: `8,185`; `filtered`: `6,722`) | `filtered.rds`: **`326,247`** 细胞 | `Figshare 33719593` `run.log` 实测 `genome_wide_filtered.rds` 中 **仅 `66` 个细胞匹配对照模式**（`<2,000`，触发伪对照）；`binA/B/C` 含全量 NTC 但总数亦未达常规大规模基线 | `c>=20`: **`1,541` / `1,233`**；`filtered`: **`1,158` / `958`** | 即使仅作参考（不适用）：<br>`filtered (1,158)`: `231.6–463.2` \| `972.3` \| `113.0` | `genome_wide_filtered.rds`: `3.61 GB`; `binA/B/C`: `10.23 GB` | **淘汰 / 排除**（图 3 全基因组筛选与 `vimentin` 筛选均在 R9 已用 **`K562`** 上进行；且 `genome_wide_filtered.rds` 实质缺乏真实 NTC 对照） |
| **`GSE343369`** (`UCSF Study 2`, MorPhiC R4.0) | [GEO GSE343369](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE343369) (`2026-08-12`) | **`KOLF2.2J`** iPSC/EB（与 R9 `nourreddine_kolf21j` 同源亚克隆，已排除） | **CRISPRi**（低 MOI，双 guide/载体） | `229` 规范化靶基因 (`perturbed_gene_list.csv.gz`) | 4 GEO 样本，16 块 `.h5ad` | 含 NTC | **`28`** / **`20`** | 按 `u6_rep=28`：<br>`20–40%`: **`5.6–11.2`**<br>`A549 (83.96%)`: **`23.5`**<br>`HEK293 (9.76%)`: **`2.7`** | `GSE343369_perturbed_gene_list.csv.gz`（小表）+ 多块 `.h5ad` | **淘汰**（同源 `KOLF2` 已用系，且 `u6_rep = 28` 在所有校准下均 `<50`） |
| **`Figshare 31276981` / `GSE320250`** (Bowness & Velten `2026-08-26`) | [Figshare 31276981](https://api.figshare.com/v2/articles/31276981) / [GEO GSE320250](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320250) | **小鼠原代 HSC**（`Mus musculus`, bone marrow `Lin- cKit+`） | CRISPRi (`CROP-seq`) | `520` 小鼠基因 | `111,003` + `19,702` 细胞（`.RDS`） | 含对照 | N/A（非人类） | N/A（物种排除） | `10.59 GB` + `977.5 MB` | **淘汰**（GEO 与正文确认为 **小鼠 `Mus musculus`**，非人类细胞） |

---

## 四、全平台 **`2026-11-22` 之后** 更新与常见陷阱系统清查记录

本轮对各指定门户截至 **`2026-11-23`** 的新增公开及近期 2026 年底条目逐一核实如下：

1. **IGVF Portal（`api.data.igvf.org`）**：
   - 全库 `AnalysisSet` 共 `3,877` 条、`MeasurementSet` 共 `7,938` 条；**最新释放时间戳分别为 `2026-11-06T20:52:10Z` 与 `2026-11-06T21:00:06Z`**，**`2026-11-22` 之后新增公开数为 `0`**。
   - 核查全库所有已释放 Perturb 类 `AnalysisSet`（`Perturb-seq` 26 组、`CC-Perturb-seq` 24 组、`TAP-seq` 9 组）：
     - `TeloHAEC`（1+15）、`GM25256 (WTC-11)`（13+9+4）、`K562`（2+1）、`H1-hESC`（1）、`HCASMC-hTERT`（2）、`CD8+`（2）：均为已用系、靶向可及染色质调控元件或 CRISPRa / 靶向 TAP-seq；
     - **`H9`**（5 组 `Perturb-seq`：`IGVFDS6924DJAZ`, `2266YDVM`, `2873EPIU`, `9825FOBP`, `4762GRIW`，诱导心肌或神经元）与 **`HFF-1`**（1 组 `IGVFDS3405NGXF`，`MOI=1`）：经样本描述核实全部为 ***"targeting accessible genome regions in many genomic loci"***（非编码调控元件，非基因级 CRISPRi/KO）；
     - **`Hs27`**（1 组 `IGVFDS2001NDKP`，`2026-08-24`）：经样本描述核实为 **`CRISPR activation`（`Sp-dCas9-VPH` CRISPRa）**，属于模态排除。
     - **`HUES8`**（97 组 Perturb `MeasurementSet`s）：`derived_analysis_sets = None`，无矩阵、且 `MOI = 10–12` 针对开放染色质区（详见第三节表）。

2. **NCBI GEO（`eutils` + `acc.cgi`）**：
   - **`2026-11-22` 至 `2026-11-23`**：`"Perturb-seq"` / `"CROP-seq"` / `"CRISPRi"` / `"CRISPRko"` / `"CRISPR + single cell"` 人类新释出 **`0` 条**。
   - **`2026-11-19..2026-11-21` 临界新增排查**：
     - `GSE349301` (`2026-11-20`): 呼吸道 `Calu-3`/`HSAEC1-KT` 干扰素筛选，无低 MOI 单细胞全转录组 Perturb-seq 矩阵；
     - `GSE333994` (`2026-11-20`): `BxPC-3` 混合 `Cas12a` 生长耐药筛选（非 scRNA-seq）；
     - `GSE277332` / `GSE277333` (`2026-11-20/21`): `K562`（已用系）图像 CRISPR 筛选配套 RNA-seq；
     - `GSE354700` / `GSE354698` (`2026-11-19`): `TF-iMG` 小胶质细胞，仅 **19 个转录因子阵列筛选（arrayed Perturb-seq）**，`R12_DATA_SCREEN.md` 已核查淘汰。
   - **2026 其他高误报 GEO 条目核实（均在元数据核验后排除）**：
     - `GSE295214` (`Public May 05 2026`, `85` GSMs, `Pan-GW-hESC`): GSM 样本属性核实 `cell line: H1 hESC`（R9 已用，排除）；
     - `GSE299590` (`Public Jun 11 2026`, `WND238` 黑色素瘤 `235` 表观调控子): GSM9041534 核实 `!Sample_molecule_ch1 = genomic DNA`（gDNA `MAGeCK` 混合测序，并非单细胞 scRNA-seq！）；
     - `GSE253842` (`Public Jun 15 2026`, `SH-SY5Y` `KRAB-dCas9`): 关联论文 `PMID 42374457` (`PMC13621727`) 实测仅做 **7 个基因** 的神经元 CRISPRi（远 `<50`）；
     - `GSE352893` (`Nov 02 2026`): `WTC11` + `K562`（均已用系），且标称 **`MOI 5 / 7`**；
     - `GSE349414` (`Nov 16 2026`, `PEP-seq`): **`HEK293T-PE7`**（已用系 + 引导编辑 `pegRNA` 非 CRISPRi/KO）；
     - `GSE347104` (`Oct 05 2026`): **`hCASMC`**（已用 `HCASMC`）阵列非编码增强子 CRISPRi；
     - `GSE349107` / `GSE349277` (`Oct 12 2026`): `WTC11` hiPSC（已用系）及 1 基因 (`MEIS2`) bulk RNA-seq；
     - `GSE342892` (`Oct 14 2026`): `K562`（已用系）配对 `pegRNA` 删除 100 cCREs；
     - `GSE279287` (`Feb 10 2026`, `786-O`/`Caki-1`): **高 MOI** CROP-seq 解析 GWAS `rs28684409` 等增强子；
     - `GSE349847` / `GSE283614` (`Oct 2026`): `hESC` 剂量化 `OCT4`/`NANOG` 少量增强子 Perturb-seq。

3. **ArrayExpress / BioStudies（`collection=arrayexpress`）**：
   - 全库最新释放日期 **`2026-11-22`**（`E-MTAB-17812`，脐带血 RRBS）；**`2026-11-22` 之后 `0` 条新增，2026 全年亦无人类低 MOI Perturb-seq/CROP-seq 新增**。`E-MTAB-14567`（`MCF-7` 仅 50 个基因、`2025-11-19`）与 `E-MTAB-14004`（酿酒酵母 `S. cerevisiae`）及 `E-MTAB-13064/13066`（`SUM159PT` 无 CRISPR 的谱系条形码追踪）、`E-MTAB-13204`（2023 年 NK 共培养有限靶点 CROP-seq）均不合条件。

4. **CELLxGENE（`api.cellxgene.cziscience.com`）**：
   - 全库 `408` 个公开合集中，匹配 `perturb`/`crop-seq`/`crispri`/`crisprko` 的合集**最新修订停留在 `2026-09-28`**（`Brain vascular single-cell multi-omics`），**`2026-11-22` 之后 `0` 更新**，无任何新增可用人类 CRISPRi/KO Perturb-seq 数据集。

5. **Zenodo & `scPerturb`**：
   - `scPerturb` 官方档案 (`Zenodo 13350497`) 仍停留于 **`2024-08-26` (`v1.4`)**，无更新。
   - Zenodo 全库 **`>= 2026-11-22`** 搜索 `Perturb-seq`/`CROP-seq`/`CRISPRi`/`CRISPRko`/`h5ad`：**`0` 条真实人类单细胞 CRISPR 数据集**；11 月上中旬新上传 `24047293` / `23470845` 为 **`CD4+ T`**（R9 已用系），`21725795` / `21282294` 为 **`K562`**（R9 已用系），`18557106`（`HT29`）为 `MAGeCK` 混合 gDNA 筛选而非 scRNA-seq。

6. **Figshare（`api.figshare.com`）**：
   - 按发表时间倒序检索 `Perturb-seq`/`CRISPRi`/`CROP-seq`：**`>= 2026-11-22` 无新人类 Perturb-seq 数据**（`2026-11-21` `35082220` 为十字花科油料作物 *Thlaspi arvense* CRISPR；`2026-11-13` `34705101` 为增强子兼容性）。
   - `34002315`（HKU `2026-10-13` 统一重处理 12 个 CRISPRi 屏幕）`README.md` 逐项核实：全部 12 个屏幕均出自 **`K562`, `RPE1`, `Jurkat`, `HepG2`, `X-Atlas/Orion (HCT116, HEK293T)`, `HipSci iPSC`, `H1 hESC`, `VCC`, `KOLF2.1J`** —— **100% 落在排除集合内**。
   - `30898151` (`OASIS`, `2026-10-13`) 实测对应 `KOLF2.1J_Perturbation_Cell_Atlas`（已用 `KOLF2.1J` 且为可溶性配体拴系筛选，非 CRISPRi/KO）；`29557997` (`SYM10X1`) 为 **CRISPR-Cas13**（仅 15 个靶点）。

7. **bioRxiv / medRxiv / EuropePMC / HuggingFace（`2026-11-18..2026-11-24`）**：
   - 完整翻页扫描 `bioRxiv` 全部 `899` 篇与 `medRxiv` 全部 `271` 篇新增预印本：`medRxiv` 命中 **`0` 篇**；`bioRxiv` 命中 `8` 篇，均为基准评测或重分析（如 `10.64898/2026.11.18.767605` `ImmunoTwin` 重分析已用 `CD4+ T` 与巨噬细胞屏幕）、已排除的 `TF-iMG` 阵列屏幕（`GSE354700`）、`FULL-seq` 多组学方法或细菌/植物论文，**无新的公开可用低 MOI 人类 CRISPRi/CRISPRko Perturb-seq 独立细胞系**。
   - HuggingFace（`syntensor/cellsim-crispri-screens` 为合成评测任务；`Boom5426/iPSC-CRISPRi-PerturbSeq` 为已排除的 Feng 2026 `HipSci iPSC` 衍生效应表；`weililab/CRISPRI_50TF` 仅 `50` 个 TF，总靶点数即使全落入 `u6_rep` 也仅 `50`，无法通过叠加 R9 重合率、`half_valid` 与 `AD/Energy` 筛选后仍 `>=50` 的硬门槛）。

---

## 五、核心缺口与给后续轮次的客观建议

1. **为什么看似存在多个“未用系”实则无一达标**：
   - 一类（如 **`IGVF HUES8 / H9 / HFF-1`**、**`GSE279287`**）属于 **高 MOI 调控元件/增强子筛选**，无法满足低 MOI 单 sgRNA 基因级干预契约；
   - 一类（如 **`Perturb-ME` `A375`**）虽为单 sgRNA CRISPRko 且名义上有 221 模型靶点（`u6_rep = 77`），但 **90% 通道经 FACS 极端 5% 蛋白表达分选**，会系统性扭曲下游转录组协方差与 DE 估计，且按 `20–40%` 或 `HEK293` 比率仍远低于 `50`；
   - 一类（如 **`GSE301119`** 原代巨噬细胞、`u6_rep` 极低或 **`28,446` 细胞仅 `33.7%` 单 guide**）测序深度与单 sgRNA 计数远低于 R9/R10/R12 Addendum-B 稳定性需求；
   - 唯一尚存理论空间的 **`GSE327644`**（原代人巨噬细胞 `SLICeVLP`, `u6_rep = 108`）必须假设 **A549 级别的极高存活率（`84%` → `90.7`）** 方能跨过 `50`，而按常规 **`20–40%`（`21.6–43.2`）** 或 R12 真实浅层样本 **`HEK293`（`9.76%` → `10.5`）** 均会失败，且单供体、`Ctrl+LPS` 混样与未公开 MOI 使风险过高。

> **最终建议**：基于截至 **`2026-11-23`** 的全面排查，**不存在合格的公开第二独立细胞系**；如需继续推进独立确认，仅应在已知 R12 结果的基础上如实报告当前 `K_conf = 1 (a549_ko)` 的跨系泛化边界，或等待真正大规模低 MOI 单 sgRNA 编码基因 Perturb-seq 数据（如未曾使用过的新细胞系达到 `>= 300` 靶点且 `u6_rep >= 200`、多通道深测序、无极端表型 FACS 门控、含充足真实 NTC）公开发布后再行立项；若后续拟对 `GSE327644` 做探索性深度核实，则必须预先向人类/主控说明其极高失败概率（`10.5–43.2` vs `90.7`）且不承诺能升格为确认线。
