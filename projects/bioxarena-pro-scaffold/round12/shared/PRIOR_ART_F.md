# PRIOR_ART_F：Round 12 既有工作核查报告（A-T5）

- 核查日期：2026-11-22
- 核查范围：针对 `R12_第十二轮行动指令.md` 第 4 节第 4 项要求的两条文献既有工作检索问题，通过公开文献接口（bioRxiv API、Crossref、Semantic Scholar、EuropePMC、NCBI E-utilities GDS、arXiv API、Zenodo API、GitHub raw 文件）及已下载到本地缓存目录的全文与代码仓进行一一手核验。
- 数据与合规边界：未执行任何 `git` 命令；未访问、下载或读取 `X-Atlas/Orion`、`VCC`、`GSE345058` 或 `XuCao2023` 的任何表达矩阵；遇到 NCBI GEO 网页前端 Google reCAPTCHA 与 DBLP Anubis 验证页时已立即停止，未作任何绕过，仅采用官方 API；全文不含禁止词。

---

## 一、核心结论总览

| 编号 | 核查问题 | 三档结论 | 一句话依据 |
| --- | --- | --- | --- |
| (1) | `Liu et al. 2026`（`GSE345058`，bioRxiv `10.64898/2026.06.01.728087`）或任何其他已发表/预印本工作是否使用该 **A549** 数据集做过跨细胞系转录组扰动预测 | **没有找到** | 原始研究与配套代码仓仅在 A549 内部比较 scRNA-seq（`sVAEplus`）与光学成像表型分辨率；Semantic Scholar 全部 4 篇施引文献、Crossref/EuropePMC/arXiv/bioRxiv 及本地库均未见使用 `GSE345058` 做跨细胞系转录组扰动预测 |
| (2) | 此前是否有研究在扰动预测评估中做过“目标系自身基线 vs 错配细胞系基线”的配对对比（对应本轮每靶点配对比较 `FCT − Mis_s`） | **部分做过** | `CanPert/CASTER`（`10.64898/2026.11.15.766815`）实现并计算了 `CASTER`（目标对照输入）vs `CASTER_wrongctrl`（错配他系对照输入）的逐靶点配对 Wilcoxon 检验；但该设计锁定把预测 shift 加回目标自身对照表达（不改变基线表达水平）、检验对象是神经模型输入消融而非独立基线对照（`FCT − Mis_s`）、每靶点系仅替换 1 个循环邻居而非全源系 Intersection-Union Test，且不使用 `r_orth_cal` |

---

## 二、问题 (1) 核查：`Liu et al. 2026`（`GSE345058`）与 A549 跨细胞系转录组扰动预测

**判断：没有找到**

### 1. 原始研究与随附代码库核验

- **论文与持久标识符**：
  - 标题：*A multimodal perturbation atlas defines the phenotypic resolution of cellular morphology*
  - 作者：Chad Liu, Alexander Hillsley, Madhurya Sekhar, ..., Manuel D. Leonetti（Biohub, San Francisco）
  - bioRxiv DOI：`10.64898/2026.06.01.728087`（bioRxiv API 核实：v1 发布于 `2026-06-02`，v2 发布于 `2026-08-25`；EuropePMC `PPR1244083`；OpenAlex `W7163223070`）
  - GEO 登录号：`GSE345058`（NCBI E-utilities GDS `esummary` `uid=200345058` 核实：`Homo sapiens`，公开日期 `2026/08/26`，含 CROP-seq 1,000 gene 批次样本）
  - 随附代码与汇总数据：GitHub `czbiohub-sf/ops-paper-analysis`、`czbiohub-sf/cyclops-monorepo`；Zenodo `10.5281/zenodo.21970895`（`altair_v2.zip`，precomputed figure tables 与基因面板参考数据）；交互门户 `https://biohub.ai/ops-explorer`
- **一手文本与代码核验要点**：
  1. 论文摘要与 GEO 摘要均明确指出，该图谱构建了 **A549 细胞中 1,000 个混合 CRISPR 敲除扰动**的多模态图谱（42 个活细胞与 13 个固定荧光标记、同一活细胞的无标记定量相差成像及 scRNA-seq，总计约 6,500 万单细胞记录），核心结论为：在同等试剂成本下，定量相差成像的表型分辨力高于荧光成像和 scRNA-seq。
  2. 直接检查 `czbiohub-sf/ops-paper-analysis` 的 `README.md` 与分析脚本目录（`notebooks/figure_1/` 至 `figure_5/`、`notebooks/target_gene_selection/`）：
     - `notebooks/figure_5/rna_ops_metrics.ipynb`、`mAP_complex_KO_scatter.ipynb`、`joint_heatmap.ipynb` 与 `rna_image_confusion_matrix.ipynb` 将 A549 CROP-seq 的 `sVAEplus` 潜在表征（`cropseq_svaeplus.h5ad`）与同一 A549 细胞系的光学成像特征（`ops_phase_only.h5ad` 等）在**单一细胞系（A549）内部**比较表型区分度 mAP、EBI 蛋白复合物一致性 mAP、Leiden 聚类混淆矩阵及 Pearson 相关热图。
     - `notebooks/target_gene_selection/gene_panel_selection.ipynb` 读取公共 `Replogle_PerturbSeq`（K562 与 RPE1）数据仅用于在**实验建库前筛选 1,000 个目标基因面板**，不存在任何从或向 A549（`GSE345058`）进行的跨细胞系转录组扰动预测建模或评估。

### 2. 全部已知引用文献与多数据库检索核验

- **引用检索计数**：Semantic Scholar 记录 `10.64898/2026.06.01.728087` 总被引 **4 篇**（Crossref `is-referenced-by-count: 3`），4 篇施引文献摘要与引用上下文已全部核查：
  1. **Li et al. (2026)**, *"When label-free morphology is sufficient for targeted cellular measurements"* — bioRxiv DOI `10.64898/2026.09.22.753654`（v1 `2026-09-28`, v2 `2026-10-13`; EPMC `PPR1330852`）。引用语境为分析该 A549 光学扰动图谱（*"nearly 10 million same-cell phase–reporter observations spanning 1,000 gene knock-outs, 52 fluorescent reporters and 73 A549 screens"*）及原代肝细胞亮场化合物筛选，研究无标记成像预测荧光报告蛋白的“测量充分性”（`MorphoSuff`）；**未使用 `GSE345058` 的 CROP-seq scRNA-seq 表达数据，亦未做跨细胞系转录组扰动预测**（Semantic Scholar 显示其尚无施引文献）。
  2. **Zheng, Hong, & Huang (2026)**, *"Virtual experiments bridge sequence and microscopy with generative models"* — bioRxiv DOI `10.64898/2026.09.13.751243`。提出蛋白序列—荧光显微图像双向生成模型 `CELL-FM`，仅在引言将 `Liu et al. 2026` 作为光学混合筛选（OPS）背景文献引用（refs 15, 16）。
  3. **Kirby et al. (2026)**, *"A one-week automated genome-wide optical pooled screen using OttoSeq"* — *Genome Biology* DOI `10.1186/s13059-026-04172-1`（`PMID: 42393721`, `PMCID: PMC13621683`）。仅在引言讨论 OPS 实验流程耗时与分析脚本复用时引用（refs 3, 7）。
  4. **Montjean et al. (2026)**, *"Evaluation of a Novel Fluorescent Marker for Continuous Live Imaging of Preimplantation Embryos"* — *J. Dev. Biol.* DOI `10.3390/jdb14030039`（`PMID: 42783547`, `PMCID: PMC13607945`）。小鼠胚胎实时荧光染料评估，与单细胞转录组扰动预测无关。
- **跨库检索与易混文献排除（截至 2026-11-22）**：
  - **Crossref**：全文检索 `"GSE345058"` 得 **0** 命中。
  - **EuropePMC**：检索 `"GSE345058"` 得 **0** 命中；检索 `"10.64898/2026.06.01.728087"` 或 `"ops-explorer"` 仅命中 `PPR1244083`（`Liu et al. 2026` 本身）。检索 `"A549" AND "perturbation prediction" AND ("cross-cell" OR "cross-line" OR "CROP-seq") AND PUB_YEAR:2026` 仅命中 2 篇：
    - **Li et al. (2026)**, *"A systematic comparison of single-cell perturbation response prediction models"*（*Sci. Adv.* `12(37):eaed3414`, DOI `10.1126/sciadv.aed3414`, `PMCID: PMC13557097`）：纳入早期多细胞系药物/CRISPR 筛选与 Virtual Cell Challenge；核查其本地 XML (`li_2026_pmc13557097.xml`) 可知，其引用的 Leonetti 文献是 2023 年 SARS-CoV-2 宿主因子 Perturb-seq 论文（`10.1038/s41467-023-41788-4`, `PMCID: PMC10558542`），而非 `Liu et al. 2026` / `GSE345058`。
    - **Ma et al. (2026)**, *"Enhancing cross-context generalization in drug perturbation prediction with a multimodal conditional diffusion framework"*（*Bioinformatics* `42(7):btag482`, DOI `10.1093/bioinformatics/btag482`, `PMCID: PMC13364675`）：研究跨细胞系**药物**扰动扩散生成，并非 `GSE345058` CRISPR 敲除图谱。
  - **arXiv 与本地预印本库排查**：`export.arxiv.org` 检索 `"GSE345058"` 与 `728087 + A549` 均为 **0** 命中；`arxiv_perturb_nov2026.json`（465 篇）、`arxiv_qbio_nov2026.json`（465 篇）及 `biorxiv_all_nov16_21.json`（1,572 篇）均未见使用 `GSE345058`。需明确排查一篇容易混淆的新近预印本：**Gao (2026)**, *"The Delta Gap: Cellular World Models for Perturbation Response Prediction Are Far Harder Than Global Correlation Implies"*（`arXiv:2611.10658v1`, `2026-11-07`），其 4.9 节声称测试从 `"Replogle K562"` 零样本迁移至 `"Norman A549"`（8 个过表达条件、741 个共享基因），但其全文表格与参考文献 `[2]` 明确写为 **Norman et al. (2019), *Nature* 571, 613–620**（公开 CRISPRa 过表达数据集），并非 `Liu et al. (2026)` / `GSE345058`。本地精读的 `CanPert`、`Molina & Zhang`、`HyperMap`、`PIE`、`X-Cell`、`CGC`、`SBB`、`scPertEval`、`PertEval`、`Miller et al. (2026)`、`Systema`、`OCCAM`（`arXiv:2611.21959`）等亦均未触及 `GSE345058`。

---

## 三、问题 (2) 核查：“目标系自身基线 vs 错配细胞系基线”配对对比

**判断：部分做过**

直接将不经拟合的**反事实基线**构造在目标系自身基线 `c_{t,g}` 与错配他系基线 `c_{s,g}` 之上，并在校准指标 `r_orth_cal` 上对全部源系 `s ≠ t` 执行每靶点配对比较（`FCT − Mis_s`）与 Intersection-Union Test（`min_s Δ_s`），此前文献**尚未做过**；但在神经跨细胞系预测器的输入对照消融、跨个体配对破坏控制与目标系对照基线设定上，已有高度相关的部分先例，必须在立论中精确界定边界。

### 1. 最核心的“部分做过”先例：`CanPert / CASTER`（Kawsar et al., 2026）

- **出处**：Kawsar, M., Dip, S. A., Paul, B. K., Noman, A. A., Zhang, L., & Moni, M. A. (2026-11-17), *"Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line"* — bioRxiv DOI `10.64898/2026.11.15.766815`；官方代码与结果表已本地核查。
- **它实际做了什么（构成“部分做过”的直接证据）**：
  1. 在 `src/canpert/train.py`（`wrong_ctrl_predictor`）、`experiments.py`、`methods_text.py`（`"Controls for per-line effects (iii) Wrong-line control cells"`）、`results_v1.py`、`results_text.py`（Figure 3B、Table S10、v5 全细胞控制分析）中明确定义了 **`CASTER_wrongctrl`（"CASTER, wrong-line controls"）**：将训练好的跨细胞系 `CASTER` 模型输入的对照细胞集替换为**另一个错配细胞系的对照细胞**（循环顺序 `k562 → hepg2 → jurkat → rpe1 → k562`，即 `other = LINES[(LINES.index(target) + 1) % len(LINES)]`）。
  2. 在 `src/canpert/review_report.py`（`controls_summary` 与 `_paired`）中，对 `CASTER`（目标系自身对照细胞输入）vs `CASTER_wrongctrl`（错配细胞系对照细胞输入）执行了**逐扰动靶点配对的双向 Wilcoxon 符号秩检验**，既汇总全部测试敲除（`_paired(wc, "CASTER", "CASTER_wrongctrl", m)`），又分细胞系检验（`wrongctrl_per_line`），并在预计算输出表给出定量落差：
     - `figure_tables/figures/fig6_data_E.csv`（对应 Figure 3B）：每靶点 `pearson_delta` 配对差值的分系中位数（`CASTER − CASTER_wrongctrl`）在 4 个细胞系中均为正值（`k562: +0.0215`, `hepg2: +0.0613`, `jurkat: +0.0287`, `rpe1: +0.2209`）。
     - `figure_tables/tables/t1_fold_means.csv`（seeds 1–2 `pearson_delta` 均大幅下滑：`hepg2` `0.326/0.319 → 0.187/0.181`；`jurkat` `0.182/0.203 → 0.141/0.133`；`k562` `0.197/0.199 → 0.150/0.152`；`rpe1` `0.385/0.390 → 0.136/0.134`）。
     - `figure_tables/tables/calibrated_normalised.csv`（`t1` 经技术重复标定的 `pearson_deltapert_top100`）：`hepg2` `0.638 → 0.169`、`jurkat` `0.373 → 0.120`、`k562` `0.496 → 0.225`、`rpe1` `0.581 → 0.149`。
- **与 Round 12 `FCT − Mis_s` 设计的四点实质边界差异**：
  1. **输出仍然强制锚定目标系自身基线，刻意不改变基线表达水平**：`wrong_ctrl_predictor`（`train.py:393–404`）提取的是错配输入下模型预测的**差值偏移** `shift = base(p, wc).mean(0) - wc.mean(0)`，随后再将其加回**目标系自身的对照细胞**（`return ctrl + shift[None, :]`）。代码注释与方法学原文均明示：*"so only the model's use of the control set changes, not the baseline expression level"*。即它仅剥离神经网络对“对照细胞特征”的条件化利用，而**并未用错配细胞系基线表达 `c_{s,g}` 替换目标系基线表达 `c_{t,g}` 来比较不同基线锚点**。
  2. **对比对象是“模型与自身错配输入消融”，而非“基线与错配基线”**：`CanPert` 中的基线（`PerturbMean`、`LineOffset`、`ContextMean`、`Linear`）本身均未做 `wrong-line` 替换；它未对比两组非学习型/可解释基线构造（如本轮的 `c_{t,g}·FC` vs `c_{s,g}·FC`）。
  3. **仅替换单一循环邻居，非全源系配对 Intersection-Union Test**：每靶细胞系仅测 1 个指定 `other` 细胞系，未对全部 `s ≠ t` 逐一构造 `Mis_s` 并检验最坏情形下限 `min_s (FCT − Mis_s) > 0`。
  4. **指标体系不同**：使用 `pearson_delta`、`discrimination`、`de_auprc` 与标定 `pearson_deltapert_top100` 等，不涉及对跨系共享偏移正交校准后的相关性度量 `r_orth_cal`。

### 2. 其他相关方法学与基准文献核验（对照定位）

1. **`CGC / IGC`（Huang, Wang, Lu / IGC Virtual Cell）**：
   - 论文 DOI：`10.64898/2026.09.30.755756`（CGC：*The recoverable resolution of cellular perturbation-response prediction*）与同系列先行工作 `10.64898/2026.08.21.746243`（IGC）。
   - 在本地代码 `cgc-counterfactual-resolution-main/.../lea_historical_analysis/scripts/level2_lea_controls.py` 中（基于 `GSE207049` 340 例配对 LCL 个体的 EtOH 基线与地塞米松 DEX 响应），作者设置了名为 **`baseline_response_mismatch`** 的阴性对照：在**训练阶段**将训练样本的基线表达打乱错配给其他样本的扰动残差（`x_train = baseline[train][permutation]`, `y_train = real_y_train`），再用固定测试集 `baseline[test]` 评估 Ridge / PCA-Ridge / MLP 模型的总样本指标。这是**训练集配对破坏对照**，而非推理/基线构造阶段比较“目标自身基线 vs 错配细胞系基线”的每靶点配对对比。同一代码仓的 `cgc_tcell_1c/adjudicate.py` 计算的 `matched` vs `mismatched` 余弦差则是**跨扰动匹配 vs 错配扰动**的检索特异性，与细胞系基线错配无关。
2. **Fan et al. (*"Zero-Shot Cross-Cell-Line Perturbation Prediction Needs Controls"* / *"Cross-Cell-Line Perturbation Prediction Needs Controls"*)**：
   - 作者与出处：Xingyu Fan, Jinghao Wang, Chang-Yu Hsieh, Chunbin Gu, Pheng-Ann Heng；NeurIPS 2026 Evaluations and Datasets Track（通过 Web 搜索定位 OpenReview 与 NeurIPS 官方会议索引摘要核实；直接 API 查询 OpenReview 返回 HTTP 403，未读取全文 PDF）。
   - 该工作针对 `Replogle Perturb-seq` 跨细胞系零样本迁移（`K562`, `RPE1`, `HepG2` → 留出 `Jurkat`）指出，总分布指标易被目标系基线细胞状态方差虚高抬升，因而要求基准测试必须配备两类基线对照：**`target-control-copy`**（直接复制未扰动的目标系对照细胞）与 **`shifted-control` 基线**（以源系平均偏移平移目标系自身对照细胞）。其核心比较是**深度学习/单细胞基础模型/流匹配模型 vs 目标系自身对照衍生的简单基线**，据目前可查公开摘要，并未建立“错配细胞系基线”以进行“目标系自身基线 vs 错配细胞系基线”的配对对照检验。
3. **同方向近期细胞系条件化 / 配对评估工作排查（均未发现对应设计）**：
   - **He et al. (2026-11-17)**, *"Cell-state conditioning predicts the expression response to an unseen stimulus but not stimulus-dependent genetic effects"* — bioRxiv DOI `10.64898/2026.11.16.767041`（`biorxiv_all_nov16_21.json` 一手核实）：评估第二代基因组序列模型对干扰素-β 刺激 vs 静息单核细胞单细胞快照的动态条件化预测及 eQTL 特异性，聚焦序列模型条件输入与静息基线限制，非跨细胞系扰动预测的基线错配对比。
   - **Li et al. (2026, `MorphoSuff`)** — `10.64898/2026.09.22.753654`：使用群内同细胞配对打乱控制证实无标记显微图像含超出扰动平均的信息，非细胞系基线对比。
   - **Rethinking Perturbation Prediction Baselines** — Sun, Zhu, & Chen (ICLR 2026 Workshops, OpenReview `ooPrpVP9wP`，Web 搜索核实)：评估基因多数表决与 kNN 在多细胞系未过滤数据上的分类表现，未涉及错配细胞系基线对比。
   - **Molina & Zhang (2026)** (`10.64898/2026.07.24.740459`)、**HyperMap** (`10.64898/2026.04.23.720505`)、**PIE** (`10.64898/2026.10.02.756297`)、**X-Cell** (`10.64898/2026.03.18.712807`)、**Cupaioli & Merelli (2026)** (`10.64898/2026.11.03.763968`, `fmaudit`)、**Miller et al. (2026)**（*Nat. Biotechnol.* 校准指标下模型 vs 同系均值配对 Wilcoxon/t 检验）、**SBB / CellSimBench**（模型 vs 同系 `control_mean` 配对 t 检验）、**OCCAM** (`arXiv:2611.21959`)、**Can We Trust In-Distribution Success?** (`arXiv:2608.00152v3`)、**The Delta Gap** (`arXiv:2611.10658v1`)、**cell-eval2**、**scPertEval** 与 **PertEval**：已核查代码与全文文本，均未包含“目标系自身基线 vs 错配细胞系基线”的每靶点配对比较。

---

## 四、对 Round 12 检验设计定位的文献启示

1. **数据集新颖性（问题 1：`没有找到`）**：使用 `Liu et al. (2026)` 的 `GSE345058` A549 CROP-seq 单细胞扰动转录组数据作为跨细胞系扰动预测（及其反事实基线特异度验证）的独立外部确认系，在现有已发表及预印本文献中无先例可查。
2. **检验设计边界与声明措辞（问题 2：`部分做过`）**：写作时不可泛称“此前从未做过错配细胞系对照检验”，因为 `CanPert/CASTER`（`10.64898/2026.11.15.766815`）已对神经网络模型做过“正确目标对照输入 vs 1 个错配他系对照输入”的逐靶点配对 Wilcoxon 检验（并明确发现替换对照会显著降分）；同时可准确论证 Round 12 `FCT − Mis_s` 与之的本质差别在于：**直接比较反事实基线表达锚点本身（保留 vs 错配基线 `c_{b,g}`，而不强制加回目标系对照均值）、在目标—错配基线之间做配对比较、覆盖全部源细胞系并采用 Intersection-Union Test，且在正交校准指标 `r_orth_cal` 上检验基线信号是否真正来自目标系特异性表达基线而非全局/共用趋势**。
