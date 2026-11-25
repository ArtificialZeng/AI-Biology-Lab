# Round 16 先验文献与官方实现核查报告（`PRIOR_ART_I`）

- **智能体**：`A-T5`（A队 · 文献与官方指标）
- **服务协议**：[`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) / [`PREREG_I_AUDIT_R1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PREREG_I_AUDIT_R1.md)
- **核查红线遵守声明**：全程仅使用公开 API / 公开网页 / 开源源码与本地文献缓存只读检索；未下载任何模型权重，未下载 `>100 MB` 文件，未安装新包，未执行 `git`，未绕过反爬/Cloudflare，**未打开任何本机已保存的 `ST-*-Replogle/.../adata_pred.h5ad` 或 `adata_real.h5ad`**。

---

## 核心结论速览（Q1–Q5）

| 问题 | 核心结论 | 对 Round 16 的关键影响 |
|---|---|---|
| **Q1（Control-swap / 消融先例）** | **未发现**在 `Arc STATE` 上固定权重与扰动输入、仅跨细胞系替换目标系输入 `ctrl_cell_emb`（matched vs. mismatched controls）的对照置换实验；最接近的机制审计是 `Zhang et al. 2026`（`arXiv:2611.05849`）对解码器内部 $L_0$ SVD 方向跨细胞系解码器的置换及 `control+batch` / `pert+batch` ridge 探针、`CellCanon`（E.2 移除自身 DMSO 背景）与 `Rezaei Jafari et al. 2026`（对合成 `scGPT+STATE` 消融扰动通道）。 | 支持 Round 16 将「输入目标系控制池跨系置换（control-swap / scramble）」作为 `Arc STATE` 上未见先例的机制检验；撰写讨论时必须明确与 `arXiv:2611.05849`（$L_0$ 内部激活置换 + ridge 探针）及 `OCOO-T`（`S1` vs. `S2–S7`）的边界区分。 |
| **Q2（共享轴移除与 `r_orth_cal`）** | `Systema`（`10.1038/s41587-025-02777-8`）仅将参考原点从对照质心平移为所有训练扰动质心均值 $O_{\mathrm{pert}}$（单向量平移），**并非**向 $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}^\perp$ 的正交补投影；近邻方法 `AmbiModBench` 做分数减法对照、`PerturbBridge` 在自编码器结构中保留响应正交补、`Bison` 做药物对比评分，均**未见**与 `r_orth_cal` / `r_dorth_cal` 相同的评估度量。 | 支持 `r_orth_cal` 作为主指标的先例新颖性；须与 `Systema`（质心平移）及 `AmbiModBench`（配对均值分数差）明确对照区分。 |
| **Q3（Split-half 后选基因偏倚与 `w(μ)=μ/(1+μ)`）** | `Nicol et al. 2026 v2`（`10.64898/2026.05.07.723486`）解决的是**共享对照引起的伪相关**（通过控制组拆半）；`Miller et al. 2026` 将细胞拆分为 `GT` 与 `TD` 但评价用 DEG 仍在 `GT` 本身计算；`Occam` 与 `cell-eval2` 用分半作为噪声上限。**未发现**针对「split-half 选基因之后残存的真实幅度赢家诅咒（winner's curse）」或扰动评估权重 `w(μ)=μ/(1+μ)` 的推导或应用（仅 `Booeshaghi & Pachter 2021` 与 `Ahlmann-Eltze & Huber 2023` 分析过 `log(1+x)` 的方差与 Taylor 小均值线性渐近性）。 | 支持 Round 16 提出残余真幅度选择偏倚诊断及 `w(μ)=μ/(1+μ)` 加权的先例空缺判定；须注意与 `Booeshaghi & Pachter 2021` 的 Taylor 展开工作区分范围。 |
| **Q4（`SE-600M` 与 `ST-*-Replogle` 暴露与协变量）** | **重要风险预警**：不仅 `Replogle K562/RPE1`（`PRJNA831566`, `2022-04-26`；`SE-600M/config.yaml` 显式列出 `rpe1_top5000_variable.h5ad`）与 `Nadig HepG2/Jurkat`（`GSE264667`, `2024-05-05`），**`Feng iPSC`（`PRJEB81502`, ENA `first_public = 2024-11-08`）同样早于 `scBaseCount` v1（`2025-02-27`/`2025-03-04`）与 `SE-600M`（`2025-06-23`）**；而 `CD4 T`（`GSE314342` `2026-02-06`）、`DLD1`（`GSE337988` `2026-07-10`）与 `A549`（`GSE345058` `2026-08-26`）严格晚于 `SE-600M` 训练时间。代码层面 `cell_type_key`（`cell_line`）仅在 `_infer.py` 用于分组抽取对照池，**不进入 `StateTransitionPerturbationModel.forward`**；所有 8 个 `ST-*-Replogle` 模型均配置 `batch_encoder: true`、`use_batch_token: false`、`batch_col: gem_group`。 | **必须修正对 Group 2 时间线先验的判断**：在私有训练集 CSV（`esm_basecount_tahoe_cellxgene_train_filtered.csv`）`未能读取` 的前提下，**`Feng iPSC` 不能按“绝对晚于 `SE-600M` 预训练截止”处理**，应按 ENA `2024-11-08` 标为“早于 `SE-600M` 潜在暴露”；只有 `CD4 T`、`DLD1`、`A549`（及使用 10x Flex 的 `VCC 2025 H1`）具有明确时间或技术栈隔离优势。 |
| **Q5（`Illumina/CellularAI` 仅用公开数据运行情况）** | **未发现**任何人仅用公开数据端到端复现或运行 `Illumina/CellularAI`（`Apache-2.0`）。仓库明确声明不包含从原始训练数据生成训练汇总包（summary bundle）的脚本，完整运行依赖受访问控制的 `Illumina Billion Cell Atlas` 衍生汇总包；GitHub 公开库仅内置 30 基因合成烟雾测试与 `vcc2025` 评分框架，截止核查时仅 1 条祝贺 Issue、0 个 PR，外部无复现记录。 | 支持不对 `CellularAI` 做跨系训练复现假设；与 VCC 2026 第 2 名 `flowstate`（纯公开 222 套数据独立实现，但未见独立外部仓库）形成对照。 |

---

## 受阻与未能读取的一手文献/清单透明披露

依据工作规程，凡遇到 Cloudflare / HTTP 403 / HTTP 429 或访问控制壁垒的文献正文与内部清单，本报告均明确标注 **`未能读取`**，不作记忆补全：

1. **`Arc STATE` 论文全文**：`bioRxiv` 全文（`https://www.biorxiv.org/content/10.1101/2025.06.26.661135v2.full`）返回 `HTTP 403 Forbidden`（Cloudflare），Elsevier API（`https://api.elsevier.com/content/article/doi/10.1016/j.cell.2026.07.052`）返回 `HTTP 429 Too Many Requests`——**`未能读取`正文**；已核实 PubMed/bioRxiv 官方摘要、Arc 官方发布博客、`ArcInstitute/state` 开源代码及 HuggingFace 配置文件。
2. **`Fan et al. 2026`（OpenReview）**：OpenReview API（`https://api2.openreview.net/notes?id=lnpzb1HvwX`）返回 `HTTP 403 Forbidden`——**`未能读取`一手 PDF/API 正文**；仅见项目内部二级审查摘要 [`PRIOR_ART_H_CONTROLLER_FAN2026.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/status/PRIOR_ART_H_CONTROLLER_FAN2026.md)（非一手文献）。
3. **`Vollenweider & Buehlmann 2026`（`SBB`）全文**：`bioRxiv` 全文页面（`https://www.biorxiv.org/content/10.64898/2026.04.20.719650v1.full`）返回 `HTTP 403 Forbidden`——**`未能读取`全文**；已通过 bioRxiv API 核实标题、DOI 与完整摘要。
4. **`CellFlow`（`Klein et al. 2025`）全文**：`10.1101/2025.04.11.648220` 仅通过 bioRxiv API 核实摘要与元数据，**`未能读取`全文**。
5. **`SE-600M` 内部训练集清单**：`config.yaml` 所指 `/large_storage/ctc/userspace/aadduri/data/auxillary/esm_basecount_tahoe_cellxgene_train_filtered.csv` 及 `arc-ctc-scbasecamp` 私有接口未公开（GCS 匿名访问受限）——**`未能读取`单样本 CSV 列表**。
6. **`Illumina/CellularAI` 专有汇总包与 GitHub REST API**：专有训练汇总包受访问控制——**`未能读取`**；GitHub REST API 未认证请求返回 `HTTP 403 rate limit exceeded`，仓库 Issue/PR 计数通过 GitHub 公开 HTML（`HTTP 200`）核实。

---

## Q1：Control-swap / Scramble / Control-ablation 先例及近 24 个月扰动模型与基线文献核查

### 1.1 直接针对「控制细胞跨系置换 / 打乱 / 移除」的核查结果

**结论**：在已核查的论文与开源基准中，**未发现**有人在训练好的 `Arc STATE`（或其他虚拟细胞模型）上固定模型权重与扰动输入、专门把目标细胞系的输入对照池 `ctrl_cell_emb` 替换为另一细胞系或打乱来源对照池（matched vs. mismatched control-swap）。最接近的工作分别在**解码器内部 $L_0$ 激活方向**、**辅助探针输入**、**自身对照编码方式**或**独立无对照基线模型**上做过相关消融：

1. **`Arc STATE`（Adduri et al., *bioRxiv* 2025 / *Cell* 2026）** — 正文 **`未能读取`**（`HTTP 403`/`429`）；已核查摘要、博客与源码
   - 来源：<https://doi.org/10.1016/j.cell.2026.07.052>（PubMed PMID `42673963`；预印本 <https://doi.org/10.1101/2025.06.26.661135>；博客 <https://arcinstitute.org/news/blog/state>；仓库 <https://github.com/ArcInstitute/state>）
   - 原文引用（*Cell* 2026 摘要）：
     > "Its cell embeddings trained on observational data from 167 million cells enable the identification of strong perturbations in cellular contexts where no perturbations were observed during training."
   - 源码核查边界：`state/tx/models/decoder_only.py` 直接以目标扰动潜在表示 `batch["pert_cell_emb"]` 输入基因解码器（oracle 潜在重构基线，而非去除控制输入的扰动推理）；`state_transition.py` 与 `_infer.py` 中未含跨细胞系对照置换检验。

2. **`Zhang et al. 2026`（NeurIPS 2026 录用预印本，ST 低秩解码器审计）** — 最接近的机制探针与内部激活置换
   - 来源：<https://arxiv.org/abs/2611.05849>
   - 原文引用（Table 7 表题）：
     > "Change in mean PDC when the top L0 direction of each perturbation is replaced by the one the same encoder gives it in a decoder trained on another cell line."
   - 关键区分：原文指出解码器“the decoder receives 64 control cells from the same batch (gem group) as context for each perturbation”，在 Replogle 四系上审计 `SE-600M+ST` / `scGPT+ST` / `Tahoe-x1+ST`，在 **C.5 做了线性 ridge 探针输入对比**（“A ridge on the control and batch inputs alone reaches PDC 0.27–0.46”，也拟合了 `perturbation+batch` 单独输入和打乱扰动标签的探针），并在 **Table 7 将解码器第 0 层注意力 $L_0$ 的主奇异方向跨不同细胞系训练出的解码器互换**；但**并未把同一个训练好的 `ST` 解码器的输入 `ctrl_cell_emb` 替换成另一细胞系对照细胞来检验全模型预测稳定性**。

3. **`Mao et al. 2026`（`VCBench`）** — 含独立 `DecoderOnly` 无对照基线，不含 `STATE` 控制置换
   - 来源：<https://arxiv.org/abs/2604.27646>
   - 原文引用（Section 4.3）：
     > "DecoderOnly predicts perturbed expression solely from perturbation and covariate identities without conditioning on control cell expression, serving as a lower bound that simulates mode collapse."
   - 关键区分：在统一基准中纳入 `STATE` 与独立基线 `DecoderOnly`（仅吃扰动和协变量标签而不读对照表达），**未**在 `STATE` 上执行控制池置换或控制端消融。

4. **`Rezaei Jafari et al. 2026`（归因诊断 `scGPT+STATE`）** — 消融扰动输入，未置换控制细胞系
   - 来源：<https://arxiv.org/abs/2611.02335>
   - 原文引用（Section H.2）：
     > "We either mask a gene embedding or ablate the perturbation embedding, providing black-box measures of how strongly the model depends on each input."
   - 关键区分：在 128 基因线性 SDE 合成数据上对关节微调的 `scGPT+STATE` 消融**扰动嵌入通道**（导致 Pearson-$\Delta$ 从 `0.506` 降至 `0.035`）或遮蔽输入基因，未做真实数据跨系对照置换。

5. **`OCOO-T`（Jiang et al. 2026）** — 目标系均值对照表征 vs. 分类细胞系嵌入
   - 来源：<https://arxiv.org/abs/2606.12838>
   - 原文引用（Section 4.4）：
     > "Overall, these results suggest a trade-off between applicability and peak benchmark performance: mean control-cell profile conditioning offers a practical route for representing unseen or fine-grained cellular contexts, whereas learned cell-line embeddings can provide stronger perturbation-specific discrimination and differential-expression recovery in seen contexts."
   - 关键区分：其 `S1`（显式细胞系嵌入）vs. `S2–S7`（目标系自身 1–64 个对照均值）比较的是条件编码方式，而非错配细胞系对照置换。

6. **`CellCanon`（Liang et al. 2026）** — 仅在自有架构内做 DMSO 背景消融
   - 来源：<https://arxiv.org/abs/2611.05890>
   - 原文引用（Appendix E.2）：
     > "All-gene Pearson decreases by 0.0265 and Spearman by 0.0091; the top-100 decreases are 0.0241 and 0.0293."
   - 关键区分：仅在 CellCanon 自身保留响应图谱的前提下测试移除自身实测 DMSO 背景特征（SciPlex3 药物响应），未测试 `STATE` 或跨细胞系错配对照替换。

7. **`CASTER`（Kawsar et al. 2026）** — 跨四系自适应转移，未做对照替换
   - 来源：<https://doi.org/10.64898/2026.11.15.766815>（代码：<https://github.com/CanPert> / 本机 `CanPert` 缓存已核查）
   - 原文引用（bioRxiv 摘要）：
     > "CASTER led only at 2.5% screened (0.62 versus 0.47) and fell below the transferred mean without a screen."

8. **`Fan et al. 2026`（OpenReview `lnpzb1HvwX`）** — **`未能读取`一手正文**（`HTTP 403`）
   - 来源：<https://openreview.net/forum?id=lnpzb1HvwX>
   - 核查说明：OpenReview API/全文受 `HTTP 403` 阻断；参考项目内部二级记录 [`PRIOR_ART_H_CONTROLLER_FAN2026.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round15/status/PRIOR_ART_H_CONTROLLER_FAN2026.md)，其涉及在目标系自身对照上的平移对照基线与 `STATE`/`CellFlow` 对比，无跨系错配对照置换；因一手页面未读取，此处不引一手原文。

### 1.2 近 24 个月（2024.10–2026.11）核心基准、竞赛总结与新模型核查表

以下逐篇确认其**方法/发现**、**一手来源**、**一句原文引用**，并确认**均未包含**对 `Arc STATE` 或同类模型做跨系目标对照池置换（control-swap）的消融实验：

- **`VCC 2025 Wrap-Up`（Arc Institute, 2025-12）**
  - 来源：<https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up>
  - 原文引用：
    > "The winning approaches combined deep learning with classical statistical features, suggesting pure end-to-end learning is yet to solve this problem."

- **`VCC 2026 Wrap-Up`（Arc Institute, 2026-11-19；含 `CellularAI`、`flowstate`、`X-Cell Transfer`）**
  - 来源：<https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up>
  - 原文引用：
    > "Together, these observations indicate that current models recover roughly what a cross-line average would predict, whereas the part of the response that is specific to a cell line, which is precisely what a virtual cell model would need in order to improve on an experiment in another model system, remains largely missing."

- **`Dobin et al. 2026`（`cell-eval2` / `vcc2026` 指标体系）**
  - 来源：<https://doi.org/10.64898/2026.11.17.767214>（规范文档：`vcc2026-metrics-brief.md`）
  - 原文引用（`vcc2026-metrics-brief.md`）：
    > "The reference's control cells are the origin on both sides: the predicted effect is measured against them, and they are the reference group of the differential-expression test for the prediction as well as for the reference."

- **`Luan et al. 2026`（VCC 2025 回顾性 KOLF2.1J 跨系统计迁移）**
  - 来源：<https://doi.org/10.64898/2026.10.20.760514>
  - 原文引用（bioRxiv 摘要）：
    > "The method centers an independent KOLF2.1J source on known training targets, transfers its response signs and magnitude ranks, and uses isotonic regression to render predicted expression means."

- **`Denvis 2026`（虚拟细胞经验对比与基线差距）**
  - 来源：<https://doi.org/10.21203/rs.3.rs-10434123/v1>
  - 原文引用（摘要）：
    > "First, on unseen single-gene perturbations no method dramatically beats a well-tuned additive/linear baseline: the best model (State) improves over a cell-mean baseline by 26% versus 19% for the linear baseline and 22% for GEARS, and several foundation-model predictors fall below the linear baseline—echoing recent critical benchmarks."

- **`Ahlmann-Eltze et al. 2025`（*Nature Methods*）**
  - 来源：<https://doi.org/10.1038/s41592-025-02772-6>
  - 原文引用（结论段）：
    > "The analysis presented in this Brief Communication shows that, despite their complexity, current deep learning models do not outperform linear baselines in predicting gene perturbation effects, thus emphasizing the importance of further method development and thorough evaluation."

- **`Csendes et al. 2025`（*BMC Genomics*）**
  - 来源：<https://doi.org/10.1186/s12864-025-11600-2>
  - 原文引用（摘要）：
    > "Surprisingly, we found that even the simplest baseline model—taking the mean of training examples—outperformed scGPT and scFoundation."

- **`PertEval-scFM`（Wenteler et al. 2024）**
  - 来源：<https://doi.org/10.1101/2024.10.02.616248>
  - 原文引用（摘要）：
    > "Our results show that scFM embeddings do not provide consistent improvements over baseline models, especially under distribution shift."

- **`Miller / Mejia et al. 2026`（*Nature Biotechnology*）**
  - 来源：<https://doi.org/10.1038/s41587-026-03307-w>
  - 原文引用（DRF 标定定义）：
    > "In a miscalibrated scenario, DRF values approach zero, reflecting limited metric sensitivity to detect improvements over the negative control."

- **`Li et al. 2026`（*Science Advances* 25 数据集系统基准）**
  - 来源：<https://doi.org/10.1126/sciadv.aed3414>
  - 原文引用（摘要）：
    > "We present a large-scale benchmark of 13 representative methods and baselines across 25 datasets spanning diverse perturbation modalities and species, including two primary immune-cell drug-response resources."

- **`Vollenweider & Buehlmann 2026`（`SBB`）** — 正文 **`未能读取`**（`HTTP 403`）；已读 bioRxiv API 摘要
  - 来源：<https://doi.org/10.64898/2026.04.20.719650>
  - 原文引用（bioRxiv 摘要）：
    > "Applying these principles across seven transcriptomic perturbation datasets, including single and double perturbations, we demonstrate that complex deep learning methods, including foundation models, often still fail to meaningfully surpass simple linear baselines, and that substantial room for improvement remains even where they do."

- **`Shoeibi & Yousefi 2026`（锁盲迁移与测序深度缠结）**
  - 来源：<https://arxiv.org/abs/2608.00152>
  - 原文引用（Section 2）：
    > "A recent line of work documents that deep and foundation-model perturbation predictors frequently fail to outperform simple baselines."

- **`Panagopoulos et al. 2026`（细胞行重复导致的 Wilcoxon 显著性膨胀陷阱）**
  - 来源：<https://arxiv.org/abs/2610.26267>
  - 原文引用（摘要）：
    > "Taking this pitfall into account in our benchmark comparison, we find that a random forest using Gene Ontology features achieves the best overall results."

- **`AmbiModBench`（Huang et al. 2026）**
  - 来源：<https://arxiv.org/abs/2609.32527>
  - 原文引用（摘要）：
    > "First, absolute metrics cannot separate target-specific predictions from a shared background response."

- **`Bison`（Liu et al. 2026）**
  - 来源：<https://arxiv.org/abs/2609.32467>
  - 原文引用（摘要）：
    > "This comparison reveals that high overall response agreement can coexist with weak prediction of drug-specific differences, despite reproducible signals across repeated measurements."

- **`S1Pert`（Sheng & Wang 2026）** — 做了目标排列诊断（Table IV），未做控制细胞置换
  - 来源：<https://arxiv.org/abs/2611.17705>
  - 原文引用（Section 1）：
    > "Across four cellular benchmarks, S1Pert consistently outperforms state-of-the-art baselines in directional decision accuracy, within-gene target discrimination, and online inference efficiency."

- **`PerturbJEPA`（Pang et al. 2026）** — 验证对比坍缩（对照-处理差缩小），未做 `STATE` 控制置换
  - 来源：<https://arxiv.org/abs/2611.20575>
  - 原文引用（摘要）：
    > "Controlled diagnostics show that treatment–control differences can shrink even when representations remain globally dispersed."

- **`X-Cell`（Wang et al. 2026）** — 依托自身多模态先验与测试时自监督适应，未做 `STATE` 对照置换
  - 来源：<https://doi.org/10.64898/2026.03.18.712807>
  - 原文引用（bioRxiv 摘要）：
    > "X-Cell-Ultra demonstrates zero-shot generalization to novel biological contexts, including unseen iPSC-derived melanocyte progenitors and primary human CD4+ T cells from multiple donors, and outperforms all baselines after self-supervised test-time adaptation."

- **`CellFlow`（Klein et al. 2025）** — 全文未读，仅核实 bioRxiv API 摘要
  - 来源：<https://doi.org/10.1101/2025.04.11.648220>
  - 原文引用（bioRxiv 摘要）：
    > "Here, we present CellFlow, a flexible framework based on flow matching that can model single cell phenotypes induced by complex perturbations."

---

## Q2：共享响应轴移除（`Systema` 与 $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}^\perp$ 正交投影）

1. **`Systema` 的精确数学定义（非正交补投影，而是单质心均值平移）**
   - 来源：Viñas Torné et al., *Nature Biotechnology* (2025), <https://doi.org/10.1038/s41587-025-02777-8>（PMC13271886）
   - 原文引用（Online Methods / Systema reference）：
     > "In Systema, we consider the same set of evaluation metrics using the average O pert of all perturbation-specific centroids O(X) as reference, which effectively emphasizes perturbation-specific effects."
   - 机制对比：`Systema` 将原本以未扰动对照质心为原点的偏移向量 $x_p - O_{\mathrm{ctrl}}$ 替换为以训练扰动平均质心为原点的向量 $x_p - O_{\mathrm{pert}}$（即减去单个均值向量）；它**不会**移除任意预测向量在 $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}$（常数全局偏移 + 源上下文均方响应方向 + 测试上下文平均扰动基线）上的各向分量缩放（例如把源上下文效应直接按标量放大/缩小仍可能残留在非正交空间内）。

2. **最接近的评估/建模近邻对比**
   - **`AmbiModBench`（Huang et al. 2026）**：<https://arxiv.org/abs/2609.32527>（在指标得分层面减去训练集均值基线得分 $d_t(S)=S(\hat{y}_t, y_t)-S(c, y_t)$，属于指标标量减法对比，不是向量空间正交补投影）。
     > "In this work, we address the three deficiencies in the same order and (1) replace the absolute score with a paired training-mean contrast, in which the reference is fitted on training targets and scored against the same held-out response, reported with target-cluster uncertainty; (2) calibrate readouts by gene-coordinate permutation and independent cell replicates, separating gene assignment from amplitude; and (3) formulate training coverage through a test-to-training embedding radius and an ambiguity modulus, with a conditional bound and a fixed-test subsampling experiment."
   - **`PerturbBridge`（Liu et al. 2026）**：<https://doi.org/10.64898/2026.09.24.754259>（在其响应锚定自编码器 `RA-AE` 内部用奇异响应轴的正交补来编码互补细胞状态变化，属模型表示学习而非评估指标）。
     > "A response-anchored autoencoder constructs a compact cellular representation that preserves dominant transcriptional response directions identified from training perturbations while encoding complementary cellular variation."
   - **`Bison`（Liu et al. 2026）**：<https://arxiv.org/abs/2609.32467>（单独评估 `Drug contrast` 与 `pooled vs. plate-matched controls` 分数，见 Appendix C.5）。
     > "Thus, the overall ranking persists when the scoring reference changes, while drug contrast separately evaluates compound-dependent differences."

3. **Q2 结论**：**未发现**在单细胞扰动评测中显式将测试/预测表达增量投影到 $\operatorname{span}\{\mathbf{1}, \bar{\Delta}^{\mathrm{src}}_p, \bar{b}_S\}^\perp$ 的先例指标（即 [`PREREG_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/prereg/PREREG_I.md) 所定义的 `r_orth_cal` 与 `r_dorth_cal`）。

---

## Q3：Split-half 后选基因偏倚与 `log1p` 饱和权重 $w(\mu)=\mu/(1+\mu)$

1. **`Nicol et al. 2026 v2`（`10.64898/2026.05.07.723486`, `PMC13192888`）的精准边界**
   - 来源：<https://doi.org/10.64898/2026.05.07.723486>
   - 原文引用（摘要）：
     > "Reanalysis of published datasets using a simple control-splitting procedure that removes this bias leads to a substantial reduction in performance previously attributed to biological signal."
   - 机制边界：该文针对**真实与预测均减去同一控制组估计均值 $\hat{\mu}_0$ 所引发的共享随机变量伪相关**，提出将控制组独立拆半以获得无共享噪声偏差的 $\hat{\rho}_\Delta^{\mathrm{split}}$；**并未**讨论「用独立 Half A 筛选 TOP DEGs 后在 Half B 上计算权重/指标时，高真实幅度基因（或与筛选正相关的系统信号）仍因选择概率单调递增而产生的残余真幅度选择偏倚」。

2. **分半评估与 DEG 选择的其他现有文献**
   - **`Miller / Mejia et al. 2026`（*Nature Biotechnology*）**：<https://doi.org/10.1038/s41587-026-03307-w> —— 将细胞单次随机拆为 `GT`（真值）与 `TD`（预测器训练）两半，**但在同一份子采样的 `GT` 上原地完成标准化、HVG 与 DEG 筛选及打分**，未隔离选基因与评估噪声。
     > "For each target cell count n ∈ {8, 16, 32, 64, 128, 256}, we sampled n /2 cells from each half and then independently repeated normalization, HVG selection and DEG computation on the subsampled data."
   - **`Occam`（Bendidi et al. 2026）**：<https://arxiv.org/abs/2611.21959> —— 将细胞与对照按批次分层拆分为 `A` 与 `B` 两半，以 `A` 作参考、“比较 `B` 与 `A` 测实验可重复性天花板”。
     > "Within each test perturbation, we divide measured cells into disjoint, equally sized halves A and B, stratified by batch; controls are split separately."
   - **`cell-eval2`（Dobin et al. 2026）**：<https://doi.org/10.64898/2026.11.17.767214>（`vcc2026-metrics-brief.md`）—— 同样将每组与对照细胞对半分 5 次（`n_splits: 5`），仅用于定义实验重复天花板 $r$ 以标定 $(u-b)/(r-b)$。

3. **`log1p` 饱和与权重 $w(\mu)=\mu/(1+\mu)$ 文献**
   - **经典预处理/方差理论**：`Booeshaghi & Pachter 2021`（*Bioinformatics*, <https://doi.org/10.1093/bioinformatics/btab085>）和 `Ahlmann-Eltze & Huber 2023`（*Nature Methods*, <https://doi.org/10.1038/s41592-023-01814-1>）系统分析了 `log(1 + X)` 对小均值基因方差/期望的截断与 Taylor 渐近行为，其中 `Booeshaghi & Pachter 2021` 明确指出对 `log(1 + \lambda)` 在 $\lambda$ 小时斜率趋近 `1`、大时趋近 `log(\lambda)`；但二者讨论的是**数据预处理/变换偏差**，而非在扰动预测误差/响应排序度量中对 `log1p` 差值显式乘以饱和权因子 $w(\mu)=\mu/(1+\mu)$。
     > "While this is approximately λ when λ is large, it is close to 1 when λ is small ( de L’Hospital, 1768 )." （Booeshaghi & Pachter 2021）
     > "We compare their strengths and weaknesses and find that the latter three have appealing theoretical properties; however, in benchmarks using simulated and real-world data, it turns out that a rather simple approach, namely, the logarithm with a pseudo-count followed by principal-component analysis, performs as well or better than the more sophisticated alternatives." （Ahlmann-Eltze & Huber 2023）
   - **术语消歧（`Bison`, `arXiv:2609.32467`）**：仅在量化 FSQ 自编码器预训练损失中出现 `"coordinate-saturation penalty"`（用于防止离散编码坐标饱和），与表达均值权 $w(\mu)=\mu/(1+\mu)$ 无关。
   - **Q3 结论**：除 `Nicol et al. 2026 v2` 之外，**未发现**其他工作专门讨论 split-half 选基因后的残余真实幅度偏倚或 $w(\mu)=\mu/(1+\mu)$ 权重函数。

---

## Q4：`SE-600M`、`ST-*-Replogle` 架构与 7 套目标数据集暴露审计

### 4.1 模型配方、公开配置文件与源码级协变量机制核查

1. **`arcinstitute/SE-600M` 配置与数据规模口径说明**
   - 来源：<https://huggingface.co/arcinstitute/SE-600M/blob/main/config.yaml>（仓库创建 `2025-06-23`；`se600m_epoch15/16.ckpt` 上传于 `2025-07-09`，`model.safetensors` 重命名于 `2025-07-15`，2026-02-24 仅重上传 safetensors）；训练入口源码：`state/emb/train/trainer.py`。
   - 原文引用（*Cell* 2026 摘要 & `SE-600M/config.yaml`）：
     > "Its cell embeddings trained on observational data from 167 million cells enable the identification of strong perturbations in cellular contexts where no perturbations were observed during training." （*Cell* 2026 摘要，DOI `10.1016/j.cell.2026.07.052`）
   - `config.yaml` 与源码实证（**解释关键数值差异**）：
     - `dataset.current: scbasecamp-cellxgene-tahoe-filtered`，标注 `num_datasets: 14420`，训练/验证 CSV 指向 `/large_storage/ctc/userspace/aadduri/data/auxillary/esm_basecount_tahoe_cellxgene_train_filtered.csv`（内部 CSV **`未能读取`**）；
     - `dataset.num_cells: 36238464` 在 `state/emb/train/trainer.py` 第 41–42 行仅用作 `TOTAL_N_CELL = cfg.dataset.num_cells` 来计算每个 epoch 的步数（`warmup_steps` 调度），而 *Cell* 2026 摘要/Arc 博客声明 `SE-600M` 实际在约 `167 million`（博客称 `"nearly 170 million"`）观测单细胞上训练；
     - **`config.yaml` 的 `dataset.overrides` 与 `validations.diff_exp` 显式写入** `/large_storage/ctc/datasets/vci/validation/rpe1_top5000_variable.h5ad` 及 `/large_storage/ctc/datasets/vci/validation/replogle_perturbation.h5ad`（`validations.*.enable: false`）。
   - `scBaseCount` 来源与口径：Youngblut et al., *bioRxiv* `2025-02-27`/`2025-03-04` (<https://doi.org/10.1101/2025.02.27.640494>)；*Cell* 2026 (<https://doi.org/10.1016/j.cell.2026.08.025>)。
     > "Built by directly mining all 10x Genomics datasets from SRA, scBaseCount is the largest freely accessible public repository of single-cell gene expression data, comprising over 502 million cells across 27 organisms and 75 tissues, offering an unbiased view of the composition of data within SRA."

2. **`ST-SE-Replogle` 与 `ST-HVG-Replogle` 源码与 8 个 `zeroshot/*` 配置核查**
   - 来源：<https://huggingface.co/arcinstitute/ST-SE-Replogle> / <https://huggingface.co/arcinstitute/ST-HVG-Replogle> 与源码 `state/_cli/_tx/_infer.py`、`state/tx/models/state_transition.py`。
   - 原文引用（HF 仓库模型卡 `README.md`）：
     > "This repository contains 4 fewshot and 4 zeroshot canonical models, with their saved configurations, loading files and available CellEval v1 outputs."
   - 全部 8 个 canonical `zeroshot/*` 配置（`hepg2/jurkat/k562/rpe1` × `X_state/X_hvg`）直接核实字段完全一致：`basal_mapping_strategy: random`、`batch_col: gem_group`、`pert_col: gene`、`cell_type_key: cell_line`、`cell_set_len: 32`、`predict_residual: true`、`batch_encoder: true`、`use_batch_token: false`。
   - 推理源码机制（`_infer.py` 第 663–697, 739–754, 822–825 行）：
     - **`cell_type_key`（`cell_line`）仅用于 `_infer.py` 分组采样对照池**，**从不进入** `StateTransitionPerturbationModel.forward` 网络权重或嵌入层；
     - `batch_encoder: true` 且 `use_batch_token: false` 时，若输入 `.h5ad` 存在 `batch_col`（`gem_group`），未在训练 `batch_onehot_map` 中见过的字符串**回退到索引 `0`** 并加上该索引对应的批嵌入向量；若输入 `.h5ad` 不存在该 `batch_col`，则置 `uses_batch_encoder = False` 并跳过 `bi = None`。

### 4.2 Round 16 目标 7 套数据集与 `SE-600M` / `ST-*-Replogle` 时间线与重合审计

| 目标数据集 / 上下文 | 一手库与文献来源（DOI / Accession） | 公开释放日期（GEO / ENA） | 对 `SE-600M`（`2025-06` 预训练）暴露状态 | 对 `ST-*-Replogle` 监督微调暴露状态 | 原文一句话核证引用 |
|---|---|---|---|---|---|
| **Replogle K562 / RPE1** | <https://doi.org/10.1016/j.cell.2022.05.013>；ENA [`PRJNA831566`](https://www.ebi.ac.uk/ena/browser/view/PRJNA831566) | ENA `first_public`: **`2022-04-26`**（Cell 2022-07） | **明确早于 `SE-600M`**；且 `SE-600M/config.yaml` 显式出现 `rpe1_top5000_variable.h5ad` | **库内训练集之一**（按各自 `*_zeroshot.toml` 留出测试上下文） | "We use transcriptional phenotypes to predict the function of poorly characterized genes, uncovering new regulators of ribosome biogenesis (including CCDC86, ZNF236, and SPATA5L1), transcription (C7orf26), and mitochondrial respiration (TMEM242)." |
| **Nadig HepG2 / Jurkat** | <https://doi.org/10.1038/s41588-025-02169-3>；GEO [`GSE264667`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667) / ENA `PRJNA1100571` | GEO `Public on`: **`May 05 2024`**；ENA `first_public`: `2024-06-23` | **严格早于 `scBaseCount` v1（`2025-02-27`）与 `SE-600M`（`2025-06-23`）**；**潜在暴露**（私有 CSV `未能读取`） | **库内训练/评测集之一**（按 `*_zeroshot.toml` 留出对应测试系） | "This repository contains data from two large-scale pooled CRISPR screens with single-cell RNA sequencing readout (ie Perturb-seq) experiments in Jurkat and HepG2 cells." |
| **Feng iPSC**（多供体 iPSC） | <https://doi.org/10.1016/j.xgen.2025.101076>（PMC12903452）；ENA [`PRJEB81502`](https://www.ebi.ac.uk/ena/browser/view/PRJEB81502) (`ERP165335`) | **ENA `first_public`: `2024-11-08`**（论文发表 2026-02） | **⚠ 风险警报：原始测序早于 `scBaseCount` v1 与 `SE-600M` 达 ~3.6–7.5 个月！虽具体 run 是否入 `14,420` 过滤集因 CSV `未能读取` 尚不能坐实，但不能标为“时间绝对隔离”** | **未见暴露**（`ST-*-Replogle` 仅含 Replogle+Nadig 4 系） | "Here, we present a genome-scale CRISPR interference perturbation map with scRNA-seq readout across many genetic backgrounds in human pluripotent cells." |
| **VCC 2025 H1**（10x Flex） | <https://doi.org/10.1016/j.cell.2025.06.008>；博客：<https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge> | 挑战公布 `2025-06`（博客 `2025-08-22`） | 与 `STATE`/`SE-600M` 同期发布；**探针定量 10x Flex 需不同于 `scRecounter` 的处理管线**（且 Cell 2026 正文 **`未能读取`**） | **未在 `ST-*-Replogle` 监督训练** | "We chose 10x Genomics Flex chemistry because it consistently produced the highest-quality data in our hands, both in pilot comparisons and in the full-scale screen." |
| **Zhu CD4+ T**（`cd4t_rest` / `cd4t_stim48`） | <https://doi.org/10.1016/j.cell.2026.08.002>（PMC13614707）；GEO [`GSE314342`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE314342) / ENA `PRJNA1359008` (`SRP643211`) | GEO `Public on`: **`Feb 06 2026`**；ENA `first_public`: `2026-02-09` | **严格晚于 `SE-600M`（`2025-07` 模型权重定型）无预训练暴露** | **无曝光** | "Here, we perturbed all expressed genes across 22 million primary human CD4+ T cells from four donors and developed a probe-based perturb-seq platform to measure the transcriptome effects in cells at rest and after stimulation." |
| **DLD1 Perturb-seq** | GEO [`GSE337988`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE337988) / ENA `PRJNA1492064` | GEO `Public on`: **`Jul 10 2026`**；ENA `first_public`: `2026-07-12` | **严格晚于 `SE-600M`，无预训练暴露** | **无曝光** | "CRISPRi Perturb-seq was performed in DLD-1 cells expressing inducible Zim3-dCas9 at multiple multiplicities of infection (MOI)." |
| **A549 CROP-seq** | GEO [`GSE345058`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE345058) / ENA `PRJNA1483600` | GEO `Public on`: **`Aug 26 2026`**；ENA `first_public`: `2026-07-03` | **严格晚于 `SE-600M`，无预训练暴露** | **无曝光** | "Here, we present a multimodal perturbation atlas of 1,000 pooled CRISPR knockouts in A549 cells, profiled by fluorescence microscopy (42 live, 13 fixed markers), label-free quantitative phase imaging of the same live cells (at single timepoints), and single-cell RNA sequencing (scRNA-seq)." |

*(注：同期核查 `GSE295214` / `PRJNA1252779` 为 Pan et al. hESC 工作，ENA `first_public = 2026-06-19`，非 VCC 2025 H1 10x Flex 数据。)*

---

## Q5：开源 `Illumina/CellularAI` 仅用公开数据端到端运行先例核查

**结论**：**未发现任何外部团队或已发表文章仅使用公开数据端到端成功运行 `Illumina/CellularAI`**（GitHub：<https://github.com/Illumina/CellularAI>，`Apache-2.0`）。

1. **开源代码本身刻意切除了训练汇总生成链路（`docs/input-files.md` 与 `README.md`）**
   - 来源：<https://github.com/Illumina/CellularAI/blob/main/docs/input-files.md>（& <https://github.com/Illumina/CellularAI/blob/main/README.md>）
   - 原文引用（`docs/input-files.md` 第 5 行）：
     > "Raw training-data preprocessing and summary generation are not included."
   - 一手机制核实（`README.md` 第 9、58–76 行）：仓库明言其核心预测依赖基于 **`Illumina Billion Cell Atlas`**（受访问控制，**`未能读取`**）与 iPSC / CD4+ T 公开数据预计算好的 training summary bundle，要求向 `BusinessDevelopment@illumina.com` 申请获取该汇总包；无汇总包时仓库只可运行 30 基因虚构合成数据烟雾测试（`examples/synthetic/run_example.py`）和基于已生成预测矩阵的 `vcc2025` 评分脚本。

2. **外部社区复现、GitHub 检索与 VCC 2026 对比证据**
   - GitHub 全库搜索 `"CellularAI"` 仅返回 5 个相关名称仓库（除 `Illumina/CellularAI` 外均为无关项目/早期原型），且 GitHub 页面 HTML（`HTTP 200`；REST API 因未认证 `HTTP 403 rate limit exceeded` **`未能读取`**）显示截至核查仅 `1` 个非技术性祝贺 Issue、`0` 个 Pull Request；Google/arXiv/bioRxiv/EuropePMC 均无外部运行汇报。
   - **重要对照（VCC 2026 第 2 名 `flowstate`）**：根据 Arc `VCC 2026 Wrap-Up`（<https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up>），第 2 名 `flowstate` 使用独立自建管线在 222 套公开数据集上完成纯公开数据迁移（非运行 Illumina 开源 `CellularAI` 仓库）：
     > "flowstate was trained exclusively on public data, namely 222 perturbation-response datasets from 128 studies that cover 86 cellular contexts."
   - **Illumina 官方介绍博客**：<https://www.illumina.com/science/genomics-research/articles/introducing-the-cellularai-algorithm--predicting-single-cell-res.html>（与 VCC 2026 第一名说明一致，见 <https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up>）：
     > "The CellularAI algorithm is a statistical transfer framework that first predicts summary statistics for the new cell line and then builds cells that follow them."
