# PRIOR_ART_D: R10-E 文献三问核查报告（目标线基线消融、`log1p`/几何平均压缩增益、与“情境特异”可被无情境基线复现）

- **审计智能体**：`A-T5`（Team A · 文献与官方指标）
- **核查时间**：2026-11-22T13:08+08:00
- **任务依据**：[`R10E_探索性消融与合稿指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R10E_%E6%8E%A2%E7%B4%A2%E6%80%A7%E6%B6%88%E8%9E%8D%E4%B8%8E%E5%90%88%E7%A8%BF%E6%8C%87%E4%BB%A4.md)
- **前置关联文档**：[`PRIOR_ART_A.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_A.md) · [`PRIOR_ART_B.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_B.md) · [`PRIOR_ART_C.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round9/p2/vcc_linespecific/audit/LIT/PRIOR_ART_C.md)
- **合规边界**：严格仅做离线源码审计与 API / 公开文献核查（OpenAlex、bioRxiv API、EuropePMC、arXiv HTML/API、GitHub 开源仓库）；未接触任何未披露或受限表达数据，未做任何对外提交或账号注册。

---

## 执行摘要：三问三档结论总览表

| 问题编号 | 核心核查问题 | 三档结论 | 精准范围与关键边界说明 |
|---|---|---|---|
| **Q1** | `CellularAI`、`AtlasShift` 或其他 VCC 2025/2026 方案中，是否有人做过“去掉目标线基线”或“用常数基线替换”的消融？ | **「没有找到」** | **VCC 2025/2026 全部公开方案与基准中均未见此消融**：`CellularAI` 固定依赖目标对照基因比例 `self.fc * exp2(l2)`，其 `Baseline anchor` 仅比较均等平均扰动响应，`Fig. 7` 仅做 `5→50` 细胞系数量与组件累加；`AtlasShift` 在 `desired_mean()` 中强制使用目标 `control_probability`，仓库零消融脚本；`X-Cell`、`Arc State`、`PRiMeFlow`、`Luan et al. (2026)` 亦未做常数基线消融。**最近邻（非 VCC、边界不同）**：`CASTER`（[10.64898/2026.11.15.766815](https://doi.org/10.64898/2026.11.15.766815)）在非 VCC 4 线部分筛查任务 `T1` 中做过编码器输入端错线对照（`CASTER_wrongctrl`）与洗牌标签对照（`CASTER_label`），但评估仍把预测 `shift` 加回**真实目标对照细胞**（`ctrl + shift`），并未在 `log1p(c·FC)−log1p(c)` 前向映射内部替换目标基线 `c` 为常数。 |
| **Q2** | 是否有人指出：仅在 fold change 空间做平均（几何平均）或 `log1p` 非线性压缩本身就能在扰动预测基准上制造增益？ | **「没有找到」**<br>（紧邻数学/评测前作仅构成部分先验基础） | 文献中**无人**提出或证明：仅凭跨系 fold-change 几何平均与 `log1p(c·FC)−log1p(c)` 正向压缩即可在“常数基线 / 无目标基线信息”条件下制造扰动预测（尤其是去振幅后 `r_orth_cal`）的虚假增益（即主控 `E1` 发现的机制）。<br>作为紧密相邻但机理不同的前作，仅有四类部分结论：① `Booeshaghi & Pachter (2021)` 与 `Kviman et al. (2025/2026)` 指出对原始计数取 `mean(log1p(x))` 会在低表达基因处扭曲差异表达 LFC 估计并推高假阳性；② [`ArcInstitute/cell-eval2`](https://github.com/ArcInstitute/cell-eval2)（`#258/#260/#264`）证明单细胞 `lognorm` 的 Jensen 凹性缺失对零离散预测施加的是**负向惩罚/质量亏空**（42.8% 亏空、`expr_mse` 偏移 `+0.8..+1.4`），迫使 `X-Cell` 改用外部参照与重校准，并区分了 `mean_calc='geometric'` vs `'arithmetic'`；③ `Liu et al. (2025, arXiv:2511.16954)` 指出 $\ell_1/\ell_2$-PDS 对全局线性倍乘 `c` 及 `log1p` 预处理范数敏感；④ `Panagopoulos et al. (2026)` 指出重复单向量行构造细胞会虚高 MWU 显著性。以上均**未**触及 `log1p(c·FC)−log1p(c)` + 几何平均如何在常数基线下凭空复现 `86%–112%` 真实 FCT 得分。 |
| **Q3** | 是否有人已经证明：声称“细胞情境特异（cell-context-specific）”的扰动预测模型，其预测收益实际上可以被“不含情境信息”的基线复现？ | **「已有人做过」** | **已有三项直接实证研究**从不同维度证明多项“情境特异 / 跨细胞系零样本”模型收益主要来自情境共享均值而非情境特异信号：① **`CGC` (`Huang et al. 2026`, [10.64898/2026.09.30.755756](https://doi.org/10.64898/2026.09.30.755756))**：直接分解发现供体预测 `Ridge/Bilinear/MLP` 的 `shared_prediction_energy_fraction = 0.9718–0.999998`（留一供体特异 `R² ≤ 0`），且审计 `Arc State` 官方零样本预测发现相比无情境 `Pert Mean`，常规 `Delta-Pearson` 虽高（Parse `0.489 vs 0.259`；Tahoe `0.578 vs 0.219`），情境特异算子恢复 `g` 却跨零或显著为负（Parse `+0.015 [-0.615, +0.726]`；Tahoe `ST-SE -0.497 [-1.334, -0.084]`）；② **`Molina & Zhang (2026)` ([10.64898/2026.07.24.740459](https://doi.org/10.64898/2026.07.24.740459))**：在 4 细胞系上将零样本输出拆解为 $\mu+\alpha_c+\beta_p+\gamma_{c,p}$，所有 8 个模型（`STATE`、`MORPH`、`BioLORD`、`scGen`、`Ridge`、`MLP`）在零样本下均只恢复情境共享扰动分量 $\beta_p$（`0.510–0.788`），情境交互分量 $\gamma_{c,p}$ 全为近零或负值（`STATE γ = -0.123..+0.134`, 残差 $\epsilon = -0.0008..+0.0226$；`MORPH γ = -0.291..-0.143`），且情境无关 `MLP/Ridge` 胜过深度模型；③ **`CASTER` (`Kawsar et al. 2026`, [10.64898/2026.11.15.766815](https://doi.org/10.64898/2026.11.15.766815))**：零样本 `Task T2` 下 `CASTER`（`0.1714`, discrimination `0.2179`）显著**输给**无目标情境信息的跨线简单平均 `PerturbMean`（`0.1981`, `0.3682`）。同时，`Systema`、`Miller 2026`、`arXiv:2608.00152` 及早期 `scGPT/GEARS` 批判文献多测**同细胞系内未见扰动**或**线性尺度/采样深度混杂**，`AMM-SimWMag` 则依赖**目标对照幅度校准 $\hat{s}_h$** 而非无情境基线——**仍留出 Team A R10-E 在 VCC/FCT 上“常数基线即可复现 `r_orth_cal`”的空白**。 |

---

## 一、Q1 详查：`CellularAI`、`AtlasShift` 与 VCC 2025/2026 方案中的“去目标线基线 / 常数基线”消融

### 1.1 结论：**「没有找到」**（VCC 2025/2026 体系内零先例）

我们逐一代码级审计了 VCC 2025/2026 关键公开仓库与论文，所有方案在计算表达变化或生成细胞时均将目标细胞系未扰动对照（或合成对照）作为硬性前置条件，**无一团队测试过“移除目标对照基因表达特异性”或“将目标对照基因向量替换为全基因常数向量”后评估指标（如 `r_orth_cal` / 传统相关）是否坍塌**：

1. **Illumina `CellularAI`（VCC 2026 #1）**：
   - **代码与公式审计**（[GitHub `Illumina/CellularAI`](https://github.com/Illumina/CellularAI)）：在 `src/cellularai/prediction/synthesis.py:30-47` 中，`Engine` 强制读取目标对照计数矩阵并计算每基因对照比例 `self.fc = sum(X_ctrl)/sum`，并在 `design_frac(l2)` 中执行 `f = self.fc * exp2(l2)`；此外 `src/cellularai/prediction/plan.py` 的 `similarity_weights()`（L30–41）、`power_score()`（L131–140）以及 `src/cellularai/prediction/totals.py` 的 `make_spread()`（L52–57）全部显式读入目标对照 `fc` / `mean_cpm`。
   - **评估与消融审计**：`docs/evaluation.md:37-38` 与 `src/cellularai/evaluation/data.py` 定义的唯一对比锚点是 `Baseline anchor`（跨靶点均值）与 `Replicate anchor`（半分割重复实验）；检索整个仓库 `src/`、`docs/`、`configs/` 与 `tests/`，**不存在任何将 `self.fc` 替换为常数向量或跨线平均对照的消融配置或测试**。
   - **官方博客与待发论文**（[Illumina Blog](https://www.illumina.com/company/news-center/feature-articles/cellular-ai-virtual-cell-challenge-2026.html)）：博文 `Figure 7` 展示的消融仅比较“细胞系数量从 5 增至 50”及逐步加入标准化/加权/残差校正模块（相对于跨细胞系均等平均基线达到 2 倍增益），并明确声明完整方法学论文要到 **2027 年初（"early 2027"）** 才发布。

2. **`AtlasShift`（`kaipengm2/Virtual-Cell-Challenge-2026`，VCC 2026 #82）**：
   - **代码与公式审计**（[GitHub `kaipengm2/Virtual-Cell-Challenge-2026`](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026)）：审阅完整 `README.md`、`model.py`（L1–444）与 `predict.py`（L1–228）。其中 `load_xatlas`（`model.py:27-28`）提供 `matched_control_probability` 与 `global_control_probability` 两个分支，但该开关仅作用于**训练源数据 X-Atlas** 的伪计数收缩先验（`prior = counts + prior_counts * ctrl`）。到了目标细胞系预测端 `desired_mean()`（`model.py:112-130`），其 `log2fc` 与 `bulk_delta` 两个空间强制执行 `control_probability * np.exp2(effect)` 与 `np.expm1(np.maximum(np.log1p(50000 * control_probability) + effect, 0))`，`control_probability` 恒为目标细胞系实测对照概率分布，**全仓无任何消融参数或常数基线测试**。

3. **其他 VCC 2025/2026 公开方案与对比基准**：
   - **Xaira `X-Cell`（VCC 2026 #3）**：[Xaira-Therapeutics/X-Cell](https://github.com/Xaira-Therapeutics/X-Cell) 文档（`docs/model.md`）显示其在训练期将对照细胞随机 25%/50%/75%/100% 基因位置替换为扰动真值并在推理端以目标未扰动对照群精细扩散演化；[Arc VCC 2026 Wrap-Up](https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up) 补充其参赛提交用参考 CRISPRi 效应按对照相似度加权并用 `X-Cell` 门控缩放。未发现去目标对照或常数基线消融。
   - **Arc `State`（VCC 2025/2026 官方基线模型）**：[Adduri et al. 2025](https://doi.org/10.1101/2025.06.26.661135) 与 [ArcInstitute/state](https://github.com/ArcInstitute/state) 中的 `PerturbMeanPerturbationModel`（`perturb_mean.py`）和 `ContextMeanPerturbationModel`（`context_mean.py`）分别计算 `ctrl_mean[c] + offset[(c, p)]` 与训练集内该细胞类型的平均扰动表达；无目标基线置常数消融。
   - **Altos Labs `PRiMeFlow`**：[`arXiv:2604.13986v2`](https://arxiv.org/abs/2604.13986)（评估于 PerturBench 与 VCC 2025 H1 私有测试集，第 4.4 节及附录 B.3）在推理期同时以非对称无分类器引导权重生成合成对照（`CFG=5`）与扰动细胞（`CFG=20`），消融实验比较 U-Net vs MLP、PCA/OT vs 基因空间及微调范围（200/250/300 靶点），未做常数基线消融。
   - **Luan et al. (2026)**（[bioRxiv `10.64898/2026.10.20.760514`](https://doi.org/10.64898/2026.10.20.760514)，回顾性评估 VCC 2025 `TEST100`）：将外部 `KOLF2.1J` 响应中心化对齐到 H1 训练靶点后用等渗回归映射到 H1 表达均值；同样不含去目标基线或常数基线消融。
   - **`ds4cabs/VirtualCellChallenge2025`**：仅封装 `State` 流程文档，无相关消融。
   - **FlowState（VCC 2026 #2）**：在 OpenAlex 与 arXiv 检索 `FlowState Virtual Cell Challenge` **无独立论文或公开消融报告**（截至 2026-11-22 仅有 Arc VCC 2026 Wrap-Up 的方法摘要）。

### 1.2 必须写明的“非 VCC 紧邻消融”（边界差异对照）

虽然 VCC 2025/2026 体系内无相关消融，但在非 VCC 的 4 细胞系（K562/RPE1/HepG2/Jurkat）部分筛查跨系迁移研究 **`CASTER`（[Kawsar et al. 2026, bioRxiv `10.64898/2026.11.15.766815`](https://doi.org/10.64898/2026.11.15.766815) / [`CanPert` 源码](file:///Users/mac/.gemini/antigravity/brain/c4d0971d-8452-497f-b5d2-466ba007691c/scratch/lit/caster_repo/CanPert/src/canpert/train.py#L224-L404)）** 中，作者为了回应审稿人对模型是否真正利用目标对照细胞的质疑，实现了三项对照替换消融（`Task T1` 含 25%–30% 目标支持扰动）：
1. **`CASTER_meanprofile`**（`ctrl_mode="mean"`）：把一组异质对照细胞替换为目标均值向量 `ix.ctrl_mean` 的重复拷贝，去掉单细胞离散度，仅保留目标系平均表达；
2. **`CASTER_label`**（`ctrl_mode="label"`，`line_label()` `train.py:224-230`）：对 `ix.ctrl_mean` 做**固定基因排列置乱**，彻底摧毁真实逐基因基线表达值、仅作为细胞系身份标签；在全细胞校正指标 `pearson_deltapert_top100` 上，`CASTER_label` 的得分几乎等同于 `CASTER_meanprofile` 与完整 `CASTER`（`results_v1.py:988-990`）；
3. **`CASTER_wrongctrl`**（`wrong_ctrl_predictor()` `train.py:393-404`）：将**另一个细胞系**的对照细胞 `wc` 喂入模型编码器计算 `shift = base(p, wc).mean(0) - wc.mean(0)`，`pearson_delta` 从 `0.2759` 降至 `0.1516`（`discrimination` 从 `0.4236` 降至 `0.1805`，见 `t1_fold_means.csv`）。

**与 Team A R10-E `E1` 的关键结构差异**：`CanPert` 的 `wrong_ctrl_predictor` 与 `caster_predictor` 在返回最终预测时均严格写回真实目标对照细胞（`return ctrl + shift[None, :]` / `return base + (out_ - ctrl)`），其源码注释直言 *“so only the model's use of the control set changes, not the baseline expression level”*（`train.py:395-396`）。它检验的是**神经网络条件编码器 `f_ctx` 是否读取目标对照特征**，而不是检查 `log1p(c · exp2(L)) − log1p(c)` 这种乘法-对数映射在目标基线 `c` 本身被置为全基因常数后是否仍能凭几何平均与凹函数压缩得分。此外，[Molina & Zhang (2026)](https://doi.org/10.64898/2026.07.24.740459) 在补充表 `s4a_cross_features` 中测试过直接用 `Target ctrl` 线性回归跨系响应（仅 K562 `-0.045`, RPE1 `+0.083`, HepG2 `+0.041`, Jurkat `+0.001`），同样属于**加性特征消融**而非前向折叠比例压缩的常数基线消融。

---

## 二、Q2 详查：fold-change 几何平均或 `log1p` 非线性压缩单独制造扰动预测增益

### 2.1 结论：**「没有找到」**（针对扰动预测基准制造增益的完整命题；但存在 4 类紧邻数学/评测前作，须明确边界）

经过对 OpenAlex、bioRxiv、EuropePMC、PubMed、arXiv（2018–2026）及 `cell-eval2` / `scPertEval` 代码库的系统检索：**尚未发现任何已发表研究或预印本指出“在目标对照基线退化（甚至为全基因常数）的情况下，仅靠在 fold-change 空间做几何平均并使用 `log1p(c · FC) − log1p(c)` 非线性前向压缩，即可在扰动预测基准（或去振幅后指标 `r_orth_cal`）上制造显著正向增益”。**

不过，文献中**已有多篇方法学与评测工程工作揭示了 `log1p` 与均值顺序失真的底层数学性质**，在论文综述与讨论中必须精准区分“前人已指出的 DE/方差/尺度现象”与“Team A 在 `R10-E` 发现的新伪影机制”：

| 文献 / 代码实证 | 核心发现（已核验原文） | 与 R10-E (`E1`) 命题的本质区别 |
|---|---|---|
| **Booeshaghi & Pachter (2021)**<br>*Bioinformatics* ([`10.1093/bioinformatics/btab085`](https://doi.org/10.1093/bioinformatics/btab085)) | 通过泰勒展开从理论上证明：先逐细胞做 `log1p(x)` 再取跨细胞算术平均 `mean(log1p(x))`，再求两组差值，会在低表达区（$\lambda \approx 1$ 附近偏离 $\log 2$ 最大）**系统性扭曲真实 fold change 估计**，其本质是可交换期望与凹函数带来的**差异表达（DE）估计偏差**。 | 讨论的是**观测量两样本 DE 统计量（LFC）估算误差**，未讨论跨细胞系扰动响应迁移，未分析 `log1p(c·exp2(L)) − log1p(c)` 如何在预测基准上与真值产生相关性或虚高去振幅指标。 |
| **Kviman et al. (2025/2026, v2 2026-11-11)**<br>*"LN's t-test"* ([`10.1101/2025.03.12.642799`](https://doi.org/10.1101/2025.03.12.642799)) | 证明单细胞软件默认的基于 `log1p` 的 t 检验在检验聚类 DE 时因 LFC 估计有偏而产生**严重膨胀的假阳性率**，提出渐近无偏 LFC 检验量。 | 针对单细胞聚类标志基因检验的 I 类错误控制；不涉及跨细胞系扰动响应模型评估，亦未讨论 fold-change 空间几何平均对扰动基准分数的贡献。 |
| **Arc `cell-eval2` (`v0.15–v0.22`)**<br>[`docs/metrics.md` §1.2](https://github.com/ArcInstitute/cell-eval2) · `src/cell_eval2/norm.py:40-75` · `CHANGELOG.md` | ① **Jensen 凹性亏空（`#258/#260/#264/#288`）**：实测逐细胞均值 `mean_c log1p(s·x/L)` 会使重构总质量 $\sum_g \operatorname{expm1}(b_{p,g})/s$ 相对 `bulk_lognorm` 亏空 **42.8%**（`0.572` @ $s=10^6$），同一组预测在两种空间下的 `expr_mse_unbiased_capped_norm` 相差 **`+0.8` 至 `+1.4`**。正因 `Xaira X-Cell` 按逐细胞 `lognorm` 均值训练而官方用 `bulk_lognorm` 评测，其在 VCC 2026 必须回退到参考效应用校准重算（[VCC 2026 Wrap-Up](https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up)）。<br>② **`mean_calc: geometric vs arithmetic`（`CHANGELOG.md`）**：明确指出以 `expm1(mean(log1p))` 定义的几何平均 LFC 依赖 `target_sum` 缩放，而先线性归一再取算术平均满足缩放不变性。 | 注意**作用方向与机理完全不同**：`cell-eval2` 的 Jensen 凹性亏空是**单细胞细胞间异质性压缩了对零离散模型（或未加分散度模型）的惩罚（penalty）**，即模型因未模拟细胞内方差而在 `lognorm` 指标下受损；而 R10-E `E1` 发现的是即使使用全基因常数基线 `c`，**跨源细胞系 fold change 几何平均 + `log1p(c·FC)−log1p(c)` 正向变换本身会在目标基因相对表达压缩曲线下创造显著正相关得分**（相当于利用 `log1p` 非线性抑制高/低区间比值权重、复现 86%–112% 的真实 FCT `r_orth_cal`）。 |
| **Liu et al. (2025/2026)**<br>[`arXiv:2511.16954`](https://arxiv.org/abs/2511.16954) | 针对 VCC 2025 使用的 Perturbation Discrimination Score（PDS）证明：**$\ell_1/\ell_2$-PDS 对预测向量的全局倍乘常数 $c\hat{\Delta}_i$ 高度敏感**，且发现对 VCC 训练原始计数采用 `per-10k + log1p` 与采用 `median` 归一化的真实扰动向量虽然余弦相似度高，但 $\ell_1/\ell_2$ 范数差异巨大；促使 VCC 2026 改用 `pds_cosine`。 | 分析的是**线性整体幅度标量 $c$ 和评价距离范数选择（$\ell_1/\ell_2$ vs 余弦）**对排行指标的影响，未探讨“去振幅后依然保留的余弦/相关增益是否仅来源于几何平均 + `log1p` 非线性压缩”。 |
| **Panagopoulos et al. (2026-10)**<br>[`arXiv:2610.26267`](https://arxiv.org/abs/2610.26267) & **Ahlmann-Eltze & Huber (2023)** ([`10.1038/s41592-023-01814-1`](https://doi.org/10.1038/s41592-023-01814-1)) | `2610.26267` 证明把单个伪批量预测向量重复拷贝成多细胞行会放大 Mann-Whitney U 检验的标准化 z 值（67.1% 检验对假阳性通过），导致 Precision@K 随截断长度 K 改变 **86.7%** 榜首。`Ahlmann-Eltze 2023` 评测 `log1p` 与残差/潜在变换对下游降维聚类的影响。 | 二者均未讨论扰动预测任务中 fold-change 几何平均与 `log1p` 压缩如何在缺失目标基线信息的情况下制造正向预测相关性。 |

---

## 三、Q3 详查：是否已有人证明“细胞情境特异”扰动预测收益可被“无情境信息”基线复现

### 3.1 结论：**「已有人做过」**（多项 2025–2026 独立研究直接证明“零样本/跨情境表面高分主要由跨情境共享均值解释，真实情境特异分量近零或被无情境基线打败”）

针对指令明确点名的各路基准与论文（并补充 VCC 2026 官方复盘与同期大型基准），我们核验后的逐项判定如下：

#### A. 直接正面证明“情境特异收益可被无情境基线（或共享效应基线）复现/击败”的核心研究（3 项）

1. **`CGC` (`Huang, Wang, Li & Wilson, 2026-10-01`, [`10.64898/2026.09.30.755756`](https://doi.org/10.64898/2026.09.30.755756), [GitHub `cgc-counterfactual-resolution`](https://github.com/yongqih/cgc-counterfactual-resolution)) — 直接定罪供体预测与 `Arc State` 官方零样本模型**：
   - **T 细胞供体预测分解（`source_data/Figure_1.csv`）**：在 Feng 人类 T 细胞 Perturb-seq 留一供体测试中，`Ridge`、`Bilinear` 与 `MLP` 模型在传统供体预测 Pearson 上高达 **`0.987–0.991`**；但拆解后发现 **`shared_prediction_energy_fraction`（跨供体共享预测能量占比）高达 `0.9718`（Ridge）、`0.99969`（Bilinear）与 `0.999998`（MLP）**，严格留一供体特异适配 `R²` 仅为 **`-0.00389` / `-7.57e-5` / `-2.23e-7`**（`operator_alpha ≈ 0.0065–0.0067`）。即 97.2%–99.9998% 的模型输出本质上是无供体特异性的共享规则。
   - **对 `Arc State` 官方零样本跨情境预测的审计（`Extended_Data_1_cytokine.csv` & `Extended_Data_1_drug.csv`）**：CGC 将 `Arc State`（`ST-SE` / `ST-HVG`）在 Parse 10M 细胞因子（5 免疫亚群）与 Tahoe-100M 药物（5 癌细胞系）上的官方零样本保留情境预测，直接与**留一情境无情境均值基线 `Pert Mean`** 对比：
     - **Parse 细胞因子**：`State` 在常规 `Delta-Pearson` 上远高于 `Pert Mean`（`0.4893` vs `0.2586`），但相对于 `Pert Mean` 与 44 个精确置换零分布（context derangements）校准后的情境特异算子恢复分数 `operator_g = +0.0153`（`95% CI [-0.6146, +0.7261]`），**统计上不能拒绝零**，其中 `Plasmablast` 情境甚至低至 **`-0.7525`**。
     - **Tahoe-100M 药物**：`State` 常规 `Delta-Pearson` 为 **`0.5782` vs `Pert Mean` `0.2192`**，但校准后的 **`Pooled ST-SE operator_g = -0.4968`（`95% CI [-1.3342, -0.0836]`，仅 `1/5` 情境为正，中位 `g = -0.7828`）**，根本原因是其预测方差能量膨胀（`kappa = 2.2877`）。
     - **加性主效应 vs 交互效应分解（`Figure_2.csv`）**：在 Tahoe `50 contexts × 93 drugs` 矩阵上，`Additive main effect`（仅情境均值+扰动均值叠加，不含交互）恢复 `g = 0.3104`（`[-0.092, +0.440]`），加入 `Low-rank interaction` 毫无提升（`g = 0.3053`，差值 `-0.0051 [-0.066, +0.031]`）；即使在观测覆盖度达 **99.98%（4,649 / 4,650 格子）** 的全掩码测试中，`Context-specific g` 也只有 `0.1107`（同期 `Full response g = 0.6945`）。

2. **`Molina & Zhang (2026-07-27)` ([`10.64898/2026.07.24.740459`](https://doi.org/10.64898/2026.07.24.740459), [GitHub `xinyizhanglab/perturbation-decomposition`](https://github.com/xinyizhanglab/perturbation-decomposition)) — 直接分解 `STATE`、`MORPH`、`BioLORD`、`scGen` 等 8 个模型在 4 细胞系零样本迁移中的分量来源**：
   - 将 K562/RPE1/HepG2/Jurkat 响应分解为四分量 ANOVA $\delta_{c,p} = \mu + \alpha_c + \beta_p + \gamma_{c,p}$（`fig1f_anova_results.json` 显示保守扰动分量 $\beta_p$ 占总信号 `29.4%`，去模板后占 `35.0%`；系×扰动交互 $\gamma_{c,p}$ 占 `23.5%` / `27.7%`）。
   - **独立成分对齐直接实证（`final/components/results.json` & `fig4_crosscl_decomp_results.json`）**：在无目标扰动（`source_only` 零样本跨系）条件下，全部 8 个模型都强烈对齐**不含系别特异性的保守扰动主效 $\beta_p$**：`Ridge 3→1` `0.613`、`MLP 3→1` `0.788`、`MORPH` `0.678`、`STATE` `0.510`、`BioLORD` `0.603`、`scGen` `0.529`；而对**细胞系特异交互 $\gamma_{c,p}$ 的相关性几乎全为零或负值**：
     - `STATE (source_only)`: `γ = +0.134 (K562), +0.020 (RPE1), -0.123 (HepG2), -0.019 (Jurkat)`；去模板残差 $\epsilon$ 相关仅 `0.0087, -0.0008, +0.0066, +0.0226`；
     - `BioLORD (source_only)`: `γ = -0.098, -0.058, -0.123, +0.005`；
     - `scGen (source_only)`: `γ = -0.186, -0.224, -0.121, -0.051`；
     - `MORPH (controls / source_only)`: `γ = -0.291` 至 `-0.133`；去模板残差 $\epsilon$ 为 `-0.0075` 至 `+0.0188`。
   - 此外，在 3→1 跨系总迁移精度上（`fig3_b_saved.json`），不使用目标对照自适应的 DepMap 简单模型 `MLP 3→1`（`0.226–0.339`）与 `Ridge 3→1`（`0.119–0.287`）均显著超越 `STATE`（`0.028–0.086`）、`BioLORD`（`0.039–0.087`）和 `scGen`（`0.077–0.103`），且往 DepMap 输入额外加入 `Target ctrl` 表达在 `s4a_cross_features` 中几乎零提升（K562 `0.0960→0.0961`, RPE1 `0.1713→0.1714`, HepG2 `0.2386→0.2383`, Jurkat `0.1804→0.1811`）。

3. **`CASTER` (`Kawsar et al., 2026-11-17`, [`10.64898/2026.11.15.766815`](https://doi.org/10.64898/2026.11.15.766815), `CanPert/figure_tables/tables/t2_fold_means.csv`) — 零样本自适应输给纯跨线无情境均值基线**：
   - 在 **零样本未见细胞系任务 `Task T2`**（四系留一跨线）中，作者设计的“利用目标对照集上下文注意力缩放跨线先验”模型 `CASTER` 在所有关键指标上均**低于无任何目标情境自适应信息的 `PerturbMean`**（其他 3 线同一扰动 $\Delta$ 的等权平均）：`pearson_delta` **`0.1714` vs `0.1981`**，`discrimination` **`0.2179` vs `0.3682`**，`de_auprc` **`0.1061` vs `0.1353`**，`lfc_spearman` **`0.3401` vs `0.3974`**（文中事前假设检验 `H2` 明确报告 *`Not supported; PerturbMean significantly better on all three`*，`results_text.py:272-273`）。在含 30% 目标支持扰动的 `Task T1` 中，完全看不到目标对照细胞的 `Linear` 基准（`0.2821 / 0.4258`）同样持平甚至微超 `CASTER`（`0.2759 / 0.4236`）。

#### B. 对指令点名其余各文献与重要近期基准的逐一边界核验（防混淆、防误引）

| 文献 / 框架 | 真实考察的任务设定与基准结论 | 是否“证明细胞情境特异预测收益可被无情境基线复现”？ |
|---|---|---|
| **Arc Institute `VCC 2026 Wrap-Up`** (2026-11-19, [`arcinstitute.org`](https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up)) | 官方验证线实测：Top-10 零样本模型在**不同细胞系间对同一靶点给出的预测余弦相似度高达 `0.77–0.89`**（真实跨系观测仅 `0.55`），行特异斜率仅 **`0.02–0.13`**（重复实验 `1.02`），且与目标实测的一致度（`0.38–0.52`）反而不如“直接拿另一个实测细胞系当预测”（`0.55`），明确指出：*“current models recover roughly what a cross-line average would predict, whereas the part of the response that is specific to a cell line... remains largely missing”*。 | **强烈佐证 Q3（群体级竞赛实证）**：虽未公布每支前十队伍单独剥离细胞系特异基线的消融得分表，但从官方诊断指标上定性定量证实 Top-10 方案集体退化为跨线均值搬运。 |
| **`Systema`** (`Viñas Torné et al. 2025/2026`, *Nat Biotechnol*, [`10.1038/s41587-025-02777-8`](https://doi.org/10.1038/s41587-025-02777-8), `PMC13271886`) | 检验 **同一数据集内未见单/双基因扰动**（`scGPT`、`GEARS` 等），证明传统 `Pearson(Δ)` 虚高主要源于“所有扰动相对于对照存在统一系统性偏移”，致使无扰动特异信息的 `perturbed mean` / `matching mean` 基线就能复现高相关（`scGPT r=0.91`、`GEARS r=0.95`）。 | **不直接回答跨情境问题（针对单系内系统偏移，非跨系特异增益）**：`Systema` 在每个数据集独立评测并构建去偏移指标，未对比跨细胞系情境特异声称与跨系无情境基线。 |
| **`Miller et al. (2026)`** (*Nat Biotechnol*, [`10.1038/s41587-026-03307-w`](https://doi.org/10.1038/s41587-026-03307-w)) | 专注**单系内未见基因扰动与组合扰动**，用 `WMSE` 与动态范围分数 `DRF` 重新校正确认“经过恰当校准的指标上深度模型可显著拉开与 `mean baseline` 的差距”；作者在局限段逐字明确声明：**“leaving the unseen context task (predicting responses in unobserved cell types, donors or conditions) unaddressed”**。 | **明确排除**：作者本人已书面声明**未覆盖**未见细胞情境任务。 |
| **`AMM-SimWMag`** (`Sun et al. 2026`, *Genes*, [`10.3390/genes17070816`](https://doi.org/10.3390/genes17070816)) | 在跨细胞类型/跨物种迁移（`Kang IFN-β`、`Haber`、`Hagai`）中证明：**零神经网络（非深度方法）即可打平并超越 `scGen`、`scPRAM`、`biolord`、`CPA`**。**关键陷阱与精细边界**：纯无情境平均 `LinearDelta`（`mean_other`，`R²-DEG = 0.810`）其实**落后于** `scGen`（`0.868`）；只有额外引入**基于目标对照 ISG 特征的幅度标量校准 $\hat{s}_h$**（升至 `0.872`，匹配 `scGen`）及**目标对照相似度加权方向 `simw_delta`**（升至 `0.924`）后才打平/反超深度模型。 | **需严格区分！** 证明的是“神经自编码器在跨类型迁移上的收益可被**利用目标对照统计量做线性校准的封闭式方法**复现”，而非“不含任何目标情境信息的无情境基线即可复现”。 |
| **`arXiv:2608.00152v3`** (`Shoeibi & Yousefi 2026`, [`arXiv:2608.00152`](https://arxiv.org/abs/2608.00152)) | 审计冻结 `Geneformer` 在 VCC 2025 H1 hESC 与 Replogle（`K562/RPE1/Jurkat/HepG2`）上的表现：发现虽在分发内获得正向 `Delta R² = +0.1645`，跨系零样本迁移却完全坍塌为负相关（`Spearman ρ = -0.139 / -0.267`），且分发内敲低幅度预测高度纠缠细胞采样深度 `n_cells`（`ρ = +0.705`；纯细胞数即可解释 `R² = +0.4325` vs 真实均值变化增量仅 `+0.0017`）。 | **揭示跨系崩溃与测序计数混杂，但边界不同**：研究的是标量敲低幅度估计从内部分发到跨系的失效与采样深度混杂，未比较跨系向量模型与“不含情境信息基线”的收益复现。 |
| **`scGPT / GEARS` 早期简单基线批判文献组**<br>（`Ahlmann-Eltze 2025` [`10.1038/s41592-025-02772-6`](https://doi.org/10.1038/s41592-025-02772-6)、`Csendes 2025` [`PMC12016270`](https://pmc.ncbi.nlm.nih.gov/articles/PMC12016270/)、`Li 2025` [`PMC12202205`](https://pmc.ncbi.nlm.nih.gov/articles/PMC12202205/)、`PertEval-scFM` [`10.1101/2024.10.02.616248`](https://doi.org/10.1101/2024.10.02.616248)、`Kernfeld 2025` [`10.1186/s13059-025-03840-y`](https://doi.org/10.1186/s13059-025-03840-y)、`Kedzierska 2025` [`10.1186/s13059-025-03574-x`](https://doi.org/10.1186/s13059-025-03574-x)、`Bendidi 2024` [`arXiv:2410.13956`](https://arxiv.org/abs/2410.13956)） | 系统证明 `scGPT`、`scFoundation`、`Geneformer`、`UCE`、`scBERT`、`GEARS`、`CPA` 在**未见单基因或双基因扰动**上不优于甚至输给简单 `mean` / `additive` / `PCA` / `linear` 基线；其中 `Ahlmann-Eltze 2025 Extended Data Fig. 9` 对 K562→RPE1 训练了跨系线性模型，但目的是比较预训练源（扰动图谱 vs 单细胞观察图谱）并说明跨系误差随基线差异上升，而非拆解情境特异收益。 | **主要批判“扰动特异预测”（同情境未见扰动）**：它们证明深度学习的收益不如训练均线/加性线；除少数附带实验外，主体并非针对“跨细胞情境特异迁移声称”。 |
| **同期补充重要跨情境基准**：<br>① **`PerturBench` (`Wu et al. 2024/2025`, [`arXiv:2408.10609`](https://arxiv.org/abs/2408.10609))** & **`PRiMeFlow` (`Yan et al. 2026`, [`arXiv:2604.13986`](https://arxiv.org/abs/2604.13986))**<br>② **`Li et al. (2026)`** (*Science Advances*, [`10.1126/sciadv.aed3414`](https://doi.org/10.1126/sciadv.aed3414), `PMC13557097`)<br>③ **`Dibaeinia et al. (2026)`** ([`10.64898/2026.02.04.703804`](https://doi.org/10.64898/2026.02.04.703804), `PMC12919078`)<br>④ **`Ghajari et al. (2026)`** ([`10.64898/2026.11.04.764088`](https://doi.org/10.64898/2026.11.04.764088)) | ① `PerturBench` 协变量迁移中，不含交互项的 `Latent Additive` (`x + f(p, cov)`) 与不输入表达谱的 `Decoder-Only` 在大数据集（`Jiang24`）平均响应指标上匹敌或超越 `CPA*`/`SAMS-VAE*`/`PRiMeFlow`（见 `PRiMeFlow Table 6`: `Decoder`/`LA logFC RMSE = 0.018` vs `PRiMeFlow 0.019`、`CPA 0.020`，但注意其接收 `cov_one_hot` 协变量加性嵌入，属于已见协变量组合缺失问题）；<br>② `Li et al. 2026` 在跨细胞类型 `Task 3` 显式设置了无情境跨系平均基线 `Baseline3_PerturbMean`（Eqs. 8–9）并进行 `B_c + G_c + S_p + C_{c,p}` 分量分解，得出无方法能稳定跨细胞类型泛化、表现主要由跨系生物学一致性主导；<br>③ `Dibaeinia et al. 2026` 以“仅按扰动平均（context-free `P`）”为基线评测 `scLDM` 供体×时间点迁移，指出尽管 `Corr-Delta ≥ 0.7` 但 `DEG-F1 < 0.1`（召回仅 `~0.09`）；<br>④ `Ghajari et al. 2026` 揭示免疫供体筛查中大量所谓“供体特异信号”源于与空白孔共享的孔板技术偏移。 | **共同强化 Q3**：完整覆盖了协变量加性近似、跨系 `PerturbMean` 显式基线、以及伪装成生物学特异性的技术孔偏移；**但以上工作均未分析在 VCC 2026 / `CellularAI` 式 `FCT` 管道中，去除了目标对照基线的常数前向映射如何仅靠几何平均与 `log1p` 压缩便在 `r_orth_cal` 上复现“表面目标系特异增益”**——而这正是 Team A `R10-E` 实证的核心缺口与独特贡献。 |

---

## 四、对 Team A `R10-E` 合稿与方法论定位的直接启示

1. **Q1 + Q2 构成 R10-E 极其罕见且未被人占先的原始诊断机制**：
   - 在 VCC 2026 最具代表性的公开参考方法（`CellularAI`、`AtlasShift`）以及所有公开发表的 `log1p`/评测反思论文（`cell-eval2`、`Booeshaghi 2021`、`Kviman 2026`、`arXiv:2511.16954`）中，**从未有人测试并揭示“目标对照常数替换（E1）+ 跨源几何平均 FC + `log1p(c·FC)−log1p(c)` 正向非线性压缩”本身就会在去振幅后相关指标 `r_orth_cal` 上人工制造出相当于真实目标基线模型 `86%–112%` 的正向增益**（`DLD-1 +0.157 vs +0.140`、`HepG2 +0.226 vs +0.261`、`H1 +0.308 vs +0.350`）。
   - 该结果与 `cell-eval2` 的发现（per-cell `lognorm` Jensen 凹性带来的是惩罚而非增益）正好互补，构成了对 `log1p` + 倍数缩放（fold-change transfer）在零样本跨系预测评测中产生虚假表面进展的独立机理贡献。

2. **Q3 文献定位与严谨引用要求（合并手稿时必须遵守）**：
   - 论文中**绝不能声称**“首次发现声称具有细胞情境特异性的扰动预测模型其实只恢复了情境无关均值（或被无情境基线击败）”——因为 [`CGC` (`Huang et al. 2026`)](https://doi.org/10.64898/2026.09.30.755756)、[`Molina & Zhang (2026)`](https://doi.org/10.64898/2026.07.24.740459) 与 [`CASTER` (`Kawsar et al. 2026`)](https://doi.org/10.64898/2026.11.15.766815) 已先后对 `Arc State`、`MORPH`、`BioLORD`、`scGen` 与 `CASTER` 零样本模式给出了严格定量证明；同时 [Arc VCC 2026 Wrap-Up](https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up) 已对 Top-10 参赛方案给出跨系平均坍塌结论。
   - **Team A 的精准、牢不可破的新颖定位应聚焦于**：在前人证明“模型跨系收益只来自情境共享平均值 $\beta_p$/`Pert Mean`”的基础上，进一步揭示**即使在明确剥离全系振幅（$\beta_t$ 移除后 `r_orth_cal`）的细粒度目标系残差检验上，当前跨系 fold-change 缩放转移（FCT）看似“仅存的少数目标系特异优势（H-C2 未通过 `0/2` 特异性，而 FCT 此前呈正）”仍几乎完全由常数基线下的几何平均 + `log1p` 压缩伪影（`E1`）驱动**，从而把“虚假情境特异增益”的拆解从**加性均值（`Pert Mean` / $\beta_p$）层面**深入到了**非加性对数折变压缩层面**，并为后续 `E6/E7` 验证提供了精确的阴性对照与机理分解准则。

---

# R11 文献补充核查（Q4 & Q5）

- **补充核查时间**：2026-11-22T13:30+08:00
- **任务依据**：[`R11_第十一轮行动指令.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/directives/R11_%E7%AC%AC%E5%8D%81%E4%B8%80%E8%BD%AE%E8%A1%8C%E5%8A%A8%E6%8C%87%E4%BB%A4.md) §A-T5

## R11 执行摘要：Q4 / Q5 三档结论总览表

| 问题编号 | 核心核查问题 | 三档结论 | 精准范围与关键边界说明 |
|---|---|---|---|
| **Q4** | 是否有人指出：在扰动预测评测中“按真实观测响应选基因”（如按真实 `|Δ|` 取前 K 个基因、按真实 DEG 评测）会造成选择偏倚，让无信息变换也得高分？ | **「部分做过（写清做到哪一步）」** | **两类直接邻近工作已揭示“依赖真值选基因”造成的得分虚高，但主流仍提倡 DEG/Top-K 评测，且无人点出 `E1`/置换在无噪声/跨系基线下的 `resp_top50` vs `all_nontarget` 偏倚**：<br>① **`Nicol, Shivakumar & Irizarry (2026)`** ([`10.64898/2026.05.07.723486`](https://doi.org/10.64898/2026.05.07.723486), [GitHub `phillipnicol/systema`](https://github.com/phillipnicol/systema) `Figure 2` & `ctl-mean-topk.ipynb`) 在模拟与 10 个 Perturb-seq 数据集上直接证明：当按 `|True − Control|` 选 `top k = 20` 基因并共享对照样本误差时，随机噪声相关从 `r = 0.528` 放大至 **`r = 0.869`**，零变化 `Null Baseline`（`ctl-mean`）最高刷到 **`+0.408`**；但该工作归因于**同数据集有限样本共享对照噪声与同批选基因耦合**，一旦做 split-half（`Split DEG + control`）即归零（`-0.044..+0.047`）。<br>② [`cell-eval2`](https://github.com/ArcInstitute/cell-eval2)（`#351`）仅修了**表达量存在性门槛** `perturbed > 5 CPM OR control > 5 CPM` 泄露 LFC 正号的问题（盲目抬高低于 5 CPM 基因曾赚得 `+0.372..+0.506` 方向保真），而在 DE 指标中**继续使用**参考显著 DEG 集合 $S_p^{\text{real}}$ 与 `top_k`。**相反地**，`Miller 2026`（明确倡导 `top 100 DEGs`/`WMSE`）、`SBB`（明确推广 `top20` DEG 指标与 `wMSE`）、`Systema`（报告 `PearsonΔ20` 但归因系统偏移）、`scPerturBench`（`DEG=100`）、`arXiv:2608.00152`（train HVGs）均未警告 `resp_top50` 会相对于 `all_nontarget` 抬高**不含目标基线信息但非常数的几何平均 `log1p` 压缩变换（`E1`）与标签置换**（`H-D4`）。 |
| **Q5** | 是否有人报告过：跨细胞类型的扰动响应在 fold-change 空间比在 `log1p` 差值空间更守恒，或者线特异的 `log1p` 差值主要来自基线表达差异？ | **「部分做过（写清做到哪一步）」**<br>（子命题 A：**「没有找到」**；<br>子命题 B：**「部分做过」**） | **子命题 A（fold-change 空间比 `log1p` 差值更守恒）：「没有找到」**。`CellularAI`、`AtlasShift` 与 `Nadig 2025`/`sigRecon` 在实践中选择在 `log2FC` 空间工作，而 `Molina & Zhang 2026`、`CASTER`、`Li 2026`、`scPerturBench` 在 `log1p(CP10k)` 差值空间评测；**未见任何公开论文或预印本在同一批多细胞系扰动数据上直接对比两种表示的跨系保守性并报道倍数变化高于 `log1p` 差值**。<br>**子命题 B（线特异 `log1p` 差值主要来自基线表达差异）：「部分做过」**。已有四类局部相关报道：① `Molina & Zhang (2026)` 证明 `log1p` 分解中的**扰动无关**系别模板 $\alpha_c$ 有 `72.3%–88.8%` 位于前 5 个未扰动对照 PC 内，但对真正依赖靶点的系别交互 $\gamma_{c,p}$ 反而指出目标对照几乎无预测力；② `AMM-SimWMag (2026)` 在免疫刺激数据证明跨类型**整体标量振幅**差异可由目标基线 ISG 表达解释（`R²-DEG 0.810→0.872`）；③ `Ahlmann-Eltze et al. (2025)` `Extended Data Fig. 9` 指出 K562→RPE1 转移误差随两系未扰动基线表达差异增大而增大；④ `Booeshaghi & Pachter (2021)` 推导了低基线 $\lambda$ 处相同倍数变化被 `log1p` 压缩的公式。**但无人提出或证实：逐基因线特异 `log1p` 差值（乃至剥离全局振幅 $\beta_t$ 后的 `r_orth_cal` 胜过 `E1` 之差）主要由“目标线实测基线 $c_{t,g}$ 对相对保守的跨系 fold-change 经 `log1p(c_{t,g}·FC) − log1p(c_{t,g})` 非线性缩放”产生。** |

---

## 五、Q4 详查：按真实观测响应选基因（`top-K` / real DEGs）导致的选择偏倚与无信息变换虚高

### 5.1 结论：**「部分做过（写清做到哪一步）」**

在扰动预测评测工程与统计学批判文献中，**已有两条直接相邻但机理边界不同的线索**揭示了“按真实观测端选基因”会抬高无信息基线；但另一方面，**绝大部分主流基准非但没有警惕 `resp_top50` 类指标相对全基因 `all_nontarget` 的选择偏倚，反而在 2025–2026 年强力推行 `top-K` / DEG 加权评测**。具体逐项核验如下：

#### A. 直接涉及“选基因依赖于观测答案会奖励无信息基线”的 2 项先导工作（精确区分到哪一步）

1. **`Nicol, Shivakumar & Irizarry (2026-05-07)`**（*Spurious correlation inflates performance in single-cell perturbation prediction*, [`10.64898/2026.05.07.723486`](https://doi.org/10.64898/2026.05.07.723486), [GitHub `phillipnicol/systema`](https://github.com/phillipnicol/systema)）：
   - **做到哪一步**：该文的主体贡献是证明共享同一批对照均值 $C$ 计算 `(True − C)` 与 `(Pred − C)` 会在有限对照样本下产生虚假相关性；但尤为关键的是，其 **Figure 2（`spurious_correlation_paper/figure_2.qmd`）** 与 **对照均线测试（`split_control/ctl-mean-topk.ipynb`）** 专门考察了 **`top k = 20` Pearson $\Delta$**：
     - 在 Figure 2a 模拟中（独立正态 `True`, `Predicted`, `Control`，`n = 1000`），按 `abs(True - Control)` 取 **`top k = 20`** 后，虚假相关从全基因 `r = 0.528` 进一步**暴涨至 `r = 0.869`**；
     - 在 10 个 Perturb-seq 数据集实测中（`ctl-mean-topk.ipynb`），代表**纯粹无信息零变化的控制均值基线（`Null Baseline: ctl-mean`）** 在“同样本选择 top-20 DEGs + 共享对照”（`Shared DEG + control`）下取得显著正向相关：**`TianCRISPRi2021 = +0.4082`、`TianCRISPRa2021 = +0.3226`、`XuKinetics2024 = +0.1549`、`FrangiehControlSingle2021 = +0.1176`**（而 `ReplogleK562 = +0.0657`, `ReplogleRPE1 = +0.0590`）。一旦改为独立样本对半挑选 DEGs 并拆分对照细胞（`Split DEG + control`），`ctl-mean` 即塌缩回 **`-0.0441` 至 `+0.0471`**（零附近）。
   - **尚未覆盖到 `R11`（`H-D4`）的本质区别**：Nicol 等剖析并解决的仍是**同一批单细胞有限样本误差耦合**（共享随机噪声 $-C$ 被 `|True − C|` 尾选放大 + 选基因半与估计半未独立）。而在 Team A 的 `R10/R11` 设定（`H-D4`）中：即便对照样本独立/足够大且脱离共享随机噪声，`E1`（`常数基线 × 跨系几何平均 FC` 经 `log1p` 压缩）是**确定性、跨系聚合的非均匀向量**（并非 `ctl-mean` 零变化，也不共享目标单细胞采样误差 $-C$）；按目标实测 `|Δ|` 选 `resp_top50` 仍会利用极值尾选择（range restriction / tail amplification）把 `E1` 与保守强响应轴的微小耦合及置换极值大幅抬高，从而在 `E1(resp_top50) − E1(all_nontarget)` 上制造出 Nicol 的 split-control/split-DEG 方案无法消除的虚假高分。

2. **`ArcInstitute/cell-eval2`（`v0.22.0`）与 VCC 2026 DE 指标**（[GitHub `ArcInstitute/cell-eval2`](https://github.com/ArcInstitute/cell-eval2), `docs/vcc2026_metrics/vcc2026-metrics.md:352-635`, `src/cell_eval2/de_compute.py:1212-1243`, `docs/metrics.md:788-808, 1108-1135, 1674-1735`）：
   - **做到哪一步**：
     - **`#351` 存在性泄漏修复**：指出旧版 `cell-eval` 用 `(perturbed > 5 CPM) OR (control > 5 CPM)` 筛选检测基因违背 **"the surviving gene set must not depend on the answer"**（`vcc2026-metrics.md:362`）——对对照 `≤ 5 CPM` 但扰动 `> 5 CPM` 的基因，`100.0000%` 的真实 `log2FC > 0`（每上下文 `26,373 / 33,969 / 26,839` 条）；一个仅在对照细胞里对这些低于 5 CPM 基因加常数量、完全不看任何靶点信息的盲目基线，竟能凭此拿下 `de_wilcoxon_direction_fidelity_yield_raw` 的 **`+0.3722 / +0.5059 / +0.3651`** 及总分 **`+0.1057`**（`de_compute.py:1222-1235`）。修复版改用**仅由对照 `> 5 CPM` 决定的测试宇宙**。
     - **防常向量作弊与固定门控**：对可加选的 `delta_pearson_degs`（`docs/metrics.md:788-808`）及 VCC 计分指标 `de_wilcoxon_lfc_nmae`（`vcc2026-metrics.md:579-620`），强制要求 `min_gate_size >= 10`（避免 2 个基因相关恒为 $\pm 1$）、要求预测不能全常（全常相关性罚 `-1`，`#92/#391`）、用 `S_p^real` 参考集固定评估队列防缺失值刷分（`#89`），并实测通用均值基线在 `de_wilcoxon_lfc_nmae` 上仅得 **`1.0009–1.0017`**（与预测全零 `1.0` 无技能点相当），同时针对 `overlap_top50..top500` 在 §5 用 $2\times 2$ MCC 校正大基因集下的随机交集膨胀。
   - **尚未覆盖到的空白**：`cell-eval2` 本身仍然**常规依赖**参考真值 DEG 集合 $S_p^{\text{real}}$（`de_wilcoxon_lfc_nmae`、`de_wilcoxon_direction_reach_raw`、`delta_pearson_degs`）与按真值 `|lfc|` 排序的 `top_k`（`de_wilcoxon_overlap_top50` 等），**未分析或警示**：当测试基线是类似 `E1` 这种“不含目标基线信息但保留几何平均方向排序”的**非均匀变换**、或置换零分布时，仅在真值高响应子集（`resp_top50`）评估连续相关性与方向对齐会因尾部极值选取比在 `all_nontarget` 上系统性地给出更高分数。

#### B. 指令点名核心基准与近期文献的对照核查（绝大多数甚至相反主张“应该仅在 Top-K / DEG 上评测”）

| 文献 / Benchmark | 对“按真实观测 DEG / Top-K 选基因评测”的真实立场 | 是否指出 Q4 选择偏倚（奖励无信息变换）？ |
|---|---|---|
| **`Systema`** (`Viñas Torné et al. 2025/2026`, *Nat Biotechnol*, [`10.1038/s41587-025-02777-8`](https://doi.org/10.1038/s41587-025-02777-8)) | 同时用 `PearsonΔ`（全基因）和 `PearsonΔ20`（按观测真值 `|ΔX|` 取前 20 DEGs）评估，发现无扰动特异信息的 `perturbed mean` 在 `PearsonΔ20` 上同样接近或击败深层模型（仅 `Frangieh18` 除外）。 | **否**：作者将其归因于**全数据集系统偏移（systematic perturbation-vs-control shift）与扰动面板选偏**，并在校正后的 `Systema` 指标中也包含 top-20 版本，**未**比较 `top20` 相对于 `all_nontarget` 对常数/置换变换本身的放大作用。 |
| **`Miller et al. (2026)`** (*Nat Biotechnol*, [`10.1038/s41587-026-03307-w`](https://doi.org/10.1038/s41587-026-03307-w)) & **`Mejia et al. (2025)`** ([`arXiv:2506.22641`](https://arxiv.org/abs/2506.22641)) | **完全相反（大力提倡）**：主张按每扰动 **top 100 DEGs**（`Extended Data Fig. 5`）或按 vs-Rest DEG 显著性加权（**`WMSE`** / $R_w^2(\Delta)$）能去除海量非响应基因的“信号稀释”、显著提高指标区分真/假对照的 `DRF` 校准度；仅注明计算 DEG 掩码/权重时应排除充当正对照的半分割技术重复细胞以防泄漏。 | **否**：不仅未警告选择偏倚奖励无信息变换，反而将其视为比全基因更优越的黄金评测子空间。 |
| **`SBB`** (`Vollenweider & Bühlmann 2026`, [`10.64898/2026.04.20.719650`](https://doi.org/10.64898/2026.04.20.719650)) | **明确提倡**：主推 `top_n_deg_*`（如 `top20_mse`、`top20_pearson_dp`）与 `wMSE`，展示 `Frangieh '21` 上全基因 `MSE` 的 `BDS ≈ 0%` 而 `wMSE` 升至 **`92%`**（`discrimination_analysis.md:335-445`）。 | **否**：未报道 `top20` 掩码会在置换或 `E1` 基线上抬高虚高分数。 |
| **`scPerturBench`** ([`10.1101/2024.12.20.629581`](https://doi.org/10.1101/2024.12.20.629581) / `Wei et al. 2026`, *Nat Methods*, [`10.1038/s41592-025-02980-0`](https://doi.org/10.1038/s41592-025-02980-0)) | 代码库（`calPerformance.py:33-56, 98, 135`, `calPerformance_genetic.py:25-46, 98, 147`）机械地同时跑 `DEG=100`（`DEG_hvg5000.pkl` 真值前 100 DEG）与 `5000` HVG 全集并对比 `trainMean.py`。 | **否**：仅作双尺度报表，无对“按观测响应选 `DEG=100` 引发选择偏倚”的任何机制批评。 |
| **`arXiv:2608.00152v3`** (`Shoeibi & Yousefi 2026`, [`arXiv:2608.00152`](https://arxiv.org/abs/2608.00152)) | 其锁定基准评估标量敲低幅度 `AD`，特征空间固定取训练集挑选的 **`d = 2,000` HVGs**，不按观测响应挑选评价基因；文中只引述 `Mejia 2025` 用 DEG 加权损失抗模态坍塌。 | **否**：未讨论按观测 `|Δ|` 选 top-K/DEG 基因的评测选择偏倚。 |
| **`scPertEval`** (`Schäfer et al. 2026`, [`10.64898/2026.07.23.740433`](https://doi.org/10.64898/2026.07.23.740433)) | 在架构上将空间明确划分为“**全数据集固定且独立于扰动的对照高表达面板** `heg_1000`（Ahlmann-Eltze 2025 方案）”与“**每扰动真值响应子集** `top_k` / `degs_padj`（Miller 2026 方案）”（`catalog.py:33-51`），并探讨细胞数影响 DE 检验力；且演示 `de_overlap_k=10` 在 `Wessels23` 上 `BDS = 0.382 < 0.5`。 | **边缘铺垫但未点破核心**：阐明了两组空间“是否依赖每扰动真值”，却仍按 `DRF`/`BDS` 把 `top_k`/`degs` 列为标准协议，未检验 `top_k` 对常数基线/置换的虚假抬高。 |
| **`Panagopoulos et al. (2026)`** ([`arXiv:2610.26267`](https://arxiv.org/abs/2610.26267)) | 指出复制单向量构造成多细胞矩阵会使 **模型预测侧** 出现 `67.1%` 的假显著 DEG 对，导致 `Precision@N` 排行随截断阈值 `N` 翻转 **`86.7%`**；并明言数值精度（如 `MSE`/`Pearson`）在不同 $N$ 下相对稳定（§4.1.1）。 | **否（批判对象不同）**：讨论的是“**预测端**伪重复膨胀导致的 `Precision@N` 不稳定”，而非“**真值端**按 `|Δ|` 选基因让 `E1`/置换等无信息变换得高分”。 |
| **`CGC`** ([`10.64898/2026.09.30.755756`](https://doi.org/10.64898/2026.09.30.755756), `Extended Data Fig. 2d` & `run.py`) | 在训练输出端检验支持度选基因 vs 全基因，发现 `support_minus_full = -0.0355 [-0.0666, -0.0044]`（选基因反而丧失下游全响应效用）。 | **角度不同**：讨论的是**模型输出支撑裁剪损失下游效用**，非“评测阶段按真实 `|Δ|` 截断造成无信息基线偏高”。 |
| **同样本“double dipping”经典统计工作**（`ClusterDE` [`10.1101/2023.07.21.550107`](https://doi.org/10.1101/2023.07.21.550107)、`countsplit`） | 针对“同一批细胞既聚类又做聚类间 DE 检验”会制造假阳性 DEG。 | **相邻统计原理，不同任务**：仅涉无监督聚类后 DE 检验的 I 类错误，未涉跨系/跨扰动预测基准的 `resp_top50` vs `all_nontarget` 评估。 |

---

## 六、Q5 详查：跨细胞系扰动在 `fold-change` vs `log1p` 差值空间的保守性及基线表达驱动假说

### 6.1 结论：**「部分做过（写清做到哪一步）」**（严格拆分：子命题 A「`fold-change` 空间比 `log1p` 差值更守恒」＝ **「没有找到」**；子命题 B「线特异 `log1p` 差值主要来自基线表达差异」＝ **「部分做过」**）

#### A. 子命题 5a（`fold-change` 空间跨细胞类型是否比 `log1p` 差值空间更守恒）：**「没有找到」**
在检查了 Perturb-seq 跨系图谱、药物迁移、以及 VCC 2025/2026 方法与基准文献后：**未发现任何已发表研究或预印本直接设计实验对比“同一批多细胞系扰动响应在 fold-change（比值/`log2FC`）空间 vs `log1p` 差值空间的跨系保守性相关/一致性”，并报告 fold-change 空间更守恒**。两类表示目前在实际工程中各走一边、缺少交叉定量对决：
- **直接在 `log2FC` / 相对比例空间做跨情境比较或迁移的代表**：
  - Illumina [`CellularAI`](https://github.com/Illumina/CellularAI)（`docs/method.md`）直接在居中 `log2FC` 空间对跨系摘要进行随机效应汇集，再通过 `self.fc * exp2(l2)` 映射到目标对照比例；完整方法学论文延至 2027 年初，目前公开源码与[官博](https://www.illumina.com/company/news-center/feature-articles/cellular-ai-virtual-cell-challenge-2026.html)均**未报告**“为什么在 `log2FC` 空间而不直接在 `log1p` 差值空间汇集更守恒”的定量消融；
  - [`AtlasShift`](https://github.com/kaipengm2/Virtual-Cell-Challenge-2026)（`model.py:89-130, 221-227`, `predict.py:196`）按经验权重同时混合 `space="log2fc"`（权重 `0.6`）与 `space="bulk_delta"`（`log1p(50000·p) − log1p(50000·c)`，权重 `0.3`），且实现了 `log1p(baseline·exp2(cd4)) − log1p(baseline)` 的跨空间映射，但未发布二者的保守性比较；
  - [`Nadig et al. (2025)`](https://doi.org/10.1038/s41588-025-02169-3)（*Nat Genet*, `PMC13063516`）使用 `TRADE` 测量 K562/RPE1/Jurkat/HepG2 间 `log2FoldChange` 的转录组级相关（中位数 `0.40–0.75`，指出未校准点估计中位相关仅 `0.16`），与 [`sigRecon`](https://doi.org/10.64898/2026.08.14.744937) 一样全程采用 `log2FoldChange`，从未与 `log1p` 差值空间对比。
- **直接在 `log1p` 差值空间做跨情境分解与评测的代表**：
  - [`Molina & Zhang (2026)`](https://doi.org/10.64898/2026.07.24.740459)、[`CASTER`](https://doi.org/10.64898/2026.11.15.766815)、[`Li et al. (2026)`](https://doi.org/10.1126/sciadv.aed3414)、[`scPerturBench`](https://doi.org/10.1038/s41592-025-02980-0) 全程直接在 `log1p(CP10k)` 差值 $\Delta$ 空间测量跨系相关性与泛化，同样**没有**在倍数变化空间再测一次跨系保守性进行横向比较。

#### B. 子命题 5b（线特异 `log1p` 差值是否主要源于基线表达差异）：**「部分做过（写清做到哪一步）」**
文献中确实存在 **4 项局部涉及“基线表达差异影响跨系 `log1p` 响应/迁移差异”的实证与数学工作**，但每一项都有严格的限定范围，且**无一揭示 R10-E / R11 的非线性折叠映射机制**：

1. **`Molina & Zhang (2026)`（[`10.64898/2026.07.24.740459`](https://doi.org/10.64898/2026.07.24.740459)）— 对全局模板 $\alpha_c$ 成立，但对靶点特异交互 $\gamma_{c,p}$ 恰恰给出否定结论**：
   - 作者在 `log1p` 四分量分解 $\delta_{c,p} = \mu + \alpha_c + \beta_p + \gamma_{c,p}$ 中把 $\alpha_c$ 定义为 *cell-line-specific baseline (template)*（占全信号 `27.8%`，各系 `23.4%–64.6%`），并在 `fig1e_ctrl_space_results.json` 中证明 **`α_c` 方差的 `72.3%–88.8%` 落在未扰动对照表达的前 5 个主成分内**（前 100 PC 内达 `91.8%–98.3%`；而特异残差在 100 PC 内仅 `34.6%–39.3%`）。
   - **关键边界**：这里的 $\alpha_c$ 是**所有扰动共享的固定细胞系偏移**（在去振幅 / 去模板后完全消失）；而在讨论随靶点变化的真正线特异部分 $\gamma_{c,p}$ 时，作者反而发现往跨系线性模型加入 `Target ctrl` 特征**几乎零增量**（`s4a_cross_features` `-0.045..+0.083`），说明加性线性映射下的基线表达并不能解释 $\gamma_{c,p}$——**Molina 并未测试非线性的 $\log1p(c_{t,g}\cdot\mathrm{FC}_{p,g}) - \log1p(c_{t,g})$ 映射如何在基线 $c_{t,g}$ 驱动下直接生成依赖扰动 $p$ 的系别差值残差**。

2. **`AMM-SimWMag`（`Sun et al. 2026`, *Genes*, [`10.3390/genes17070816`](https://doi.org/10.3390/genes17070816)）— 仅针对免疫刺激的单一细胞类型振幅标量**：
   - 在 `Kang IFN-β`、`Haber`、`Hagai` 等细胞因子/病原体刺激跨类型迁移中（`log1p(CP10k)`），指出单一无缩放平均 `LinearDelta` 低估强响应类型、高估弱响应类型，而仅用**目标细胞类型未扰动对照的基线 ISG 表达均值**推导一个**全基因共享的标量放大系数 $\hat{s}_h$**，即可将 `R²-DEG` 从 `0.810` 提升至 `0.872`（逼近先知标量 `0.877`）。
   - **关键边界**：$\hat{s}_h$ 是**全基因同一个乘数**（在去除了目标系全振幅 $\beta_t$ 的 `r_orth_cal` 评测中会彻底消去），并未研究 CRISPRi 多靶点下**逐基因基线表达 $c_{t,g}$ 对 fold-change 进行非线性 `log1p` 压缩所造成的基因级分布差异**。

3. **`Ahlmann-Eltze et al. (2025)`（*Nat Methods*, [`10.1038/s41592-025-02772-6`](https://doi.org/10.1038/s41592-025-02772-6), `Extended Data Fig. 9`）— 基线差异与跨系预测误差的经验相关性**：
   - 作者在 `K562 ↔ RPE1` 跨系预训练线性模型分析中专设 **`Extended Data Fig. 9`**：*"Analysis how differential expression between K562 and RPE1 effects prediction accuracy of transfer learning"*，定量给出读出基因预测误差（panel b）及按扰动 `Pearson Δ`（panel c）均随 K562 与 RPE1 间**未扰动基线表达差异**的增大而显著变差，正文中明确表述 *"The predictions were more accurate for genes that were more similar between K562 and RPE1"*。
   - **关键边界**：报告的是跨系预测错误率随基线表达差异上升的经验关系，**未区分**这种误差下降是源自真实调控重连，还是仅仅因为 `log1p` 空间里相同的相对倍数变化在两端基线不同时体现出不同的绝对 `log1p` 增量。

4. **`Booeshaghi & Pachter (2021)`（*Bioinformatics*, [`10.1093/bioinformatics/btab085`](https://doi.org/10.1093/bioinformatics/btab085)）与 `CASTER`（[`10.64898/2026.11.15.766815`](https://doi.org/10.64898/2026.11.15.766815)）**：
   - 前者在单细胞数学公式上说明了 $\log(1+\mathrm{FC}\cdot \lambda)-\log(1+\lambda)$ 仅在大基线 $\lambda\gg 1$ 时趋近 $\log(\mathrm{FC})$、小基线 $\lambda$ 处随 $\lambda$ 缩小（`btab085 Fig. 1c`），但未联系跨细胞系扰动响应差异；后者虽在合成模拟（`results_text.py:208-212`）验证当情境特异分量显式依赖基因基线时注意力网络可从对照恢复它，但在真实 4 系零样本 `Task T2` 中同样发现目标对照无法预报先验残差（反而输给 `PerturbMean` `0.1714 vs 0.1981`）。

---

## 七、对 Team A `R11` 前瞻检验（`PREREG_D` / `H-D1..H-D4`）与合稿叙事的直接启示

1. **Q4 为 `H-D4` 与主掩码切换（`resp_top50 → all_nontarget`）提供了清晰且不被抢先的方法学定位**：
   - 虽然 `Nicol et al. (2026)` 已证明**同数据集共享对照误差 + 同批选取 `top-20` 基因**会把随机噪音和零变化 `ctl-mean` 相关推高至 `+0.408..+0.869`（且 `cell-eval2 #351` 封堵了 `>5 CPM` 存在性门控泄漏），但业界权威基准（`Miller 2026`、`SBB`、`scPerturBench`、`Systema`、`cell-eval2 delta_pearson_degs`）目前仍然普遍认可甚至推荐在**真实高响应/显著 DEG 子集**上评测相关性。
   - Team A 在 `R10` 事后发现并于 `R11` **严格前瞻检验**的 **`H-D4`**（即即便在跨细胞系 `r_orth_cal` 指标下，没有任何目标对照基因特异性的**确定性常数基线变换 `E1`** 与标签置换，仍会在 `resp_top50` 上拿到虚假正向得分、而在切换至 `all_nontarget` 后塌缩回 0 附近），属于对现有基于观测响应子集评测体系的重要盲区补充；引用 `Nicol et al. (2026)` 与 `cell-eval2 #351` 作为先导并在边界上明确区分（样本共享误差 vs 确定性非均匀变换在极值尾上的对齐放大），论据将无懈可击。

2. **Q5 为 `H-D1..H-D3`（全基因口径下 `FCT > E1`）精准划定了生物学可解释性边界**：
   - 在文献中**无人**证实“跨系扰动在 `fold-change` 空间比在 `log1p` 差值空间更守恒”，亦**无人**报道过“跨系逐基因特异 `log1p` 差值主要由保守倍数变化经目标线实测基线 $c_{t,g}$ 压缩产生”。因此，若 `R11` 在全新前瞻细胞系上同时证实了 `H-D1..H-D3`（`FCT` 在 `all_nontarget` 下显著为正、越过错配与置换零线、且显著击败常数基线 `E1`）与 `E6/E7` 检验，Team A 将同时确立两项文献尚属空白的贡献：**① 揭示 `FCT` 对目标线真实基线表达 $c_{t,g}$ 的非线性 `log1p(c_{t,g}·FC) − log1p(c_{t,g})` 缩放确实携带跨系可迁移的非加性信号（超越 `Molina & Zhang` 与 `CASTER` 在加性/黑盒特征下“目标基线无助于 $\gamma_{c,p}$”的悲观结论）；② 揭示该表面增益同时伴随的评估伪影边界（`E1` 在 `resp_top50` 下的选择偏倚与在 `all_nontarget` 下的剥离）**。
