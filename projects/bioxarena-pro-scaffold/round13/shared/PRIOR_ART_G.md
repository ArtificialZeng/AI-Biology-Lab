# R13 文献检索与 Prior-Art 审计报告（PRIOR_ART_G v2）

- **智能体**：A-T5（文献与官方指标）
- **修订版本**：v2（2026-11-23 依据主控第 1 轮核验 [`PRIOR_ART_G_REVIEW.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/status/PRIOR_ART_G_REVIEW.md) 修订；旧版已备份为 [`PRIOR_ART_G_v1.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round13/shared/PRIOR_ART_G_v1.md)）
- **检索截止日**：`2026-11-23`（覆盖截至当日之 `bioRxiv`、`arXiv`、`EuropePMC`、`PubMed`、`Crossref`、`OpenReview` 及公开基准代码仓）
- **检索对象**：针对 [`R12_VERDICT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/R12_VERDICT.md) 两条核心主张（**Claim 1**：`FCT` 对全源错配基线 `Mis_s` 的跨系与跨模态 `r_orth_cal` IUT 检验；**Claim 2**：`E1_ConstMed` 在分半响应幅度筛选基因集 `resp_top50_split` 下于 `r_orth_cal` 维持正偏置）进行投稿前先占评估与诚实重叠度判定（`完全 / 部分 / 无`；不写任何对外宣传文字）。

---

## 一、主源可读性与信息边界声明（`未能读取` 严格标注）

1. **Nicol, Shivakumar & Irizarry (2026)**（`doi:10.64898/2026.05.07.723486`）：
   - **v1 全文 (`2026-05-12`)**：通过 `PMC13192888.1` XML 全文读取完毕（仅含 Figure 1 与全基因拆分对照分析，无 top-k DEG 分析）。
   - **v2 预印本正文/补充说明 (`2026-08-11`)**：`bioRxiv` HTML 与 JATS XML 均返回 `HTTP 403 Forbidden`，**`未能读取`**（仅 Crossref / EuropePMC API 摘要与元数据可读）。
   - **v2 官方配套代码与输出摘要（`nicol_systema.zip`, `github.com/phillipnicol/systema`）**：完整读取 `spurious_correlation_paper/figure_1.qmd`, `figure_2.qmd`, `figure_2_recompute.py`, `figure_2_files/figure-html/unnamed-chunk-{7,8}-1.png`, `split_control/eval_split_control_multifold.ipynb`, `split_control/eval_split_control_topk.ipynb`, `split_control/ctl-mean-topk.ipynb`, `evaluation/README.md`；另注 `results/*.csv` 原始逐条件文件未随仓提交，相关汇总均读自已执行并保存输出单元格的 notebook。
2. **Heidari, Karimpour, Srivatsa & Montazeri (2026)**（`doi:10.64898/2026.02.14.705879`）：
   - **v1 预印本全文 (`2026-02-17`)**：`bioRxiv` HTML 与 JATS XML 返回 `HTTP 403 Forbidden`，且未在 NCBI PMC / arXiv 收录，正文与论文结果表 **`未能读取`**（仅 EuropePMC `PPR1223517` 与 Crossref API 元数据及摘要可读）。
   - **官方配套代码仓（`github.com/LCBB-UT/scPerturb-eval`, `main` 分支）**：完整读取 `README.md`, `sc_proj/run_pipeline.py`, `sc_proj/utils/assessment.py`, `sc_proj/utils/common.py`, `sc_proj/models/{baseline_logfc.py,noperturb.py,random_all.py}`, `sc_proj/configs/data/*.ini`, `sc_proj/configs/models/*.ini`。
3. **我方内部定稿基准**：[`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md)（公式与数据集出处）及 [`R12_VERDICT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/R12_VERDICT.md)（数值结果与范围限定）。其中 **A549 测试集（`a549_ko`）出处为 Liu et al. (2026), `GSE345058`, `doi:10.64898/2026.06.01.728087`**（`CRISPRko`, `N=445`, `|G_eval|=6,121`），并非 "Norman A549"。

---

## 二、核心结论摘要（依据主源核验更新）

1. **Claim 2 新颖性判定（诚实下调为「特定条件下之窄增量现象，机制未知且存在显著 prior-art 边界重叠」，评级：`部分`）**：
   - **先对 `E1_ConstMed` 与机制边界正本清源**：依据 [`PREREG_E.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/prereg/PREREG_E.md) L152，
     $$\hat{\Delta}^{E1}_{t,p,g} = \log(1 + c_{\text{med}} \cdot \overline{FC}_{p,g}) - \log(1 + c_{\text{med}}), \quad c_{\text{med}} = \mathrm{median}_g \max(\mu_{t,g}, 0)$$
     `E1_ConstMed` **包含源细胞系对同一扰动 $p$ 的平均倍数变化 $\overline{FC}_{p,g}$**，因此它**不是**「与扰动靶点无关」的空预测器，而是一个**「不使用目标细胞系自身对照基线谱 $\mu_{t,g}$ 的有信息迁移预测器」**。它与 Nicol 等人的零信息对照 `ctl-mean`（控制细胞分半均值）并不构成直接对等参照；`ctl-mean` 在分半选 DEG 下归零（`[-0.044, +0.047]`）既不与 Claim 2 冲突，也不能用来反衬 Claim 2 的独特性。
   - **机制状态诚实声明**：依据 [`R12_VERDICT.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round12/status/R12_VERDICT.md) §5.6，**Claim 2 现象的底层机制目前未知（unexplained）**；第 11 轮基于真实 $\mu$ 的无噪声模拟在全基因掩码（`all_nontarget`）下于 **`8/12` 条细胞系上符号预测错误**（模拟预测为正而实测为负）。`log1p` Jacobian 权重与基线丰度–幅度耦合**仅是未获证实的候选假设之一**，不能作为已确立机理解释。
   - **Nicol et al. (2026) v2 对「有信息但不针对目标/特征简化的基线」的实测核查与对 Claim 2 的新颖性影响**：
     - 从 v2 代码 `eval_split_control_topk.ipynb`（Cell 5）与 `eval_split_control_multifold.ipynb`（Cell 4）核实：**Nicol v2 已经正式提出了 `Split DEG + control`（Half A 独立选 top-20 DEG，Half B + 拆分控制评估），并直接报告了有信息、但忽略目标具体扰动身份的 `nonctl-mean`（Perturbed Mean）在该评估下在多个数据集上仍保持明显正相关**（如 `Adamson 0.7323`, `Norman 0.4993`, `ReplogleRPE1 0.4919`, `ReplogleK562 0.2322`, `TianCRISPRa 0.2123`, `XuKinetics 0.1409`, `FrangiehCoculture 0.1301`；`gears` 亦呈相似走势）。
     - 精密区别在于：（a）在 Nicol 的常规控制中心化 `PearsonΔ` 下，`nonctl-mean` 在 `Split DEG + control top-20` 与 `all genes`（`corr_split`: `0.7089`, `0.5313`, `0.4634`, `0.1888`, `0.2206`, `0.0865`, `0.1094`）之间**并无系统性大幅抬升**（10 数据集均差仅 `+0.0043`，范围为 `-0.0536` 至 `+0.0544`，且 `figure_2.qmd` 只在同一 `k=20` 下对比 `Shared` 与 `Split` 而未拿 `top-k` 对比 `all`）；（b）而我方 Claim 2 是在 **Systema 式 `Transfer_Mean` 正交去均值残差空间（`r_orth_cal`）** 下，`E1_ConstMed`（省略目标线基因特异基线，仅用标量 `c_med` + 源线 $\overline{FC}_{p,g}$）在全基因上归零（`-0.0201`），但在 `resp_top50_split` 上回升至 `+0.2565 [+0.2195, +0.2981]`（且来源端独立选基因 `H-E2b` 同样为 `+0.2766 [+0.2412, +0.3164]`）。
     - **诚实降级结论**：由于（i）`Split DEG + control` 协议及“有信息简化基线在 `Split DEG + control` 上可保留显著正分”已由 Nicol v2（代码与 Fig 2b）先占，（ii）`Heidari et al. (2026)` 已系统测试 `CrossSplit`、`*_vs_topdeg` 与 `baseline_logfc` / `rigorous_subset_pearson_delta`，且该两篇的 **bioRxiv 完整正文均因 `HTTP 403` 而 `未能读取`**，（iii）我方自身仅在单个高扰动数据集（Liu et al. 2026 `GSE345058` A549 CRISPRko）的 `r_orth_cal` 空间观测到此现象且**机制未知**——**Claim 2 的新颖性必须诚实下调：它不再宜作为普适性评估原理的核心发现，而应降级收窄为「在 `r_orth_cal` 正交残差空间与单数据集（GSE345058 A549）上观察到的待解实证细节：即使剔除目标线特异基线并采用 `split-half`（及来源端 `H-E2b`）选基因，`E1_ConstMed` 仍从全基因近零回复至 `~+0.26–0.28`」。**

2. **Claim 1（`2026-11-15` 至 `2026-11-23` 增量检索 + `Kawsar et al. 2026` 与 `Heidari et al. 2026` 边界核查）**：
   - `2026-11-15` 至 `2026-11-23` 新增检索中，**无任何文献**测试跨细胞系与跨扰动模态（`CRISPRi -> CRISPRko`）的“纯未扰动基线生物变换（`FCT`）对全部错配源基线（`Mis_s`）在 `r_orth_cal` 下的配对 Intersection-Union Test（`min_s (FCT - Mis_s) > 0`）”。
   - **Heidari et al. (2026)** `scPerturb-eval/sc_proj/models/baseline_logfc.py` 确实实现了一个概念上相近的 `logFC` 乘法迁移基线 `predict_baseline_logfc`（用源细胞系加权 `avg_log2FC` 乘以目标未扰动细胞计数 `ood_adata.X` 再行 CPM+log1p），但其仅在单模态 OOD 细胞类型/物种替换中评测，未做全源错配 `Mis_s` 的 IUT 检验、未做 `r_orth_cal`、亦未测 `CRISPRi -> CRISPRko`。
   - **Kawsar et al. (2026, `CanPert / CASTER`)** `CASTER_wrongctrl` 仅用于部分筛选场景 `t1`/`t1_grouped` 的神经网络控制输入单次循环替换消融，零样本 `t2` 中完全未运行（且 `CASTER` 本身劣于 `PerturbMean`），与 Claim 1 为 **`部分`** 重叠，不构成完全先占。

---

## 三、Claim 2 专项审计：两篇关键基准主源核查与领域对照

### 3.1 Nicol, Shivakumar & Irizarry (2026) v1/v2 主源逐项核查（含 `nonctl-mean` 与 `gears`）

- **标识符与主源分工**：`doi:10.64898/2026.05.07.723486`
  - **v1 (`2026-05-12`)**：`PMC13192888.1` XML（可读）。推导共享控制抽样噪声协方差 $\sigma^2 c_0$ 对 `PearsonΔ` 的上偏（Assumptions 1–3, Propositions 1–2），仅测试全基因下的 `nonctl-mean`、`ctl-mean` 与 `GEARS`（Figure 1）。
  - **v2 (`2026-08-11`) 预印本正文**：`bioRxiv` HTML/JATS XML 返回 `HTTP 403`（**`未能读取`**）。
  - **v2 代码与汇总输出（`nicol_systema.zip`, `github.com/phillipnicol/systema`）**：完整读取：
    - 全基因：`split_control/eval_split_control_multifold.ipynb`（Cell 4 输出，3 seeds × 100 次拆半控制对，**标准未正交 `PearsonΔ`**，含全预测条件）。
    - Top-20 DEGs：`split_control/eval_split_control_topk.ipynb`（Cell 5 输出，及 `spurious_correlation_paper/figure_2.qmd` + `figure_2_recompute.py`，限定 `1-gene` 且 `>=3` 细胞条件；`stratified_half_split` 将细胞对半分为 `selection_adata` Half A 选 `top_non_dropout_de_20`，与 `estimation_adata` Half B + 巢式拆半控制计算 `corr_split_topk`）。
    - Null 基线 Top-20：`split_control/ctl-mean-topk.ipynb`（Cell 5 输出）及绘图结果 `spurious_correlation_paper/figure_2_files/figure-html/unnamed-chunk-7-1.png`、`unnamed-chunk-8-1.png`。另注 `eval_split_control.ipynb` 虽定义过 `corr_all_allpert`/`corr_20de_allpert` 列但该笔记本仅含进度输出、未附汇总数值表。

- **10 数据集实测逐项数值对照（直接摘自 `eval_split_control_multifold.ipynb` Cell 4、`eval_split_control_topk.ipynb` Cell 5 及 `ctl-mean-topk.ipynb` Cell 5）**：

| 数据集 | `nonctl-mean` 全基因 `corr_all` / `corr_split` | `nonctl-mean` Top-20 `corr_topk` / `corr_split_topk` | `Δ_nonctl` (`split_topk − all_split`) | `gears` 全基因 `corr_all` / `corr_split` | `gears` Top-20 `corr_topk` / `corr_split_topk` | `ctl-mean` Top-20 `corr_shared_topk` / `corr_split_topk` |
|---|---|---|---|---|---|---|
| `Adamson2016` | `0.703510` / **`0.708907`** | `0.730102` / **`0.732299`** | `+0.023392` | `0.672909` / `0.668710` | `0.728977` / `0.730768` | `-0.001773` / `-0.008998` *(Fig 2 PNG ~`+0.11`/`+0.066`)* |
| `Norman2019` | `0.535293` / **`0.531299`** *(含2-gene)* | `0.520382` / **`0.499326`** *(1-gene)* | `-0.031973` | `0.526735` / `0.509418` | `0.575652` / `0.551945` | `-0.008885` / `-0.040212` *(Fig 2 PNG ~`+0.038`/`+0.031`)* |
| `ReplogleRPE1_v2` | `0.490658` / **`0.463411`** | `0.587265` / **`0.491871`** | `+0.028460` | `0.237314` / `0.196953` | `0.416549` / `0.327457` | `0.059049` / `-0.002293` |
| `ReplogleK562_gwps` | `0.265219` / **`0.188764`** | `0.308956` / **`0.232247`** | `+0.043483` | `0.096187` / `0.048291` | `0.122876` / `0.073066` | `0.065688` / `-0.004059` |
| `TianCRISPRa2021` | `0.500208` / **`0.220634`** | `0.695956` / **`0.212316`** | `-0.008318` | `0.296783` / `0.135980` | `0.561045` / `0.179153` | `0.322631` / `0.005242` *(Fig 2 PNG ~`+0.26`/`-0.04`)* |
| `TianCRISPRi2021` | `0.512875` / **`0.020516`** | `0.747156` / **`0.028830`** | `+0.008314` | `0.247289` / `0.039298` | `0.452512` / `0.041243` | `0.408225` / `-0.030718` |
| `XuKinetics2024` | `0.321754` / **`0.086467`** | `0.487223` / **`0.140881`** | `+0.054414` | `0.129382` / `0.061270` | `0.244837` / `0.148121` | `0.154863` / `-0.006349` |
| `FrangiehCocultureSingle2021` | `0.147906` / **`0.109382`** | `0.134603` / **`0.130103`** | `+0.020721` | `0.039099` / `0.036252` | `0.027646` / `-0.036682` | `0.030036` / `-0.044059` |
| `FrangiehControlSingle2021` | `0.151942` / **`0.124094`** | `-0.140260` / **`0.070477`** | `-0.053617` | `0.039196` / `0.042139` | `0.011711` / `-0.057360` | `0.117597` / `0.047121` |
| `FrangiehIfnSingle2021` | `0.121614` / **`0.092796`** | `-0.032583` / **`0.050984`** | `-0.041812` | `0.050425` / `0.047195` | `0.034640` / `-0.071639` | `0.077953` / `-0.006056` |

*(注：图 `unnamed-chunk-7-1.png` 的 Null Baseline 来自 `figure_2_recompute.py` / `results/avg_corrs_ctl_mean_topk.csv`，而 `ctl-mean-topk.ipynb` 为另一次单次随机种子输出，数值在 `Adamson`、`TianCRISPRa`、`XuKinetics` 等处略有差异，上表已并行标出；而 `nonctl-mean` 与 `gears` 两笔记本数值与 `unnamed-chunk-7-1.png` 散点高度吻合。)*

- **代码逻辑与作者意图核查（`figure_2.qmd`, `unnamed-chunk-8-1.png`, `PMC13192888.1`）**：
  1. 在 `figure_2.qmd` Panel a（`unnamed-chunk-8-1.png`）中，作者用纯随机高斯零信号证明：若共用控制细胞并再取 top-20 `|True - Control|`，相关性会从全基因 `0.528` 暴涨至 `0.869`；Panel b 对比固定 `k=20` 下的 `Shared DEG + control` 与 `Split DEG + control`，文字说明仅写 *"keeps differential-expression (DE) selection independent from metric estimation for the split estimator"*。
  2. **作者是否报告「有信息基线在分半 top-k 仍明显高于全基因」及如何解释**：从仓库代码与表格本身看，**`nonctl-mean` 在 `Split DEG + control` Top-20 下确实在多数数据集保留高正分（0.13–0.73），但相比自身全基因 `corr_split` 并无明显系统性跃升（`6/10` 微增 `+0.008..+0.054`，`4/10` 反而下降 `-0.008..-0.054`；均值 `0.2589` vs `0.2546`）**，因为在未去系统均值的 `PearsonΔ` 下，`nonctl-mean` 在全基因上本身就已带有共享扰动均值信号（`0.19–0.71`）。代码仓内未对 `nonctl-mean` 做 `corr_split_topk` vs `corr_split` 跨尺度对比；**至于 v2 论文正文中是否对此另作讨论，因 bioRxiv v2 全文 `HTTP 403`，`未能读取`。**
  3. **与我方 Claim 2 的精确对比**：`nonctl-mean` 是“忽略目标具体扰动 $p$、使用目标训练集所有受扰细胞均值 $\hat{\mu}_{\text{pert}}$”、且在常规未正交 `PearsonΔ` 下评估；而我方 `E1_ConstMed` 是“保留源线特定扰动平均倍比 $\overline{FC}_{p,g}$、但不用目标线基因基线 $\mu_{t,g}$”、在 `r_orth_cal` 下评估。虽然对象与度量化去均值空间不同，但鉴于 Nicol v2 已公开建立了 `Split DEG + control` 并在 Fig 2b 展示过信息基线在拆分选基因下的持续正分，Claim 2 已不再具有“发现一般性 split-half DEG 失效规律”的广义新颖性，应诚实降级并严格限域。
- **重叠度判定**：**`部分`**（评估方法学与拆分 DEG + control 框架先占）。

---

### 3.2 Heidari, Karimpour, Srivatsa & Montazeri (2026) 主源逐项核查

- **标识符与可读来源**：
  - 论文：*"Evaluating Single-Cell Perturbation Response Models Is Far from Straightforward"*，`doi:10.64898/2026.02.14.705879`（bioRxiv `2026-02-17`）。EuropePMC (`PPR1223517`) 与 Crossref 元数据/摘要可读；**bioRxiv 预印本全文 PDF/HTML 返回 `HTTP 403`，`未能读取`**。
  - 官方代码仓：`github.com/LCBB-UT/scPerturb-eval`（完整读取 `README.md`, `run_pipeline.py`, `utils/assessment.py`, `utils/common.py`, `models/{baseline_logfc.py,noperturb.py,random_all.py}`, `configs/`）。
- **逐项核实 5 个问题（全部引自代码仓，标注全文不可读部分）**：
  1. **评测范围与是否做 `top-DEG / top-k`**：**是**。覆盖 4 个数据集配置（`kang`, `datlinger`, `sciplex3`, `hagai`），在 `OOD`（未见细胞类型/物种，仅保留其对照）与 `PID`（含 20% 受扰样本）设定下系统评测。`run_pipeline.py` 与 `assessment.py` 明确包含：
     - `*_vs_topdeg` 曲线族（`spearman_vs_topdeg`, `norm_spearman_vs_topdeg`, `mmd_vs_topdeg`, `wasserstein_vs_topdeg`），按 MAST `p_val_adj` 排序依次评测 `num_top_genes = 2, 4, ..., 2048` 及全部基因；另有 `rigorous_spearman_topdeg`, `rigorous_mixing_index_{louvain,kmeans}_topdeg`, `mixing_kmeans_vs_topgenes`, `rigorous_deg_f1`。
     - `rigorous_subset_*` 基因分层族（含 `rigorous_subset_pearson_delta`, `rigorous_subset_pearson_r2`, `rigorous_subset_spearman`, `rigorous_subset_logfc_correlation`, `rigorous_subset_mmd`, `Edistance`, `c2st_distance`, `knn_loss`, `wasserstein`, `mixing_index`, `deg_f1`），通过 `define_trivial(true_adata, ctrl_adata, true_logfc_df)` 将 MAST `p_val <= 0.05` 显著基因分为 **`trivial`**（在对照或受扰中全零，或非零表达细胞占比失衡 `>10×`）与 **`non_trivial`** DEGs，另单独列 **`non_significant`** 基因。
     - 另外，`run_pipeline.py` L386/L1058 显式实现了 `# pearson_delta_brbic` -> **`rigorous_pearson_delta_pert`**（以非 OOD 源条件均值作为去偏参考），并与 `rigorous_pearson_delta`（以目标对照均值作为参考）对照。
  2. **是否用独立细胞选基因？**：**除 `rigorous_deg_f1` 外，代码中 `top-DEG` 与 `trivial/non_trivial` 基因集筛选未做独立细胞拆分**。核实 `run_pipeline.py`（L2074–2400，及 `*_vs_topdeg`, `semi_rigorous_nontrivial_*`）：`true_logfc_df = seurat_deg(ctrl_true_adata, ...)` 在进入 `CrossSplit`（`split_0..split_9` 将目标受扰细胞 `true_adata` 50/50 划分为 `train`/`test` 评估真值内在一致上限与模型得分）**之前**，即已在全量 OOD `ctrl_true_adata` 上一次性算出；随后循环用同一固定基因榜在 10 个 `split` 上打分。这不同于 Nicol v2 的 `stratified_half_split` 独立选基因与我方 `SA-8` `resp_top50_split`。
  3. **报告了哪些基线？**：代码 `sc_proj/models/` 与 `configs/models/` 包含 **3 个核心非深度学习基线** + 10 个深度模型（`ae, cae, cpa, scGPT, scGPT_pretrained, scPRAM, scVIDR, Squidiff, trVAE, trVAE_tf`）+ 评测内部参考：
     - **`noperturb`**（`models/noperturb.py`）：直接输出目标 OOD 细胞群对照细胞；
     - **`random_all`**（`models/random_all.py`）：从训练集**所有非目标细胞系受扰细胞**中随机抽取相同细胞数作为预测（即跨细胞系错配扰动池基线）；
     - **`baseline_logfc`**（`models/baseline_logfc.py`）：按源细胞系细胞比例对 MAST `avg_log2FC` 加权求平均，以 $\exp(\ln 2 \cdot \overline{\text{log2FC}})$ **乘以目标 OOD 未扰动对照细胞计数 `ood_adata.X`**，再做 `calculate_cpm`（`1e4`）+ `log1p`（属于保留目标对照基线的倍数变化迁移基线）；
     - `assessment.py` 中另有 CrossSplit `true`（目标细胞内部对半上限）与部分函数可选的 `best`（事后挑选最佳单个源细胞系）。
  4. **是否报告某个基线在 top 基因上被抬高？**：**论文正文中的具体对比数值与排序结果因 bioRxiv 全文 `HTTP 403` 而 `未能读取`**。从公开摘要（EuropePMC `PPR1223517`）及代码设计仅能确认：其揭示相关系数与分布距离受尺度、稀疏性与维度影响，并系统在 `top-DEG`（`*_vs_topdeg`）及 `trivial / non_trivial / non_significant` 子集上对比 `noperturb`、`random_all`、`baseline_logfc` 与深度学习模型；不可推测其正文中具体哪个基线得分高低或抬升幅度。
- **重叠度判定**：
  - 与 **Claim 1**：**`无`/`部分`**（`baseline_logfc` 与 `random_all` 思路涉足跨背景倍比迁移及随机错配对照，但无全源错配基线集合 `Mis_s` 的精确配对 `IUT`、`r_orth_cal` 正交去偏及跨 `CRISPRi -> CRISPRko` 模态验证）。
  - 与 **Claim 2**：**`部分`**（共享 `CrossSplit` 分半框架、按 `top-DEG` 随 $k$ 扫描、`trivial` vs `non_trivial` DEG 分解及 `rigorous_pearson_delta_pert` 去系统均值的想法；但缺少独立拆分选基因、缺少去掉目标基因基线的 `E1_ConstMed` 在 `r_orth_cal` 上的检验，且全文实测细节 `未能读取`）。

---

### 3.3 Claim 2 全领域对照文献清单（保留已核实基石文献并校正 A549 出处与 E1 表述）

*(注：本节保留 v1 已核实的全部 `bioRxiv`/`arXiv`/`Nature`/`OpenReview` 对照工作；以下条目已按主控审查要求去除「`E1` 与靶点无关」、「机制已确立为 `log1p` Jacobian」及「Norman A549」错误表述。)*

| 文献 / 工作 | 精确标识符 | 作者做了什么（已核实） | 未做什么（与当前边界后 Claim 2 差异） | 重叠度 |
|---|---|---|---|---|
| **Nicol, Shivakumar & Irizarry (2026)** | `doi:10.64898/2026.05.07.723486` (`PMC13192888.1` v1 + GitHub `phillipnicol/systema` v2；v2 bioRxiv 正文 `未能读取`) | v1 推导共享控制抽样噪声 `σ² c₀`；v2 代码建立 `Split DEG + control`（Half A 选 top-20，Half B + 拆半控制评估），令 `ctl-mean` 归零，并实测 `nonctl-mean` 仍获 `+0.13..+0.73`（与全基因接近）。 | 使用标准 `PearsonΔ` 而非正交校准 `r_orth_cal`；未检验省略目标线基线的 `E1_ConstMed`（全基因 `-0.020` -> `resp_top50_split +0.257` / `H-E2b +0.277`）。 | **`部分`** |
| **Heidari et al. (2026, `scPerturb-eval`)** | `doi:10.64898/2026.02.14.705879` + `LCBB-UT/scPerturb-eval`（bioRxiv 全文 `未能读取`） | 建立 CrossSplit、`*_vs_topdeg` 曲线、`define_trivial` 三分基因子集、`baseline_logfc` / `random_all` / `noperturb` 与 `rigorous_pearson_delta_pert`。 | `true_logfc_df` 取全量 OOD 样本而非独立拆半选基因；未测试不含目标基线的 `E1_ConstMed`（全文数值 `未能读取`）。 | **`部分`** |
| **Fukutani et al. (2026)** (`npj Syst Biol Appl`) | `doi:10.1038/s41540-026-00722-2` (`PMC13287743`) | 解析推导标准化中共享不稳定参照 `Z, W`（`X-Z`, `(X-Z)/W`）对独立随机变量诱导虚高相关性。 | 纯共享随机参照噪声效应；不涉及 `top-k` 响应幅度选基因或正交校准后 `E1_ConstMed`。 | **`部分`** |
| **Liang et al. (2026, `SafeDelta`)** | `doi:10.64898/2026.10.05.756906` (`bioRxiv`, `2026-10-13`) | 针对 few-shot 细胞系微调，证明共享控制组在反比例测序深度下携带控制噪声虚高内积/PCC；提出双半控制门控 `SafeDelta`。 | 同样研究控制抽样噪声携带（常数向量 PCC 设为 `NaN` 弃权），无 `resp_top50_split` 下 `E1_ConstMed` 升高现象。 | **`部分`** |
| **Huang (2026)** | `arXiv:2610.32962v2` (`2026-10-23` / v2 `10-26`) | 证明在混合 CRISPR 筛选中随机拆半共享同一转导克隆而导致克隆假重复（`G σ_a² / K_eff`），虚高 `cell-eval` 上限。 | 机理为细胞克隆假重复及静态全局基因筛选；未研究单扰动 `resp_top50_split` 下 `E1_ConstMed` 在 `r_orth_cal` 下的正偏移。 | **`部分`** |
| **Miller, Mejia et al. (2026)** (`Nat Biotechnol`) | `doi:10.1038/s41587-026-03307-w` (`arXiv:2506.22641v1`) | 指出未校准 `MSE` 与控制偏移指标因信号稀释偏袒均值基线；提出 DEG 加权 `WMSE` 与 `$R^2_w(\Delta)$` 使均值基线跌至零/负。 | 其发现为 DEG 加权压低均值基线；未考察 `r_orth_cal` 下 `E1_ConstMed` 在 `resp_top50_split` 反而呈现 `+0.2565` 正残留。 | **`部分`** |
| **Vollenweider & Bühlmann (2026, `SBB`)** | `doi:10.64898/2026.04.20.719650` (`bioRxiv`, v2 `2026-05-27`) | 提出 Signal–Bounds–Baselines 原则、`GlobalShiftModel`（`X_0 + b`）、DEG 加权/PDS 过滤指标与 BDS/DRF 校准。 | 未考察 `split-half` 选基因（`resp_top50_split`）下省略目标基线的 `E1_ConstMed` 于 `r_orth_cal` 中的残余得分。 | **`部分`** |
| **Schäfer et al. (2026, `scPertEval`)** | `doi:10.64898/2026.07.23.740433` (`bioRxiv`, `2026-07-27`) | 系统比较 `pearson_pert_top_k`、`pearson_ctrl_top_k`、`heg_k` 等空间对 DRF/BDS 排序的影响。 | 其 `top_k` 取自同侧 GT（无独立 `split-half`），未报告 `E1_ConstMed` 在 `resp_top50_split` 于 `r_orth_cal` 中保持 `~+0.26` 正值。 | **`部分`** |
| **Viñas Torné et al. (2026, `Systema`)** | `doi:10.1038/s41587-025-02777-8` (`Nat Biotechnol`) | 识别扰动与控制间的系统性共享变异并提出以 `all_perturbed_mean` 为参考的 `pearson_pert`。 | 关注全转录组去系统偏移；未发现 `r_orth_cal` 正交去均值后 `E1_ConstMed` 仍在 `resp_top50_split` 上显著回升。 | **`部分`** |
| **Huang, Yang, Yu, Chen & Li (2026, `AmbiModBench`)** | `arXiv:2609.32527v1` (`2026-09-26`) | 批评绝对指标混淆共享背景响应，采用 Top-100 / Top-2% / Top-20（`|y_t|`）方向一致率并与训练集均值作差；做基因置换实验（Fig 3c）。 | 将训练集均值作为减项参考，未审计 `E1_ConstMed` 在 `r_orth_cal` + `resp_top50_split` 下是否呈现正值。 | **`部分`** |
| **Panagopoulos, Omidi, Tong & Pang (2026)** | `arXiv:2610.26267v1` (`2026-10-19`; OpenReview withdrawn `1YbFu1Ad0l`) | 揭示 Varying top-N DEG 时 Precision 排名 `86.7%` 翻转，因重复均值向量抬高 Mann-Whitney 检验统计量（`67.1%` 基因对过阈）。 | 揭示细胞行数假重复膨胀显著基因集，与单向量 `E1_ConstMed` 在 `r_orth_cal` 下正值无关。 | **`部分`** |
| **Gao (2026, `The Delta Gap`)** | `arXiv:2611.10658v1` (`2026-11-07`) | 量化全局 Pearson 与 delta Pearson 鸿沟（GEARS `0.977` vs `0.076`），加测 top-20（`|δ_obs|`）Pearson（训练平均 delta 在 K562 达 `0.267`）。 | 未行 `r_orth_cal` 正交去共享响应，亦未测 `E1_ConstMed` 在 `split-half` 选基因下的行为。（注：该文所涉 A549 为其自身引用，与我方 R12 的 Liu et al. 2026 `GSE345058` 无关。） | **`部分`** |
| **Arc Virtual Cell Initiative & Goodarzi (2026, `cell-eval2`)** | `doi:10.1016/j.cell.2026.08.004` (`Cell`, `2026-08-26`) | 规定 DE 指标先由对照均值 `> 5 CPM` 单侧门控后行 BH 检验与方向/LFC 评分。 | 仅基于对照均值设下限门控，未研究 `E1_ConstMed` 与 `resp_top50_split`。 | **`部分`** |
| **Ahlmann-Eltze et al. (2025) / Wong et al. (2025) / Csendes et al. (2025)** | `doi:10.1038/s41592-025-02772-6`; `doi:10.1093/bioinformatics/btaf317`; `doi:10.1186/s12864-025-11600-2` | 在千余高表达基因（`heg_k`）或基准集上展示简单均值/线性基线匹敌深度学习模型。 | 使用固定 `HEG` 或标准 delta 空间，不涉及 `r_orth_cal` 下 `E1_ConstMed` 的 `resp_top50_split` 现象。 | **`部分`** |
| **Sun et al. (2026) / Kalygina et al. (2026) / ICLR 2027 Anon.** | `OpenReview:EZ2KZwM1mR`; `OpenReview:U5Gq071xcI`; `OpenReview:DVGHFaScK9` | 分别揭示数据筛选对众数基线虚高、Baseline Saturation、化学扰动免费基线三轴诊断。 | 均不针对 `E1_ConstMed` 在 `r_orth_cal` 下的 `resp_top50_split` 实测残留。 | **`部分`/`无`** |
| **David & Peng (`ReliFlow`, 2026) / Raza et al. (2026) / Sanchez et al. (2026)** | `arXiv:2611.23999v1`; `doi:10.64898/2026.11.19.767931`; `doi:10.64898/2026.10.02.756269` | 训练集内方差加权 / 分库级伪重复 / 患者队列聚类后 DE 检验。 | 与 Perturb-seq `E1_ConstMed` 无关。 | **`无`** |
| **表达式转换与选择后推断基础文献**：Booeshaghi & Pachter (2021), Ahlmann-Eltze & Huber (2023), Townes et al. (2019), Kriegeskorte et al. (2009), Neufeld et al. (2022), Song et al. (2026) | `doi:10.1093/bioinformatics/btab085`; `doi:10.1038/s41592-023-01814-1`; `doi:10.1186/s13059-019-1861-6`; `doi:10.1038/nn.2303`; `doi:10.1093/biostatistics/kxac047`; `doi:10.1177/15578666251383562` | 分别奠定 `log1p` 表达均值相关尺度压制及选择后推断之理论基础。 | `log1p` Jacobian 与基线–幅度耦合对 Claim 2 的解释至今仍是**未获证实的假说之一**（R11 真实 $\mu$ 无噪声模拟在 `8/12` 线符号错误）；上述工作亦均未考察 `E1_ConstMed` + `r_orth_cal`。 | **`部分`**（假说背景，非确立机制） |

---

## 四、Claim 1 审计：`Kawsar et al. 2026` 精准边界与 `2026-11-15` 后增量检索

### 4.1 Kawsar et al. (2026, `CanPert / CASTER`) 全源码与结果核查

- **文献标识**：Kawsar, Dip, Paul, Noman, Zhang & Moni (2026), *"Transferring CRISPRi knockdown responses to new cell lines requires screening part of each line"*, `doi:10.64898/2026.11.15.766815`（bioRxiv v1 `2026-11-17`；代码仓 `CanPert` `CITATION.cff v5.1.0`, `date-released: 2026-11-15`）。
- **已核实内容（`experiments.py`, `train.py`, `results_v1.py`, `t1_fold_means.csv`, `t2_fold_means.csv`, `H2_tests.csv`, `calibrated_metrics.csv`）**：
  1. 覆盖 4 个人类细胞系 CRISPRi 筛选（`K562, RPE1, HepG2, Jurkat`，6,546 共享基因上取 top 2,000 HVGs）。
  2. 仅在**部分筛选适配任务 `t1`（含 `t1_grouped`，目标细胞系含 25% 训练扰动、5% 验证、70% 测试）** 中测试 `CASTER_wrongctrl`（学习曲线 `lc` 中显式关闭 `wrongctrl=False`；**零样本未见细胞系任务 `t2` 完全不含 `CASTER_wrongctrl`**）。
  3. `wrong_ctrl_predictor`（`train.py`）使用邻系 `LINES[(idx+1)%4]` 对照细胞输入已微调网络 `CASTER`，求出增量 `shift = base(p, wc).mean(0) - wc.mean(0)` 后加回目标系真实对照均值上（`ctrl + shift`，仅改变模型控制输入条件化、保持目标基准表达）。
  4. 在 `t1` 中 `CASTER` (`0.2759`) 优于单邻系 `CASTER_wrongctrl` (`0.1516`)，但无任何控制输入的纯线性基线 `Linear` 反达 `0.2821`（校准表为 `0.2350` vs `CASTER 0.2105`）。
  5. **零样本未见细胞系 `t2`（`H2_tests.csv`）**：`CASTER` (`mean = 0.1709`) **显著劣于简单训练集均值 `PerturbMean` (`0.1967`)**（Holm 单侧 `p = 1.0`, median diff `-0.0142`；`de_auprc`: `0.1074` vs `0.1359`）。
- **与 Claim 1 的边界差异**：Kawsar et al. 测试的是经 25% 靶向扰动微调的神经网络控制输入条件化消融，且其零样本迁移彻底败于 `PerturbMean`；完全未研究纯未扰动基线生物变换 `FCT` 对全源错配基线集合 `Mis_s` 在 `r_orth_cal` 下的 Intersection-Union Test（`min_s (FCT - Mis_s) > 0`），也未验证跨扰动模态（我方从 CRISPRi 跨至 Liu et al. 2026 `GSE345058` `A549` CRISPRko 的 `6/6` IUT 成立，最小边际 `+0.0796 [+0.0693, +0.0900]`，`nadig_jurkat`）。
- **重叠度判定**：**`部分`**。

### 4.2 `2026-11-15` 至 `2026-11-23` 新增公开文献检索结果

| 文献 / 工作 | 精确标识符与发布日 | 核心内容（已核实） | 与 Claim 1 重叠度 |
|---|---|---|---|
| **Sharma et al. (2026)** | `doi:10.64898/2026.10.20.760537` v2 (`2026-11-18`, bioRxiv) | 用 guide abundance 代理单细胞扰动强度，指出 `STATE` / `Cell2Sentence` 在 `HCT116/HEK293T` 下检索不能区分扰动身份。 | **`无`** |
| **Wei et al. (2026, `UPM`)** | `doi:10.64898/2026.10.25.761625` v2 (`2026-11-18`, bioRxiv) | 构建 110M 细胞、1B+ 参数统一 VAE 预测器，在 Tahoe-100M 与细胞因子基准上对比 `State`/`Stack`。 | **`无`** |
| **Cupaioli & Merelli (2026, `fmaudit`)** | `doi:10.64898/2026.11.03.763968` v3 (`2026-11-19`, bioRxiv) | 测试 5 个 scFM 与 `STATE` 在 K562、类器官与鼠脑皮层，发现基因敲除不优于表达匹配对照基因删除及训练均值。 | **`无`** |
| **He et al. (2026)** | `doi:10.64898/2026.11.16.767041` v1 (`2026-11-17`, bioRxiv) | 用二维序列–功能模型条件化 IFN-β 单细胞状态预测响应，并检出比例统计量按大值项门控失真。 | **`无`** |
| **Bendidi et al. (2026, `Occam`)** | `arXiv:2611.21959v1` (`2026-11-16`) | 去噪 MLP 跨细胞系/双扰动预测，采用独立 A/B 控制组；基线含 `Transfer` 与 `Mean perturbation`。 | **`无`**（纯基线变换 IUT 方面）；**`部分`**（基准对照） |
| **David & Peng (2026, `ReliFlow`)** | `arXiv:2611.23999v1` (`2026-11-18`) | 训练集内 cross-fitting 估计单扰动可靠性用于组合扰动流匹配。 | **`无`** |
| **Tarikere (2026)** | `doi:10.21203/rs.3.rs-11768707/v1` (`2026-11-18`, Research Square) | Adamson 单数据集（82 扰动）上分析 OOD 预测器自适应选择与 Embedding Coverage 相关性。 | **`无`** |
| **Liu et al. (`PerturbBridge`) / Hutchins et al. (`CIRCE`) / Raza et al.** | `doi:10.64898/2026.09.24.754259` v2 (`11-16`); `doi:10.64898/2026.11.19.767979` (`11-20`); `doi:10.64898/2026.11.19.767931` (`11-20`) | 分别研究潜在 Schrödinger 桥传输、组合信号强度的玻尔兹曼监督分级、单细胞耐药分库噪音阈值。 | **`无`** |
| **更早补充确认**：Huang et al. (`CGC`, `2026-10-01`)；Zhang et al. (`2026-11-04`, 低秩瓶颈) | `doi:10.64898/2026.09.30.755756`; `arXiv:2611.05849v1` | 分析跨上下文特异性压缩（CGC）；证明 STATE 解码器呈 4–32 维低秩瓶颈且受训练最近邻距离限制。 | **`部分`/`无`** |

---

## 五、给主控决策的文献与写作用语规范（无对外文字）

1. **严格限定 `E1_ConstMed` 定义**：必须表述为「**含源细胞系平均倍数变化 $\overline{FC}_{p,g}$、但不使用目标细胞系自身对照基线谱的迁移预测器**」，严禁再写成「与扰动靶点无关」。
2. **严格限定 Claim 2 新颖性与机制边界**：
   - Nicol et al. (2026) v2（代码与 Fig 2b，因 bioRxiv v2 正文 `403` 故仅可据源码和图 2b 表述）已经提出 `Split DEG + control` 且其中有信息基线 `nonctl-mean` 在多个数据集 top-20 上仍保持高正分（与全基因 `corr_split` 相近，显著高于坍缩至 ~0 的零信息 `ctl-mean`）；
   - Heidari et al. (2026)（bioRxiv 全文 `403 未能读取`，源码可见）已经涵盖 `CrossSplit`、`*_vs_topdeg`、`trivial` vs `non_trivial` 分层、`rigorous_pearson_delta_pert` 与 `baseline_logfc`；
   - 因此，**Claim 2 应诚实降级为特定条件下的次要补充发现**：仅限于在去除共享均值的 `r_orth_cal` 空间、单数据集 Liu et al. (2026) `GSE345058`（`A549` CRISPRko，非 "Norman A549"）上观察到省略目标基线的 `E1_ConstMed` 在全基因归零（`-0.0201`）后仍可在 `resp_top50_split`（`+0.2565`）与 `resp_top50_source`（`+0.2766`）呈现显著正相关，且**机制至今未知**（R11 真实 $\mu$ 无噪声模拟 `8/12` 线符号错误；`log1p` Jacobian 与丰度耦合只作为未证实假说之一提及）。
3. **Claim 1 继续保持主发现地位**：Kawsar et al. (2026) 与 Heidari et al. (2026, `baseline_logfc`) 均仅为 **`部分`** 重叠，截至 `2026-11-23` 无文献提出并实证通过基于未扰动对照基线的 `FCT` 对全源错配基线 `Mis_s` 的跨系与跨模态（`CRISPRi -> Liu et al. 2026 GSE345058 A549 CRISPRko`）`6/6` IUT 检验（`r_orth_cal`）。
