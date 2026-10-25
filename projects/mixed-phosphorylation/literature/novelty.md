# 文献与新颖性审查：混合机制磷酸化振荡的高激酶极限

审查日期：2026-10-26（Asia/Shanghai）。审查者：本任务的独立文献子任务。

## 本次结论与可用范围

在本次检索并实际检查的材料中，未找到一篇可核验的后续完整论文，已经给出 Conradi–Mincheva–Shiu (2019) §6、Conjecture 6.2 的证明、反例，或该固定 Table 1 模型高激酶总量极限下周期/振幅的首阶校正。然而，2023 年已有直接研究同一鲁棒性问题的官方海报题目，2025 年又有同一研究者与原作者的直接相关会议摘要。因此当前的新颖性状态必须写为 **status-uncertain：存在明确研究重叠，尚未取得其结果全文**。不能写“无人研究”“首次研究浓度鲁棒性”或“已确认尚未解决”。

原作者 2019 年已发现振荡及其对较大激酶总量的近似不敏感；复现这些图不是新发现。候选增量应限于模型特有的严格高浓度极限、可验证的有限浓度误差/校正、量化失效边界，或具有明确量词的反例。普通参数扫描、通用周期轨道灵敏度公式、在未经认证的双曲周期轨道假设下援引一般持久性定理，本身不足以构成已经确认的新结果。

本文件是有边界的文献筛查，既非数学裁决，也非完整优先权清查。没有联系作者、会议组织者或任何外部人员。

## 1. 冻结源及语义边界

### 原始论文

- C. Conradi, M. Mincheva, A. Shiu, *Emergence of Oscillations in a Mixed-Mechanism Phosphorylation System*, Bulletin of Mathematical Biology 81 (2019), 1829–1852, DOI [10.1007/s11538-019-00580-6](https://doi.org/10.1007/s11538-019-00580-6).
- [arXiv:1809.02886v2](https://arxiv.org/abs/1809.02886v2)，2019-02-01；[HTML](https://arxiv.org/html/1809.02886v2)。HTML 全文可读，实际检查 §6、Conjectures 6.1/6.2、Figure 5 与其后的灵敏度建议。
- [正式出版版本的 NSF 存档 PDF](https://par.nsf.gov/servlets/purl/10090542)，可检索到相同猜想及其上下文。作者个人站点的 [mixed.pdf](https://people.tamu.edu/~annejls/papers/mixed.pdf) 是较早版本，不以其日期替代最新 arXiv 版本。

**已知基线**：Table 1 固定速率；Figure 5(c) 的 P=8、S=40，K=100、1000、10000 是数值例子。§6 的猜想仍保留对 P、S 的存在量词，不能用这个单一例子的有限 K 差异冒充驳倒整个存在命题。

**关键措辞**：正文使用“when projected to x5”；并建议分别研究周期与振幅。仅就字面而言，下列三个对象不等价：

1. 单坐标集合投影：闭轨道 Γ_K 投影为区间 [min x5,max x5]。
2. 周期和振幅：T(K)、A(K)，可能相同而波形不同；也可能投影区间相同而周期不同。
3. 物理时间的输出波形：存在相位平移 φ_K，使 x5,K(t+φ_K) 对所有 K 完全相同。该说法通常同时约束最小周期、极值与时间分布。

本次对原文的解释是：生物计时语境和紧接的灵敏度建议支持作者关注 2/3，但形式陈述没有给等价关系，因此不能静默把 1 强化为 3。若计算只发现周期改变，尚未反驳最弱的集合读法；若发现经过可靠误差控制的极值改变，也仅针对所考察参数和轨道。应在 problem.md/后续成果中分别声明解释层级。

## 2. 直接重叠的后续公开记录

### 2023 SBHD 海报：同一激酶总量鲁棒性题目

- [SBHD2023 official poster session programme](https://www.sbhdberlin.org/pdf/SBHD2023_Poster_Session_Program.pdf)，官方 PDF 3 页；实际打开并用 pdftotext 复核第 3 页 / printed p.52，Poster ID 48。
- 题目：*Analysis of the robustness of oscillations in a mixed mechanism double phosphorylation system with respect to total kinase concentrations*。
- 作者/报告者：Kim-Joao Grasse。
- **能证明的事实**：2023 年已公开报告直接相关的研究方向。
- **不能证明的事实**：该材料只有题目与姓名，未给定理、模型版本、方程、周期渐近式、证明或数据。不能据此断言猜想已解决，也不能忽略它宣称方向首次。
- 补充检索了完整题目、作者姓名、kinase/robustness/thesis/dissertation 等组合，未取得海报正文或对应学位论文。未联系作者。

### 2025 BIOMATH 摘要：原作者参与的相关工作

- K.-J. Grasse, C. Conradi, *On the robustness of oscillations in a mixed mechanism phosphorylation system against perturbations in the total amount of the processive enzyme*，BIOMATH 2025，June 15–20, 2025, Sofia。
- [官方单页摘要 PDF](https://biomath.math.bas.bg/biomath/index.php/bmcs/article/download/1787/1594/7828)；[官方摘要集，printed p.53](https://www.biomath.bg/2025/BIOMATH2025_Book_of_Abstracts.pdf)。实际读完整单页。
- 摘要描述其工作对象为反向的混合机制：distributive phosphorylation / processive dephosphorylation，考察 processive phosphatase 总量的扰动。
- 摘要将 2019 年关于激酶鲁棒性的工作概括为已经证明，但没有提供新证明，只引用 2019 年论文。该概括不足以把原文猜想升级为定理；也可能只是以“robustness”泛指已观察的持续振荡。
- 摘要没有给出模型的高浓度极限系统、有限浓度周期/极值校正，也没有给其研究方向的最终定理。
- [BIOMATH2025 official detailed programme](https://biomath.bg/2025/BIOMATH2025_DetailedProgram_Final.pdf)确认该报告列在 June 20、11:00–11:30。日期来自节目而非搜索抓取时间。

### 2025 ASIM 会议项目

- [ASIM Göttingen 2025 official programme](https://www.asim-gi.org/fileadmin/user_upload_asim/ASIM_OA/FG_SUG/WS_G%C3%B6ttingen_2025_Program_final.pdf)，March 26–28, 2025。
- Grasse 与 Conradi 的题目为 *Differences in robustness in cyclic distributive and mixed-mechanism phosphorylation networks*。
- 访问：官方 PDF 的搜索索引文字可读，未获得讲稿/结果全文。仅用作“已有相关公开工作”的证据，不用作具体数学结论来源。

## 3. 最近的完整 primary papers：哪些结果不能混用

| 来源与实际访问 | 已建立/研究对象 | 对本题的意义与局限 |
|---|---|---|
| Suwanmajo & Krishnan 2015, [DOI 10.1098/rsif.2014.1405](https://doi.org/10.1098/rsif.2014.1405)，[PMC全文](https://pmc.ncbi.nlm.nih.gov/articles/PMC4590494/)。实际检查 §3.1.1、§A.2/A.4。 | 混合机制中的振荡、酶截留机制、稳态剂量响应、改变额外反应的数值影响。 | 已有关于高 K 的稳态饱和值和自由中间态比例的分析；这不能直接替代周期轨道及其有限 K 校正。所谓 asymptotic amplitude 在稳态剂量响应上下文出现，不能自动解释为振荡峰谷差。 |
| Tung 2018, [arXiv:1712.03594](https://arxiv.org/abs/1712.03594)，[DOI 10.1016/j.mbs.2018.10.008](https://doi.org/10.1016/j.mbs.2018.10.008)。访问期刊摘要/原始预印本记录。 | 通过Dulac类判据排除双位点系统 Michaelis–Menten 近似的振荡，包括混合机制。 | 小酶/底物过量的标准 MM 约化不等于固定 P,S、K→∞。不得以“MM无振荡”否定高 K 全质量作用律系统的振荡，也不得把两个极限混写。 |
| Rao 2018, [DOI 10.1016/j.mbs.2018.05.001](https://doi.org/10.1016/j.mbs.2018.05.001)。期刊摘要与假设可读。 | 混合机制的 MM 近似唯一平衡；双位点情形全局稳定性。 | 模型方向与缩放假设均须逐项比较；不能直接作为当前高 K 极限定理。 |
| Rao & Gasparyan 2024, [官方作者机构记录](https://biblio.ugent.be/publication/01J8HN58SQB6S7R6T6JKTNTVGJ)，[accepted manuscript](https://backoffice.biblio.ugent.be/download/01J8HN58SQB6S7R6T6JKTNTVGJ/01J8HNKD93491CGN7GK96DWW73)，[DOI 10.1016/j.automatica.2024.111837](https://doi.org/10.1016/j.automatica.2024.111837)。实际读摘要及引言的 QSSA 假设说明。 | 一类封闭 MM 网络的唯一性及全局渐近稳定性，涵盖过程式/混合机制的约化网络。 | 进一步确认“MM稳定”是另一模型/另一极限；没有在检查部分提供当前固定率高 K 周期修正。未逐页审计该文全部证明。 |
| Conradi & Mincheva 2024, [DOI 10.1007/s00285-024-02114-8](https://doi.org/10.1007/s00285-024-02114-8)。全文HTML可读；实际检查引言、研究网络、§4.2/§4.3及参考文献。 | 循环分配式双磷酸化中的Hopf机制；与顺序分配式双稳态的比较。 | 研究网络与本项目不同。引用2019作为混合机制振荡基线，没有在检查内容中解决Conj6.2的高 K 鲁棒性。 |
| Conradi, Obatake, Shiu, Tang, *Dynamics of ERK regulation in the processive limit*, J. Math. Biol. 82, 32 (2021)，[作者PDF](https://people.tamu.edu/~annejls/papers/ERK-2.pdf)。实际检查摘要、§6.1/§7。 | 研究混合ERK模型随着“过程性概率”接近1时的多稳态/Hopf。 | 此文也有“Conjecture6.2”，但讨论最小双稳态ERK子网的正平衡点最多3个，**不是**本项目2019论文的同编号猜想。检索编号不能替代比较命题。过程性极限也不是K总量极限。 |
| Dajani & Mincheva 2025, [官方PDF](https://biomath.math.bas.bg/biomath/index.php/biomath/article/download/j.biomath.2025.03.285/pdf/7751)，[DOI 10.55630/j.biomath.2025.03.285](https://doi.org/10.55630/j.biomath.2025.03.285)。实际检查摘要、引言与研究网络。 | 三位点混合机制的多稳态条件。 | 不是当前双位点周期轨道高 K 分析。 |

## 4. 灵敏度方法已有，可能的新结果必须更具体

- Ingalls 2004, *Autonomously oscillating biochemical systems: parametric sensitivity of extrema and period*, DOI [10.1049/sb:20045005](https://doi.org/10.1049/sb:20045005)。本次核对书目信息，未获取并逐页审计其全文；它是2019原文直接建议的方法来源之一。
- Ingalls, Mincheva, Roussel 2017, *Parametric Sensitivity Analysis of Oscillatory Delay Systems with an Application to Gene Regulation*, DOI [10.1007/s11538-017-0298-x](https://doi.org/10.1007/s11538-017-0298-x)，[PubMed官方记录](https://pubmed.ncbi.nlm.nih.gov/28608044/)。实际读摘要：区分周期/极值灵敏度与会受相位漂移影响的直接轨道灵敏度。
- 标准的一阶周期变化公式、相位规范、常微分方程光滑参数依赖与双曲周期轨道持久性不是本项目新方法。模型特定的非零系数、消去后的显式简化、可认证的误差阶与参数范围，才可能产生可辨识增量。

## 5. 推荐在成果中采用的声明

可用：

> 我们研究2019年提出的浓度鲁棒性问题。后续存在直接相关会议研究记录；在本次查阅的公开完整材料中，尚未定位与本文具体定理等价的证明。本文的贡献范围是[填入经过审核的具体命题]，并不把一般鲁棒性方向主张为首次。

禁止在当前证据下使用：

- 已彻底解决Conjecture6.2（除非证明其完整量词与解释）。
- 第一个研究混合磷酸化高激酶鲁棒性的工作。
- 只因图形略有不同或周期变化，就宣布2019原猜想错误。
- 已发现真实细胞内的新振荡机制（本项目是给定模型的数学性质）。
- 搜索未返回结果，所以证明了新颖性。

## 6. 检索覆盖与未覆盖事项

已检索：论文全名、DOI、Conjecture6.2与kinase/oscillations组合；mixed-mechanism/mixed mechanism与large kinase、arbitrarily large、asymptotic、singular perturbation、Fenichel、period、amplitude、robustness的组合；Grasse、Conradi、Mincheva、Shiu的相关作者/题目；2019–2026后续线索。阅读了Shiu个人论文目录及相关原作者论文。搜索引擎返回的日期标签未用作正式发表日期。

未完成：不是对所有引用文献的穷尽citation-graph爬取；未取得2023海报正文、2025演讲讲稿、可能存在但未索引的学位论文/手稿；未通过作者沟通确认未发表工作；未逐页重验所有邻近论文证明；尚未把本项目最终定理与每一篇相关工作的完整数学表述逐项比对（最终定理尚未冻结）。

下一次最有价值的文献测试：拿到最终候选定理后，针对其精确变量缩放、误差阶、周期/振幅系数、参数范围重新检索与比较；若后续能合法获得Grasse相关报告正文，应优先审查。当前这些缺口应限制优先权措辞，不妨碍继续独立研究与验证。
