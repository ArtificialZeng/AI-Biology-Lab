# AI Biology Lab：本轮结果与可核查材料

更新：2026年10月26日。AI辅助推导、编码和审查；没有外部人类同行评审或新的生物实验，优先权尚未确认。

## 主要数学结果

对Conradi–Mincheva–Shiu2019年的混合机制双磷酸化模型，固定原文Table 1速率、P=8、S=40，本轮给出严格正吸引周期轨道的区间证书，并证明其对所有充分大有限激酶总量K持续存在于对应守恒类内。

在模型原有时间单位下，局部分支满足T(K)=T∞+T1/K+O(K⁻²)，其中T∞∈[4.69209598,4.69209603]，T1∈[7.65771717,7.68396573]。非零系数说明周期趋稳仍有有限浓度偏差；没有算出明确K0或余项常数。

[完整中文研究简报](../projects/mixed-phosphorylation/validated-cycle/docs/BRIEF.zh-CN.md) · [证明和复现](../projects/mixed-phosphorylation/validated-cycle/README.md) · [可用于展示的图](../projects/mixed-phosphorylation/validated-cycle/figures/period_correction.png)。图中灰色有限K点是浮点模拟，蓝色带是严格认证的渐近系数范围，不是有限K误差带。

原论文已展示数值振荡，2023/2025会议工作直接覆盖质性鲁棒性。本轮尚未在所查资料中找到这个固定实例的严格非零周期系数证书，但详细后续资料不完整，不能宣称首次或完整破解猜想。周期变化也不等于反驳标量集合投影的所有读法。

## 两个补充审计结果

- [Protein Metric Auditor](../projects/protein-metric-auditor/docs/BRIEF.zh-CN.md)：完整重放217实验、95模型配置、约246万条记录，量化NDCG副指标的行顺序敏感性。最多两名汇总位移、前五不变，未推翻Spearman主榜；并列处理已在已有文献说明。
- [MEMOTE方向审计](../projects/memote-direction-auditor/README.md)：对质量和电荷守恒的最小网络，仅将反应改成数学等价写法就可令原版能量检查产生假阳性或假阴性；真实完整ATP入口和候选两行修复均已复跑。没有证明真实Yeast-GEM受影响，也未获维护者确认。

可用报道角度：个人AI实验室公开了一套从提出问题到严格核验的研究过程，并为一个生化振荡模型给出了可重跑的周期偏差定律。对外报道前，仍应由领域专家核对新颖性。本文不是已发送的新闻稿，没有联系媒体。
