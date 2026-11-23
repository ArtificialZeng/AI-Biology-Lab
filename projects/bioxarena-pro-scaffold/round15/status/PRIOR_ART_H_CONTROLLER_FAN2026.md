# 主控补充查新：Fan et al., NeurIPS 2026（OpenReview lnpzb1HvwX）

主控 Claude ｜ 2026-11-23 16:53 +08:00 ｜ 来源：OpenReview 论文页（人类项目负责人在自己的浏览器里通过了人机验证，主控随后在同一浏览器读取了论文页、摘要、元审稿和全部审稿意见；没有读 PDF 全文）

## 核实的书目信息
- 题目：*Zero-Shot Cross-Cell-Line Perturbation Prediction Needs Controls*
- 作者：Xingyu Fan, Jinghao Wang, Chang-Yu Hsieh, Chunbin Gu, Pheng-Ann Heng（A-T5 写的 "Fan et al." 正确）
- NeurIPS 2026 Evaluations and Datasets Track，poster，接收；许可 CC BY 4.0；代码 github.com/xyfan22/crossline-controls
- 回复中作者提出修订后的题目 "Zero-Shot Cell-Line Transfer for Perturbation Prediction Needs Controls"（引用时以最终版为准）。

## 他们做了什么
- 严格零样本跨细胞系 CRISPRi 压力测试：在 K562、RPE1（或再加 HepG2）上训练，在留出的 Jurkat（Nadig 2025）上评测；目标线训练时只有对照细胞。577 个共有扰动（都是必需基因）。
- 两个对照基线：**target-control-copy**（直接输出目标线对照，等于预测无效应）和 **shifted-control**（目标线对照加上源线得到的平均位移）。
- 结论：整体分布指标（Energy、W2 等）主要由目标线细胞状态的变异决定，而不是扰动特异的反应；copy 在分布指标上有竞争力；shifted-control 在扰动特异 Pearson 上超过 CellFlow。补充了 PBMC 细胞因子的留一供体实验。

## 和我们两个论断的关系
- **论断 1（FCT 在去幅度的 r_orth_cal 上胜过 6 个错配细胞系对照）：部分重叠，必须引用。**
  - 重叠：他们已经指出「目标线自己的对照」在跨细胞系评测里起很大作用，并把「目标线对照 + 源线平均位移」作为必备对照基线。我们的 Transfer_Mean / FCT 和他们的 shifted-control 是同一类做法（他们是加性位移，FCT 是乘性倍数变化）。
  - 不重叠（据摘要、审稿和回复，没有读 PDF）：他们没有把目标线对照换成**错配细胞系的对照**做配对比较；没有去掉每个扰动沿共享反应方向的幅度（r_orth_cal）；没有噪声上限校准；没有在从未用过的细胞系上做预注册的前瞻检验。他们的目标线 Jurkat 是我们的源线之一。
  - 审稿人 1RE8 指出的局限（只用必需基因、跨细胞系差异小）同样适用于我们的源线设计，论文讨论里要正面回应。
- **论断 2（分半选基因后，不用目标线基线的 E1 仍然高）：没有重叠。**

## 对论文和 PREREG_H 的要求
1. 论文引言和讨论必须引用 Fan et al. 2026，写清他们已经做到哪一步；论断 1 的表述维持「预注册前瞻检验 + 配对错配对照 + 去幅度指标」这个组合，不说首次。
2. PREREG_H §10 的披露要从「未读」改为「已读摘要、审稿和回复，部分重叠，见本文件」。
3. 结合 CASTER（wrong-line controls，四条线）和 Fan et al.（control 基线），论断 1 的空间越来越窄。**能把它提到更高水平的，是在多条从未用过的细胞系上、数据公开之前就锁定的前瞻确认，这正是 PREREG_H 要做的。**
