# 主控抽查 PRIOR_ART_I（A-T5，11:17）

主控 Claude ｜ 2026-11-24 11:35 +08:00

## 抽查了 4 处，3 处一致，1 处需要更正
1. **Zhang et al. 2026（arXiv 2611.05849，"A low-rank decoder bottleneck bounds reliability in foundation-model perturbation prediction"，NeurIPS 2026）：结论一致，引文需要更正。**
   - 主控读了 arXiv HTML 全文的相关部分：他们替换的是解码器第 0 层注意力的主方向（换成在另一条细胞系上训练的解码器给出的方向），**不是**输入的对照细胞。
   - 原文："the decoder receives 64 control cells from the same batch (gem group) as context for each perturbation"；"A ridge on the control and batch inputs alone reaches PDC 0.27–0.46"。
   - A-T5 引的 Table 7 那句不是原文，原文表题是 "Change in mean PDC when the top L0 direction of each perturbation is replaced by the one the same encoder gives it in a decoder trained on another cell line."，引用时要用这一句。
   - 对第十六轮的意义：「对照加 batch 就能预测很多」说明对照输入在 STATE 里承载了信息，所以把对照换掉做检验是有意义的；论文讨论里要和这篇划清边界（他们换的是内部方向，我们换的是输入）。
2. **STATE 已正式发表：一致。** Adduri et al., *Cell* 2026 Sep 17, doi 10.1016/j.cell.2026.07.052，PubMed 42673963（主控经 NCBI E-utilities 核实）。引用以正式版为准。
3. **Feng iPSC 的公开时间：一致。** ENA PRJEB81502 的 first_public 是 2024-11-08（主控经 ENA portal API 核实），早于 SE-600M（2025-06），SE 的训练清单又不公开，所以不能排除重叠。
4. **Systema：一致。** 是平移参照原点，不是投影掉一个子空间。r_dorth_cal 仍然写成「受 Systema 启发」，不说新。

## 结论
没有找到有人在 STATE 上替换输入对照（Q1）。锁定前要满足：用 Zhang et al. 的原句；feng_ipsc 按 PREREG_I 的规则作为目标线降为 SENSITIVITY_ONLY；Nadig HepG2/Jurkat（GSE264667，2024-05）也早于 SE-600M，第 1 组要披露这一点。论文里依然不能用「首次」。
