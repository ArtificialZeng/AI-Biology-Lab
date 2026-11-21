# A-T5 文献、官方指标定义与竞争监测（写 p2/vcc_linespecific/audit/LIT/）
1. 找到 Arc VCC 2026 结项报告中「线特异斜率」的**精确定义**，包括 Alexander Dobin 的评测指标预印本和 cell-eval2 文档，写成可执行的公式说明 `LINE_SPECIFIC_METRIC.md`。A-T2 和 A-T3 都要用，**最先交付**。
2. 读 CIPHER v2 方法部分，给出预测公式、单参数拟合方式和数据预处理，写 `CIPHER_METHOD.md`（引用到具体章节）。
3. 读 VCC 2026 前 3 名方案（BioInsight/Illumina CellularAI、FlowState、Xaira X-team）的公开说明，以及 Kawsar CASTER、Huang CGC、HyperMap，判断有没有人已经用目标线共波动预测线特异部分。写 `PRIOR_ART_A.md`。
4. 之后每 2 小时检索一次 bioRxiv/arXiv 新文献（共波动、线特异、cross-cell-line、VCC 2026），有撞车风险立即写入 `p2/status/A-T5.txt` 并标「ALERT」。
