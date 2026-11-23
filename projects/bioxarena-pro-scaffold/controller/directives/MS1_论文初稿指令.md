# 论文初稿指令（主控 Claude，2026-11-23 10:30 +08:00）

用户已同意：**现在按单线确认的结果写完整论文初稿。** 执行人：A-T4。

## 文件和发布规则
- 写在 `paper/MANUSCRIPT_DRAFT_v1.md`（英文，期刊论文格式，Markdown，公式用 LaTeX）。
- **不提交 git，不推送，不放进任何公开位置。** 仓库是公开的；初稿要先经用户审核。main 代理的定时同步也必须跳过 `paper/`（主控另行通知 main）。
- 不写新闻稿、推文或任何对外文字。人类角色只称 "Human project lead Zijian Zeng"，不用 PI。

## 依据（每个数字都必须能在这些文件里找到，写在 Markdown 注释 `<!-- src: 文件 -->` 里）
- `round12/status/R12_VERDICT.md`（确认性核心）、`round12/results/*holdout_r12*`、`round12/prereg/PREREG_E.md`、`round12/prereg/DEVIATIONS_R12.md`、两份裁定。
- `round13/status/R13_VERDICT.md`、`round13/prereg/PREREG_F.md`、`round13/status/RULING_R13_READER.md`。
- `round13/status/PRIOR_ART_CONCLUSION_R13.md`（含补充）、`round13/shared/PRIOR_ART_G.md` v2。
- `round10/paper/OUTLINE.md`（v3，结构和措辞约束）以及第 9–11 轮的判定文件。

## 结构
1. Title（用大纲候选 1 或 2，不用含 "establish" 的候选 3）。
2. Abstract（≤ 250 词）：必须写明只有一条确认性细胞系（K_conf = 1，推广性有限）；两个论断的假设来自第 9–11 轮的事后观察；论断 1 是特异性不是精度（Source_Grand_Mean 0.405 > FCT 0.190）；论断 2 是方法学观察；第十三轮第二条线门槛失败；第十四轮独立 A549 数据集同样门槛失败（32 < 50）。
3. Introduction：零样本扰动预测的评估问题；正确与错配对照的已有工作（必须引用 Kawsar et al. 2026 CASTER，不声称首次）；选基因偏倚的已有工作（Nicol et al. 2026、Heidari et al. 2026、Miller et al. 2026、Systema）。
4. Results：
   - 4.1 第 9–11 轮：怎样得到这两个假设（简短；写明是事后观察）。
   - 4.2 第十二轮确认（A549，N = 445）：全部 H-E1a/b/c、H-E2a/b 数字和 97.5% CI；SA-Cas9、SA-Src8；HEK293 只作描述；三方复核一致性。
   - 4.3 FCT 不是最好的预测器：同一指标下的全部预测器表。
   - 4.4 第十三轮：第二条细胞系的事先登记尝试在门槛失败（3 < 50）及原因；N = 3/5 的数字不写成复现或未复现。
   - 4.5 第十四轮（10:30 已知门槛结果）：独立实验室的 A549 数据集 GSE337804（PREREG_G）在门槛失败：评测靶点 32 < 50，K_conf = 0，不做确认性判定（`round14/status/CONFIRMATORY_LINES_R14.txt`、`round14/shared/R14_QC_COUNTS.md`）。主控观察到的原因（描述性）：Cas9 持续表达、第 16 天取样，必需基因敲除细胞大量丢失（例如 CMTR1、EIF4G1、GMPPB、RIOK2、SPCS3 每个只有 11–33 个细胞），而源线是必需基因筛选。描述性指标「[R14 DESCRIPTIVE PENDING]」，等 R14_VERDICT.md 后填写；同样不能写成复现或未复现。
5. Discussion：论断 1 的含义和限度；论断 2 作为评估建议（检验靶线特异性要用全部非靶基因）；单线局限；对 VCC 类评测的意义；未来需要什么数据（覆盖必需基因的大文库、新细胞系）。
6. Limitations（单独一节）：K_conf = 1；同一扰动方式外推；偏离 D1–D4 与 HEK293 曾被列为确认性线后退出（SUPERSEDED_v1）；r_orth_cal 指标的范围；论断 2 机制未知；两篇先前工作的全文尚未人工读过（作为内部提醒，用 `<!-- TODO -->`，不写进正文）。
7. Methods：数据集（Liu 2026 GSE345058、Xu 2023、Jung 2026 GSE327644、第 9 轮源线）；预处理与 QC；FCT、错配基线、E1 的精确定义；r_orth_cal；掩码；Addendum-B 靶点准入；IUT 加 Bonferroni、bootstrap、置换；预注册和冻结流程；三方独立复核；双 AI 分工。
8. Data and code availability（预注册提交哈希、冻结提交哈希、数据登录号）。
9. AI involvement and author contributions disclosure（主控 Claude；执行代理 Antigravity 的 A-T1–A-T5 分工；人类项目负责人）。
10. References：**只列已经核实过的文献**（DOI 或 PMID 必须来自依据文件）；任何没有核实的条目标 `[VERIFY]`，不要编造。

## 写作要求
- 不夸大：不用 "establish"、"prove"、"first"；不把论断 2 写成主要发现；不暗示能推广到其他细胞系。
- 不要捏造数字或文献；拿不准就写 `[CHECK]`。
- 图表：先写图注和表格（数据来自结果文件），图片之后再做。
- 完成后在 `round14/status/A-T4.txt` 最后单独一行写 `MS1_DRAFT_DONE`。
