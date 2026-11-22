# 主控审稿：第九轮预印本初稿 R1

审稿人：Claude（主控） ｜ 2026-11-22 ｜ 对象：`paper/PREPRINT_DRAFT.md`（sha256 aebcae7f…7a93）

**总体：** 结构完整，「未确认」写得清楚，探索性部分隔离得当，局限诚实。下面 B1–B4 必须改，S1–S7 也请照改。

## 阻断项
**B1 预注册文字和冻结代码不一致，必须披露。**
- PREREG_A v4 §4.2 写的是 P⊥ = I − q1q1ᵀ，q1 是源线扰动向量的第一主成分（秩 1）。
- 冻结代码 1ca9f60 的 `src/predictors.py::fit_shared_subspace` 实际去掉的是 3 维子空间 span{1, 源线总均值, 残差 PC1}。初稿 2.3 节第 4 条写的是代码版本，却没有指出它和预注册文字不一致。
- 要求：
  - A-T2 在 `vcc_linespecific/DEVIATIONS.md` 记录这一条。
  - A-T2 按预注册文字（秩 1 的 q1）补跑敏感性分析 S2b，A-T3 独立复算。
  - 论文方法部分新增「与预注册的偏差」小节，写明：这个差异在锁定前就存在，是审稿时发现的；判定以冻结代码为准；给出 S2b 的结果。
  - H-A1（斜率 > 0.13）在任何投影下都不通过，所以 A 组的总判定不会变。但必须照写。
- 同时说明：A-T3 的「独立实现」复刻了同一个 3 维投影，两套实现在这一点上并不独立。写进局限。

**B2 作者分工写错了。** FCT 的探索性特异性检验是 A-T2 做的（见 `status/A-T2.txt` 和 `results/EXPLORATORY_FCT.md`），不是 A-T4。第 4 节第一段和第 7 节都要改。

**B3 措辞过强。**
- 标题里的 "Rejects"、正文的 "rejected every hypothesis"、"confirmatory rejection" 都改成 "did not confirm" 或 "NOT CONFIRMED"。预注册的阈值没达到，不等于假设被证伪。
- "multi-laboratory evaluation" 删掉（数据来自不同实验室，但评测不是多个实验室做的）。
- "unblinded … cohort" 改为 "held-out cohort"。
- 摘要末句 "establish … baselines required before …" 改为 "suggest that amplitude orthogonalization and matched-control nulls should be applied before …"。

**B4 摘要超长。** 现在约 450 词，压到 250 词以内。必须保留：
- 两组判定都是 NOT CONFIRMED，分别失败在哪几条；
- GlobalScale 说明（单个标量的原始相关范围，去幅度后为 0）；
- R14 的措辞；
- 确认性线只有 3 条，且 H1 只有 80 个靶点；
- 三方独立复算一致。

数字放主表，不要全堆进摘要。

## 建议修改
- **S1 独立性的局限。** 加一条局限：A-T2 和 A-T3 都是同一家族模型（Antigravity / Gemini）的智能体。主控 Claude 用另写的代码独立复算了 B 组 tmvs 在 3 条确认性线上的 r_orth_cal（HepG2 +0.0419，Jurkat +0.0947，H1 −0.0772，一致），但没有复算 A 组。
- **S2 机制推测要标明是推测。** 5.1 节第 1 条 "Mechanistically, multiplying … captures shared perturbation-variance alignment …" 是推测，改成 "One possible explanation is …"。
- **S3 选择偏倚。** 5.1 节第 4 条引用 resp_top50 的 L2 = 14.1%–18.9% 时，要注明 resp_top50 是按观测 |Δ| 选的基因，有选择偏倚（S6），并优先引用 all_nontarget 的数字。
- **S4 先验性措辞。** "no prior perturbation study extended … or demonstrated …" 改为 "in our bounded literature search (PRIOR_ART_A–C), we found no study that …"。
- **S5 第 4 节末句语病。** "remains untested pre-registeredly" 改为 "has not been tested in a pre-registered design"。
- **S6 披露角色。** 第 5 行改为：执行智能体是 Antigravity（Google Gemini 模型），主控是 Claude（Anthropic），负责审核预注册、核验结果和独立复算；人类负责人 Zijian Zeng 没有生物学背景。
- **S7 数据可得性。** 第 7 节加上仓库链接 https://github.com/ArtificialZeng/AI-Biology-Lab 和各数据集的原始来源链接。

改完后更新 `status/PAPER_DRAFT_READY.txt`（新的 sha256），写一行 `R1_REVISED`。暂不生成 PDF。
