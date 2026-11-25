# 主控审核 R1：PROPOSAL_R16 与 PREREG_I v1-DRAFT

主控 Claude ｜ 2026-11-24 10:51 +08:00 ｜ 对象：`round15/status/PROPOSAL_R16_FOUNDATION_AUDIT.md`（10:41）、`round16/prereg/PREREG_I.md` v1-DRAFT（10:47，主线 Antigravity 会话 19f7514f）

## 0. 授权范围（主控核实过）
- 人类项目负责人 10:43 在 Antigravity 界面里给主线的原话：「批准按 PROPOSAL_R16_FOUNDATION_AUDIT.md 先行起草 PREREG_I 与测试脚手架」。**批准的只是起草和用合成数据的测试脚手架。**
- **尚未批准：** 下载任何模型权重、安装 `arc-state`、克隆 CellularAI、在任何真实细胞上运行任何模型。STATE 许可证写明「下载、复制、执行即视为接受许可」，所以下载这件事要由人类决定；主控会单独去问。
- 方向本身值得做：Arc 2026-11-19 的结项文写到，参赛模型预测出的细胞系特异成分很小（斜率 0.02–0.13，重复实验是 1.02），并说 "current models recover roughly what a cross-line average would predict"。用错配对照替换去检验开放权重模型是否真的用到目标线对照，是对这句话的直接因果检验。主控目前没有找到有人对 STATE 做过这个替换（检索范围见 §2），但查新还没有做完。

## 1. 必须改的地方（v2 之前）

**I1 先查新，后锁定。** 草稿没有查新记录。A-T5 写 `round16/status/PRIOR_ART_I.md`，主控抽查原文之后才能锁定。在那之前不写「首次」「first」，提案里的 "first pre-registered quantitative proof" 要删掉。

**I2 主指标换回 r_orth_cal。** r_dorth_cal 是在看到 a549_ko 上 Source_Grand_Mean 0.405 > FCT 0.190 之后才设计的，而且和 Systema（Viñas Torné et al., Nat Biotechnol 2025, doi 10.1038/s41587-025-02777-8，用全部扰动细胞的质心作参照来分离扰动特异效应）很接近。
- H-I1、H-I2 的主指标用 r_orth_cal，和 PREREG_E–H 保持一致。
- r_dorth_cal 作为事先登记的次要指标，写成「按 Systema 的思路去掉共享扰动轴」，引用 Systema，不当作新指标，也不叫定理。

**I3 「定理 1」按构造就成立。** 投影掉某个子空间，落在这个子空间里的预测器得分必然是 0。所以它改叫「性质（按构造成立）」；H-I3a 是单元测试，不是假设，挪到 `test_unit`。

**I4 「定理 2」只是假设。** 草稿没有给出证明，只有合成数据。第十一轮的「Theorem 1」就是因为在 8/12 条线上符号错了才被撤掉的。
- 改名为「机制假设 M2」，H-I4 降为次要分析，只描述，不判定。
- 合成模拟只能算支持性证据。

**I5 主要论断只留两个，主要模型只留一个。** 主要模型是 STATE_SE（模型卡写明零样本下 X_state 更好）；STATE_HVG 是次要模型。
- 共同主要论断：H-I1（匹配对照与错配对照比较）、H-I2（模型与 FCT 比较）。
- 两个论断做 Bonferroni 校正，用 97.5% CI，和以前几轮一样。

**I6 H-I1 必须能支持「模型没有用到目标线对照」这个结论。** 草稿只检验了「匹配 > 错配」。如果 STATE 不通过，那只是「没有证据」，不能说明「没有用到」。
- 事先写死一个等价界 δ，用 TOST：配对差的 95% CI（两论断 Bonferroni 校正后）完全落在 (−δ, δ) 里，才算 `EQUIVALENT`。
- 主控建议 δ = 0.04（r_orth_cal 单位，是第十二轮最弱的 FCT−Mis 配对差 +0.080 的一半）。如果起草方有别的取值，要在锁定前写清依据。
- 每条线、每条错配线分成四类：`MATCHED_BETTER` / `EQUIVALENT` / `INCONCLUSIVE` / `MISMATCHED_BETTER`。

**I7 CellularAI_Pub 不能当作「VCC 2026 冠军算法」。** README 写明完整运行需要 "the separately supplied, access-controlled summary bundle"。用 r9_6 自建摘要表是一种没有验证过的改造。
- 降为次要分析，名字改成「CellularAI 开源代码在公开源线上的重拟合」。
- 不改的话就删掉；FCT 已经属于同一类方法。

**I8 事实更正。** 以下是主控 10:45 查询 HF API 的结果：
- zeroshot `final.ckpt` 是 202 MB（ST-SE）和 164 MB（ST-HVG），不是 ~525 MB；510 MB 的是 fewshot。
- SE-600M 的 `model.safetensors` 是 2.86 GB，另有 `protein_embeddings.pt` 0.41 GB。
- 三个仓库都没有门禁（gated=False）。
- X-Cell 仓库只有 README，与草稿一致。

**I9 Arc 已公开零样本预测结果。** 仓库里有 `zeroshot/<line>/eval_final.ckpt/adata_pred.h5ad`（ST-SE 每条线 3.15 GB）。
- 锁定之后用来做适配器核对：我们用匹配对照推断出的结果，至少在一条第 1 组细胞系上要复现 Arc 存档的预测。
- **锁定之前不准打开**，因为它会提前揭示第 1 组匹配对照这一臂的结果。

**I10 门槛顺序要改。** 权重和词表不是数据，看了不会知道结果。所以在得到人类批准之后：
1. 先下载权重，只检查结构：输出基因表、扰动 one-hot 词表、是否需要 batch/cell-type 协变量、输入格式。
2. 再锁定 PREREG_I，把确切的 G_eval、每条线的可预测扰动集、协变量的处理方式全部写死，删掉「冻结前可调整门槛」这类条款。
3. CODE_FREEZE_R16 之前，不对任何真实细胞做推断；用合成输入通过真实权重做冒烟测试是允许的。

**I11 泄漏和披露。**
- (a) 查清 SE-600M 预训练用了哪些数据；只要第 2 组某个目标数据集的细胞在里面，这条线就降为 `SENSITIVITY_ONLY` 并披露。
- (b) 第 1 组的错配对照来自检查点自己的训练线，要披露。
- (c) 第 2 组细胞系过不过门槛、FCT 结果是多少，我们都已经知道，所以只能写「预注册的回顾性模型审计」，不能写「前瞻确认」，也不和第十二轮的论断合并。
- (d) a549_ko 是 CRISPRko，STATE 是在 CRISPRi 数据上训练的。

**I12 VCC 2026 伴随分析的时间条件。**
- 只有当 `PREREG_I_COMMIT` 早于 2026/ 前缀出现第一个对象的时间（以 `round15/status/VCC2026_RELEASE_WATCH.log` 为准），这部分才算前瞻。
- 2025 年的先例是结项后 10 天发布，所以本轮目标是 **2026-11-27 前锁定**。
- 不改 PREREG_H 和第十五轮的冻结文件；伴随分析写成独立文件。

**I13 可行性也要写清楚。**
- SE-600M 要用原始计数和它自己的基因表，列出每条线原始 h5ad 是否在本地（只查计数，不看结果）。
- 估算在 MPS 上的推断时间。

## 2. 主控已做的检索（供 A-T5 继续）
- Arc 结项文：已读（斜率 0.02–0.13；领奖方案都改造目标线自己的对照细胞；没有做对照替换）。
- arXiv 2604.27646（Mao et al., 2026, in-the-wild benchmark）和 2611.02335（Rezaei Jafari et al., attribution diagnostics）：只读了摘要，都没有提到对照替换，需要读全文。
- ST-SE-Replogle 模型卡：没有写训练数据、协变量和 SE 预训练数据。
- STATE 论文全文（bioRxiv 2025.06.26.661135v2）：主控访问返回 403，没有读到。
- Systema：读了摘要和 EPFL 页面，与 I2 相关。

## 3. 对 v2 的复核（主控，2026-11-24 11:10 +08:00）
v2（sha256 549b3743…，11:04）逐条处理了 I1–I13，主控抽查了对应条文，确认已改。还剩一个设计缺口：
**I14 错配对照池必须是本地已有原始计数的细胞系。**
- 问题：v2 中第 2 组的错配池写的是「r9_6 全部 6 条线」，但 SE-600M 需要每条错配线的原始对照细胞。
  - replogle_k562_essential、replogle_rpe1 和 vcc2025_h1 都不在本地，共约 45 GB。
  - nourreddine_kolf21j 在 v2 里完全没有提到。
- 改法：主分析的错配池只用原始计数已在本地（或同步盘上）的线，不为主分析新增数据下载。
  - 第 2 组：S_mis(t) = {nadig_hepg2, nadig_jurkat} ∪（第 2 组中除 t 以外、过了门槛的线）；kolf21j 若有本地原始计数也加入，没有就写明排除。
  - 分两类报告：「检查点见过的线」（hepg2、jurkat）和「没见过的线」（其他第 2 组线）。H-I1 的线级判定要求两类都成立。
  - 第 1 组主分析只用 nadig_hepg2 和 nadig_jurkat；k562、rpe1 和 vcc2025_h1 作为可选扩展，需要人类另外批准下载，没有批准就不做，在论文里如实写。
- FCT 作为 H-I2 的对照，仍然用 r9_6 的伪 bulk 数据（已在本地），这一点不变。

## 4. 对 v3 的复核（主控，2026-11-24 11:15 +08:00）
- I14 已关闭：主分析只用本地已有原始计数的 6 条线（第 1 组 hepg2、jurkat；第 2 组 a549_ko、dld1、cd4t_rest、feng_ipsc），不新增数据下载；kolf21j 因为本地没有原始计数被排除；k562、rpe1、vcc2025_h1 列为需要另外批准的扩展。
- 锁定前补一句披露：第 2 组用 4 个检查点的平均，hepg2 和 jurkat 分别各被其中 3 个见过，所以「见过」的意思是「4 个检查点中有 3 个见过」。
- 下一步：等 A-T3 独立审核（PREREG_I_AT3_AUDIT_DONE）、A-T5 查新（PRIOR_ART_I_DONE），以及人类对下载权重的决定。
