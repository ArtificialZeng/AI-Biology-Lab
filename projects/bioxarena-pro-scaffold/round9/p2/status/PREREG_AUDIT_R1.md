# 主控核验 PREREG_A.md（第 1 轮）：退回修改，暂不锁定

核验人：Claude（主控） ｜ 时间：2026-11-21 21:10 +08:00
核验对象：`p2/prereg/PREREG_A.md`（sha256 16c30a97…852afa）

**结论：不通过。** 改完下列 R1–R11，A-T2 和 A-T3 在 Replogle 开发对上跑通新增对照并补好单元测试后，重新提交 `COMMIT_REQUEST_A-T3.txt`。在主控写出 `PREREG_COMMIT.txt` 之前，规则 3 照旧：不得在留出细胞系上计算任何指标。

## 已核实、无需修改
- 参考文献 1–15：Crossref 和 arXiv API 逐条核对，标题、第一作者、作者人数、年份都对得上。
- 零样本边界、λ 和 P⊥ 只用源细胞系、split-half 噪声天花板、联合 bootstrap、「未确认」硬规则：写法正确。
- 开发阶段如实写明 H-A1 的斜率只有 0.005–0.009，没有事后放大：保留。

## 必须修改（阻断项）

**R1 措辞。** 删掉 "Principal Investigator" 和 "Human PI"，改成 "Human lead: Zijian Zeng"。全文不出现 PI。

**R2 确认性细胞系只用真正留出的线。** 配置（log1p、shrink=0、mask_self、Hybrid_ProjShared）是看着 K562 和 RPE1 的结果选出来的。这两条线的 6 线留一结果照常报告，但不计入判定。

**R10 最低靶点数。** 进入确认性判定的线，主配置下可评测靶点必须 ≥ 50。按元数据（不是结果），nourreddine_kolf21j 只有 18 个，降为次要线，单独报告。
- 确认性线因此是 3 条：nadig_hepg2、nadig_jurkat、vcc2025_h1。
- 阈值按卡片比例取整，只能更严：H-A1、H-A2、H-A3 要求 3/3 通过（⌈0.8×3⌉=3）；H-A4 要求 0/3 条线 CI 严格小于 0（⌊0.2×3⌋=0）。
- 第 6 节的 N_pass 和 N_inf 公式按这个改。

**R3 错配协方差不能只挑一条源线。** 现在 `run_eval.py` 用 `src_inputs[t_idx % len]` 选 s'，等于任意挑一条（例如 KOLF 配到 K562，而不是和它最像的 H1）。改成：
- 对全部 5 条源线 s' 都算 `NegCtrl_HybridProj_MismatchedCov(s')`。
- 目标线 t 上 H-A3 通过的条件：`matched_ci95_lo > max_{s'} mismatch_ci95_hi(s')`。
- 5 个错配结果全部报告。

**R4 区分「共波动结构」和「目标线边际方差」。** 新增两个预注册的次要对照，不进入四项判定，但结论措辞受它们约束：
- (a) `NegCtrl_HybridProj_SrcCorrTgtVar`：C = D_t^{1/2} · R̄_S · D_t^{1/2}，其中 D_t 是目标线对照的基因方差，R̄_S 是源线相关矩阵的平均。也就是目标线方差加通用相关结构。
- (b) `Transfer_Mean_VarScaled`：Δ̂_{p,t,g} = Δ̄src_{p,g} · sqrt(C_t,gg / mean_s C_s,gg)，只用目标线边际方差。
- 写死的措辞规则：
  - 只有当 Hybrid_ProjShared 在全部确认性线上，ls_slope 和 cal_r 都以不重叠 CI 胜过 (a)，结论里才可以说「细胞系特异的共波动结构」。否则只能说「目标线对照统计量」。
  - 如果 cal_r 不胜过 (b)，必须写明增益可由边际方差解释。

**R5 H-A4 不留选择余地。** 写明非劣性必须在每一种口径上同时成立：ctrl_* 和 cent_*，A-T2 和 A-T3 的 disc 口径都算，取 AND。主表只报 A-T2 口径（自身靶点置零、ctrl_*），其余进附表。

**R6 斜率口径的局限。**
- VCC 结项报告和 cell-eval2 都没有公开「0.02–0.13」这个斜率的估计式（A-T5 已核实）。所以 H-A1 是在和一个估计式未经核实的数字比较，原文必须写明。
- OLS 斜率（预测对观测）等于 r·sd_pred/sd_obs，依赖尺度；MSE 最优的预测，斜率约等于 r²。
- 每条线额外报告 r、sd_pred/sd_obs、Deming 斜率和「观测对预测」斜率，作为次要结果。

**R7 和 CIPHER v2 的边界写清楚。**
- CIPHER v2 摘要（bioRxiv 10.1101/2025.06.27.661814 v2）已经写了：去掉基因间协方差、只留方差，性能下降数倍；相关结构能在同一细胞类型的不同研究之间迁移。
- 所以「对照协方差有信息」不是我们的新发现。预注册第 1 节要说明，我们检验的只是：跨细胞类型零样本时，它能否恢复线特异残差。
- A-T5 把 Luan 等 2026（10.64898/2026.10.20.760514，KOLF2.1J→H1 中心化迁移，用的数据和我们的两条线相同）补进 `PRIOR_ART_A.md` 的审计表。

**R8 如实披露。** nadig_hepg2 和 nadig_jurkat 在第 7、8 轮被别的跨系方法（S5 等）分析过。仓库检索确认，第 9 轮之前没有任何代码用过对照协方差。第 2 节写明这两点。

**R9 前瞻部分更新。**
- 人类负责人已于 2026-11-21 19:28（+08:00）发邮件给 help@virtualcellchallenge.org，询问赛后能否用 controls-final 做研究、能否先公开预测哈希。目前在等回复。
- 智能体仍然不得下载受限数据。

**R11 增加一项描述性指标。** 每条线报告 r_cal² 及其 CI，含义是「恢复的可解释线特异方差比例」，供论文和对外报道如实引用。

## 执行分工
- **A-T3**：改 PREREG_A.md。`audit_runner.py` 用独立代码实现 R3 和 R4，补单元测试。重新发 COMMIT_REQUEST。
- **A-T2**：在 `src/` 实现 R3（全部 5 个错配）、R4（a、b）、R6 的附加指标。只在 K562↔RPE1 开发对上运行，补单元测试，结果写进 `results/*_dev_*`。dev 上源线只有一条，R̄_S 就是那条源线的相关矩阵。
- **A-T5**：R7 的 Luan 审计。
- **主线执行者**：等 A-T3 重新提交后，只提交 `p2/prereg/PREREG_A.md`、`p2/.gitignore` 和 `status/`。不要提交 `shared/lines/`（已由主控加进 `p2/.gitignore`，6.8 GB）。push 后把 commit 哈希写进 `p2/status/PREREG_LOCK_PUSHED.txt`。`PREREG_COMMIT.txt` 由主控核验远端 commit 后再写。

## 补充（21:05）：第 1 版已经提交
- 主线执行者在 20:56:39 提交并推送了第 1 版（commit 222cf75），早于这份核验送达。当时还没有任何留出评测，PREREG_COMMIT.txt 也没写，所以规则 3 没有被违反。
- 第 2 版开头加「§0 修订记录」，写明三件事：
  - 第 1 版 = 222cf75（2026-11-21 20:56:39 +08:00），只做过 K562↔RPE1 开发评测；
  - 第 2 版按主控核验 R1–R11 修改，逐条列出改了什么；
  - 判定只以第 2 版为准。
- 第 1 版不要删，也不要改写 git 历史。
