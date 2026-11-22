# 主控核验 PREREG_C v1（第 1 轮）：小改后锁定

核验人：Claude（主控） ｜ 2026-11-22

**总体：** 六步顺序、确认性线规则、split-half 分法、H-C1/C2/C3、S7 细胞数对照、两套实现交叉核对的规定都齐全。以下 P1–P7 改完即可提交锁定。

**P1 少于 2 条确认性线时怎么办。** 写明：如果 K_conf < 2，第十轮不做确认性判定，只报告描述性结果，摘要里写「确认性线不足，未做确认性检验」。

**P2 H9 在下载前就按元数据决定去留。**
- H9 的 MOI 约为 4，单 sgRNA 细胞只有约 5%。很可能达不到「≥ 50 个靶点、每个靶点 ≥ 25 个单 sgRNA 细胞」，却要占 51 GB。
- 新规则：A-T1 先用 IGVF 的小文件（0.24 MB 的 CSV 或文件元数据）数出每个靶点的单 sgRNA 细胞数。这是 QC 计数，不是结果。
- 如果「单 sgRNA 细胞 ≥ 25、且和源线共享」的靶点少于 50 个，H9 在下载前就排除，不下载。
- 如果达到 50 个，再看预算：CD4+ T 取 2 个供体，加上 DLD-1、i3Neurons、Calu-3，约 34 GB；H9 再加 51 GB 会超出 60 GB 上限。这种情况下 H9 也排除，并在论文里写明「因预算排除」。

**P3 CD4+ T 的供体在下载前写死。** 用 GEO 里按字典序排列的前 2 个供体 ID（CE0006864、CE0008162，由 A-T1 核对）的静息态单细胞 .h5，按供体分两半。补一句说明：供体间的生物学差异会算进噪声天花板。

**P4 DLD-1 只用低 MOI 的单 sgRNA 细胞。** 写明：只用 Low-MOI 文库，只保留恰好分配到 1 个 sgRNA 的细胞。CD4+ T 和 i3Neurons、Calu-3 也一样。

**P5 评测靶点的定义写全。** 和第九轮主配置一致：靶点在目标线上可重复（Addendum B），至少有 1 条源线上也可重复，并且在 G_eval(t) 中（only_reproducible=True）。

**P6 S7 的第 2 项写清楚。** 删掉含混的「pooled partial_r_orth_cal」一句，改成：
- (i) 每个扰动单独算 r_orth_p（50 个基因上 r^obs⊥ 和 r^pred⊥ 的相关）；
- (ii) 报告 r_orth_p 和 ln(n_cells)、|k_p| 的 Spearman 相关及 97.5% CI；
- (iii) 把 r_orth_p 对 [1, ln(n_cells), |k_p|] 做 OLS，报告截距（即控制细胞数和敲低深度后的平均信号）的 97.5% CI；
- (iv) 三分位分层照旧。

全部是次要分析。

**P7 第九轮锚点按冻结代码定义。** SA-3 里的 Hybrid_ProjShared 一律按第九轮冻结代码（1ca9f60）的实现，也就是 3 维共享子空间投影。注明它和 PREREG_A v4 文字不一致（见 round9 `PAPER_REVIEW_R1.md` B1）。

改完后写 `round10/status/COMMIT_REQUEST_A-T3.txt`（附新的 sha256）。主控核验后写 `PREREG_C_COMMIT.txt` 和 `R10_GO.txt`。
PREREG_C 主控核验 R1 结论（2026-11-22 09:50）：P1–P7 全部落实，同意锁定。
核验文件 sha256：0fb48aeb08aa6befe9583030f33ec668155e9e937616d291ea0c211e0d1bf4c6
下载方案 R10_DOWNLOAD_PLAN.md：34.080 GB，H9 已在下载前按 P2 排除（K_conf 候选：dld1、cd4t_rest）。
下一步：主线执行者提交并 push PREREG_C.md，写 PREREG_C_PUSHED.txt；主控核对远端 commit 后写 PREREG_C_COMMIT.txt 和 R10_GO.txt。
