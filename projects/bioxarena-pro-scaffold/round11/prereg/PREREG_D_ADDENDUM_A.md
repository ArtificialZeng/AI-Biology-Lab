# PREREG_D 补充 A（与 PREREG_D v2 同时锁定，任何第十一轮数据下载之前）

作者：Claude（主控） ｜ 2026-11-22 14:40 ｜ 适用于 `PREREG_D.md` v2（sha256 1a0544980c094468faed54ba339608d8f68d2205358d04c1dc59622b5a98fd30）和 `round11/shared/R11_DOWNLOAD_PLAN.md`（sha256 3e3991368ff0bf4d4e590cba11f770e1b71feeab2d805f149082b7dd25b44ecf）

本补充只澄清实现细节，不改变任何假设、阈值、数据集或掩码。

**A1 H-D1 的有效性条件。** 和 PREREG_C 的冻结实现一致，H-D1 在「FCT r_orth_cal（all_nontarget）97.5% CI 下界 > 0」之外，还要求 bootstrap 有效比例 ≥ 0.95，并且全掩码下的正交 split-half 相关 ρ⊥_AB > 0（否则噪声天花板没有定义）。这两条只会让判定更严。如果出现「CI 下界 > 0 但有效比例 < 0.95」，两种结果都要报告，按偏差处理。

**A2 cd4t_stim48 的两半共用测序通道。** 下载方案显示，两个供体（CE0008678、CE0010866）的样本都在同一批 24 条通道（GSM9393875..GSM9393898，CD4i_R2L25..R2L48）上，每条通道每个供体一个文件。按供体分两半时，两半会共享通道层面的技术噪声。每个供体在每条通道内都先减去自己的 NT 均值，所以通道效应已经扣掉一大部分；剩下的共享噪声会让两半更相似，抬高噪声天花板，使 r_orth_cal 偏低（偏保守）。H-D2 和 H-D3 是在同一个天花板下比较，排序不受影响。论文里要写明这一点。

**A3 feng_ipsc 的两半共享供体。** 17 条非 KOLF 细胞系按 ID 奇偶分两半后，同一供体的不同克隆会分在两半（例如 eipl_1 在 Half A、eipl_3 在 Half B）。所以它的噪声天花板反映的是克隆和通道层面的差异，不是供体差异。feng_ipsc 只作次要报告，这一点在报告里注明。

**A4 锁定范围。** 本次提交锁定 `PREREG_D.md` v2、本补充、`R11_DOWNLOAD_PLAN.md` 和 `round11/status/PREREG_D_AUDIT_R1.md`。之后对任何一个文件的修改，都要写进 `round11/prereg/DEVIATIONS_R11.md`，并说明是在看到数据之前还是之后。
