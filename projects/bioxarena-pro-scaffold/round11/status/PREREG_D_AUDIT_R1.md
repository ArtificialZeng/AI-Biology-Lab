# 主控核验 PREREG_D v1（第 1 轮）：改完即可锁定

核验人：Claude（主控） ｜ 2026-11-22 13:32 ｜ 对象：`round11/prereg/PREREG_D.md`（sha256 87e428b4…9f1f）

**总体：** 主掩码改为 all_nontarget、三条主假设、H-D4 方法学假设、两种掩码下都算置换、功效预估和「不得事后放宽」的声明都很好。主控核对了冻结代码：`metrics_r10.py::_extract_masked_idx_per_pert` 里 all_nontarget 就是 G_eval(t) 去掉自身靶基因后的全部有限基因，和 §3.1 一致。以下 Q1–Q8 改完再发 COMMIT_REQUEST。

## 必改
**Q1 置信水平改为 97.5%（和第九、十轮一致）。**
- 主要论断要求 H-D1、H-D2、H-D3 在全部确认性线上**同时**通过，这是交集-并集检验（intersection–union test）：每一项在 α 水平下检验，整体错误率就不超过 α，三项之间不需要再做 Bonferroni 校正。
- 本轮需要校正的是「论断」的个数：主要论断和 H-D4 方法学论断，共 2 个，所以用 Bonferroni 2 项，双侧 97.5% CI（[1.25%, 98.75%]），和 PREREG_C 完全相同。H-D_util 沿用 97.5%。
- 在 §4.1 写明这个理由。§6 功效预估按 z = 2.2414 重算。

**Q2 在下载之前写死候选线、文件和 split-half 分法。** 不能留到 R11_GO.txt 再定。依据 `round10/shared/R11_DATA_SCREEN.md`：
- **确认性候选（3 条）：**
  1. `telohaec`：TeloHAEC_2285（IGVFDS7943IIWZ / GSE210681，静息态，1.47 GB）；
  2. `hcasmc`：HCASMC-hTERT 971 基因（IGVFDS6248NDUW，`.hdf5` IGVFFI5999NXAF，42.89 GB）；
  3. `cd4t_stim48`：GSE314342 Stim48hr，供体 CE0008678 和 CE0010866（27.39 GB）。必须写明：这和第十轮 cd4t_rest 是同一实验室、同一细胞谱系，独立性较弱；但供体和细胞状态（TCR 刺激 48 小时）都是新的。
- **次要（不进 K_conf）：** `feng_ipsc`：Feng 2026 TargetedScreen，排除 kolf_2、kolf_3 后的 43 条 HipSci iPSC（4.48 GB）。它和第九轮源线 KOLF2.1J 是同一种细胞类型，所以预先写一个描述性问题 **H-D5**：在 all_nontarget 下，用 KOLF2.1J 对照均值做错配基线时，得分是否接近 FCT。这回答「FCT 的特异性是细胞系层面还是细胞类型层面」。
- **split-half 分法（逐条写死）：** telohaec 按排序后的 10x 文库 / 通道 ID 奇偶分；hcasmc 按 Parse 子文库 ID 奇偶分；cd4t_stim48 按供体分，按字典序 Half A = CE0008678，Half B = CE0010866；feng_ipsc 按排序后的细胞系 ID 奇偶分。如果某个数据集没有可用的单元元数据，就按 seed = 42 分层随机对半分，并在锁定前就在文中标明会这样做。
- **单 sgRNA 规则（P4）：** 逐条写明每个数据集如何判定「恰好 1 个靶向 sgRNA」，以及 NT / safe-targeting 对照怎么用（TeloHAEC 有 NTC 和 safe-targeting 两类，预先写明只用 NTC 还是两者都用；建议只用 NTC，和前几轮一致）。

**Q3 下载预算。** 主控把第十一轮的下载上限定为 **80 GB**（四个数据集合计 76.23 GB）；磁盘空余 2.4 TB。如果 hcasmc 的 `.hdf5` 在下载后、计算任何指标之前就证明无法按块读取（内存或格式问题），按「可行性排除」处理，在 QC 报告里写明原因。这必须发生在 CONFIRMATORY_LINES_R11.txt 之前，不能看过任何结果再排除。

**Q4 删掉 §2.2 里「ideally prioritizing…」「if … ever screened as a fallback」这类留有余地的话。** 候选池就是 Q2 的 3 + 1 条，不再增减。

## 建议修改
**Q5 H-D2 里错配对照的 CI 也用 97.5%**，和 FCT 的 CI 同一水平（随 Q1 一起改）。

**Q6 新数据和旧数据不重叠。** 写明：feng_ipsc 排除 kolf_2 和 kolf_3（KOLF2.1J 源自 kolf_2）；cd4t_stim48 的两个供体不同于第十轮的 CE0006864 和 CE0008162。A-T1 下载后先核对，再做 QC。

**Q7 摘要必须报告的内容。** 不论结果如何，摘要都要给出：每条确认性线上 H-D1、H-D2、H-D3 的判定；H-D4 的判定；K_conf 的数值；cd4t_stim48 的独立性限制（如果它进入 K_conf）。

**Q8 H-D4 的地位。** 写明 H-D4 是 2 个确认性论断之一（方法学论断），单独判定，同样要求在全部确认性线上通过。

改完后更新 `round11/status/COMMIT_REQUEST_A-T3.txt`（新的 sha256）。主控核验通过后，由主线执行者提交推送，主控写 `PREREG_D_COMMIT.txt` 和 `R11_GO.txt`。

## 查新补充后的两条（13:33 追加）
**Q9 H-D4 要引用并对标 Nicol, Shivakumar & Irizarry (2026)。** PRIOR_ART_D 的 Q4 结论是「部分做过」：他们已经证明，按真实 |Δ| 选前 k 个基因会因为共享对照噪声而虚高分数，改成用一半数据选基因、另一半评测后虚高消失。所以：
- §0 和 §4.3 写明 H-D4 是在他们的结论上做延伸：检验这种选择偏倚是否也奖励「不含目标线基线、但非常数的压缩变换」（E1）和标签置换，并且发生在去幅度后的跨细胞系指标里。不得声称首次发现选基因偏倚。
- 新增描述性次要分析 **SA-8：分半选基因掩码 `resp_top50_split`**。用 Half A 的 |Δ| 选前 50 个基因、在 Half B 上评测，再反过来，两个方向取平均；对 FCT、E1 和置换分别报告。这直接检验 resp_top50 的虚高是否就是 Nicol 等人说的机制。只作描述，不进入任何判定。

**Q10 主要论断的生物学表述。** Q5 查新显示「跨细胞类型的扰动响应在倍数变化空间比在 log1p 差值空间更守恒」没有找到先例。§0 的论断措辞改为：「在全基因口径下，用目标线自身基线把跨线平均倍数变化映射回 log1p 空间（FCT），能预测一部分去幅度后的线特异响应，并且胜过错配基线、置换和常数基线。」不要写成「线特异响应主要由基线决定」，因为 r_orth_cal 只有 0.05–0.19，解释的方差很小，要在局限里如实写出。
