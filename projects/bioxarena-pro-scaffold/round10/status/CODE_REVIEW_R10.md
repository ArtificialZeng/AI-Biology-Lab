# 主控代码审查：第十轮 C 组管线（A-T2 `round10/src/`）和独立复算（A-T3 `round10/audit/`）

审查人：Claude（主控） ｜ 2026-11-22 11:55 ｜ 依据：PREREG_C（dc4b4a5，sha256 0fb48aeb…）

**结论：通过，可以冻结。** 没有阻断项。下面两条是从第九轮冻结代码继承下来的命名和严格程度差异，记录在案，不改代码。

## 核对过的内容
1. **FCT 公式**（`fct_pipeline.py::compute_source_transfer_and_fc_bar`、`fct_from_fc_bar`）：和 PREREG_C §3.1 / PREREG_A v4 §4.6 一致，源线 μ_s 下限 1e-8，几何平均倍数变化，log1p(μ_t·FC̄) − log1p(μ_t)。目标线 μ_t 的下限写成 0（文字是 1e-8），两者差别约 1e-8·(FC̄−1)，可忽略。
2. **6 个错配对照均值**：每条源线一个，用源线的 μ_s′ 替换 μ_t，FC̄ 不变。
3. **100 次置换**：`np.random.default_rng(k).permutation(G)`，k = 0..99，只置换 μ_t，FC̄ 不变。
4. **H-C2**（A-T2 `adjudicate_r10.py`、A-T3 `run_r10_audit.py`）：FCT 的 97.5% CI 下界 > max(6 个错配对照 CI 上界, 置换分布 97.5 百分位)，在全部 K_conf 条线上都要通过。两套实现一致。
5. **P1**：K_conf < 2 时不做确认性判定，摘要用规定措辞。两套实现都有。
6. **冻结闸门**：两套代码都要求 `PREREG_C_COMMIT.txt` 和 `CODE_FREEZE_R10.txt` 存在且非空，才能读第十轮的线。
7. **主控独立复算（自写代码 `controller/checks/r10_fct_check.py`，只用第九轮的线做留一线测试）：**
   - H1：FCT r_orth_cal = +0.3505，置换 97.5 百分位 = +0.2445，和 A-T2 完全一致。
   - HepG2：FCT = +0.2611（N = 799），A-T2 为 +0.2614（N = 797，HGNC 去重后少 2 个靶点）。错配对照上界 0.1317、置换阈值 0.1696，两边一致。

## 记录在案的两条差异（不改代码）
- **D-R10-1 命名。** PREREG_C §4.2 写的 H-C3 六种口径是 `ctrl_disc_ownmasked`、`ctrl_disc_panel`……，代码沿用第九轮冻结代码的名字：`ctrl_disc_raw` 就是 own-masked（每行只遮自己的靶基因），`ctrl_disc_raw_panel` 就是 panel（遮掉全部评测靶基因）；`cent_*` 同理。定义相同，只是名字不同。论文里统一用预注册的名字，并注明对应关系。
- **D-R10-2 H-C1 多了两个保护条件。** A-T2 的 H-C1 除了「CI 下界 > 0」，还要求 bootstrap 有效比例 ≥ 0.95、正交 split-half 相关 ρ_AB⊥ > 0。这沿用了第九轮冻结代码，只会让判定更严。如果出现「CI 下界 > 0 但有效比例 < 0.95」的情况，两种结果都要报告，并按偏差处理。

## 冻结范围
`round10/src/`、`round10/audit/`、`round10/shared/scripts/`（含 `prepare_r10_lines.py`），以及 `controller/checks/r10_fct_check.py`。冻结之后改动任何一个文件都要写进 DEVIATIONS。
