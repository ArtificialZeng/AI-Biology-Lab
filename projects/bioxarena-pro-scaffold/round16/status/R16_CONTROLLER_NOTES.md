# 第十六轮主控记录

## 2026-11-25 23:07 A-T3 独立实现（合成数据）的核验
- 主控用 round16/.venv 重跑了 `round16.audit.test_audit_r16`：6/6 通过，返回码 0，用时 12 秒。test_6 让真实的 ST-HVG-Replogle 权重（zeroshot/hepg2）跑合成的 Poisson 计数，没有用任何真实细胞。
- **盲态说明：** A-T3 的函数名（build_mismatched_strata_r16、classify_pair_tost_r16）和主线代码一样。原因是主控 11-24 让 A-T3 审核 PREREG_I 时，要求它读主线的 adjudicate_r16.py 和 test_unit_r16.py。所以 A-T3 的判定逻辑不能算完全独立于主线；推断适配器、对照抽样和指标部分是 A-T3 自己写的。
- **补救：** 主控只根据 PREREG_I 正文，独立写判定和配对比较的核对代码（controller/checks/r16_check.py），作为三方复核的第三方；以后不再让任何一方读另一方的代码。
- 结构线索（来自 A-T3 测试打印的模型结构）：ST-HVG 有 `batch_encoder: Embedding(56, 128)`，也就是 56 个 gem group 批次的嵌入。对没见过的线怎么设 batch 输入（A-T3 实现了 `bypass_zero` 和 `index_zero` 两种）是 Gate 2 的关键决定，必须在锁定前写死，并且匹配臂和错配臂完全相同。

## 2026-11-25 23:09 主控独立判定代码 controller/checks/r16_check_adjudicate.py
- 只根据 PREREG_I v4 §3 正文编写，不读 round16/src 和 round16/audit；28/28 分支测试通过。
- **发现 PREREG_I 的一处歧义（v5 必须改）：** §3.1 第 4 条规定「不是 MATCHED_BETTER、不是 MISMATCHED_BETTER，且 within_delta」就判 EQUIVALENT；第 5 条的括号却说 `ci975_lo == 0.0` 这种平局判 INCONCLUSIVE。当 ci975_lo 恰好等于 0、而 95% CI 落在 ±0.04 里面时，两条互相矛盾。
- 主控的处理：取保守的读法，97.5% CI 的任一端恰好等于 0 都判 INCONCLUSIVE。v5 要把这一条写成单独的一行，主线和 A-T3 的代码都要照此实现，并各加一个对应的单元测试。

## 2026-11-25 23:27 A-T5 第二个模型盘点（PRIOR_ART_I_ADD1）的抽查
- A-T5 的结论：真正满足「权重公开、免登录、以对照细胞为输入、没有细胞系 ID、能在 M3 Pro 上推断」的非 STATE 模型很少。主推 Occam，备选是 scGPT 微调权重和 TxPert。
- 主控核实了以下几项：
  - **Occam 论文**（arXiv 2611.21959，Bendidi, Kenyon-Dean, Denton, Genovesio, Noutahi，2026-11-16）读了 HTML 全文，没有任何替换、打乱或去掉输入对照的实验，所以第十六轮的新颖性不受影响。最接近的是它的 Transfer 基线（"Transfer adds the same perturbation's effect measured in other training cell lines to a control profile from the evaluation cell line"）。在没见过的线上，这个基线的 Pearson Δ 比 Occam 高（K562 0.271 对 0.212，H1 0.232 对 0.107），和我们的问题方向一致。
  - **Occam 权重：** HF bendidiihab/Occam，gated=False，仓库共 42.9 GB；只需要 checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth（2.55 GB）加上 3 个 rhj_fewshot 特征文件（0.38 GB），合计 2.93 GB。代码许可是 MIT（GitHub IhabBendidi/Occam 的 LICENSE）；权重仓库没有单独写许可。
  - **Occam 的训练线**是 RPE1、HepG2 和 Jurkat（RHJ），所以第 1 组的 hepg2 和 jurkat 不能作为 Occam 的目标线，Occam 只评第 2 组；hepg2 和 jurkat 在 Occam 那里属于「训练时见过」的错配来源。
  - **scGPT 微调权重**（matthewshu/scgpt-replogle-base-ft，207.5 MB）是第三方上传的，HF 上没有写许可（scGPT 代码本身是 MIT）；仓库里另有 6.3 GB 的预测文件，不下载。
- 主控建议：第二个模型用 Occam，写成 PREREG_I 的次要模型，或者在任何 Occam 推断之前锁定一份增补。下载要人类批准。

## 2026-11-25 23:52 Gate 2 复核（R16_GATE2_SCHEMA.md，23:49）
- 主控直接读权重文件核对：G_model 2,000；pert 2,024；batch_encoder 56×128（ST-SE 行范数 9.78–12.41，中位 11.09；ST-HVG 8.81–12.48）；ST-HVG 的 gene_decoder 在 4 折之间完全相同（没有训练过），ST-SE 的各折不同。都和报告一致。
- **重要发现，待真实数据检验：** 在合成输入上，batch 嵌入（范数约 11）远大于对照细胞编码（0.4–1.2）；只换 batch 编号对输出的改变是只换对照的 3.6–5.7 倍。对没见过的线，官方代码一律用编号 0。真实细胞上是否也是这样，由 §4.4 的二乘二拆解来检验；在那之前不能写成结论。
- 主控裁定（已发给主线做 v5）：
  - 采纳全部 TBD_GATE2 常数、平局规则和 §4.4；
  - N_eval^R16 写成确定的数字，允许读 npz 里的可重复性标记来计数；
  - ST-HVG 主分析用官方预处理，1e4 归一化作为敏感性；
  - 新增 gem_group 编号与训练线的对应关系（描述性）；
  - 同时准备带硬门的真实数据运行脚本 run_eval_r16.py。

## 2026-11-25 23:58 A-T1 新数据集筛选（R16_CANDIDATES.md）
- 结论：没有任何从没用过细胞类型的公开数据集，能满足「与必需基因重叠 ≥ 50」以及原始 UMI、guide 和 NTC 都公开。
- 第 1 名 GSE295214（H1 hESC，CRISPRi，基因组规模，10x 5'，2,971 个靶点，其中 1,563 个在必需基因库里，8.34 GB）：主控在 NCBI E-utilities 上核实了题目、物种、公开日期（2026-05-05）和 85 个 H5 样本。H1 与 R9 的 vcc2025_h1 是同一细胞系，但 STATE 的检查点没见过 H1，数据晚于 SE-600M，本项目也从没分析过。
- 决定：作为第 3 组前瞻留出线写进 PREREG_I v5，下载之前锁定；下载要人类批准。
- 第 2 名 Pfizer TeloHAEC（Zenodo 18792213，CC BY 4.0，49.6 GiB 单个文件，291 个必需基因靶点）：暂不采用（文件太大，而且 R11 的 TeloHAEC 没过门槛）。

## 2026-11-26 00:15 A-T3 真实数据运行脚本（run_audit_r16.py）的核验
- 主控重跑 test_audit_r16：8/8 通过（61 秒）。
- 主控没有用「真的启动一次真实运行」的办法去测试硬门：万一门是坏的，就会在锁定前碰到真实细胞。所以主控改为读代码：run_real_r16 第一步调用 assert_freeze_gates_r16，在读取任何数据之前检查 PREREG_I_COMMIT.txt、CODE_FREEZE_R16.txt、R16_GO.txt 三个文件是否都存在、都非空，并要求 R16_GO 比 R16_GATE2_SCHEMA.md 新。逻辑正确。
- 补充防护：这三个门文件只能由主控创建。监控脚本现在只要发现其中任何一个出现，就立刻报警。
