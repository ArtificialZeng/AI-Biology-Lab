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

## 2026-11-26 00:28 锁定，以及主控复算代码 r16_check.py
- **PREREG_I v5 已锁定：** 提交 c1fc7c6，2026-11-26 00:24:20 推送到 origin/main，核实 sha256 43d242e4…；锁定前一刻（00:23:53）2026/ 前缀仍是空的。锁定记录 56a9ce9。git 配置的代理 127.0.0.1:7892 没有在运行，所以只这两次推送临时绕开了代理（-c http.proxy=），没有改用户的 git 配置。
- **A-T3 锁定前终审：** 没有必须改的项目，7 项补记由主控裁定，写进数据前增补 A。
- **controller/checks/r16_check.py：** 指标代码原样沿用冻结的 r15_check.py，判定调用 r16_check_adjudicate.py；从导出的 observed.npz 和 preds_*.npz 复算。合成数据测试结果：
  - 用到对照 → LINE_MATCHED_SPECIFIC；
  - 没有信号 → LINE_EQUIVALENT_NULL_SIGNAL；
  - 有信号但换对照只加噪声 → LINE_MATCHED_SPECIFIC_SUB_DELTA（配对差约 +0.0015）。
- 第三种情形说明需要一个参照，所以增补 A 加了同线重抽参照臂（第 8 条）、SUB_DELTA 的解读规则（第 9 条）和 EXACT_INVARIANCE 边界情况（第 10 条）。

## 2026-11-26 01:01 A-T3 对齐 v5 之后的复跑：发现一个不稳定的测试
- A-T3 报告 12 个测试都通过；主控复跑时 test_11 在 22 次里失败了 2 次（G4 对称检查的估计值是 NaN）。另外 20 次通过，包括 8 个固定和 12 个随机的 PYTHONHASHSEED，以及单线程和默认线程。
- 两次失败都发生在 A-T3 正在运行的时候，推测是测试读写了共享目录，或者合成夹具落在 NaN 边界上。
- 要求（冻结前的硬门）：测试只写临时目录；两套实现都要有跨进程复现性测试（CPU 上逐位相同，MPS 上 ≤ 1e-6）；真实运行不能靠遍历 set 来决定顺序。
- 两套实现都已确认对模型调用了 .eval()（排除了 dropout 带来的随机性）。

## 2026-11-26 01:31 复现性：SE-600M 嵌入的随机打破平局
- A-T3 第 7 阶段找到 test_11 不稳定的原因：夹具的两半是互相独立的噪声，ρ⊥ 在 0 附近，导致 NaN；另外它自己的嵌入代码用了不稳定的排序。修好之后，主控在两个 PYTHONHASHSEED 下各跑一次完整测试，都是 13/13 通过（第三次在后台跑）。
- 更重要的一点：官方 SE-600M 的 collator 故意用 torch.randperm 随机打破未表达基因之间的平局，所以嵌入取决于 torch 的全局随机状态。主线照官方方式调用并设了种子；A-T3 自己重写了这一步，平局的打破方式和官方不同。
- 裁定（增补 A 第 11 条）：两方都调用官方 collator，统一种子、细胞顺序、基因列顺序和 batch_size=16，并导出每个对照池的嵌入，先单独比对嵌入（≤ 1e-5）。

## 2026-11-26 01:46 主控直接比对两方的 SE-600M 嵌入（增补 A 第 11 条落实后）
- 主控自己写的脚本：从 a549_ko 的 genes.txt 取 4,000 个真实基因符号，生成 32 个细胞的合成 Poisson 计数（种子 424242）；主线的 encode_controls_se600m_r16 和 A-T3 的 CleanRoomSE600MEncoder.encode_cells 都在 CPU 上运行，batch_size=16，种子 20261601。
- 结果：两方都映射到 3,132 个基因；嵌入形状都是 (32, 2058)；**逐位相同**（最大绝对差 0.0）。主线自己重跑一次，结果也逐位相同。
- 结论：嵌入这一步在 CPU 上逐位可复现、两方逐位一致。三方复核里预测差那一层（≤ 1e-4）剩下的差异只可能来自 MPS 运算或我们自己的代码。

## 2026-11-26 02:15 冻结完成；等外置硬盘
- 增补 A 锁定 4ecfc2c；代码冻结 848bdac（FREEZE_MANIFEST_R16 共 45 项，和 origin 核对 45/45）；CODE_FREEZE_R16.txt 3c357c7。
- 阻碍：外置硬盘「雷雳传输」没有挂载，Nadig HepG2/Jurkat 的原始文件只有那一份（Spotlight 没有找到别的副本）。已请人类接上；另一个办法是从 GEO 重新下载约 14 GB，需要人类批准。
- A-T4 的 paper/MS2_METHODS.md 主控已抽查：锁定的数字都对，没有结果，没有优先权的说法，没有进 git。

## 2026-11-26 02:40 监控误报：别的项目里出现了大文件
- watch_r16d 报警：projects/project4_b1_faers_algorithmic_failure/data/processed/master_demo_dedup.parquet（211 MB，02:38），没有被 git 忽略。
- 这是同一个仓库里的另一个项目（FAERS，02:19 有自己的预注册锁定文件 PREREG_B1_LOCK.json）。最近 12 个 Antigravity 会话里都没有提到它，所以不是本项目的代理做的，和第十六轮无关。主控不动它。
- 主控提交时只 add 明确列出的第十六轮路径，所以不会把它带进 git；已告诉人类这个文件没有被忽略。监控的大文件检查已改为只看 bioxarena-pro-scaffold 目录。

## 2026-11-26 15:00 — Gate 5 released (controller)
- 14:41 the sync share 雷雳传输 came back (it is SMB //192.168.3.2, not a USB drive). At the human lead's request, data now lives on the internal disk: Nadig raw files copied to /Users/mac/.cache/r9_p2_raw/nadig2025 and verified against the GEO sha256 sidecars (14:45); the rest of the share's external_data folder mirrored to /Users/mac/.cache/sync_external_data (55 files, 11 GB).
- Only the two Nadig raw h5ad files were read by Round 16 from the share; every other input was already local (round9/p2, round10-12 shared lines and raw). Both runners accept --raw-nadig-dir, so no frozen code changed.
- FREEZE_MANIFEST_R16_DATA.txt: 183 files, 82.4 GB, written and then re-hashed (VERIFY OK). feng_ipsc and a549_ko line files equal FREEZE_MANIFEST_R11 / R12 (5/5 each). Commit 6f6c303.
- Re-verified before GO: FREEZE_MANIFEST_R16 45/45, R16_MODEL_CHECKSUMS 77/77, no round16/results, no runner process. R16_GO.txt committed 7106a3b.
- 14:55 GO sent to A-T3 (controller/msgs/r16_GO_AT3.txt). A-T3 runner started 14:55 with the exact dispatched command (pid 60770).
- The archived nadig_hepg2 adata_pred.h5ad (3,151,145,504 bytes, LFS sha256 505fcfd3...ee22) is fetched by the controller via controller/checks/r16_fetch_arc_archive.sh into /Users/mac/.cache/r16_arc_archive (outside the repo). The main-line GO (controller/msgs/r16_GO_MAIN.txt) goes out after that file verifies, because the main run includes the Addendum A section 3 check.
- Run-phase watcher controller/watchers/watch_r16e.sh. A bash-3 marker-count bug was fixed before it was relied on.
- 15:04 Arc archive verified (size + LFS sha256; R16_ARC_ARCHIVE_MANIFEST.txt). Single-connection direct speed was about 740 KB/s (the local proxy is not running), so the fetch used 12 parallel byte-range segments, about 11 MB/s. 15:05 main-line GO sent (controller/msgs/r16_GO_MAIN.txt).
- 15:44 main line recorded R16_MAIN_DEVIATION: cd4t_rest lanes use two gene orders (13 lanes genomic, 34 Ensembl-ID; same 18,130 unique genes). The frozen main-line exact-order check raised before any inference. Controller ruling: DEVIATION_R16_D1 (reorder by name onto parts[0], as Round 10 and A-T3 do). Controller test R16_D1_TEST_PASS; main-line unit suite 18/18 OK. FREEZE_MANIFEST_R16_D1 (47 entries) supersedes FREEZE_MANIFEST_R16 for execution. 15:54 rerun ordered (controller/msgs/r16_D1_MAIN.txt). A-T3 unaffected and still running.
- 15:58 A-T3 recorded R16_AT3_DEVIATION (feng_ipsc Cell_ID format). The main line's feng reader has the same gap, so a HOLD went to the main line at 16:01 (it stopped its D1 rerun at 16:02:35, before inference). The controller then dry-ran raw-control extraction (no model) for all 6 lines with both parties' code. This found D3 (main line, Nadig legacy categoricals, 0 NTC), D4 (A-T3 Nadig row-id barcodes) and D5 (main line, Arc archive positional gene index; the section 3 check would have halted). DEVIATION_R16_D2_D5.md records them. The controller's first D2 message misnamed the feng barcode form; it was corrected before any run. 17:25 dry run: all 6 lines exactly identical between parties. FREEZE_MANIFEST_R16_D2 (51 entries) supersedes D1. 17:32 both parties released (controller/msgs/r16_GO2_*.txt).
