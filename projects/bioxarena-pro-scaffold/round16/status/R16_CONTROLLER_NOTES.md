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
- 20:37 A-T3 export crashed on JSON serialization of numpy values (no metric computed). 20:47 A-T3 patched its own run_audit_r16.py (DEVIATION_R16_D6, controller-reviewed; tests 13/13 and self-check OK) and reran pass 1. FREEZE_MANIFEST_R16_D3 (52 entries) supersedes D2. The main line's Groups 1-2 rerun has continued since about 17:35.

## 2026-11-26 22:07–22:15 — main line Groups 1–2 done; controller recompute on main-line arrays

- `R16_MAIN_REAL_DONE` at 22:07:30. The main line then started the sensitivity cohort (feng_ipsc target, pid 18295). A-T3 pass 1 is still running (pid 86771).
- Addendum A §3 adapter check (`round16/results/arc_adapter_check_nadig_hepg2.json`): 639 perturbations, 1,843 genes, median per-perturbation Pearson 0.9749, pooled 0.9733. PASS (min 0.462, q25 0.955).
- Controller recompute: `controller/checks/r16_check.py`, which never imports round16/src or round16/audit, run on `round16/results/arrays/<line>/` for all five lines with the preregistered strata (PREREG_I §3: seen = other Nadig line(s); unseen = Group 2 ∪ feng_ipsc minus target). n_boot 4000, seed 20261601. Outputs are in `controller/checks/r16_results/main/`.
- `controller/checks/r16_compare_main.py` compared every estimate, CI bound and valid_frac, the per-source categories, and the H-I1 and H-I2 calls against `round16/results/summary_r16.json`: **R16_MAIN_VS_CONTROLLER_MATCH**, worst |diff| 5.6e-16 (tolerance 1e-5).
- Control barcodes (`round16/results/control_barcodes_r16.json`): ordered barcodes and batch labels equal the controller dry run for all six lines, and the within-line resample equals the dry run for the five target lines. feng_ipsc has no within-line resample in Groups 1–2 (`eligible: false`; it is not a Groups 1–2 target).
- What this check does not cover: it verifies the metric and adjudication layer on the main line's own arrays. Whether those arrays (embeddings, predicted Δ, observed Δ, FCT) are right is checked by the comparison with A-T3, which is still pending. No main-line result has been passed to A-T3.

## 2026-11-26 22:15–22:35 — staged embedding check fails on 4 lines → DEVIATION_R16_D7 (A-T3 SE-600M gene mapping)

- `controller/checks/r16_compare_embeddings.py`: main `round16/results/embeddings/x_state_*.npz` vs A-T3 `round16/audit/results/se600m_embeddings_r16.npz` (written 20:35 by the run that later failed on D6; same SE code as the current run).
  - Ordered barcodes are equal on all 11 pools.
  - X_state is bitwise equal (max |diff| 0) on a549_ko and cd4t_rest, primary and resample.
  - It differs on dld1 (0.049; resample 0.042), feng_ipsc (0.052), nadig_hepg2 (0.044; 0.048) and nadig_jurkat (0.094; 0.095). The tolerance is 1e-5.
- These four are exactly the lines whose raw gene lists contain duplicated symbols (dry run: dld1 and feng_ipsc 10 each, each Nadig line 1).
- `controller/checks/r16_d7_mapping_test.py` (no model) builds the arc-state 0.11.1 `FilteredGenesCounts` map both ways on the dry-run gene lists:
  - Main line follows Addendum A §11 item 3 (`var_names_make_unique()` after setting raw symbols): 0 duplicate global IDs.
  - A-T3's `encode_with_loaded_se_model` sets `adata.var_names = gene_syms` without make_unique. The library maps by var_names, so duplicated symbols map twice: 8 duplicate global IDs on dld1 and feng_ipsc (TBCE, CYB561D2, MATR3, HSPA14, GOLGA8M, GGT1, ARMCX5-GPRASP2, TMSB15B), 1 on each Nadig line (HSPA14).
  - On a549_ko and cd4t_rest the two maps are identical.
- An MPS rerun of A-T3's encoder (`r16_d7_embed_test.py`) was stopped unfinished at 22:32, because it was queued behind both live runs on the GPU. The mapping test explains the pattern without it. The bitwise confirmation will come from the post-fix comparison.
- Ruling (`controller/msgs/r16_D7_AT3.txt`, sent 22:35): HOLD A-T3 pass 1. The fix is one line, `adata.var_names_make_unique()`, in `audit_r16_state_adapter.py`, with tests and a self-check. Then pass 1 reruns from scratch, then pass 2. No change to the preregistration, Addendum A, the main line or the data. The main line is compliant, and its metrics and status were not disclosed.
- Watcher relaunched (MAN_EXCLUDE / STRAY_EXCLUDE for audit_r16_state_adapter.py only; extra marker R16_AT3_D7_PATCHED).

## 2026-11-27 00:11–00:17 — main-line sensitivity verified; A-T3 rerun hung on MPS and was restarted

- `R16_MAIN_SENSITIVITY_DONE` at 00:11 (feng_ipsc as a SENSITIVITY_ONLY target, N 88, G 1,575).
  - Controller recompute: `controller/checks/r16_check.py` on `round16/results/sensitivity/arrays/feng_ipsc`. Strata follow the main summary: seen = nadig_hepg2, nadig_jurkat; unseen = a549_ko, cd4t_rest, dld1. n_boot 4000, seed 20261601. Output is `controller/checks/r16_results/main_sensitivity/feng_ipsc.json`.
  - `r16_compare_main.py` against `round16/results/sensitivity/summary_r16.json`: **R16_MAIN_VS_CONTROLLER_MATCH**, worst |diff| 3.3e-16.
  - As for Groups 1–2, this covers the metric layer only; the arrays are checked against A-T3's pass 2.
- A-T3 pass-1 rerun (pid 34240, started 22:43) wrote its control barcode files at 23:55:45. After that it gained under 2 s of CPU time in about 20 minutes.
  - Its main thread was in a single `-[_MTLCommandBuffer waitUntilCompleted]` (`MPSStream::commitAndWait` from `Tensor.cpu()`) in every stack sample.
  - A separate MPS probe (20 matmuls of 1024×1024) finished in 0.14 s, so the GPU was usable. The same kind of hang hit the controller's own `r16_d7_embed_test.py` at 22:32. Both hangs happened while another MPS process was running; this time the main-line sensitivity run was on the GPU until 00:11.
  - A-T3's earlier pass 1 (20:35) needed about 2 minutes for all 11 SE-600M pools.
- 00:15: operational restart sent (`controller/msgs/r16_D7b_AT3_restart.txt`).
  - Kill pid 34240, re-verify FREEZE_MANIFEST_R16_D4, and relaunch the identical pass-1 command with only `PYTHONUNBUFFERED=1` added, which changes no computation.
  - If it hangs again, A-T3 reports R16_AT3_DEVIATION with a stack sample and changes no code.
- A-T3 relaunched at 00:16 (pid 71974). No code, data or protocol change. Rule until A-T3's pass 1 finishes: no other MPS job runs.

## 2026-11-27 00:20–00:52 — Addendum B (Occam) draft and prior-art refresh: review R1
- 00:20 task sent to the main line (`controller/msgs/r16_ADDB_MAIN.txt`): CPU only, no downloads, no GPU, files only under `round16_addB/`.
- 00:46 `R16_ADDB_DRAFT_READY`: `round16_addB/PREREG_I_ADDENDUM_B_DRAFT.md` (6117709a…) and `round16_addB/PRIOR_ART_I_ADD2.md` (9362e83c…).
- Review R1 (`round16_addB/ADDB_CONTROLLER_REVIEW_R1.md`): the design is sound. Changes are needed before lock:
  - P1: three close works are missing. Molina & Zhang 2026, bioRxiv 10.64898/2026.07.24.740459; CELLAUDIT, arXiv 2609.27234; MDPI Cells 10.3390/cells15181711. One characterization is unsupported ("The Construct Is Part of the Target").
  - P2: a separate venv and pinned code provenance.
  - P3: Gate 5B uses the authors' split assignments; the 621 vs 611 count needs reconciling.
  - P4: a pairwise common-axis sensitivity for the gene-coverage confound.
  - P5: the exact forward definition.
  - P6: an MPS-hang fallback to CPU.
- Molina & Zhang abstract (bioRxiv API, read by the controller): "the perturbation and cell-line specific components are high-dimensional and cannot be recovered from control expression". The full text is behind Cloudflare (403) and is not in Europe PMC (PPR1286577), so it is not read yet. This matters for how new the R16 finding is. It does not change the R16 analysis.
- Nothing about A-T3 or the cross-party check was passed to the main line.

## 2026-11-27 01:03–01:15 — Addendum B v2 reviewed and locked
- 01:03 `R16_ADDB_V2_READY`: `PREREG_I_ADDENDUM_B_DRAFT_v2.md` (537467db…) and `PRIOR_ART_I_ADD2_v2.md` (2e432d47…). P1–P6 are all addressed:
  - Molina & Zhang, COMPASS, CELLAUDIT and Stańczak & Niazi are added. The unverified OpenReview items are removed.
  - Separate venv and code provenance.
  - 621 → 611 explained (min_cells 10). Author splits are regenerated with their blake2b recipe and checked against their 20 hashes, with a tighter tolerance if they match.
  - pairwise_common_axis arm.
  - Per-cell forward, deterministic.
  - MPS stall → CPU rule.
- Controller checks:
  - scanned both files for R16 main-line numbers: none found;
  - read Molina & Zhang's `models/state/predict.py`: it scores saved predictions, with no control swap.
- Locked with two binding clarifications (`round16_addB/PREREG_I_ADDENDUM_B_LOCK.md`):
  - L1: Nadig legacy categoricals are decoded as in D3.
  - L2: the human approval must name the weights, the code archive and the interpreter install.
- Nothing downloaded. Occam waits on the human project lead.

## 2026-11-27 01:22 — Addendum B model-free adapter parts (main line)
- `R16_ADDB_ADAPTER_PREP_DONE`: `round16_addB/src/` holds six files (sha256 listed in MAIN.txt).
- Controller reran `test_unit_addb` with `round16/.venv/bin/python`, without installing anything: 9/9 OK. There are no Occam imports.
- Reported Gate 5B panel on the local HepG2 file (CPU, no model): 611 perturbations, 4,976 controls, 42,402 rows. This matches L1 and v2 §5.2.
- Not committed yet: it will be frozen together with the model-dependent parts after the human approval. Occam is still waiting on the human project lead; R16 is waiting on A-T3.

## 2026-11-27 01:47 — Correction: the 00:15 "MPS hang" was most likely slow GPU work, not a hang
- The relaunched run (pid 71974, unbuffered log) rewrote its barcode files at 01:25:55. It then logged 4 of 11 SE-600M pools by 01:46, about 5 minutes per pool, while its CPU time barely moved.
- That matches the main line's own measurement: SE-600M on MPS float32 takes 331 s per 512 cells (MAIN.txt, synthetic smoke test).
- The controller's 00:15 reasoning had two mistakes:
  - It assumed the 20:35 pass 1 needed about 2 minutes for all embeddings. That was read from the order of lines in a block-buffered log, which is not chronological.
  - It took near-zero CPU time and a single long `waitUntilCompleted` as proof of a hang. During GPU-bound work both are normal. The small MPS probe finishing quickly also proves nothing, because the GPU is time-shared.
- So pid 34240 was most likely still encoding when it was stopped, 20 minutes into the stage, about 4 pools in. The controller's own `r16_d7_embed_test.py` at 22:32 was probably also just slow, queued behind two runs. Neither event shows an MPS deadlock.
- Consequence: about 90 minutes lost to a restart. Correctness is not affected: same code, same command, all pools recomputed from scratch.
- 01:47: told A-T3 not to stop the run, withdrew the 10-minute stall rule, and said to stay on MPS (`controller/msgs/r16_D7c_AT3_keep_running.txt`). MPS and CPU embeddings differ (the main line measured cos 0.9994), so a CPU fallback would break cross-party agreement for STATE.
- Addendum B's P6 rule (10 minutes with no progress → CPU for both parties) still stands for Occam. Its encode is under a second per pool, and both parties switch together. The lock record's reasoning for it, however, came from the misdiagnosis above.
- Revised controller practice: judge a GPU stage by per-unit progress lines, set against the measured per-unit time, never by CPU time.

## 2026-11-27 02:27 — staged verification step 1 passes
- A-T3 SE-600M embeddings vs the main line: all 11 pools bitwise identical, with barcodes equal (DEVIATION_R16_D7.md itself is frozen in FREEZE_MANIFEST_R16_D4 and stays unchanged; this note is the confirmation record). Pools took about 5.4 minutes each on MPS.
- Waiting for `R16_AT3_REAL_DONE` (STATE inference on 5 lines), then step 2 (predicted Δ) and step 3 (metrics).
- Detail of the step-1 result: A-T3's `se600m_embeddings_r16.npz` was written 02:26. `r16_compare_embeddings.py` → R16_EMBEDDINGS_MATCH: 11/11 pools, ordered barcodes equal, shape 512×2058, max |diff| 0, including dld1, feng_ipsc, nadig_hepg2 and nadig_jurkat, which differed before D7. A-T3's `control_barcodes_r16.json` (rewritten 02:26 with SE-600M x_state sha256 fields added) matches the main line's barcodes and batch labels on all 6 primary pools and 5 resample pools.
- 02:30: the controller briefly appended this confirmation to DEVIATION_R16_D7.md. That file is frozen in FREEZE_MANIFEST_R16_D4, so the watcher guard fired, and the controller reverted the append with `git checkout`. D4_FREEZE_OK again.

## 2026-11-27 02:46–02:52 — A-T3 pass 1 done; primary STATE_SE verified across three parties; DEVIATION_R16_D8 (STATE_HVG, Nadig totals)
- `R16_AT3_REAL_DONE` at 02:46.
- Step 2 (`controller/checks/r16_compare_arrays.py`, run by `r16_step23.sh`):
  - STATE_SE predicted Δ is bitwise identical to the main line on all 5 lines, for matched, every mismatched source and resample.
  - Observed Δ, halves, transfer, grand mean and FCT are within 1e-5. Targets and genes are exactly equal.
- Step 3 (controller `r16_check.py` on A-T3's STATE_SE arrays vs the main line's `summary_r16.json`): **R16_MAIN_VS_CONTROLLER_MATCH**, worst |diff| 5.6e-16. That covers every estimate, CI bound, per-source category and the H-I1/H-I2 calls.
- A-T3's own STATE_SE verdict text in A-T3.txt agrees with the main line's to six decimals on all 5 lines.
- **So the primary R16 result (STATE_SE, Groups 1–2) now agrees across main line, A-T3 (blind) and controller.** It is not released to A-T3 or committed until A-T3's secondary work is done, so that verification stays blind.
- STATE_HVG (secondary): Nadig-fed arms differ by 3.7e-4 to 2.8e-3. Cause: A-T3 uses obs UMI_count for the Nadig normalization totals instead of the sum over the file's features. This is DEVIATION_R16_D8 (record in `round16/status/DEVIATION_R16_D8.md`; ruling in `controller/msgs/r16_D8_AT3.txt`). No STATE_HVG line call changes.
- A-T3's pass 2 (pid 38323) is on hold. Next: R16_AT3_D8_PATCHED → STATE_HVG-only pass-1 rerun into `round16/audit/results/hvg_d8/` → pass 2 from scratch.
- 02:56 `R16_AT3_D8_PATCHED`.
  - Controller `git diff`: exactly 1 changed line in `round16/audit/audit_r16_raw_controls.py` (`sel_umis = row_sum_umis`). New sha256 82d7abe5…, which matches A-T3's report.
  - A-T3's py_compile, test_audit_r16 (13/13) and its Nadig self-check (totals = row sums, barcodes unchanged on the hepg2 and jurkat primary and resample pools) all pass.
  - The new code manifest is FREEZE_MANIFEST_R16_D5 (54 entries; D4 plus the new hash plus the D8 record): D5_FREEZE_OK.
  - A-T3's STATE_HVG-only pass-1 rerun started at about 02:56 (pid 43630) into `round16/audit/results/hvg_d8/`.

## 2026-11-27 04:14–04:25 — A-T3 STATE_HVG rerun after D8: verified against the main line
- `R16_AT3_HVG_D8_DONE` at 04:14 (pid 43630, about 75 minutes on MPS, output in `round16/audit/results/hvg_d8/`).
- Step 2 (`r16_compare_arrays.py`, each `hvg_d8/STATE_HVG_*__<line>.npz` vs `round16/results/arrays/<line>/`): **R16_ARRAYS_MATCH on all 5 lines**, 17 checks per line, worst predicted-Δ |diff| 0. Every Nadig-fed arm that differed before D8 (matched, resample, and mismatched arms from a Nadig source) is now bitwise identical.
- Step 3b (new `controller/checks/r16_compare_at3_metrics.py`: A-T3's own metric JSONs vs the main line's `summary_r16.json`, STATE_HVG section): **R16_AT3_METRICS_MATCH**, 490 checks, worst |diff| 1.4e-17. That covers the matched, double-orthogonal matched, model − FCT, per-source paired and resample estimates with all CI bounds, every H-I1 category, stratum summary and line verdict, and the H-I2 call (derived for the main line from its 97.5% CI under the PREREG_I rule).
- Negative control: the same comparator on A-T3's pre-D8 pass-1 HVG files gives R16_AT3_METRICS_DIFF (80 of 490 numeric checks, worst 1.6e-4, all on Nadig-fed quantities; no categorical difference). So the comparator does detect the deviation.
- Correction to the D8 record (DEVIATION_R16_D8.md is frozen in D5 and stays unchanged): its "up to about 5e-5" for the pre-fix metric gaps came from the matched and model − FCT estimates only. Over all HVG quantities the largest pre-fix gap was 1.6e-4, on the nadig_hepg2 within-line resample delta. The conclusion stands: no STATE_HVG call changed, and STATE_SE was unaffected.
- **So STATE_HVG pass 1 (Groups 1–2) now agrees across main line and A-T3 (blind); the controller's own HVG recompute is implied by bitwise-identical arrays.**
- A-T3 moved the stopped pass 2 to `secondary_pre_d8/` and started pass 2 from scratch at about 04:13 (pid 77604, `--models STATE_SE STATE_HVG --include-hvg-perm100 --include-g4`, output `round16/audit/results/secondary/`). Its log `secondary_run.log` is block-buffered (still 0 bytes), so progress is judged from output files and the completion marker `R16_AT3_SECONDARY_DONE`, not from log silence.
- 04:40: pass-2 verification tools prepared (no A-T3 output read yet):
  - `controller/checks/r16_compare_at3_pass2.py` covers the G4 decomposition, the (17, 0) symmetry check, the single-index spread and the HVG Perm100 null;
  - `controller/checks/r16_pass2_check.sh` drives it together with the arrays and HVG-metric comparators.
  - Self-test on synthetic inputs built from the main-line summary: MATCH on faithful copies, DIFF on 2e-5 to 5e-5 perturbations.
  - A-T3's Perm100 uses the same seed rule as the main line (`default_rng(k).permutation(G)`, k = 0..99), so the permutation scores are compared element by element.
- 04:45: a Round 17 preregistration-draft task was dispatched to the main line (`controller/msgs/r17_DRAFT_MAIN.txt`). It is draft and model-free inventory only: no inference, no MPS, no downloads, writes only under `round17/`. The topic follows from R16 secondary diagnostics, so the message file stays uncommitted, and its content unrecorded here, until A-T3's pass 2 is verified.

## 2026-11-27 04:41 — human approval: Occam bundle and GSE295214
- Human project lead: "Approve both" (question: "Approve the Occam bundle and GSE295214 downloads?").
- Records: `round16_addB/R16_ADDB_HUMAN_APPROVAL.txt` (lock L2 scope: Minimal Inference Bundle, archive 9ddb890e, Python 3.8.19 via uv or the §6.1 fallback chain) and `round16/status/R16_GSE295214_HUMAN_APPROVAL.txt` (PREREG_I §2.4).
- Pre-download HEAD checks:
  - GSE295214_RAW.tar Content-Length 8,340,520,960;
  - best_model.pth X-Linked-Size 2,550,066,734, with X-Linked-ETag equal to the locked SHA-256 11d5ba55…;
  - all sidecars resolve at revision acd080b7.
- Downloads started by the controller at 04:42, each into a new empty folder:
  - Occam under `round16_addB/shared/models/` (git-ignored by a new root `.gitignore` rule);
  - GSE295214 under `/Users/mac/.cache/r16_gse295214/` (internal disk, outside the repository).
- Conditions: no MPS for Occam until A-T3 has finished R16 and `R16_ADDB_GO.txt` exists. Group 3 runs only after counts-only QC and its own data freeze.
- 04:50: Round 17 draft received from the main line (`R17_PREREG_DRAFT_READY`). Controller review R1 was sent (`controller/msgs/r17_REVIEW_R1_MAIN.txt`; uncommitted for the same blinding reason as the dispatch).
  - Main ruling affecting R16: gse295214_h1 also becomes an R17 confirmatory line. The R16 Group 3 GO is therefore held until PREREG_J is locked, so that no STATE inference on GSE295214 precedes either preregistration.
  - The R16 Group 3 counts-only QC is model-free and is not held.

## 2026-11-27 04:46–05:00 — downloads verified; TeloHAEC approved for Round 17
- **Occam** (`round16_addB/R16_ADDB_CHECKSUMS.sha256`):
  - best_model.pth matches the locked SHA-256 11d5ba55…;
  - resolved_config.yaml, gene_axis.csv, perturbations.csv, cohort.json, splits_manifest.json and selected_scores_by_seed_split.csv each match the Hub git blob id at acd080b7.
  - Correction to Addendum B v2 §1.3: `ef64f2ba…` is Occam's `gene_axis_sha256` content hash (it appears in splits_manifest.json, and the main line computed it from the local HepG2 file), not the byte hash of gene_axis.csv. The file's byte SHA-256 is a8eb6e4b…, and its blob id matches the Hub.
- **Occam source archive** (0253be5c…, 2,446 entries, no absolute or `..` paths, no symlinks) is not code-only. It also ships `data/`, about 251 MB:
  - two GenePT `.pkl` tables, which are outside the approved Minimal Inference Bundle;
  - five author `splits.npz` files, which §1.3 excludes so that Gate 5B regenerates splits independently.
  - These 7 members were removed from the extracted tree (`round16_addB/R16_ADDB_REMOVED_ARCHIVE_MEMBERS.sha256`). The tarball itself is kept unchanged as the hash-pinned record. The main line was told never to re-extract it or read those members.
- **GSE295214** (`round16/shared/R16_GSE295214_CHECKSUMS.sha256`):
  - the tar is 8,340,520,960 B (sha256 a925b04a…), 85 `.h5` entries, no unsafe paths;
  - extracted to `/Users/mac/.cache/r16_gse295214/raw/`, where all 85 file sizes match GEO `filelist.txt`.
- **TeloHAEC:** Human project lead: "Approve TeloHAEC" (record: `controller/approvals/R17_TELOHAEC_HUMAN_APPROVAL.txt`). The download of Zenodo 18792213 imru_full.h5ad (53.29 GB, CC BY 4.0) started at 04:55 into `/Users/mac/.cache/r17_telohaec_pfizer/`. The md5 check follows, and the file is not opened until PREREG_J v2 is locked; structure-only reads are allowed.
- 05:00: the main line was sent the two CPU-only R16 tasks (`controller/msgs/r16_APPROVED_MAIN.txt`): Group 3 counts-only QC, and the Occam venv, adapters and CPU parts of Gates 5A/5B. No MPS until A-T3 is done.
- 05:10–05:15: main-line `R16_ADDB_DEVIATION`. Occam's unchanged `occam.vcell.pipeline` imports `torch_geometric`, which Addendum B v2 §6.1 wrongly called unnecessary; Occam's `pyproject.toml:48` pins it.
  - Ruling: DEVIATION_R16_ADDB_D1 (`round16_addB/DEVIATION_R16_ADDB_D1.md`). Install `torch-geometric==2.6.1` hash-pinned (1.1 MB wheel, sha256 8faeb353…) into `round16_addB/.venv` only, with no change to any §6.1 pin, and no pyarrow or compiled PyG extensions.
  - The main line's venv is on the primary tier (Python 3.8.19, torch 2.4.1).
- 05:10: R17 v2 accepted in substance (`controller/msgs/r17_V2_ACCEPT_MAIN.txt`). The lock waits for the TeloHAEC download, then a structure-only check.
- The run-phase watcher is paused while the Group 3 QC writes large temporary shard files under `round16/shared/lines/`. Instead, the D5 manifest is checked by hand at each event (05:06: 0 failures, HEAD c4f0028). A marker waiter covers QC, gates, deviations and A-T3.
- 05:16: main-line `R16_ADDB_GATES_CPU_DONE`. Controller check:
  - Venv pins: torch 2.4.1, numpy 1.24.4, scipy 1.9.3, scikit-learn 1.1.3, pandas 2.0.3, anndata 0.9.2, h5py 3.9.0, plus torch-geometric 2.6.1 (D1). pyarrow and the PyG compiled extensions are absent. Import check from `round16_addB/.venv` OK.
  - Gate 5A items 1–2: PASS as reported. Item 1 uses `gene_axis_sha256` as a content hash, which matches. Item 2: CPU determinism and cross-process max_abs 0.
  - Gate 5B, Case A: 20/20 split hashes. The main line's five per-split pearson_delta values equal the authors' `selected_scores_by_seed_split.csv` rows (rhj_fewshot, ABL_N_ALL, seed 42, nocal) to 6 decimals: 0.392850, 0.397797, 0.395306, 0.392278, 0.392659. The mean 0.394178 differs from the target by 2.2e-7; mae and the experimental replicate are also within 3e-7.
  - **Gate 5B PASS confirmed by the controller.**
  - Gate 5A item 3 (MPS vs CPU parity) waits until A-T3 has finished R16. Then come A-T3's independent adapter and gates, the code freeze, and `R16_ADDB_GO`.
- 05:19: main-line `R16_MAIN_G3_QC_DONE`. gse295214_h1 is **PASS_GROUP3_PRIMARY**:
  - N_eval^R16 406 (n_eval with P5 reproducibility, excluding vcc2025_h1: 427), |G_eval^R9| 6,131, |G_eval^R16| 1,837;
  - per-half NTC 28,929 / 29,224;
  - 85 channels (43 / 42 by the GSM-parity half rule); 1,499,727 cells after QC; 892,356 single-guide cells.
  - Controller recomputed the sha256 of all 6 written files (equal to the main line's report). The temporary shard files are gone.
- 05:22: Group 3 data freeze `round16/status/FREEZE_MANIFEST_R16_DATA_G3.txt` (92 entries: tar, 85 .h5, 5 line files, QC record). Its raw hashes equal `R16_GSE295214_CHECKSUMS.sha256`.
  - The R16 Group 3 GO is still held for the PREREG_J lock.
  - A-T3's Group 3 run follows its R16 pass 2, blind, from the raw files.
- 06:16–06:24:
  - TeloHAEC download complete: size exact, md5 equal to Zenodo, sha256 4b22ad11….
  - Structure-only check (`controller/checks/r17/r17_telohaec_struct_check.py`): every column and label the reader resolves exists, and the mixed symbol/Ensembl var index is handled by the frozen HGNC resolver.
  - **PREREG_J locked by hash commitment, commit 9789b0b** (`round17/prereg/PREREG_J_LOCK.md`; verified on origin/main).
  - Then **R16 Group 3 GO** (`round16/status/R16_G3_GO.txt`, `controller/msgs/r16_G3_GO_MAIN.txt`): main line only, `--results-dir round16/results/group3`, MPS concurrent with A-T3 as in the original R16 GO. A-T3's blind Group 3 run comes after its pass 2 is verified.
  - The main line may also run the model-free TeloHAEC counts-only QC (Round 17).
- 06:53 (Round 17): main-line `R17_MAIN_TELOHAEC_QC_DONE`, after the lock, with the locked reader unchanged.
  - pfizer_telohaec_basal PASS_CONFIRMATORY_PRIMARY: N_eval 114, |G_eval^R9| 6,141, per-half NTC 2,414 / 3,087, flask-parity halves.
  - pfizer_telohaec_cytokine (secondary, same family): N_eval 116.
  - So the R17 confirmatory cohort is K_conf = 3: cd4t_stim48 401, gse295214_h1 406, pfizer_telohaec_basal 114.
  - The main line ran the QC with `/opt/anaconda3/bin/python3 -I` (numpy 2.2.3), not `round16/.venv` (numpy 2.4.6); PREREG_J does not fix the interpreter. The controller is rerunning the basal QC in `round16/.venv` into the scratchpad to check that the outputs do not depend on the environment (the P5 energy test uses PCA and permutation p-values).
- Main-line Group 3 run in progress (pid 33705, SE-600M encoding).

## 2026-11-27 07:20 — Addendum B: A-T3 clean-room Occam adapter passes the CPU gates
- `R16_AT3_ADDB_GATES_CPU_DONE`: Case A, 20/20 split hashes; pearson_delta 0.394178 (diff 2.19e-7), mae 0.030157, experimental 0.372174; PASS. MPS is held.
- Independence check (controller): A-T3's `round16_addB/audit/` shares only 7 of its 463 substantive code lines (1.5%) with the main line's `round16_addB/src/`. Those are generic path strings and hash helpers, and A-T3's code never references the main line's files.
- Cross-party Gate 5B parity: mean pearson_delta 0.39417821857779 (A-T3) vs 0.39417821859052 (main), |diff| 1.3e-11; per-split |diff| < 5e-11; experimental replicate identical. This is far inside the Addendum B layer-5 tolerance (1e-5).
- Remaining for Addendum B, in order: both parties' STATE Group 3 runs finish (no concurrent MPS for Occam) → Gate 5A item 3 (MPS vs CPU parity), both parties → code freeze of `round16_addB/src/` and `round16_addB/audit/` → `R16_ADDB_GO.txt`.

## 2026-11-27 07:41–07:47 — Round 17: main-line implementation frozen; GO
- Main line `R17_MAIN_IMPL_READY`: 6 files.
  - Controller rerun: unit tests 17/17 OK. Gate refusal verified: `BlindingGateError` without `R17_CODE_FREEZE.txt`.
  - Synthetic end-to-end run exits 0, and its from-disk recompute differs by 0.
- 12 interpretation points were listed. They are resolved in `round17/prereg/PREREG_J_CLARIFICATIONS_C1.md`, which binds both parties.
  - C1.1: seeds 20261701 / 20261702 per the locked text. The controller's dispatch "20261601" was a controller error and is withdrawn.
- `round17/status/R17_CODE_FREEZE.txt` (15 entries, 0 failures), `round17/status/PREREG_J_COMMIT.txt`, `round17/status/R17_GO.txt`.
  - The main line starts after its R16 Group 3 run finishes; MPS is shared with A-T3.
  - Expected runtime about 5.3 h for 11 lines.
- A-T3 dispatched to its Round 17 clean-room implementation (`controller/msgs/r17_AT3_IMPL.txt`), from PREREG_J + C1 + the array schema only. Its code is frozen separately before it runs.
- 08:09–08:12 (Round 17): A-T3 `R17_AT3_IMPL_READY`, clean-room code under `round17/audit/` (8/8 tests, synthetic end-to-end recompute exact).
  - Independence: 29 of 1,244 substantive lines (2.3%) are shared with `round17/src/`, all schema keys and paths. No references to the main line's modules.
  - The controller compared both codebases on A-T3's 12 interpretation points (no data involved) and wrote `round17/prereg/PREREG_J_CLARIFICATIONS_C2.md`.
    - C2.1 Holm and C2.2 permutation p: A-T3 adapts to the main line's frozen behaviour.
    - C2.3: mac_ctrl's undefined cross-fit is reported as INCONCLUSIVE (VALIDITY_FAIL, HALF_SIZE_LT_2) for both parties. The main line's code has no such guard; C2.3 corrects this at reporting level, and the frozen code is not touched.
    - C2.4: barcodes compared after removing party-specific prefixes.
    - C2.5: the observed-NTC profile definition.
  - Next: A-T3 patches for C2.1/C2.2 (`R17_AT3_C2_PATCHED`) → controller freeze `R17_CODE_FREEZE_AT3.txt` → `R17_GO_AT3.txt`.
- 08:20–08:27 (Round 17): A-T3 `R17_AT3_C2_PATCHED` (Holm changed per C2.1; its permutation p already followed C2.2). Hashes match its report. Controller rerun of the tests: 8/8 OK.
  - Frozen as `round17/status/R17_CODE_FREEZE_AT3.txt` (12 entries, 0 failures). `R17_GO_AT3.txt`: A-T3 starts after its R16 Group 3 run.
  - Both parties' Round 17 code is now frozen before either has touched real cells.
  - A controller subagent is writing `controller/checks/r17_check.py`: independent of both codebases, validated on both parties' synthetic exports.
- 08:36 (Round 17): controller checker `controller/checks/r17_check.py` with driver `r17_check_synthetic.sh` (self-test 42/42).
  - Written by a controller subagent from PREREG_J, C1, C2 and the schema only, with no reading of either party's code. It reuses the R16 controller metric unchanged.
  - On the synthetic exports: main line 22,182 checks and A-T3 20,371 checks, matching values within 5e-14.
  - Disagreements, all settled in `round17/prereg/PREREG_J_CLARIFICATIONS_C3.md`:
    - A-T3's p_Δ used a two-sided count, against C1.5 (C3.1). **A-T3 is on hold to patch and be re-frozen; its run had not started.**
    - The main line's float Spearman misses exact J5b ties; the report uses the controller's exact values (C3.2).
    - Holm tie order: tie-order-dependent rejections are reported as not rejected (C3.3).
    - mac_ctrl override (C3.4 = C2.3).
  - No call changed in the synthetic runs.

## 2026-11-27 09:20–09:30 — Group 3 verified across parties
- Main-line `R16_MAIN_G3_DONE` 09:20 (elapsed 10,491 s). A-T3's blind Group 3 outputs were written at 09:18; its marker is pending, as it is also patching for C3.
- Arrays: STATE_SE and STATE_HVG predicted Δ bitwise identical (worst 0).
- Controller recompute on A-T3's STATE_SE arrays vs main: R16_MAIN_VS_CONTROLLER_MATCH (5.6e-16). A-T3's HVG metrics: 109/109. Pass-2 items: 433/433. Neutral panel (controller): 33/33.
- Verdict drafted in `round16/status/R16_VERDICT_G3.md`. STATE_SE is 0.022 against FCT 0.152, FCT superior under every batch setting, with low signal.
- Not committed until A-T3's Round 17 code is re-frozen after the C3.1 patch, so that no Group 3 number is in view while A-T3 edits code.
- The main-line Round 17 run started at about 09:20 (pid 11669).
- 09:25–09:33: A-T3 `R17_AT3_C3_PATCHED`. The diff is one line (`p_delta` → one-directional helper, `min(1, 2(1+#{d≤0})/(B+1))`) plus tests. Controller rerun of the tests: OK.
  - Controller checker on a fresh A-T3 synthetic export: **R17_CHECK_MATCH** (20,371 checks, 0 fail). The main-line synthetic export has only the 10 expected report-level differences: 6 on mac_ctrl (C2.3) and 4 J5b exact ties (C3.2).
  - Re-frozen as `round17/status/R17_CODE_FREEZE_AT3_v2.txt`; the old file is marked superseded. A-T3 is cleared to start its Round 17 run (its R16 Group 3 run is done and verified).
  - Group 3 verdict and results committed now.
- 09:35 (Addendum B): main-line Gate 5A item 3. max|z_ctrl MPS−CPU| 2.34e-5 > 1e-5 (FAIL); max|Δ| 1.5e-6 ≤ 1e-4.
  - By the pre-locked §5.1.3 rule, **CPU float32 is locked for both parties**, and the CPU Gate 5B stands.
  - Occam therefore never shares MPS with the Round 17 STATE runs.
- Both parties dispatched to implement the full Addendum B evaluation:
  - the main line first writes `round16_addB/src/ADDB_ARRAY_SCHEMA.md`, mirroring the R16 export layout;
  - A-T3 works clean-room from the locked addendum plus that schema;
  - real targets only after `R16_ADDB_CODE_FREEZE*` and `R16_ADDB_GO*`.
- The Round 17 runs are both in progress: main-line pid 11669 from about 09:20, A-T3 pid 15653 from about 09:33.
- 09:40: main-line R16_ADDB_SCHEMA_READY, round16_addB/src/ADDB_ARRAY_SCHEMA.md (sha256 e610a67d79cc2061d18d8ac02443971b91db73d54d5b8ef8b5435b3783596754); A-T3 notified.
- 09:50: main-line R16_ADDB_IMPL_READY: run_eval_addb.py + tests (18/18; controller rerun of test_eval_addb 9/9 OK); synthetic smoke with exact from-disk recompute; CPU estimate ~6.4 h for 5 lines. 8 interpretation points listed; review together with A-T3's when it reports.

## 2026-11-27 09:52–09:55 — Round 17: DEVIATION_R17_D1 (a549_myllia gene axes)
- The main-line Round 17 run crashed at 09:52 in control extraction for a549_myllia: the two Seurat channel exports have different gene axes (18,711 vs 18,808; 18,480 shared). No arrays or results had been written.
- A-T3's frozen code would have used a different rule (ch_1 axis, zero-filled).
- Ruling (`round17/status/DEVIATION_R17_D1.md`): use the R10 common-key intersection, the axis the line's observed side was built on, with per-cell totals as row sums over it.
  - Both parties patch their control path, and A-T3 stops its run.
  - Each party then runs a controls-only pre-check for all 17 pools. The controller compares the two before re-freezing and restarting both runs.
