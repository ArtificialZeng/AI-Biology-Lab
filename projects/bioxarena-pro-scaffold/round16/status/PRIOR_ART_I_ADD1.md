# Round 16 可换对照公开权重扰动模型增补核查报告（`PRIOR_ART_I_ADD1`）

- **智能体**：`A-T5`（A队 · 文献与官方指标）
- **台账依据**：[`TOPJOURNAL_LEDGER.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/controller/TOPJOURNAL_LEDGER.md)（`G3`：除 `Arc STATE` 外，至少补充 1 个可接受控制细胞输入且公开权重的其他扰动预测模型，执行跨细胞系「换对照」对照测试）
- **衔接报告**：[`PRIOR_ART_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PRIOR_ART_I.md) · [`MAIN.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/MAIN.txt)
- **目标评估队列参照**：
  - Group 1：`Replogle K562 / RPE1`、`Nadig HepG2 / Jurkat`、`Feng iPSC`、`VCC 2025 H1`
  - Group 2（[`MAIN.txt`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/MAIN.txt) `PREREG_I` v4）：主验证队列 `{A549 KO, DLD1, Zhu CD4+ T resting}`；`Feng iPSC` 标为 `SENSITIVITY_ONLY`
- **本机硬件与环境实测**（只读 `sysctl` 核实）：`Apple M3 Pro`，物理统一内存 **`38,654,705,664` bytes（36 GiB）**，无 NVIDIA CUDA GPU。
- **合规红线遵守声明**：全程仅使用 HuggingFace API、GitHub raw/API、arXiv API、bioRxiv/EuropePMC API 与本地 `scratch/` 缓存只读查询；未下载任何模型权重或 `>100 MB` 文件，未安装新包，未执行 `git`，未绕过访问限制，未打开任何本机 `.h5ad` 文件。报告中所有受阻或未开源一手文件均明确标为 **`未能读取`**。

---

## 一、核心结论与可换对照候选总览表

在逐一检索与源码核验近 24 个月及经典扰动模型（共 26 个模型/分支）后，结论如下：

1. **真正满足「公开权重 + 免登录非门控 + 原生读取控制细胞 + 可脱离离散细胞系 ID 执行跨系换对照 + M3 Pro 36 GB 单机可推理」的非 `Arc STATE` 独立模型极稀缺**。多数近期顶会/竞赛模型要么**未释放权重**（`X-Cell`、`flowstate`、`OCOO-T`、`CellCanon`、`S1Pert`、`PerturbJEPA`、`PerturbBridge`、`Bison`、`CASTER`、`Wu et al. 2026 scLDM`、`scLAMBDA`、`SAMS-VAE`），要么**推理根本不读取目标对照细胞而是纯条件噪声采样**（`CZI/Biohub scLDM` `biohub/scldm_cd4`）。
2. **`ArcInstitute/PIE`（`PIE_replogle_xdataset`）虽实现了对 Replogle 四系 response label 全持留（`replogle: 0.0`），但工程上隐含 `~104 GB` 原生证据资产下载**（因 `PieDataModule` 在推理阶段回溯训练预处理分片与 11 个知识源），且其网络仅输出基因级伪批量摘要统计 `(p_de, lfc_pred, delta_p_pred)`、并同时消费 `context_text` 词元，不输出单细胞分布。
3. **本报告最终主推 1–2 个真实可执行候选**（均需提前向主控报备 `>100 MB` 权重下载许可）：
   - **主推候选 1（首选独立全转录组模型）**：**`Occam`**（`bendidiihab/Occam`，`checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth`，`2.55 GB`，**`MIT License`**）。单文件自包含、无需外挂知识库、纯 PyTorch CPU 回退、推理无离散细胞系 ID 且用 `FiLM(multiplicative_only, mean_free=True)` 直接交互 `control_expression`。
   - **主推候选 2（轻量级非 STATE 不同架构对照）**：**`scGPT` 微调扰动权重**（`matthewshu/scgpt-replogle-base-ft`，`207.5 MB`，`bowang-lab/scGPT` `MIT` 代码框架；按控制细胞表达进行 MLM 生成推理，无需 FlashAttention）**或** **`TxPert` 第三方复现权重**（`bendidiihab/Occam-competitor-checkpoints/txpert/rhj_fewshot/seed_42/checkpoint.ckpt`，`102.9 MB`，GNN+MLP 对照编码，`Recursion Non-Commercial EULA`）；若仅作 `Arc STATE` 内部「去除 `SE-600M` 编码器先验」消融，可另配 `arcinstitute/ST-HVG-Replogle`（`164.1 MB`）。

### 总览矩阵（五维准入核查）

| 分类 | 模型 / 仓库 | 权重公开与体量 | 许可 / 门控 | 训练集 vs 目标队列重合 | 吃控制细胞 & 换对照可行性 | M3 Pro 36 GB 可行性 |
|---|---|---|---|---|---|---|
| **主推 1** | `Occam` (`bendidiihab/Occam` `rhj_fewshot/seed_42`) | `2.55 GB` 单文件自包含 | `MIT` / 免登录(`gated=False`) | 微调：`HepG2, Jurkat, RPE1`；**未微调 `K562, H1` 与 Group 2**（预训练 Encoder 用 DiverseRNA） | **直接吃 `control_expression`**；无离散细胞系 ID；**原生可换对照**（按 `6,983` 交集轴对齐） | **可行**（纯 PyTorch CPU/MPS，需注意 Py3.8 锁定检查公差） |
| **主推 2A** | `scGPT` (`matthewshu/scgpt-replogle-base-ft`) | `207.5 MB`（`+vocab 1.32 MB`） | 公开 HF / 免登录(`gated=False`)；基座 `scGPT` 为 `MIT` | 微调：`K562, HepG2, Jurkat`；**未微调 `RPE1, H1` 及 Group 2**（预训练 CellXGene `2023-05`） | **直接吃对照细胞表达 `x[:, 0]`**；无细胞系 ID；**原生可换对照** | **可行**（PyTorch CPU/MPS；`flash_attn` 为可选项） |
| **主推 2B** | `TxPert` (`Occam-competitor` `txpert/rhj_fewshot`) | `102.9 MB`（另需 GO/STRING 小图表） | 公开 HF / 免登录；代码为非商 `Recursion EULA` | 重训：`Jurkat, RPE1, HepG2(fewshot)`；**未训 `K562, H1` 与 Group 2** | **直接吃基底对照 `cntr`**（`BasalMLP` + 图邻接）；无细胞系 ID；**可换对照** | **可行**（PyTorch + PyG，`102.9 MB`；需 `torch>=2.6`） |
| **同族对照** | `ST-HVG-Replogle` (`arcinstitute/ST-HVG-Replogle`) | `164.1 MB`/折（无 `SE-600M`） | Arc 非商许可 / 免登录(`gated=False`) | 每折训 3 系，留 1 系；对 `H1` 与 Group 2 零暴露 | 吃 `X_hvg` 对照池；可换对照（但仍受同构 HVG 列对齐限制，且非独立新模型族） | **可行**（`164.1 MB` PyTorch，但属 `STATE` 家族内部变体） |
| **重资产备选** | `PIE` (`arcinstitute/PIE_replogle_xdataset`) | `966.4 MB` 权重（但**原生推理触发 `~104 GB` 数据源**） | Arc 非商 / `CC BY-NC-SA 4.0` / 免登录 | **`replogle: 0.0`**（训于 Tahoe, Orion, Jiang, VCC25 H1） | 吃均值对照 `ctrl_means` **但同时进 `context_text` 文本 token**；仅输出伪批量 `(p_de, lfc, dp)` | **权重可跑，资产极重**（原生脚本强拖 `44 GB` `PIE_sources` + `59 GB` 训练集，需改造旁路） |
| **有重但偏/重训** | `PerturbDiff` (`DeepGraphLearning/PerturbDiff` / `Occam-competitor`) | 官方 `1.63–3.55 GB`；Occam 版 `1.98 GB` | `MIT` / 免登录(`gated=False`) | 官方训 Replogle/Tahoe/PBMC；Occam 版训 `Jurkat, RPE1, HepG2(fewshot)` | 吃 `X_hvg` 基底采样，但**官方配置含 `celltype_encoding: llm` 离散词典查表** | 勉强可行但较重（分布扩散采样；且细胞系编码与外部词典耦合） |
| **有重但受阻** | `CellFlow` (`Occam-competitor` / `matthewshu/cellflow-norman`) | Occam 版 `18.94 GB`；Norman 版 `3.92 GB` | 代码 `PolyForm-NC-1.0.0` / 免登录 | Occam 版训 `Jurkat, RPE1, HepG2`；Norman 仅训 `K562` | JAX OT-FM 吃源分布，**但 `sample_covariates: [cell_type]` 需目标细胞系协变量注册** | **高风险/不推荐**（`18.94 GB` pickle 逼近 36 GB 内存红线，JAX Apple CPU/Metal 兼容性差） |
| **有重但架构不符** | `CZI/Biohub scLDM` (`biohub/scldm_cd4`) | `830.5 MB` | `MIT` / 免登录(`gated=False`) | **直接训练于 `Zhu CD4+ T`**（存在目标队列直接暴露） | **否定**：`predict_step` 排除 `COUNTS` 并从 `randn` 纯离散条件采样，**不吃目标对照细胞** | 技术权重可载入但**逻辑不可做换对照** |
| **有重仅编码器** | `tahoebio/Tahoe-x1` / `tahoebio/Rhaister` | 权重公开 | `Apache-2.0` / 免登录 | Tahoe-100M / DE 汇总表 | `Tahoe-x1` 无扰动解码头；`Rhaister` 仅吃 DE 汇总统计不吃单细胞 | **不支持单细胞换对照** |
| **有代码无官方可用跨系权重 / 未开放** | `ST-HVG/SE-Parse`, `ST-HVG/SE-Tahoe`, `scGPT+ST`, `Tahoe-x1+ST`, `X-Cell`, `flowstate`, `GEARS`, `CPA`, `OCOO-T`, `CellCanon`, `S1Pert`, `PerturbJEPA`, `PerturbBridge`, `Wu et al. 2026 scLDM`, `Bison`, `CASTER`, `scLAMBDA`, `SAMS-VAE`, `BioLORD`, `chemCPA`, `scGen` | 见下文逐条核查 | — | — | 或无公开权重，或仅单品/细胞因子词汇，或强绑单系图/Linux CUDA | **不可作为即用公开权重候选** |

---

## 二、分模型详细核查（按五维标准）

### A. `ArcInstitute/STATE` 其他编码器与变体族

#### 1. `arcinstitute/ST-HVG-Replogle`（去 `SE-600M` 的高变基因空间 `STATE`）
- **来源**：<https://huggingface.co/arcinstitute/ST-HVG-Replogle> · <https://github.com/ArcInstitute/state>
- **原文引用**：
  > "State Transition (ST) models using **X_hvg** inputs on **Replogle**."
- **五维核查**：
  1. **许可、门控与体量**：`arc-institute-model-license`（非商业科研可用），`gated=False` 免登录。含 4 个 `fewshot` 与 4 个 `zeroshot` 模型（如 `zeroshot/k562/checkpoints/final.ckpt` 为 **`164,071,711` bytes ≈ `164.1 MB`**），**无需下载 `SE-600M`**（`config.yaml` 中 `embed_key: X_hvg`, `llm_name: null`）。
  2. **训练集与防泄漏**：每一折训练 Replogle/Nadig 四系中的其余三系并留出一系（如 `zeroshot/k562` 仅训 `hepg2, jurkat, rpe1`）；对 **`VCC 2025 H1`** 及 Group 2（`Zhu CD4 T`、`DLD1`、`A549`、`Feng iPSC`）零训练暴露。
  3. **输入类型**：`_infer.py` 直接读取 `adata.X` / `adata.obsm['X_hvg']` 的对照细胞矩阵 `ctrl_cell_emb`（经 `var_dims.pkl` 对齐固定高变基因维度）；`cell_type_key` 仅用于抽取对照池分组，不入 `forward`；配置含 `batch_encoder: true`, `use_batch_token: false`。
  4. **换对照可行性**：技术上可直接固定权重并置换控制细胞池；**但本质仍与 `ST-SE-Replogle` 同属 `Arc STATE` 同一架构族（不满足 `G3` 独立第三方架构要求）**，且 `X_hvg` 维度锁定在训练集指定 HVG 顺序，跨技术栈基因缺失时零填充风险高。
  5. **M3 Pro 36 GB 可行性**：**可行**。仅 `164.1 MB`，与现有 `STATE` 推理路径一致，CPU 即可快速完成推理。

#### 2. `arcinstitute/ST-HVG-Parse` & `arcinstitute/ST-SE-Parse`
- **来源**：<https://huggingface.co/arcinstitute/ST-HVG-Parse> · <https://huggingface.co/arcinstitute/ST-SE-Parse>
- **原文引用**：
  > "State Transition (ST) models using **X_hvg** inputs on **Parse**."
- **五维核查**：
  1. **许可、门控与体量**：`arc-institute-model-license` 并附带 Parse 数据合规提示（`"Certain uses of a model trained with Parse data may require a license from Parse Biosciences, Inc."`），`gated=False`。`ST-HVG-Parse` `fewshot/split_0..3/checkpoints/best.ckpt` 约 **`540.1 MB`**；`ST-SE-Parse` 约 **`578.6 MB`**（另需 `SE-600M`）。
  2. **训练集与防泄漏**：训练于 Parse 10M PBMC 细胞因子扰动数据集（`B_Intermediate_Memory` 等细胞亚型），未训于 Replogle/Nadig CRISPRi。
  3. **输入与换对照可行性**：虽读取控制细胞 `X_hvg`/`X_state`，但 **`pert_onehot_map.pt` 词表为 PBMC 细胞因子干预，而非基因敲低 CRISPRi 靶标**，无法直接应用于 Group 1/2 基因敲低评估。
  4. **M3 Pro 36 GB 可行性**：`ST-HVG-Parse` 权重本身可在 36 GB CPU 载入，但扰动空间与 CRISPRi 基因敲低不匹配。

#### 3. `arcinstitute/ST-HVG-Tahoe` & `arcinstitute/ST-SE-Tahoe`
- **来源**：<https://huggingface.co/arcinstitute/ST-HVG-Tahoe> · <https://huggingface.co/arcinstitute/ST-SE-Tahoe>
- **原文引用**：
  > "State Transition (ST) models using **X_hvg** inputs on **Tahoe**."
- **五维核查**：
  1. **许可、门控与体量**：`arc-institute-model-license`，`gated=False`。`ST-HVG-Tahoe` 为 `fewshot/.../best.ckpt` (**`1.072 GB`**) 与 `zeroshot/.../final.ckpt` (**`1.084 GB`**)；`ST-SE-Tahoe` 约 **`1.11–1.12 GB`**（另需 `SE-600M`）。
  2. **训练集与防泄漏**：训练于 Tahoe-100M 癌细胞系药物筛选（`C32, HOP62, HepG2-C3A, Hs 766T, PANC-1` 为留出测试系），未见 Replogle/Nadig CRISPRi。
  3. **输入与换对照可行性**：`_infer.py` 默认 **`pert_col: drugname_drugconc`**（小分子化合物-浓度表），无法预测遗传基因敲低 `gene` 列表。
  4. **M3 Pro 36 GB 可行性**：可加载运行但**不适用基因敲低队列**。

#### 4. `scGPT+ST` & `Tahoe-x1+ST`（`Zhang et al. 2026`, `arXiv:2611.05849`）
- **来源**：<https://arxiv.org/abs/2611.05849> · <https://github.com/HowIII/lowrank-decoder-reliability>
- **原文引用**：
  > "The scripts that load STATE checkpoints also need `pip install arc-state==0.10.4`."
- **五维核查**：
  1. **许可、门控与体量**：`HowIII/lowrank-decoder-reliability` 仅发布审计与训练脚本（README 仅引导下载官方 `SE-600M`、`ST-SE-Replogle`、`ST-SE-Tahoe`、`tahoebio/Tahoe-x1` `70m/*` 与 `scGPT_human`）。**未见作者公开发布已训练好的 `scGPT+ST` 或 `Tahoe-x1+ST` 解码器权重（`未能读取`）**。
  2. **训练与换对照/硬件**：需用户自行复训后才存在对应权重；当前无法作为「现成公开权重」直接执行推理。

#### 5. `Occam` 仓库重训 `STATE`（`bendidiihab/Occam-competitor-checkpoints/state`）
- **来源**：<https://huggingface.co/bendidiihab/Occam-competitor-checkpoints/tree/main/state>
- **原文引用**：
  > "Arc Research Institute (“**Licensor**”) is releasing certain materials collectively defined below as the “**State Model**”."
- **五维核查**：
  1. **许可、门控与体量**：遵循同一 `MODEL_LICENSE.md`（非商），`gated=False`。`state/rhj_fewshot/seed_42..46/checkpoint.ckpt` 为 **`555,646,557` bytes ≈ `555.6 MB`**（另有 `norman_combo` 版）。
  2. **训练集与防泄漏**：`rhj_fewshot`（`Jurkat + RPE1` 全量 + `HepG2` 训练折，`6,983` 共同 Ensembl 基因，`embed_key: X_hvg`, `hidden_dim: 336`, `cell_set_len: 64`）；未训 `K562`, `VCC 2025 H1` 与 Group 2。
  3. **输入、换对照与本机**：与官方 `ST-HVG-Replogle` 完全同构（仅多细胞集合长度不同），可换对照、M3 Pro 36 GB 可跑，但**仍属 `Arc STATE` 架构本身**，非独立模型家族。

#### 6. `ArcInstitute/PIE`（`PIE_replogle_xdataset` & `PIE_replogle_wdataset`）
- **来源**：<https://doi.org/10.64898/2026.10.02.756297> · <https://huggingface.co/arcinstitute/PIE_replogle_xdataset> · <https://github.com/ArcInstitute/pie>
- **原文引用**：
  > "Target-context control expression is an inference input, so zero-shot here means no target response-label training or fine-tuning."
- **五维核查**：
  1. **许可、门控与体量**：模型权重采用 `pie-model-license`（非商），代码采用 `CC BY-NC-SA 4.0`，`gated=False`。单文件权重 `PIE_replogle_xdataset/best_auprc.ckpt` 为 **`966,412,013` bytes ≈ `966.4 MB`**（另 `PIE_replogle_wdataset` 四折各约 `950 MB`；注意配套的 `arcinstitute/PIE_eval_datasets` 返回 `HTTP 401 Unauthorized` **`未能读取`**）。**关键大文件预警**：审查 `src/pie/predict.py` 与 `src/pie/data/datamodule.py` 源码发现，`predict_loaded` 会无条件用 `evidence_source=train_data` 初始化 `PieDataModule` 并调用 `setup_stats()`，经 `snapshot_download` 拉取全套训练集预处理分片（`PIE_tahoe100m` `45.14 GB` + `PIE_x_atlas_orion` `11.40 GB` + `PIE_replogle_nadig_essential` `1.98 GB` + `PIE_jiang` `0.41 GB` + `PIE_arc_vcc_25` `0.13 GB` ≈ `59.06 GB`）及全部 11 个知识源 `arcinstitute/PIE_sources`（合计 **`44.00 GB`**，仅 `ncbi_text/embeddings.npy` 就占 `43.18 GB`），**原生推理总拉取量高达 `~104 GB`**。
  2. **训练集与防泄漏**：`PIE_replogle_xdataset` `config.yaml` 设定 `dataset_weights.replogle: 0.0`，仅训练于 `Tahoe (0.68)`、`Orion (0.30, HCT116/HEK293T)`、`Jiang (0.01)` 与 `ARC VCC 25 (0.01)`；因此对 **`Replogle K562 / RPE1`、`Nadig HepG2 / Jurkat` 及 Group 2 (`Zhu CD4 T`、`DLD1`、`A549`、`Feng iPSC`) 的响应标签均保持 0% 暴露**（注：`VCC 2025 H1` 在训练集内）。
  3. **输入与换对照可行性**：架构为 Perceiver IO（非 `STATE` 家族）。模型解码端读取背景均值 **`ctrl_means`** 经 `evidence_ctrl_enc` 与 `evidence_*` 融合；但编码端 `_encode` 仍同时向 KV 序列注入 **`context_text`（Cellosaurus 细胞背景文本嵌入）** 及训练库汇总证据，且输出为基因级伪批量摘要统计 `(p_de, lfc_pred, delta_p_pred)` 而非逐细胞表达矩阵。若做纯对照端换对照，必须在编码器侧对 `context_text` 掩蔽或固定为无源载体，且需绕开原生的训练证据全量下载。
  4. **M3 Pro 36 GB 可行性**：模型本体参数规模小（`d_model=768`, `4+6+4` 层, `inference_chunk_size=2048`），加载后内存可控（`<4 GB`）；**但受制于原生推理管线对 `~104 GB` 远程资产的强绑定与磁盘/下载政策制约，开箱即用性低**。

---

### B. 独立公开权重 / 近年模型与基线详查

#### 7. `Occam`（`IhabBendidi/Occam` · `bendidiihab/Occam`, `arXiv:2611.21959`）— **主推候选 1**
- **来源**：<https://arxiv.org/abs/2611.21959> · <https://github.com/IhabBendidi/Occam> · <https://huggingface.co/bendidiihab/Occam>
- **原文引用**（模型卡）：
  > "No separate knowledge download is needed for inference."
- **五维核查**：
  1. **许可、门控与体量**：代码与发布说明为 **`MIT License`**，`gated=False` 免登录。核心虚拟细胞权重 `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth` 为 **`2,550,066,734` bytes（`2.55 GB`）**（官方提供 `seed_42`–`seed_46`；另发布双基因激活版本 `norman_combo_NF32ALLA` 同为 `2.55 GB`）。单文件内已封装编码器、解码器、GO+STRING+GenePT+ESM2 扰动图模块、FiLM 模块、`44,728` 基因表与归一化元数据，**推理阶段无需再下载图谱或文本嵌入**。
  2. **训练集与防泄漏**：`rhj_fewshot`（`ABL_N_ALL`）监督微调仅使用 `Nadig HepG2`、`Nadig Jurkat` 与 `Replogle RPE1`（三系交集 `6,983` Ensembl 基因，按模型基因表零填充至 `44,728` 维；其中 `Jurkat`、`RPE1` 全扰动入训，`HepG2` 划分训练 1,657 / 验证 103 / 测试 621 扰动）。因此对 **`Replogle K562`、`VCC 2025 H1` 与 Group 2（`Zhu CD4 T`、`DLD1`、`A549`、`Feng iPSC`）零微调暴露**（官方 `official_tables.md` 已报告在 `h1_seen` 与 `k562_seen` 零微调结果）。**须披露限制**：(a) 底层冻结编码器 `occam_small_encoder_decoder.pth` 预训练自 `DiverseRNA-1.4M`（`arXiv:2605.31562`，其全部原始 GEO 样本清单未逐条对照 Group 1/2 排查）；(b) `rhj_fewshot` 在 `6,983` Ensembl 基因轴训练，将新队列全基因直接输入会偏离冻结编码器训练分布；(c) 扰动输入仅支持训练图谱涵盖的 `2,381` 个 Replogle/Nadig 基因（不支持新小分子药物或未入库基因名）。
  3. **输入类型**：`VirtualCellModel.forward(control_expression, pert_indices)` 直接读取单细胞 `control_expression`（按 `target_sum=2000.0`, `log1p` 预处理，支持 `Ensembl ID` 对齐检查 `GeneSpace.prepare`）。**不接收任何离散细胞系 ID 或批次 ID**。
  4. **换对照可行性**：**原生极佳适配**。`configs/rhj_fewshot_official.yaml` 设定 `interaction_type: film`, `film_mode: multiplicative_only`, `film_mean_free: true`，即 $\Delta z$ 显式与控制隐变量 `z_ctrl` 相乘调制后再解码 $\hat{x}=\mathrm{Decoder}(z_{\mathrm{ctrl}}+\Delta z)$；固定权重与靶基因列表，仅将 `controls` 替换为同系 vs 跨系对照细胞即可直接测出模型对控制细胞背景的真实敏感度与相关性塌缩。
  5. **M3 Pro 36 GB 可行性**：**可行**。纯 PyTorch MLP + GNN（无 `flash_attn`/CUDA 专用算子），`src/occam/api.py` 中的 `_pick_device` 在无 CUDA 时自动打印提示并回退 `"cpu"`（亦可在源码层检验 `mps`）。注意 `pyproject.toml` 为复现微米级浮点一致性固定 `requires-python = ">=3.8,<3.9"` 与 `torch==2.4.1`，在实际运行环境中可先做 `verify_checkpoint.py` 数值公差验证再行推理。

#### 8. `TxPert`（`valence-labs/TxPert`, `arXiv:2505.14919`）与 `Occam-competitor` 权重 — **主推候选 2 备选 B**
- **来源**：<https://arxiv.org/abs/2505.14919> · <https://github.com/valence-labs/TxPert> · <https://huggingface.co/bendidiihab/Occam-competitor-checkpoints/tree/main/txpert>
- **原文引用**：
  > "In particular, we present: (i) TxPert, a new state-of-the-art method that leverages multiple biological knowledge networks to predict transcriptional responses under OOD scenarios; (ii) an in-depth analysis demonstrating the impact of graphs, model architecture, and data on performance; and (iii) an expanded benchmarking framework that strengthens evaluation standards for perturbation modeling."
- **五维核查**：
  1. **许可、门控与体量**：官方代码为 `Recursion Non-Commercial End User License Agreement`（`license.pdf`）；`valence-labs/TxPert` GitHub 自身不含 checkpoint 文件，但 **`bendidiihab/Occam-competitor-checkpoints`（`gated=False` 免登录）公开了 `txpert/rhj_fewshot/seed_42..46/checkpoint.ckpt`，单折仅 `102,901,694` bytes ≈ `102.9 MB`**！（推理需配套 `valence-labs/TxPert` 仓库中体积较小的 `data/graphs/go` 与 `string` 表）。
  2. **训练集与防泄漏**：`rhj_fewshot` 版采用仅公开 `GO + STRING` 知识图谱（故意移除官方闭源的 `pxmap/txmap`），训练于 `Jurkat + RPE1` 与 `HepG2` 训练折，**未训 `K562`, `VCC 2025 H1` 及 Group 2 队列**。
  3. **输入与换对照可行性**：`gspp/models/txpert.py` 显示前向接口为 `forward(cntr, pert_idxs, p_emb)`，由 `BasalMLP` 将单细胞对照表达 `cntr` 编码为 `z_intrinsic`，再与图模块 `pert_z` 相加并送入 `decoder` 输出预测（`no_basal_model: false`），**不含离散细胞系标签输入**，技术上**可直接做跨系换对照检验**。
  4. **M3 Pro 36 GB 可行性**：**可行**。仅 `102.9 MB` 权重 + `MLP` 架构（依赖 `torch==2.6`, PyG，无自定义不可替代 CUDA 核；需注意官方 `pyproject.toml` 声明 Python `>=3.12` 与 `module load CUDA/12.4.0` 为群集默认指引）。

#### 9. `scGPT` 扰动微调（`bowang-lab/scGPT` & `matthewshu/scgpt-replogle-base-ft`）— **主推候选 2 备选 A**
- **来源**：<https://github.com/bowang-lab/scGPT> · <https://huggingface.co/matthewshu/scgpt-replogle-base-ft> · <https://huggingface.co/matthewshu/scgpt-replogle-esm-ft>
- **原文引用**（HF 模型卡）：
  > "Fine-tuned the scGPT whole-human pretrained checkpoint on this split with no additional gene prior."
- **五维核查**：
  1. **许可、门控与体量**：`bowang-lab/scGPT` 主仓库采用 **`MIT License`**；官方主仓未直接随包附带 Replogle 全基因组微调权重（原 `scGPT_human` 经 Google Drive 分发），但公开评测仓库 `mattshu0410/sc-interp` 在 HuggingFace（`gated=False` 免登录）托管了 **`matthewshu/scgpt-replogle-base-ft`**（`best_model.pt` **`207,506,241` bytes ≈ `207.5 MB`**，`vocab.json` `1.32 MB`；注 `sc-interp` GitHub 仓库本身未放独立 `LICENSE` 文件）以及增强先验版本 `matthewshu/scgpt-replogle-esm-ft`（`1.46 GB + 1.24 GB ESM2`，另存在 `scGPT-norman-ft` 约 `207 MB`）。
  2. **训练集与防泄漏**：`scgpt-replogle-base-ft` 微调自 `scGPT_human`（`cellxgene_census_human-May23-08-36-2023`，2023 年 5 月前观测语料，**严格早于 Group 2 的 `CD4 T (2026-02)`、`DLD1 (2026-07)`、`A549 (2026-08)`**），扰动监督微调仅在 `State-Replogle-Filtered` 的 **`K562 + Jurkat + HepG2`（1,383 扰动）**上进行；**对 `Replogle RPE1`、`VCC 2025 H1` 以及全部 Group 2 队列零扰动微调暴露**。
  3. **输入与换对照可行性**：源码 `scgpt/model/generation_model.py` 中 `TransformerGenerator.pred_perturb` 直接取每个样本的真实对照基因表达 `ori_gene_values = x[:, 0]` 与基因级扰动标记 `pert_flags = x[:, 1]`，经 Transformer 生成受扰表达 `mlm_output`，**无细胞系 ID 输入**。因此固定权重与靶基因标记、将对照细胞表达 `x[:, 0]` 替换为异系对照池即可执行规范「换对照」检验。
  4. **M3 Pro 36 GB 可行性**：**可行**。仅 **`207.5 MB`**；`bowang-lab/scGPT` README 及 `run_scgpt.py` 明确将 `flash_attn` 设为可选项（并显式忽略缺失警告），可在 PyTorch CPU/MPS 上直接以评估批执行（全基因模式 `max_seq_len=1536` 或按需截断，CPU 单次推理耗时较 GPU 长但内存占用 `<4 GB`）。

#### 10. `PerturbDiff`（`DeepGraphLearning/PerturbDiff`, `arXiv:2602.19685`）
- **来源**：<https://arxiv.org/abs/2602.19685> · <https://github.com/DeepGraphLearning/PerturbDiff> · <https://huggingface.co/katarinayuan/PerturbDiff_release_ckpt> · <https://huggingface.co/bendidiihab/Occam-competitor-checkpoints/tree/main/perturbdiff>
- **原文引用**（`DeepGraphLearning/PerturbDiff` `README.md`）：
  > "PerturbDiff addresses this by moving from cell-level generation to **distribution-level generation**."
- **五维核查**：
  1. **许可、门控与体量**：代码 **`MIT License`**，`gated=False`。官方 `katarinayuan/PerturbDiff_release_ckpt` 权重为 `finetuned_replogle.ckpt` (**`2.98 GB`**)、`from_scratch_replogle.ckpt` (**`1.68 GB`**)、`pretrained.ckpt` (**`3.55 GB`**)；`Occam-competitor-checkpoints` 版 `perturbdiff/rhj_fewshot/seed_42/checkpoint.ckpt` 为 **`1.98 GB`**。
  2. **训练集与防泄漏**：官方 `pretrained` 混训 CellXGene + Tahoe100M + PBMC + Replogle；`Occam-competitor` 版仅训 `Jurkat + RPE1 + HepG2(fewshot)`。
  3. **输入与换对照可行性**：配置含 `basal_mapping_strategy: random` 读取 `X_hvg` 对照，**但官方 `trixie_onehot.yaml` 同时启用了 `celltype_encoding: llm` 与 `batch_encoding: onehot` 并外挂 `new_all_emb.pkl` / `idx_to_pertemb.pkl` 词典映射**。如果跨细胞系预测需要同步提供目标 `cell_line_key` 嵌入，则对照组不纯粹是「仅改变控制细胞表达」，存在协变量混淆风险。
  4. **M3 Pro 36 GB 可行性**：扩散多步采样推理较慢且须额外重建/匹配外部嵌入词典，工程复杂度显著高于 `Occam` 与 `scGPT`。

#### 11. `CellFlow`（`theislab/CellFlow`, `10.1101/2025.04.11.648220`）
- **来源**：<https://doi.org/10.1101/2025.04.11.648220> · <https://github.com/theislab/CellFlow> · <https://huggingface.co/bendidiihab/Occam-competitor-checkpoints/tree/main/cellflow> · <https://huggingface.co/matthewshu/cellflow-norman>
- **原文引用**：
  > "CellFlow is a framework for predicting single-cell phenotypes induced by complex perturbations."
- **五维核查**：
  1. **许可、门控与体量**：框架许可为 **`PolyForm-Noncommercial-1.0.0`**。官方仅开源库未放跨系大模型权重；第三方权重有 `Occam-competitor/cellflow/rhj_fewshot/seed_42/CellFlow.pkl`（**`18,936,497,131` bytes ≈ `18.94 GB`**！）及 `matthewshu/cellflow-norman`（**`3.92 GB`**，仅 Norman 单基因激活 K562）。
  2. **训练集、输入与换对照**：OT-FM 流匹配读取源分布（对照细胞），但 `config.yaml` 配置 `sample_covariates: [cell_type]`、`split_covariates: [cell_type]` 与 PCA-256 空间，即显式将目标 **`cell_type` 协变量作为离散条件嵌入输入**，限制了跨全新未注册细胞系的纯粹控制替换。
  3. **M3 Pro 36 GB 可行性**：**极不推荐**。`18.94 GB` JAX pickle 反序列化与 OTT-JAX 图编译峰值极其接近 36 GB 系统内存上限，且 JAX 在 Apple Silicon 上缺少官方稳定加速后端。

#### 12. `CZI / Biohub scLDM`（`czi-ai/scldm`, `arXiv:2511.02986`, `biohub/scldm_cd4`）与 `Wu et al. 2026 scLDM`
- **重要同名模型区分**：文献存在两套完全独立、均名 `scLDM` 的模型：
  - **(a) CZI / Biohub `scLDM`**：<https://arxiv.org/abs/2511.02986> · <https://github.com/czi-ai/scldm> · <https://huggingface.co/biohub/scldm_cd4>
    - **原文引用**：
      > "`scLDM` is a latent-diffusion model consisting of a novel fully transformer-based VAE architecture for exchangeable data that uses a single set of fixed-size, permutation-invariant latent variables."
    - **核查结果（明确排除）**：`MIT License`，`gated=False`，`830,489,233` bytes（`830.5 MB`）。然而两项硬伤使之**完全不可用于本轮换对照**：(1) 直接基于细胞因子刺激的 **`Zhu CD4+ T`**（存在目标队列直接暴露）；(2) 源码 `LatentDiffusion.predict_step` 明确设置 `exclude_keys = (COUNTS, GENES, COUNTS_SUBSET, GENES_SUBSET, LIBRARY_SIZE)` 并在 `sample` 中以 `z = torch.randn(...)` 纯高斯先验结合离散标签生成，**推理端根本不接受控制细胞表达作为输入**。
  - **(b) `Wu et al. 2026 scLDM`**：<https://doi.org/10.1093/bioinformatics/btag727>（`PMC13661255`）· <https://github.com/samrogers1233/scLDM> · <https://doi.org/10.5281/zenodo.22658165>
    - **原文引用**：
      > "scLDM couples a VAE that maps gene-expression space into a compact latent with a conditional diffusion backbone that learns the control-to-perturbed distribution shift directly in latent space."
    - **核查结果（无预训练权重）**：`MIT License`。模型原理虽以控制细胞潜向量 $z_{\mathrm{ctrl}}$ 与 OT 配对（`gene_sample.py`），但经核验 GitHub 文件树（86 个源码文件，0 个 `.pt/.ckpt/.pth`）与 Zenodo `22658165` 记录仅含代码快照，**未公开发布预训练权重（`未能读取`）**。

#### 13. `GEARS`（`snap-stanford/GEARS`）
- **来源**：<https://github.com/snap-stanford/GEARS> · <https://huggingface.co/matthewshu/gears-norman>
- **原文引用**：
  > "GEARS is currently not designed to handle training across multiple cell types or cross-cell type transfer of predictions"
- **五维核查**：
  1. **许可与体量**：`MIT`；官方库无通用跨系权重；第三方 `matthewshu/gears-norman` 仅 `7.5 MB`（训于 Norman K562）。
  2. **输入与换对照**：架构设计明确面向单细胞系内未见扰动组合推理（取单批控制均值 `ctrl_adata` 为起点叠加图效应），官方明确声明**不支持跨细胞系迁移**，不满足 `G3` 跨系换对照要求。

#### 14. `CPA`（`theislab/cpa`）
- **来源**：<https://github.com/theislab/cpa>
- **原文引用**：
  > "`CPA` is a framework to learn the effects of perturbations at the single-cell level."
- **五维核查**：
  1. **许可与体量**：`MIT`。仅提供训练代码，**无通用跨细胞系/Replogle 公开预训练权重**。
  2. **输入与换对照**：`predict()` 虽支持传入对照组表达作反事实解码，但训练与推理需同时绑定已知离散细胞系/批次协变量映射并逐库重新拟合。

#### 15. `X-Cell`（`Xaira-Therapeutics/X-Cell`）
- **来源**：<https://huggingface.co/Xaira-Therapeutics/X-Cell>
- **原文引用**（`README.md` 与 `MODEL_CARD.md`）：
  > "**Status: Model weights and inference code coming soon.**"
- **五维核查**：`CC BY-NC-SA 4.0`，`gated=False`，HF 仓目前仅约 `0.3 MB` 文本与概览图文件（3 个文件），**无任何已发布权重或推理代码（`未能读取`）**；宣称训于 `X-Atlas/Pisces`（含 7 套 CRISPRi 屏幕），待后续开源后方可实测。

#### 16. `VCC 2026 #2 flowstate`（Khodaveisi & Bagheri）
- **来源**：<https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up>
- **原文引用**：
  > "Together, these observations indicate that current models recover roughly what a cross-line average would predict, whereas the part of the response that is specific to a cell line, which is precisely what a virtual cell model would need in order to improve on an experiment in another model system, remains largely missing."
- **五维核查**：公开摘要信息表明其基于 222 套公开数据集 across 86 contexts 训练，但截至核查时全网 GitHub、HuggingFace 与 Zenodo **均未检索到该竞赛版本的公开代码仓库与模型权重（`未能读取`）**。

#### 17. `OCOO-T`（`arXiv:2606.12838`）
- **来源**：<https://arxiv.org/abs/2606.12838>
- **原文引用**：
  > "This paper introduces OCOO-T, a minimalist flow-matching-based AIVC model for transcriptional perturbation response prediction."
- **五维核查**：截至核查时，arXiv 元数据、GitHub 与 HuggingFace **未见任何官方开源代码或公开模型权重（`未能读取`）**。其架构采用流匹配去噪 + 上下文词元（含对比 `S1` 显式细胞系嵌入与 `S2–S7` 自身对照均值条件，详见 [`PRIOR_ART_I.md`](file:///Users/mac/AI-Biology-Lab/projects/bioxarena-pro-scaffold/round16/status/PRIOR_ART_I.md)）。

#### 18. `CellCanon`（`ZhaokangLiang123/CellCanon`, `arXiv:2611.05890`, HF dataset `X1Agu/CellCanon`）
- **来源**：<https://github.com/ZhaokangLiang123/CellCanon> · <https://huggingface.co/datasets/X1Agu/CellCanon>
- **原文引用**：
  > "CellCanon combines biological knowledge and cellular observations to predict perturbation responses, explain source-supported hypotheses, and update predictions after new measurements."
- **五维核查**：GitHub 代码为 `MIT License`，但 **HuggingFace `models/X1Agu/CellCanon` 不存在（`401 Unauthorized`），`datasets/X1Agu/CellCanon` 仅公开三套 SciPlex3 处理数据包（353 MB + 119 MB + 126 KB），无预训练模型权重发布（`未能读取`）**；且安装说明书（`environments/README.md`）要求 Python 3.11 + Linux `x86_64` + NVIDIA CUDA 13.0，无法在 macOS M3 Pro 直接运行全栈。

#### 19. `S1Pert`（`liujiawen-jpg/S1Pert`, `arXiv:2611.17705`）
- **来源**：<https://arxiv.org/abs/2611.17705> · <https://github.com/liujiawen-jpg/S1Pert>
- **原文引用**：
  > "S1Pert predicts **Up / Unchanged / Down** for a query gene under a genetic intervention in one cell line."
- **五维核查**：仓库（`NOASSERTION`）公开代码并依赖 HuggingFace 通用决策大模型骨架 `Mapika/decider-4b`（`Apache-2.0`），但**未公开发布已训练好的 `S1Pert` 扰动解码器权重**；任务聚焦单细胞系内上/无变化/下调三分类（非全连续表达谱输出），且 README 明确要求 Linux 与 CUDA 12.8 环境（`decider-4b` 4B 参数在本地重训/全量端到端推理亦远超轻量化标准）。

#### 20. `PerturbJEPA`（`BGI-Hangzhou-OmicsAI/PerturbJEPA`, `arXiv:2611.20575`）
- **来源**：<https://arxiv.org/abs/2611.20575> · <https://github.com/BGI-Hangzhou-OmicsAI/PerturbJEPA>
- **原文引用**：
  > "Trained with separate drug and dose ranking objectives, the model predicts candidate responses from a common control population and matches them to the observed response without cell pairing or expression reconstruction."
- **五维核查**：`MIT License`，仓库仅含 Parse10M Tiny 复现代码与配置（51 个代码/配置，**0 个模型权重**）；且设计目标明确声明**不做基因表达重构**，仅做潜空间药物/剂量检索，无法产出基因空间扰动表达供 `r_orth_cal` / DE 指标评估。

#### 21. `PerturbBridge`（两套同名：`echQoQ/PerturbBridge` 与 `yeung-lab/PerturbBridge`）
- **来源**：<https://github.com/echQoQ/PerturbBridge>（转录组条件潜空间 Schrödinger Bridge）· <https://github.com/yeung-lab/PerturbBridge> / <https://arxiv.org/abs/2610.17066>（显微图像表型生成）
- **原文引用**（`echQoQ/PerturbBridge` README）：
  > "An NVIDIA GPU with CUDA support is required for the training pipeline."
- **五维核查**：转录组版本 `echQoQ/PerturbBridge` 为 `MIT License`，**无预训练权重（`weight_files=[]`）**且依赖 PyTorch 1.13.1 + Python 3.7.12 + NVIDIA CUDA；成像版本 `yeung-lab/PerturbBridge` 用于 BBBC021/JUMP 显微图像，无开源权重亦不面向 scRNA-seq。

#### 22. `Bison`（`arXiv:2609.32467`）
- **来源**：<https://arxiv.org/abs/2609.32467>
- **原文引用**：
  > "To exploit complementary chemical supervision while targeting these differences, we introduce Bison: a shared gene representation connects native panels, while two discrete diffusion models compose context-dependent responses with molecular deviations learned through matched drug-contrast supervision."
- **五维核查**：面向 8 套小分子化合物屏幕的未见化合物预测（非 CRISPRi 基因扰动），截至核查时未检索到独立官方公开模型权重与推理包（`未能读取`）。

#### 23. `CASTER`（`Kawsar et al. 2026`, `10.64898/2026.11.15.766815`）
- **来源**：<https://doi.org/10.64898/2026.11.15.766815>（相关 `CanPert`）
- **原文引用**：
  > "CASTER led only at 2.5% screened (0.62 versus 0.47) and fell below the transferred mean without a screen."
- **五维核查**：上下文自适应转移模型，目标场景依赖在目标细胞系抽取 ≥2.5% 扰动样本做适配，未见 HuggingFace/GitHub 独立公开可推理权重发布。

#### 24. `tahoebio/Tahoe-x1` & `tahoebio/Rhaister`
- **来源**：<https://huggingface.co/tahoebio/Tahoe-x1> · <https://huggingface.co/tahoebio/Rhaister>
- **原文引用**：
  - `Tahoe-x1`：> "Tahoe-x1 is a family of perturbation-trained single-cell foundation models with up to 3 billion parameters, developed by Tahoe Therapeutics."
  - `Rhaister`：> "Rhaister predicts perturbation responses in unseen (cell, perturbation) combinations from aggregated summary statistics."
- **五维核查**：二者均为 `Apache-2.0`、`gated=False`；但 `Tahoe-x1` 仅为表示编码器、无独立下游未见扰动预测器，而 `Rhaister` 仅读取差分分析摘要统计（`log2fc` / Mann-Whitney U 检验 p 值 / delta 均值）做高速线性推断，**均无法读取原始单细胞对照池以执行单细胞换对照实验**。

#### 25. 经典与近端方法补充核查（`scLAMBDA`, `SAMS-VAE`, `BioLORD`, `chemCPA`, `scGen`）
- **`scLAMBDA`**：<https://github.com/gefeiwang/scLAMBDA>（`GPL-3.0`，`10.1101/2024.12.04.626878`）。GitHub 仅含代码与 4 系/Norman 示例 notebook，**无预训练权重**；原文引用：
  > "Its disentangled representation learning framework enables the modeling of single-cell level perturbations by separating basal cell representations from salient representations associated with perturbation states, allowing for single-cell level generation and *in silico* prediction of perturbation responses."
  *(注：HF 检索到的 `YuhanQiu/*_scLAMBDA` 系文本大语言模型微调产物，非单细胞扰动权重。)*
- **`SAMS-VAE`**：<https://github.com/insitro/sams-vae>（NeurIPS 2023）。代码仓**未声明 GitHub 标准开源许可（`license: None`）且无预训练权重**；原文引用：
  > "The results in the paper were generated using the Linux environment."
- **`BioLORD`**：<https://github.com/nitzanlab/biolord>（`BSD-3-Clause`）。单数据集解耦框架、**无单细胞扰动预训练权重**；原文引用：
  > "A deep generative framework for disentangling known and unknown attributes in single-cell data."
  *(注：HF 上 `FremyCompany/BioLORD-2023` 等同名热门仓为生物医学概念文本嵌入模型，与单细胞扰动模型无关。)*
- **`chemCPA`**：<https://github.com/theislab/chemCPA>（`MIT`）。专为化学小分子化合物设计，**无预训练权重与通用 CRISPRi 遗传扰动能力**；原文引用：
  > "Code accompanying the [NeurIPS 2022 paper](https://neurips.cc/virtual/2022/poster/53227) ([PDF](https://openreview.net/pdf?id=vRrFVHxFiXJ))."
- **`scGen`**：<https://github.com/theislab/scgen>（`GPL-3.0`，*Nature Methods* 2019）。按单细胞系内/两条件潜在向量差做线性迁移，**无跨千靶点全局预训练权重**；原文引用：
  > "scGen is a generative model to predict single-cell perturbation response across cell types, studies and species [(Nature Methods, 2019)](https://www.nature.com/articles/s41592-019-0494-8)."

---

## 三、明确推荐方案（1–2 个可执行候选及行动前报备说明）

基于**非 `Arc STATE` 独立架构、公开权重且免登录、输入严格吃对照细胞、技术支持跨细胞系固定权重换对照、以及 Apple M3 Pro 36 GB 单机可推理**五大硬约束，向主控推荐如下 **1 主 1 副** 执行方案：

### 推荐 1（首选主推）：`Occam`（`bendidiihab/Occam` · `checkpoints/occam_seeds/rhj_fewshot/all_sources/seed_42/best_model.pth`, `2.55 GB`）
1. **选择理由**：
   - 唯一兼具 **`MIT` 许可 + 免登录 + 单文件自包含（无需下载额外知识库/预训练语料）+ 独立 MLP+GNN+FiLM 架构 + 直接 `control_expression` 乘法交互** 的成熟公开虚拟细胞权重。
   - 前向计算零依赖离散细胞系 ID，可最纯粹地执行「固定 `pert_indices`（在 2,381 靶基因交集内）与权重，比较 matched control vs. mismatched control 预测迁移表现」。
   - **清晰可评估的零微调防泄漏靶标**：该模型只在 `Nadig HepG2/Jurkat` 与 `Replogle RPE1` 上做扰动监督微调，**完全未微调过 `Replogle K562` 与 `VCC 2025 H1`**（且天然独立于 `Zhu CD4 T`、`DLD1`、`A549` 扰动标签，仅需以规范化方式标注底层编码器来自 `DiverseRNA-1.4M` 且按官方 `6,983` 共同基因轴评估）。

### 推荐 2（次推轻量架构对照）：二选一配对
- **选项 2A：`scGPT` 扰动微调权重**（`matthewshu/scgpt-replogle-base-ft`，`207.5 MB`，`bowang-lab/scGPT` `MIT` 代码架构）：
  - **优点**：Transformer 自注意力架构（与 `Occam` 的 MLP/FiLM 完全正交），体量小（`207.5 MB`），无需 `flash_attn`，直接吃对照表达向量 `x[:, 0]` 且不输入细胞系标签；微调只用 `K562+HepG2+Jurkat`（留出 **`Replogle RPE1`**、**`VCC 2025 H1`**，以及预训练截止期 `2023-05` 之后的 **`CD4 T / DLD1 / A549`**）。
  - **注意**：系社区公开研究库 `sc-interp` 托管的第三方微调权重点，`sc-interp` GitHub 自身缺少显式 `LICENSE` 文件，作为次要跨模型交叉验证参考极稳妥。
- **选项 2B（同预算最小物理体积候选）：`TxPert` 第三方复现权重**（`bendidiihab/Occam-competitor-checkpoints/txpert/rhj_fewshot/seed_42/checkpoint.ckpt`，`102.9 MB`，`Recursion Non-Commercial EULA`）：
  - **优点**：与 `Occam` 共享同一套冻结 `rhj_fewshot`（6,983 Ensembl 基因，未见 `K562` / `VCC H1`），仅 `102.9 MB`，GNN + `BasalMLP` 架构；适合与 `Occam` 做严格同分割对比（需从 `valence-labs/TxPert` 仓库拉取 GO+STRING 小型图定义文件，并配置 PyG 环境）。
- *（补充对照参考：`arcinstitute/ST-HVG-Replogle` `164.1 MB` 可视作 `STATE` 家族内部仅去除 `SE-600M` 的同构消融；而 `ArcInstitute/PIE_replogle_xdataset` 因原生推理触发 `~104 GB` 训练数据与 `PIE_sources` 资产下载且不输出单细胞分布，除非主控批准专项裁剪证据管线，否则不列入优先执行。）*

### ⚠️ 重要合规报备（致主控 Claude）
上述**所有可用预训练扰动预测权重均超过 `100 MB` 阈值**：
- `Occam`（`seed_42/best_model.pth`）：**`2.55 GB`**（`2,550,066,734` B）
- `scGPT`（`scgpt-replogle-base-ft/best_model.pt`）：**`207.5 MB`**（`207,506,241` B）
- `TxPert`（`Occam-competitor/txpert/rhj_fewshot/seed_42/checkpoint.ckpt`）：**`102.9 MB`**（`102,901,694` B）
- `ST-HVG-Replogle`（`zeroshot/k562/checkpoints/final.ckpt`）：**`164.1 MB`**（`164,071,711` B）

**当前阶段严守红线，未下载上述任何权重文件**。后续阶段一旦启动 `G3` 跨模型实机推理与「换对照」实验，请主控按工作流程审批确认对应 `>100 MB` 权重下载额度并指派对应沙箱路径。

---

## 四、未能读取的一手资料透明清单（`未能读取`）

1. `arcinstitute/PIE_eval_datasets`：HuggingFace API 返回 `HTTP 401 Unauthorized` — **`未能读取`**。
2. `X-Cell` 权重与推理代码：HF 与 GitHub 官方库明确处于 `"coming soon"` 状态（无文件可读取） — **`未能读取`**。
3. `flowstate`（VCC 2026 #2）代码与权重、`OCOO-T`（`arXiv:2606.12838`）代码与权重、`Bison`（`arXiv:2609.32467`）权重、`CASTER`（`10.64898/2026.11.15.766815`）预训练跨系权重、`scGPT+ST` / `Tahoe-x1+ST` 解码器权重（`HowIII/lowrank-decoder-reliability`）、`Wu et al. 2026 scLDM` 权重：公开通道均未见已上传文件 — **`未能读取`**。
4. `CellCanon`（`huggingface.co/models/X1Agu/CellCanon` 返回 `HTTP 401`，仅有数据集仓 `datasets/X1Agu/CellCanon`）模型权重 — **`未能读取`**。
5. `CellFlow`（`10.1101/2025.04.11.648220`）bioRxiv 正文及 `PerturbDiff` 引证摘要补充外部网页（EPMC 查询返回超时的个别分支） — 通过各自官方 GitHub README/YAML、HF API 与 `DeepGraphLearning/PerturbDiff` 仓库一手核实，未作臆测。
