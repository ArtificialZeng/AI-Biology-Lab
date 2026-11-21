# VCC 2026 前瞻押题：人类操作与合规说明（HUMAN_STEPS）

**更新日期**：2026-11-21  
**核查智能体**：A-T4（只读公开网页与公开 API；未注册账号、未接受任何条款、未向竞赛提交）

---

## 1. VCC 2026 测试集（`AC16`、`MOLM13`、`U-87MG`）对照数据与靶点列表如何获取

### 1.1 官方数据集结构与匿名化说明
- 官方页面 `https://virtualcellchallenge.org/datasets` 原文：
  > "For the 2026 Challenge we performed Perturb-seq experiments on six cell lines originating from different tissues, using CRISPRi and 10x Flex scRNA-seq. The six lines are split into two panels of three: three cell lines for the validation round, and three cell lines reserved for the final test. Each line is released only as an anonymized cell context — labelled A, B and C for validation, and D, E and F for the final test. We do not publish which cell line a context is. What you receive instead is that context’s unperturbed expression profile, which is the only thing your model needs to condition on."
- 结项公告 `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up`（2026-11-19）公布了六个细胞系身份：
  > "The validation lines, Cal-33 (head and neck squamous cell carcinoma), HeLa (cervical adenocarcinoma), and Jurkat (T-cell leukemia), drove the live leaderboard, while the test lines, AC16 (cardiomyocyte), MOLM13 (acute myeloid leukemia), and U-87MG (glioblastoma), determined the final ranking. The test panel included 300 knockdowns, 129 of which were shared with the validation panel and 171 of which were new."
- **最终测试阶段压缩包内容**（`vcc_2026_controls_final.zip`，约 **615 MiB**，2026-10-22 起对参赛者开放）：
  - `context_D.h5ad`、`context_E.h5ad`、`context_F.h5ad`：每条线 18,400 个非靶向对照细胞（46 个非靶向 sgRNA × 400 细胞，中位 ~20,000 UMIs/细胞）× 18,533 基因的原始计数矩阵（`target_gene == "non-targeting"`，且 `ntc_id` 标识对照引导序列）。
  - `gene_names.csv`：18,533 个基因 symbol（含表头 `gene_name`），规定提交矩阵的严格列顺序。
  - `pert_counts.csv`：测试集需预测的 300 个靶基因列表（`target_gene` 列）。
  - `manifest.json`：记录 `contexts`（`["D", "E", "F"]`）、`partition` 与 `panel_id`。

### 1.2 获取方式与账号注册门禁（核心阻塞）
- **需要注册账号并登录（且含 Persona 人脸自拍防重验证）**：
  - `https://virtualcellchallenge.org/datasets` 原文：
    > "Competitors download the control data once they register and are logged into the application, either from the website or with the command-line tool. The final phase bundle becomes available on October 22, 2026, from the Download Final Datasets button in the app or with `vcc datasets download controls-final`. It has the same structure for contexts D, E and F, its own `pert_counts.csv`, and downloads as `vcc_2026_controls_final.zip` (about 615 MiB), under its own name."
  - `https://virtualcellchallenge.org/rules` 原文：
    > "To register for the Challenge you must: Confirm that you are eligible to enter; Agree that you will follow all of these rules; Agree to the Terms of Use..."
    > "You must upload the final entry to the Challenge website by 11:59pm UTC on November 5, 2026."
  - `https://virtualcellchallenge.org/faq` 原文：
    > "Why do I need to provide a selfie? To keep the competition fair. A selfie check confirms you are a real, live person and also is designed to confirm no one registers more than once, which the rules do not permit."
- **截至 2026-11-21 公共存储桶实测为空**：
  - 查询公开 GCS JSON API `https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/2026` 返回空对象（`{'kind': 'storage#objects'}`，无文件与子前缀）。因此目前**无法免账号匿名下载**测试集对照与靶点列表。

---

## 2. 使用条款（Terms of Use / Rules / FAQ）对赛后研究与公开预测的约束

### 2.1 官方原文核查（限制条款 vs 鼓励科研条款）
1. **严格限制条款**（`https://virtualcellchallenge.org/terms`，Effective Date: August 4, 2026）：
   - > "The Content is protected by copyright and other intellectual property laws, and forms part of our Services, even if you download it. Your use of our Content is subject to these Terms, including the rules and restrictions in the following section."
   - > "You agree that you will not do any of the following: ... **Use your account for any purpose other than participating in the Challenge;**"
   - > "**Use, copy, reproduce, modify, translate, publish, broadcast, transmit, distribute, perform, upload, display, license, sell, create derivative works based on, or otherwise exploit for any purpose any of our Services except as expressly permitted by these Terms;**"
   - > "We may also terminate your access to our Services after the Challenge has concluded or been terminated."
   - > "Termination of your account or your withdrawal from the Challenge does not affect any rights, obligations, or licenses that by their nature are intended to survive..."
2. **规则优先顺序条款**（`https://virtualcellchallenge.org/rules`，Revised September 16, 2026）：
   - > "The Arbitration Agreement and Terms of Use form part of the binding contract governing participation in the Challenge. If there is any conflict between these rules and any other guidance about the Challenge, these rules control."
3. **FAQ 鼓励学术发表条款**（`https://virtualcellchallenge.org/faq`）：
   - > "**I'm a researcher or student. Can I publish my results?** We encourage publishing and sharing research with the scientific community, if that is allowed by your employer or academic institution."
   - > "**Will Arc own my model?** No. Arc does not obtain any rights to your model by virtue of your participation in the Challenge."

### 2.2 合规操作结论（给人类研究者的硬性建议）
- **智能体侧**：按 `P2_RULES.md` 第 6 条，智能体**绝对不注册账号、不接受条款、不点击登录、不绕过门禁**。
- **情景 A（人类尚未在 2026-11-05 截稿前注册/下载 `vcc_2026_controls_final.zip`）**：
  - 赛事已于 **2026-11-05 23:59 UTC** 截稿并于 **2026-11-19** 公布获奖者；`/terms` 明确禁止 *"Use your account for any purpose other than participating in the Challenge"*。
  - **推荐稳妥路径**：**不要现在为赛后研究而注册新账号**；直接等待 Arc 按 2025 年先例将完整 VCC 2026 数据以 **CC0 1.0** 公开发布到 `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2026/`，或先发邮件至 **`help@virtualcellchallenge.org`** 书面询问赛后仅用于科研目的获取 `controls-final` 与提前发布预测哈希的许可。
- **情景 B（人类本人已于 2026-11-05 前作为合规参赛者下载过官方 `vcc_2026_controls_final.zip`）**：
  - **严禁**公开再分发从参赛门户下载的原始文件（`context_D/E/F.h5ad`、`gene_names.csv`、`pert_counts.csv`），因其受 `/terms` "Content ... forms part of our Services, even if you download it" 约束。
  - 虽 `/faq` 鼓励科研人员发表研究结果且 Arc 不主张模型所有权，但 `/terms` 字面禁止参赛以外用途及衍生品再发布且 `/rules` 约定规则/条款优先于其他指引。因此**在公开发布基于受限 `controls-final` 生成的预测文件（或哪怕仅公开预测矩阵本身）前，建议先发邮件至 `help@virtualcellchallenge.org` 取得书面确认**；在等待确认期间可先在本地离线运行推理与格式校验。

---

## 3. 真值预计何时公开（官方原文与 2025 先例实证）

- **2026 官方是否给出明确日历日期**：**没有**。经逐字检查 `https://arcinstitute.org/news/virtual-cell-challenge-2026-wrap-up`（2026-11-19）、`https://arcinstitute.org/news/behind-the-data-virtual-cell-challenge-2026`（2026-10-01）、`https://virtualcellchallenge.org/{datasets,rules,terms,faq,evaluation}` 以及 `ArcInstitute/cell-eval2` 仓库，**官方均未写出 VCC 2026 测试集真值的具体公开日历日期**。Wrap-Up 仅在 *"Looking ahead"* 末尾写道：
  > "The Virtual Cell Challenge will return next summer. In the meantime, you can: ... Explore our open datasets in the Virtual Cell Atlas"
- **VCC 2025 公开先例（已通过 Arc 新闻时间戳与 GCS JSON API `timeCreated` 直接核实）**：
  - VCC 2025 结项公告 `https://arcinstitute.org/news/virtual-cell-challenge-2025-wrap-up` 发布日期：**2025-12-06**。
  - VCC 2025 完整数据集（`adata_Training.h5ad`、`adata_Validation.h5ad`、`adata_Test.h5ad`、`gene_names.csv`、`pert_counts_*.csv`）上传至公开桶 `gs://arc-institute-virtual-cell-atlas/virtual-cell-challenge/2025/` 的 `timeCreated` 时间：**2025-12-16（`17:08:04Z` – `17:36:22Z`）**，即**结项公告后第 10 天**；并在 `https://arcinstitute.org/tools/virtualcellatlas` 标注许可证为 **CC0 1.0**：
    > "The Atlas includes the training, validation, and held-out test perturbation datasets leveraged throughout the Challenge. ... © 2025. This work is openly licensed via CC0 1.0."
- **基于 2025 先例的推断窗口（非官方承诺）**：
  - 2026 Wrap-Up 发布于 **2026-11-19**。若沿用 2025 年"结项后约 10 天至数周内释放至 Arc Virtual Cell Atlas（CC0 1.0）"的节奏，预计公开时间窗口在 **2026 年 11 月底至 2026 年 12 月中旬**（仍存在延迟或同时释出 controls 与 GT 的可能）。

---

## 4. 人类逐步操作清单（当拿到合规输入或真值公开后）

### 步骤 1：监测公共释放 / 合规准备输入目录
- **方式 A（推荐，无需登录监测 GCS `2026/`）**：
  ```bash
  curl -sL "https://storage.googleapis.com/storage/v1/b/arc-institute-virtual-cell-atlas/o?prefix=virtual-cell-challenge/2026"
  ```
  当返回不再为空时，检查其 `items`（如是否先公开 `controls` 或直接公开完整 `test` 真值），再按 `DATA_SOURCES.md` 登记 URL/许可证/大小/sha256 后下载。
- **方式 B（仅当您赛前已合规持有 `vcc_2026_controls_final.zip` 且已确认条款允许）**：
  如果此前通过官方 `vcc-cli` 下载（仅供参考命令）：
  ```bash
  # uv tool install vcc-cli
  # vcc login --token-stdin
  # vcc datasets download controls-final -d ~/vcc-final
  unzip -o -j ~/vcc-final/vcc_2026_controls_final.zip "context_*.h5ad" "gene_names.csv" "pert_counts.csv" "manifest.json" -d ~/vcc-final
  ```
  **注意**：切勿将 `vcc_2026_controls_final.zip` 及其解压文件提交到任何公开 git 仓库或转发给他人。

### 步骤 2：运行端到端前瞻预测脚本（本地生成合规预测 + SHA-256 + 时间戳说明）
```bash
python3 p2/vcc_linespecific/src/prospective/prospective_predict.py \
  --bundle-dir ~/vcc-final \
  --source-lines-dir p2/shared/lines \
  --method Hybrid_ProjShared \
  --outdir p2/vcc_linespecific/results/prospective_vcc2026 \
  --seed 42
```
输出文件包括：
- `prediction.h5ad`（严格符合 VCC 2026 格式：360,000 × 18,533 稀疏 CSR 非负整数 raw counts，列序严格对齐 `gene_names.csv`，`obs["context"]` 为 `D/E/F`，`obs["target_gene"]` 为 300 靶点各 400 细胞，不含 `non-targeting`）
- `prediction.prep.vcc`（若安装 `zstandard` / `vcc-cli` 可直接打包并通过 `vcc prep --dry-run` 校验）
- `SHA256SUMS.txt`、`prediction_manifest.json`、`TIMESTAMP_DECLARATION.md`

### 步骤 3：可选本地 `vcc prep --dry-run` 交叉校验（若您已自装 `vcc-cli`）
```bash
vcc prep p2/vcc_linespecific/results/prospective_vcc2026/prediction.h5ad \
  -g ~/vcc-final/gene_names.csv \
  --perts ~/vcc-final/pert_counts.csv \
  --dry-run
```

### 步骤 4：真值公开前锁定先验时间戳（GitHub + 可选 OpenTimestamps）
1. 在确认不违反 Arc 条款（且不上传任何 Arc 受限原始输入文件）的前提下，将**预注册文件、模型代码、以及 `SHA256SUMS.txt`、`prediction_manifest.json`、`TIMESTAMP_DECLARATION.md`**（仅哈希与声明，不含受限数据；若条款允许公开预测结果可一并附上 Release artifact）提交并推送至**公开 GitHub 仓库**（也可附加 `ots stamp SHA256SUMS.txt` 生成比特币区块链 OpenTimestamps 证明 `.ots`）。
2. 记录公开 GitHub commit SHA 与时间戳 URL。

### 步骤 5：真值在 GCS CC0 公开后的客观本地评测（`cell-eval2`）
真值释放后（预期含每个 context 的真实 `.h5ad`），先用附带脚本把真实 `non-targeting` 对照拼接到提交版 `prediction.h5ad` 上（还原平台内部 `control_source: real` 所需的 `validate_pair` 对齐输入），再用 `cell-eval2 --preset vcc2026` 本地复现官方 6 项指标与 0–1 基准/重复锚点缩放分（详见 `CELL_EVAL2_REPRO.md`）。
