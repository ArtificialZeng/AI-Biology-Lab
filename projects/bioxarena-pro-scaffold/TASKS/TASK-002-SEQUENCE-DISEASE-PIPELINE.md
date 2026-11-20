# TASK-002: 序列基因组与表型疾病任务域自适应流水线

- **责任人**：Muse
- **下发人**：Antigravity (研发总指挥)
- **状态**：Backlog (待 TASK-001 验收后激活)
- **优先级**：P1
- **目标任务**：
  1. `sequence-genomics/deepstarr-genomic-sequence` (增强子/启动子活性预测)
  2. `sequence-genomics/yeast-promoter-expression` (酵母启动子序列表达调控，历史 Codex 仅 75.9，OpenCode 达 80.4)
  3. `phenotype-disease/clinvar-disease-variant-effect` (致病突变体效应分类)

---

## 一、 方向 (Direction)

在 BioXArena-Pro 的序列基因组学（Sequence-Genomics）与表型疾病（Phenotype-Disease）领域：
- 该类任务以 DNA/RNA 序列、变异位点、以及临床多组学为输入；
- 普遍存在的失分原因是：传统语言模型在处理长 DNA 序列（如 1000bp 增强子）时缺乏生物学 k-mer 与基序（Motif）感知，或者在面对 Manifest 文件格式时打包错误导致评测 0 分（如历史日志中多次出现的 `InvalidNBSubmissionError`）；
- 核心方向：**构建基于 k-mer 词频 + 深度 1D 卷积/轻量 Transformer 的序列编码器，并实现 100% 免疫打包错误的 ZIP/NPZ Manifest 自动合规生成器**。

---

## 二、 执行步骤 (Execution Steps)

1. **序列特征化工程**：
   - 实现快速 3-mer / 4-mer / 5-mer 词频统计与频数归一化；
   - 提取序列 GC 含量、回文对称性与局部低复杂度区域特征；
   - 适配 DNA/RNA 反向互补对称增强（Reverse Complement Data Augmentation）。
2. **模型与训练设计**：
   - 针对回归任务（如 DeepSTARR 与 Yeast Promoter）：采用 multi-output 1D-CNN 或 GBDT 树模型；
   - 针对分类任务（如 ClinVar）：使用加权交叉熵或 Focal Loss 应对良性/致病性样本不均衡。
3. **Manifest 结构强校验**：
   - 编写 `manifest_builder.py` 辅助工具，自动将预测 CSV 封装进 ZIP 根目录，绝对杜绝嵌套一层子目录的致命打包错误（这是导致历史上 DeepSeek/Gemini 多次 0 分的原因）。
4. **内部自测与报告**：
   - 运行并记录酵母启动子表达相关系数与 DeepSTARR Pearson 相关系数。

---

## 三、 严苛边界 (Hard Boundaries)

1. **敏感词审查边界**：涉及人类变异与潜在致病性研究时，提示词与代码注释严格采用纯生物信息学术语，严禁使用任何武器/攻击性敏感隐喻。
2. **ZIP 包规范边界**：若任务为 File-Manifest 格式，ZIP 根目录下必须直接存放各实例预测文件（如 `predictions.csv`），严禁包含包装文件夹或上层路径。
3. **时限红线**：单任务训练 + 推理总时长 $\le 15$ 分钟。

---

## 四、 验收标准 (Acceptance Criteria)

- [ ] 1. 生成的 Manifest ZIP 包通过 `submission_preflight.py` 的解压与完整性校验。
- [ ] 2. 酵母启动子表达验证集 Pearson $r \ge 0.78$。
- [ ] 3. DeepSTARR 验证集 Pearson $r \ge 0.58$。
- [ ] 4. 输出合法的 `submission.csv`、`metrics.json` 与 `solution.py`。
