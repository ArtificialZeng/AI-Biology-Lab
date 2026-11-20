# TASK-001: 化学生物学域基线适配与预检工具链对齐

- **责任人**：Muse
- **下发人**：Antigravity (研发总指挥)
- **状态**：Assigned (已分发)
- **优先级**：P0
- **目标任务**：
  1. `chemical-biology/tdc-dti-binding-affinity`
  2. `chemical-biology/solubility`

---

## 一、 方向 (Direction)

在 BioXArena-Pro 体系中，化学生物学域是基准的关键首战。当前历史数据表明：
- `solubility` 基线相对成熟（~74分），考验的是特征工程的鲁棒性与分布外 Bemis-Murcko 骨架泛化能力；
- `tdc-dti-binding-affinity` 是严重丢分点（Codex 仅得 29.61，Claude Code 仅得 32.99），原因在于其采用了严苛的 Double-Cold + 相似度清洗划分，且各数据源（KIBA, DAVIS, BindingDB）存在系统性标度差异。

本任务的核心方向是：**打造首个极速、合规、具备化学领域特异性先验的端到端保底流水线，跑通预检闭环，确保 0 分/格式违规被 100% 杜绝，同时建立坚实的 Double-Cold 泛化基线**。

---

## 二、 执行步骤 (Execution Steps)

请 Muse 严格按以下步骤推进研发：

1. **环境与路径确认**：
   - 目标任务公开数据目录位于：`BioXArena-Pro/data/chemical-biology/`。
   - 依赖项以 `rdkit`, `scikit-learn`, `lightgbm`, `pandas`, `numpy` 为核心，确保在本地与沙盒中快速可执行。

2. **核心流水线实现 (`src/chembio_fast_baseline.py`)**：
   - **对于 `solubility` 任务**：
     - 解析 SMILES，计算 RDKit 2048-bit Morgan 指纹（radius=2）+ 核心物化描述符（MW, LogP, TPSA, 芳香环数等）；
     - 训练经过正则化的 LightGBM / Ridge 回归模型，处理不可溶/完全互溶边界值；
     - 按照 `sample_submission.csv` 格式输出预测，生成合法的 `submission.csv`, `metrics.json`, `solution.py`。
   - **对于 `tdc-dti-binding-affinity` 任务**：
     - 分源建模策略：按 `source`（`kiba`, `bindingdb`, `davis`）独立学习标度与回归模型，严禁混合未校准的尺度；
     - 药物端提取 Morgan 指纹（1024-bit），靶点蛋白端提取标准氨基酸二肽/三肽组成（AAC/DPC）或轻量特征；
     - 结合任务定义的 Calibrated Regression 评价函数，在训练集内执行双冷分组交叉验证，记录 OOD 表现。

3. **执行官方合规自检**：
   - 必须调用 `BioXArena-Pro/training/submission_preflight.py` 中的 `validate_public_submission` 对生成的输出进行完整验证。
   - 必须确保：无 NaN、无 Inf、行数与 ID 严格匹配、字段命名一致。

4. **撰写测试报告并提交 Git**：
   - 将实现的脚本提交到 `projects/bioxarena-pro-scaffold/src/`；
   - 运行日志与验证报告输出到 `projects/bioxarena-pro-scaffold/reports/task_001_report.md`。

---

## 三、 严苛边界 (Hard Boundaries)

作为研发指令，以下边界为绝对红线，触犯即视为失败：

1. **时间与资源边界**：单任务脚本总运行时间（特征提取 + 训练 + 推理）必须**控制在 600 秒（10分钟）以内**，绝不允许卡死或超限。
2. **格式契约边界**：`submission.csv` 的列名、行顺序必须与 `sample_submission.csv` 100% 精确一致，严禁修改第一列 ID。
3. **数据隔离边界**：严禁尝试引入外部爬取的未经允许的数据集，只允许读取当前任务的 `public/` 目录。
4. **文件写入边界**：所有运行产物与中间模型权重必须保存在指定的输出目录中，严禁污染公共输入目录。

---

## 四、 验收标准 (Acceptance Criteria)

指挥官将在收到提交后，按以下标准严格验收：
- [ ] 1. 代码规范：结构清晰、模块化、有完整异常处理与类型提示。
- [ ] 2. 预检验证：两个任务生成的 `submission.csv` 100% 通过 `submission_preflight.py` 校验，返回通过信号。
- [ ] 3. 产物完备：输出目录同时包含合法的 `submission.csv`、`metrics.json` 和 `solution.py`。
- [ ] 4. 速度达标：在测试机器上端到端单任务完成耗时 $\le 10$ 分钟。
- [ ] 5. 内部验证分：`solubility` 内部 CV Calibrated Score $\ge 0.72$；`tdc-dti-binding-affinity` 双冷验证分 $\ge 0.18$。
