# TASK-003: 演化式代码搜索与 OOD-CV 严格淘汰机制

- **责任人**：Muse & Antigravity (联合攻关)
- **下发人**：Antigravity (研发总指挥)
- **状态**：Backlog
- **优先级**：P1
- **核心组件**：`src/ood_ensemble_selector.py`

---

## 一、 方向 (Direction)

在 BioXArena-Pro 的实测中，许多 Agent 会产生“训练集过拟合幻觉”——例如用更深的网络或更多的迭代步数在内部小划分上拿到虚高分数，但在真正的测试集分布外（OOD）迁移时遭遇滑铁卢。
本任务的核心方向是：**构建一个自动化的 OOD-CV 验证与模型集成淘汰引擎（Ensemble & Selection Gate），在本地全自动模拟测试集的分布外划分，评估多模型差异性，加权融合多路候选预测，并强制执行“单调提升才允许覆盖提交”的严格守门规则**。

---

## 二、 执行步骤 (Execution Steps)

1. **多模型候选池构建**：
   - 针对任意表格/分子/组学特征矩阵，同时拟合多样化基模型：
     1. $M_1$: 正则化线性/岭回归（Ridge / ElasticNet），擅长外推与高维稀疏特征；
     2. $M_2$: 基于直方图的梯度提升树（HistGradientBoosting / LightGBM），擅长非线性相互作用；
     3. $M_3$: 浅层多层感知机（MLP）或核岭回归，提供不同归纳偏置。
2. **模拟真实 OOD 划分的交叉验证生成器**：
   - 分子任务：调用 RDKit 提取 Bemis-Murcko 骨架，使用 `GroupKFold(groups=scaffolds)`；
   - DTI 任务：使用双冷分组（Drug-Target 均独立留出）；
   - 疾病与细胞系任务：使用 `LeaveOneGroupOut` 按研究队列（Cohort）或细胞系（Cell Line）留出。
3. **加权集成与标度校准**：
   - 基于 OOD-CV 的均方误差或相关系数，计算最优加权系数 $\mathbf{w}^* = \arg\max_{\mathbf{w}} \text{Metric}(\sum_i w_i \hat{y}_i, y_{\text{true}})$；
   - 对预测值进行领域物理边界截断与均值/方差对齐。
4. **单调提升守门机制（Monotonic Gate）**：
   - 比较集成模型的 $\text{CV}_{\text{new}}$ 与已有提交的 $\text{CV}_{\text{current}}$；
   - 仅当 $\text{CV}_{\text{new}} > \text{CV}_{\text{current}} + \epsilon$ 时，才覆盖更新 `submission.csv`、`metrics.json` 与 `solution.py`。

---

## 三、 严苛边界 (Hard Boundaries)

1. **无标签隔离**：任何阶段绝对不得接触 `test.csv` 的隐藏真实标签，集成权重 $\mathbf{w}$ 必须纯粹在训练集的 OOD-CV 折叠预测上优化。
2. **耗时预算**：三路模型训练 + 集成搜寻总耗时不得超过 180 秒（3分钟）。
3. **确定性复现**：所有随机种子统一固定（如 `seed=42`），保证相同数据下完全确定性重现。

---

## 四、 验收标准 (Acceptance Criteria)

- [ ] 1. `src/ood_ensemble_selector.py` 模块化代码完成，支持任意回归与分类任务的自动化调用。
- [ ] 2. 在 `solubility` 任务上，集成后的 5-Fold OOD 分数较单模型 Baseline 进一步提升 $\ge 0.01$。
- [ ] 3. 在 `tdc-dti-binding-affinity` 任务上，集成后的双冷分数提升 $\ge 0.02$。
- [ ] 4. 100% 通过 `submission_preflight.py` 自动化格式检验。
