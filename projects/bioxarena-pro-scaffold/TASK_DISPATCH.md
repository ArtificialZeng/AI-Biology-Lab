# BioXArena-Pro 协同作战任务看板 (Task Dispatch Board)

> **指挥官**：Antigravity  
> **执行员**：Muse  
> **当前战役状态**：进行中 (Phase 1 启动)  
> **同步规则**：所有任务由指挥官在此看板分配；Muse 领取任务后按任务卡逐步实施，完成后在任务卡中填写交付物并提交 git commit，通知指挥官验收。

---

## 一、 任务总览与状态流转

| 任务编号 | 任务名称 | 负责成员 | 优先级 | 当前状态 | 对应任务卡 | 验收结论 |
|:---:|:---|:---:|:---:|:---:|:---|:---:|
| **TASK-001** | 化学生物学域基线适配与预检工具链对齐 | **Muse** | **P0 (紧急)** | **Assigned (待执行)** | [`TASKS/TASK-001-CHEMBIO-PIPELINE.md`](TASKS/TASK-001-CHEMBIO-PIPELINE.md) | 待提交 |
| **TASK-002** | 任务模态路由与自适应特征生成器开发 | Muse | P1 | Backlog | 待排期 | - |
| **TASK-003** | 演化式代码搜索与 OOD-CV 严格淘汰机制 | Muse | P1 | Backlog | 待排期 | - |
| **TASK-004** | 8大任务域 64 任务全量流水线批处理与压测 | Muse | P2 | Backlog | 待排期 | - |

---

## 二、 当前正在执行的任务 (Active Task)

### 👉 [TASK-001] 化学生物学域基线适配与预检工具链对齐
- **指派人**：Antigravity
- **指派对象**：Muse
- **详细指令**：见 [`TASKS/TASK-001-CHEMBIO-PIPELINE.md`](TASKS/TASK-001-CHEMBIO-PIPELINE.md)
- **交付要求**：
  1. 完成 `chemical-biology/tdc-dti-binding-affinity` 与 `chemical-biology/solubility` 的自动化保底求解脚本。
  2. 脚本必须在 5 分钟内完成特征提取、模型训练与推理。
  3. 生成的 `submission.csv` 必须 100% 通过官方 `submission_preflight.py` 严格校验，杜绝任何 NaN/Inf/格式错误。
  4. 提交测试日志与验证报告。

---

## 三、 指挥与验收流程规范

1. **领取任务**：Muse 查阅对应 `TASKS/TASK-xxx.md`，确认明确理解【方向】、【执行步骤】与【严苛边界】。
2. **代码提交**：
   - 源码放在 `projects/bioxarena-pro-scaffold/src/` 或相应模块目录下。
   - 提交 commit 信息格式：`feat(task-xxx): <简洁描述> (by Muse)`。
3. **发起验收**：Muse 在 `ACCEPTANCE.md` 的待验收区填写本次交付物清单、运行耗时及评估指标，并在讨论区通知指挥官。
4. **指挥官审核**：指挥官在 1 小时内执行代码审计与复测，并在 `ACCEPTANCE.md` 给出【通过】或【打回修改及改进意见】。
