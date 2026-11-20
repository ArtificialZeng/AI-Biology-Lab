# BioXArena-Pro Autonomous Research Project

> **Project Lead**: Antigravity  
> **Current Stage**: **Level 2 Confirmed Breakthrough** (`xatlas-orion-perturb`: `PanPerturb-Bridge`, 10-split `+0.0485 ± 0.0085`, `2.43×` official `0.02` noise gate, $p = 3.63 \times 10^{-8}$)  
> **Primary Breakthrough Target**: `perturbation/xatlas-orion-perturb`  
> **Disciplines**: Official grader only (D1), reproducible code for every number (D2), full progressive ablation with exact test-model parity (D3), preflight ≠ score (D4), zero leakage (D5), honest staging (D6), strict D7 comparability (same dataset, same metric, same inference-time availability).

---

## 核心成果速览 (`PanPerturb-Bridge`)

- **Benchmark A（10 分割官方评分器递进消融 `S0→S7`）**：`0.0613 ± 0.0040` → **`0.1099 ± 0.0087`**（`+0.0485 ± 0.0085`, `10/10` 正增益, $p = 3.63 \times 10^{-8}$）；其中 5 个从未触碰的确认分割均值 **`0.1064 ± 0.0069`**（相对 S0 `+0.0443`）。
- **Benchmark B（10 分割超严格双细胞系联合留出 + 重拟合 Baseline & Axis 审计）**：`0.0658 ± 0.0033` → **`0.1115 ± 0.0089`**（`+0.0457 ± 0.0081`, `10/10` 正增益, $p = 3.77 \times 10^{-8}$；Untouched Held-5 `0.1097`）。
- **官方 455 行提交产物**：[outputs/perturbation_baselines/xatlas-orion-perturb/submission.csv](outputs/perturbation_baselines/xatlas-orion-perturb/submission.csv) 与 455 个 `.npz` 通过官方 `submission_preflight` 100% 校验。

## 快速导航
- **[作战总纲 (Mission Charter)](MISSION_CHARTER.md)**：分级科研标准（L0–L4）、研究方向与 D1–D7 硬性边界。
- **[任务分发主看板 (Task Dispatch)](TASK_DISPATCH.md)**：任务追踪、验收红线与状态看板。
- **[成果验收登记表 (Acceptance Ledger)](ACCEPTANCE.md)**：8 阶段递进消融、超严格重拟合审计、翻译器 5 折 CV 与阶段定级。
- **[主算法源码 (src/xatlas_panhuman_solver.py)](src/xatlas_panhuman_solver.py)**：`PanPerturb-Bridge` 统一求解器与全量审计脚本。
- **[结构化报表 (reports/)](reports/)**：`xatlas_10split_ablation.csv`、`xatlas_strict_union_refit_audit.csv`、`xatlas_translator_5fold_cv.json`、`unified_evaluation_summary.csv`。
- **[学术论文草稿 (docs/PAPER_DRAFT.md)](docs/PAPER_DRAFT.md)**：`PanPerturb-Bridge` 完整方法论、公式与严谨实验分析。
