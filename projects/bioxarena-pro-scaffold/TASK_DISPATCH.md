# BioXArena-Pro 协同作战任务看板 (Task Dispatch Board — Round 6 完成)

> **当前定级**：**L2 维持**（Codabench 官方单次私测 `Submission 1062221 = 0.0998`，属 `0.02` 噪声范围内领先；人机协作多轮迭代与 2 小时自主 agent 非同一比较类别）  
> **第六轮总报告**：[reports/ROUND6_EXECUTION_REPORT.md](reports/ROUND6_EXECUTION_REPORT.md)（`ROUND6_DONE`）

| 任务编号 | 任务名称 | 状态 | 核心产物 / 结论 |
|:---:|:---|:---:|:---|
| **T0** | 冻结 Codabench 提交包（独立 commit + SHA-256 + preflight） | **已完成** | `2a8b2e0`, `Submission 1062221 = 0.0998` |
| **T1** | 第四、五轮 P0 文档全面整改（删 p 值、删夸大术语、不利结果上主表） | **已完成（主控核验通过）** | 5 份核心文档全部合规，已复验 grep |
| **T2** | 全 5 指标主表（含原始 Δ 与 v5 空间扰动区分度 + 行级 95% CI） | **已完成** | `t2_all_5metric_summary.json`, `t2_discrim_delta_row_ci_summary.json` |
| **T3** | 行级配对自助法 95% CI（4,000 次分层重抽，xatlas + Replogle） | **已完成（主控核验通过）** | `row_bootstrap_ci_summary.json` |
| **T4** | 双数据集逐模块留一消融（LOMO + 行级 95% CI） | **已完成** | `t4_lomo_summary.json`（确认 4 跨集核心贡献，`S6` 外部强负 `−0.0153`） |
| **T5** | 本地官方 `ReactomePathways.gmt` 超几何富集（评测全集背景） | **已完成（主控核验通过）** | `t5_local_reactome_enrichment_summary.json`, `t5_enrichment_inputs.json` |
| **T6** | Replogle CRISPRi 基因组邻近 / 双向启动子假设检验 | **已完成（主控核验通过）** | `t6_replogle_crispri_proximity_summary.json`（负 z `< 10 kb` 474×；xatlas 无法排除） |
| **T7** | 可复现性加固（缓存脚本、去硬编码路径、依赖锁、SHA-256 与零差值核验） | **已完成（主控核验通过）** | `0e7d3a7`, `prepare_replogle2022_cache.py`, `requirements-lock.txt` |
