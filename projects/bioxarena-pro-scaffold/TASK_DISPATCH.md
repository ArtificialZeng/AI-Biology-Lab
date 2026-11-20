# BioXArena-Pro 协同作战任务看板 (Task Dispatch Board — Round 6)

> **当前定级**：**L2 维持**（Codabench 提交包已于 `2a8b2e0` 冻结，待人类 PI 单次上传评测）  
> **当前轮次**：第六轮（Claude 主控 · 阶段 A & 阶段 B 严格整改与验证）

| 任务编号 | 任务名称 | 状态 | 对应提交 / 脚本 |
|:---:|:---|:---:|:---|
| **T0** | 冻结 Codabench 提交包（独立 commit + SHA-256 + preflight） | **已完成** | `2a8b2e0`, `CODABENCH_FREEZE_MANIFEST.json` |
| **T1** | 第四、五轮 P0 文档全面整改（删 p 值、删夸大术语、不利结果上主表） | **已完成（待回填 T3 行级 CI）** | 5 份核心文档已通过 grep 核验 |
| **T3** | 行级配对自助法 95% CI（4,000 次分层重抽，xatlas + Replogle） | **运行中** | `src/compute_row_bootstrap_ci.py` |
| **T7** | 可复现性加固（缓存脚本、去硬编码路径、依赖锁、SHA-256 与零差值核验） | **已完成** | `0e7d3a7`, `src/prepare_replogle2022_cache.py`, `requirements-lock.txt` |
| **T5** | 本地官方 `ReactomePathways.gmt` 超几何富集（评测全集背景） | **预先声明已锁定，待运行** | `src/round6_t5_local_reactome_enrichment.py` |
| **T6** | Replogle CRISPRi 基因组邻近 / 双向启动子假设检验 | **预先声明已锁定，待运行** | `src/round6_t6_crispri_proximity_test.py` |
| **T2** | 全 5 指标主表（含独立扰动区分度百分位，xatlas + Replogle） | **预先声明已锁定，待运行** | `src/round6_t2_multi_metric_bench.py` |
| **T4** | 双数据集逐模块留一消融（LOMO + 行级 95% CI） | **预先声明已锁定，待运行** | `src/round6_t4_leave_one_out_ablation.py` |
