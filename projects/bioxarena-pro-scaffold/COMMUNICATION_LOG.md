# BioXArena-Pro 协同作战通信日志 (Communication Log)

> 本日志记录项目推进、方向调整、多轮独立审计纠偏与里程碑验证的完整历史。

---

## [2026-11-20 16:05 - 23:44] 初始部署、第三轮纠偏 (`59a7753`) 与核心模型开发 (`c479bc6`, `a9ab830`)

- 彻底废止历史 `gdsc` 0.1942（官方实测 `0.1199`）与早期夸大表述。
- 完成 `PanPerturbBridgeModel` 开发并在 xatlas 与 Replogle 2022 上进行初步验证。

---

## [2026-11-21 02:35 - 03:10] 第六轮阶段 A 完成（T0, T1, T3, T7）与 Codabench 官方单次私测归档

1. **T0 独立冻结提交（`2a8b2e0`）与官方 Codabench 单次实测**：
   - 人类 PI 授权下由主控上传冻结包一次（`Submission ID 1062221`，配置 `panperturb-bridge__human-ai-collab`，存档 `gemini/Codabench官方结果_submission1062221.html`）。
   - 官方私有测试集实测得分：**`0.0998`**（较 22 个 2 小时自主单次运行 agent 最佳 `0.0825` 高 `+0.0173`，位于官方 `0.02` 噪声阈值内，定性为「噪声范围内领先」；实测验证→测试偏移 `−0.0101`；人机协作多轮迭代与 2 小时自主 agent 非同一比较类别）。严禁据此改模型或二次提交。
2. **T7 可复现性加固（`0e7d3a7`）**：彻底清除 `src/*.py` 硬编码绝对路径，提交 `prepare_replogle2022_cache.py` 与 `requirements-lock.txt`；验证路径改造后 Seed 42 (`0.070151 / 0.126020`) 与 455 测试 `.npz` 最大差值为 `0.0`，临时目录复算 6 个 Replogle 缓存文件 SHA-256 全部逐位一致。
3. **T3 行级配对自助法 95% CI（4,000 次分层重抽）**：完成 xatlas（BenchA + BenchB）与 Replogle（含简化 GNN 逐行对比）全部 10 分割行级 CI 计算并回填至所有文档。关键发现：Replogle 上 `S5 vs AE`（10/10 CI > 0）与 `S5 vs GNN`（7/10 CI > 0）显著稳健，而 `S7 vs GNN` 仅 3/10 分割 CI > 0、`S6 vs S5` 在 9/10 分割上显著为负（`−0.0109`, `[-0.0169, -0.0050]`）。
4. **T1 文档全面纠偏**：5 份核心文档全部完成 T1 替换表与行级 CI 回填，并提交 `reports/ROUND6_PHASE_A.md`（末行 `PHASE_A_DONE`）。阶段 B（T2/T4/T5/T6）继续推进。
