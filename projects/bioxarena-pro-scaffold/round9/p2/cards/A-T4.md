# A-T4 前瞻押题准备（写 p2/vcc_linespecific/src/prospective/）
1. 查清 VCC 2026 测试集（AC16、MOLM13、U-87MG 等）的对照细胞数据和靶点列表怎么获得：是否需要注册账号、条款是否允许赛后研究使用和公开预测、真值预计何时公开（要官方原文）。**不要注册，不要接受条款。** 把人类需要做的步骤写清楚，放进 `p2/status/A-T4.txt` 和 `prospective/HUMAN_STEPS.md`。
2. 写好端到端脚本：输入对照细胞 h5ad 加靶点列表 → 输出 VCC 2026 提交格式的预测 → 计算 sha256 → 生成时间戳说明。先用 Nadig 或 Replogle 的对照数据模拟跑通。
3. 写清楚 cell-eval2 官方评测怎么在本地复现（仓库 ArcInstitute/cell-eval2），用 VCC 2025 公开数据验证评测脚本能跑。
