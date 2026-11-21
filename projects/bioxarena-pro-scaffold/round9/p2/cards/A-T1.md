# A-T1 数据工程（A 队唯一可写 p2/shared/ 的智能体）
目标：为至少 5 条人类 CRISPRi 细胞系准备统一输入。
- 每条线需要：(a) 非靶向对照的单细胞矩阵（随机抽样 ≤ 20,000 细胞即可）；(b) 每个扰动的伪 bulk Δ（相对对照）；(c) 统一的基因集合与靶点集合。
- 优先级：Replogle 2022 K562-essential、RPE1（Figshare+ 10.25452/figshare.plus.20029387）；Nadig 2025 HepG2、Jurkat（本地已有，见 `src/prepare_nadig2025_cache.py` 和 R7 缓存）；Arc VCC 2025 H1 hESC（公开数据，查许可证）；Nourreddine 2026 KOLF2.1J iPSC（10.1038/s41587-026-03199-w）。拿不到的写明原因。
- 预处理沿用 Addendum B（文库大小归一化、log1p、按对照 z-score 的选项都要保存）。不同线之间用基因 symbol 对齐，记录对齐损失。
输出：`p2/shared/lines/<line>/{ctrl_cells.npz, pseudobulk_delta.npz, genes.txt, targets.txt}`、`p2/shared/DATA_SOURCES.md`、`p2/shared/COMMON_GENES.txt`、`p2/shared/README.md`（每条线的细胞数、靶点数、和其他线的共同靶点数）。完成写 `p2/status/A-T1.txt` 的 TASK_DONE。尽量先交付 Replogle 两条线，让 A-T2 能开始开发。
