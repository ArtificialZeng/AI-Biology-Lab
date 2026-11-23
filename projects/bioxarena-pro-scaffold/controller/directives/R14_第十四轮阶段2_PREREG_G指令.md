# 第十四轮阶段 2：GSE337804 预注册（主控 Claude，2026-11-23 09:55 +08:00）

## 已确认
- A-T1 可行性报告（`round14/shared/R14_FEASIBILITY_GSE337804.md`）：可行。218 个靶点（654 条 guide）；20 条 non-targeting（`non-targeting_1..20`）可以和 20 条 gene desert safe-targeting 区分；2 个通道（ch_1、ch_2）；`u6_rep ∩ R9_COMMON = 61`，其中 15 个也是第十二轮 A549 的靶点，46 个是第十二轮没有测过的。
- 主控核对来源：Adduri et al., *Cell* 2026（State），Myllia 生成，GEO 2026-07 公开。**不属于** Virtual Cell Challenge（2026 年挑战用 Arc 自己的 CRISPRi Flex 数据）或 X-Atlas/Orion；也不是 Myllia「Echoes of Silenced Genes」挑战的答案（那个是 CRISPRi，2026-03 已结束）。

## A-T3：起草 `round14/prereg/PREREG_G.md`（PREREG_E 的增补；格式参照 PREREG_F）
1. 论断 1、论断 2、全部检验、阈值、指标、掩码、r9_6 源线、IUT 加 Bonferroni、4,000 次 bootstrap、100 次置换，**一字不改**。新种子 `20261224`。
2. 线：`a549_myllia`（唯一确认性候选线）。
3. **解读限定（写在看数据之前）：** 这是「同一细胞系（A549）、独立实验室、不同化学方法和 guide 的独立数据集复现」，**不是推广到其他细胞系**。通过门槛且两个论断都确认时，只能表述为「在两个独立的 A549 数据集上、两轮分别预注册的前瞻检验中各自确认」；不合并两轮数据。任一论断不确认，就如实写成在独立数据集上没有复现。
4. 细胞和 guide：只用作者提供的矩阵中的**整数原始计数**（X、layers 或 raw 中哪一个是整数原始 UMI，必须在看任何指标之前机械判定并记录；没有整数原始计数 → `EXCLUDED_FEASIBILITY`）；单 guide 细胞 = 作者逐细胞 guide 分配恰好一个（若只有 guide 计数，则按 PREREG_F §3.2 的备用规则）；NTC = 只有 20 条 non-targeting，safe-targeting 细胞排除；≥ 300 UMI 和 ≥ 200 个基因，不另设线粒体上限。
5. 分半和扣除：ch_1 → A 半，ch_2 → B 半；扣除单位 = 通道。
6. 门槛：与 PREREG_E §2.2 相同（n_eval_targets ≥ 50、|G_eval| ≥ 3,000、每半 NTC ≥ 15）。不过门槛 → 不做确认性判定。
7. **事先登记的描述性分析：** 把评测靶点分成「第十二轮 A549 也测过的」和「第十二轮没测过的」两组，分别报告 FCT 和 6 个配对差（只描述，不判定），用来看结果是不是只来自重叠的靶点。
8. 事先披露：在知道第十二轮、第十三轮结果之后设计；A-T1 的估计（98% 比例约 60，84% 比例约 55，保守规则 12–24）；评测靶点上限只有 61，离门槛很近；同一细胞系。
9. 下载方案：只下载 GEO 的两个 `GSE337804_A549_ch_*.RDS.h5ad`（563,808,587 B 和 608,829,330 B）和 sgRNA 库文件；记录字节数和 sha256；放在 `round14/shared/raw/gse337804/`；只用 `python3 -I` 读取。
完成后 `round14/status/A-T3.txt` 最后单独一行写 `PREREG_G_DRAFT_READY`。

## A-T2：准备读取代码（只用合成数据）
在 `round14/src/` 写 GSE337804 的读取代码，按上面第 4、5 条实现（h5ad：自动判定整数原始计数所在位置；逐细胞 guide 分配；20 条 non-targeting 为 NTC，safe-targeting 排除；靶点符号走冻结的 HGNC 表；ch_1/ch_2 分半、按通道扣除），接入冻结的 `process_r10_line`，并实现第 7 条的分组描述（评测代码可复用 round13/src 的做法，种子 20261224，`CODE_FREEZE_R14.txt` 硬门）。只用合成的 h5ad 测试，不下载真实数据，不改 round12/、round13/ 的冻结文件。完成后报告 sha256，`round14/status/A-T2.txt` 最后单独一行写 `R14_CODE_READY`。
