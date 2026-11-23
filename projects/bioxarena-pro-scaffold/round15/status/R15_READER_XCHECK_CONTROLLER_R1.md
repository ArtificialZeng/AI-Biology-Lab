# 主控读取代码交叉核对（第 1 轮）

主控 Claude ｜ 2026-11-23 17:14 +0800

**做法：** 主控自己生成合成 h5ad（A-T2、A-T3 都没见过；AnnData 磁盘格式，h5py 写成；3,800 个 R9_COMMON 基因，30 个源线可重复靶点，每条线 400 个 NTC + 30×60 个扰动细胞 + 低 UMI 和多靶点的杂细胞；4 个 Flex 通道），用 A-T2 的 `prepare_r15_line` 和 A-T3 的 `standardize_vcc_line_to_audit_line` 分别跑，**不设容差**逐个比较两边输出的全部数组。

| 情形 | 应有结果 | A-T2 | A-T3 | 数值数组（18 个） |
|---|---|---|---|---|
| 1 一个文件混三条线，有 cell_line 列和 batch（第一档），含 Jurkat | cal33、hela 读入 | 读入 | 读入 | **逐位相同** |
| 2 一个文件混三条线，只有 context 列，无 batch（第二档），给官方对照表 | 读入 | 读入 | 读入 | **逐位相同** |
| 2' 同上但不给对照表 | EXCLUDED | EXCLUDED | EXCLUDED | — |
| 3 每条线一个文件，无细胞系列，给 manifest（增补 A 规则 c） | 读入 | 读入 | **报错（见问题 2）** | — |
| 3' 同上但不给 manifest | EXCLUDED | EXCLUDED | EXCLUDED | — |
| 4 X 为对数归一化的非整数 | EXCLUDED | EXCLUDED | EXCLUDED | — |

基因集合与顺序、靶点集合与顺序在情形 1、2 都相同；delta、两个半的 delta、n_cells、ad_counts、reproducible、half_valid、energy p 值、敲低深度、对照均值全部 max|diff| = 0。

## 问题 1（A-T2，冻结前必须改）：靶点标签字段错位
- 冻结的 `process_r10_line` 要求 `pert_meta[key] = (ens_id, canon_sym, raw_sym, hgnc_status)`（`round10/shared/scripts/prepare_r10_lines.py` 第 233、1277 行）。
- A-T2 `_classify_vcc_obs_cell` 返回 `(t, canon, ensg, guide_lbl)`（`round15/src/dataset_readers_r15.py` 第 388 行），`_classify_em_guide_key` 第 434 行同样。结果：输出里 `target_ensembl_ids` 是基因符号，`target_raw_symbols` 是 ENSG，`target_hgnc_status` 是 guide 名。
- 第二位（canon_sym）是对的，所以基因、靶点和全部数值不受影响；但保存的线文件元数据是错的，以后任何用 Ensembl ID 或状态字段的分析都会出错。要求改为 `(ensg, canon, t, status)`，status 取 `hgnc.resolve` 的返回值；检查 guide 库备用路径（第 429 行 `guide_to_target`）的元组顺序；单元测试加断言（ENSG 开头、状态为 approved 等）。

## 问题 2（A-T3）：每条线一个文件时，不属于该线的文件被当成错误
- 给 A-T3 传入整个发布目录的文件时，`cal33` 因为 `hela_challenge.h5ad` 里没有 cal33 而被判 EXCLUDED_FEASIBILITY。真实发布是一个目录多个文件，这样会把所有线都排除。
- 要求：某个文件的细胞全部属于其他线（按 a/b/c 规则判定）时，该文件对这条线贡献 0 个细胞并记录在审计里，不报错；只有这条线在所有文件里都是 0 个细胞时才 EXCLUDED。A-T2 现在就是这样做的。

两处改完后，主控用同一批合成文件重跑；全部情形一致才写 CODE_FREEZE_R15。

## 第 2 轮（2026-11-23 17:18 +0800，A-T2 R15_CODE_FIX1_DONE、A-T3 R15_READER_FIX1_DONE 之后）
同一批合成文件重跑：**9/9 通过。** 6 个读入情形（情形 1、2、3 各两条线）两边全部数组和标签字段逐位相同（包括 target_ensembl_ids、target_hgnc_status、target_raw_symbols）；A-T2 的标签现在是 ENSG / approved / 原始符号。3 个排除情形两边都排除。读取部分可以冻结；评测代码的回归检验（重现第十二、十四轮已提交结果）已派给 A-T2 和 A-T3。
