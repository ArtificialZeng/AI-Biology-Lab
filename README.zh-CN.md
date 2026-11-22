# AI Biology Lab

定量发育生物学、计算生物学与数学生物学研究。各项目列明证据、适用条件和审阅状态。由 [Shawn Zeng（Zijian Zeng）](https://github.com/ArtificialZeng) 维护。

[English](README.md) · [完整项目索引](docs/PROJECT_INDEX.md) · [复现说明](docs/REPRODUCIBILITY.md) · [同行评阅指南](docs/VALIDATION.md)

## 代表成果

按科学问题的影响范围和证据强度排序。独立实证验证和精确数学证明随对应结论列出。

### 1. 胚胎细胞分裂时序

在两套独立采集的秀丽隐杆线虫成像档案中，分裂标注的选取会改变估计出的时序关联。RIKEN 的 52 个胚胎中，一项匹配对照使平均斜率上移 **0.673**，95% 置信区间为 **0.277–1.070**。采用另一种端点定义时，区间跨过零。手稿还展示了同一份粗采样记录何时容许方向相反的斜率。这是供领域同行评阅的观察性研究稿。

[项目说明](projects/embryo-division-timing/README.md) · [v7 主文](projects/embryo-division-timing/v7/manuscript_main.pdf) · [补充材料](projects/embryo-division-timing/v7/supplementary_info.pdf) · [主图](projects/embryo-division-timing/v7/figures/Fig1_optical_boundary_control.pdf)

### 2. BioXArena-Pro 单细胞评测

包含预注册留出检验、校准分析和原生基线比较。比较结果随指标及评测设置变化，预印本与审计台账同时记录正向结果和反向证据。

[预印本 v2.2](projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.2.pdf) · [审计台账](projects/bioxarena-pro-scaffold/reports/AUDIT_LEDGER.csv)

### 3. 混合机制磷酸化

在指定的已发表模型速率下，精确代数证明给出平衡点失稳的充分界，覆盖该界以上的全部有限激酶总量。另有区间认证和持续性论证，给出充分大激酶总量下的吸引周期分支，以及严格为正的一阶周期修正。原猜想的完整结论和成果优先权仍待确认。

[项目与精确证明](projects/mixed-phosphorylation/README.md) · [周期分支与周期修正](projects/mixed-phosphorylation/validated-cycle/README.md)

## 其他数学结果

| 成果 | 结论范围 |
| --- | --- |
| [泄漏树的传感器布置](projects/tree-sensor-linkages/README.md) | 给出任意已知非空泄漏集合下的精确秩与泛型局部可辨识性判据。 |
| [正则图上的固定概率界](projects/two-graph-moran-matching/README.md) | 对连通非完全正则图证明两种类型方向的精确界。 |
| [周期输入下的适应](projects/periodic-adaptation/README.md) | 给出特定可逆反应模型中的吸引周期响应和带符号的平均差恒等式。 |
| [更多证明与反例](docs/PROJECT_INDEX.md) | 各研究包收录问题、精确检查和未决情形。 |

## 软件与数据审计

[论文代码审计](projects/paper-code-auditor/README.md) · [ORD 数据审计](projects/ord-data-auditor/README.md) · [蛋白质指标审计](projects/protein-metric-auditor/README.md) · [MEMOTE 方向审计](projects/memote-direction-auditor/README.md) · [完整索引](docs/PROJECT_INDEX.md)

各项目记录复现的具体问题、验证方法及上游审核状态。

## 仓库导航

| 目录 | 内容 |
| --- | --- |
| [projects](projects/) | 研究手稿、数学证明、评测与审计 |
| [docs](docs/) | 完整索引、复现、评阅说明与研究历史 |
| [scripts](scripts/) | 磷酸化项目的复现入口 |
| [.github/workflows](.github/workflows/) | 托管验证流程 |

各项目说明区分观察性成果、数学证书和评测结果。[AI 使用与溯源](docs/AI_USE.md) · [引用信息](CITATION.cff) · [许可证](LICENSE)
