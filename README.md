# AI Biology Lab

Research in quantitative developmental biology, computational biology and mathematical biology. Each project records its evidence, assumptions and review status. Maintained by [Shawn Zeng (Zijian Zeng)](https://github.com/ArtificialZeng).

[中文](README.zh-CN.md) · [All projects](docs/PROJECT_INDEX.md) · [Reproduction](docs/REPRODUCIBILITY.md) · [Reviewer guide](docs/VALIDATION.md)

## Featured results

Scientific scope and strength of evidence determine the order. Independent empirical checks and exact proofs are identified alongside the claims they support.

### 1. Division timing in developing embryos

In public lineages from two independently acquired *Caenorhabditis elegans* imaging archives, the choice of division label changes estimated timing associations. In the 52-embryo RIKEN archive, one matched swap changes the mean slope by **+0.673** (95% CI **0.277 to 1.070**). Another endpoint definition yields an interval spanning zero. The manuscript also shows when the same coarse timing record permits opposite slopes. This is an observational manuscript prepared for specialist review.

[Project](projects/embryo-division-timing/README.md) · [Manuscript v7](projects/embryo-division-timing/v7/manuscript_main.pdf) · [Supplement](projects/embryo-division-timing/v7/supplementary_info.pdf) · [Main figure](projects/embryo-division-timing/v7/figures/Fig1_optical_boundary_control.pdf)

### 2. BioXArena-Pro single-cell benchmark

A computational benchmark with preregistered held-out analyses, calibration checks and comparisons with native baselines. The comparison depends on the metric and evaluation regime. The preprint and audit ledger report both gains and counterevidence.

[Preprint v2.2](projects/bioxarena-pro-scaffold/docs/preprint/BioXArena_xatlas_preprint_v2.2.pdf) · [Audit ledger](projects/bioxarena-pro-scaffold/reports/AUDIT_LEDGER.csv)

### 3. Mixed-mechanism phosphorylation

For specified published model rates, exact algebra certifies equilibrium instability for all finite kinase totals at or above a sufficient bound. A separate interval certificate and persistence argument give an attracting periodic branch for sufficiently large totals and a strictly positive leading correction to its period. The full source conjecture and priority remain open.

[Project and exact proof](projects/mixed-phosphorylation/README.md) · [Periodic branch and period correction](projects/mixed-phosphorylation/validated-cycle/README.md)

## Selected mathematical results

| Result | Scope |
| --- | --- |
| [Sensor placement on leaking trees](projects/tree-sensor-linkages/README.md) | Exact rank and a criterion for generic local identifiability for any nonempty known leak set. |
| [Fixation bounds on regular graphs](projects/two-graph-moran-matching/README.md) | Exact bounds in both type orientations for connected noncomplete regular graphs. |
| [Periodic adaptation](projects/periodic-adaptation/README.md) | An attracting periodic response and a signed mean-deficit identity for a published reversible motif. |
| [Other proofs and counterexamples](docs/PROJECT_INDEX.md) | Project-level statements, exact checks and unresolved cases. |

## Software and data audits

[Paper-Code Auditor](projects/paper-code-auditor/README.md) · [ORD Data Auditor](projects/ord-data-auditor/README.md) · [Protein Metric Auditor](projects/protein-metric-auditor/README.md) · [MEMOTE direction audit](projects/memote-direction-auditor/README.md) · [Full index](docs/PROJECT_INDEX.md)

Project records identify the tested failure, verification method and upstream review status.

## Repository guide

| Directory | Use |
| --- | --- |
| [projects](projects/) | Individual manuscripts, proofs, benchmarks and audit packets |
| [docs](docs/) | Full index, reproduction notes, reviewer guidance and research history |
| [scripts](scripts/) | Reproduction entry points for the phosphorylation study |
| [.github/workflows](.github/workflows/) | Hosted verification workflows |

Project documentation distinguishes observational results, mathematical certificates and benchmark comparisons. [AI use and provenance](docs/AI_USE.md) · [Citation metadata](CITATION.cff) · [License](LICENSE)
