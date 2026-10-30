# Bounded open-problem research — 26 October 2026

We screened fifteen leads and advanced two mathematically precise targets that could be checked with exact algebra on a CPU. This is a report of work performed, not a calibrated probability of finding a theorem, obtaining journal acceptance, or receiving media coverage.

## Results advanced for review

1. [Directed-cycle identifiability](../projects/directed-cycle-identifiability/README.md): a positive-integer counterexample to the output-leak endpoint in Ahmed et al.'s Conjecture 4.6, with an explicit nonzero Jacobian minor and rational inverse. A separate general tangent-kernel argument proves the strictly interior-leak inclusion. The source's main theorem is unaffected.
2. [Dual-futile-cycle symmetry slice](../projects/dual-futile-symmetry/README.md): an exact necessary-and-sufficient stability threshold, including its zero and unstable regimes, for an eight-parameter equivariant subfamily. This leaves Vassena's unrestricted sixteen-parameter Conjecture 3.1 unresolved. Priority for the symmetry corollary remains unconfirmed.

Both packages contain a mathematical statement, credited primary sources, executable exact checks, deliberately corrupted controls, and explicit limitations. A separate AI context rechecks the mathematics; that does not constitute external human review. We request upstream author assessment rather than asserting author endorsement.

## Screening decisions

The table records bounded selection decisions. “New question” means our subsidiary formulation, not a source-authored open conjecture. An old GitHub OPEN label alone does not show that current code remains unfixed.

| Lead | Primary source | Decision |
| --- | --- | --- |
| Directed-cycle singular-locus Conjecture 4.6 | [Ahmed et al., §4](https://arxiv.org/html/2412.05283v3#S4.SS1) | Exact endpoint counterexample and restricted repair advanced |
| Noncollocated catenary positive singular sets | [Ahmed et al., §5](https://arxiv.org/html/2412.05283v3#S5) | New restricted question; priority and ideal computations need another pass |
| Two-leak cycle global ambiguity / positive fibers | [Identifiability review](https://arxiv.org/html/2507.04496v1) | New restricted question; degree and positivity need separate certification |
| Minimal extra outputs for two-leak trees | [Identifiability review](https://arxiv.org/html/2507.04496v1) | Finite atlas feasible; new general graph rule not yet established |
| Support-aware loopless-FBA equivalence | [SEA formulation paper](https://arxiv.org/html/2502.00807v2) | Input-domain audit first; no claim of a new duality theorem or implementation defect |
| General dual-futile-cycle Conjecture 3.1 | [Vassena, §3](https://arxiv.org/html/2608.27081v2) | Explicitly open in checked source; no full resolution produced |
| Equivariant eight-parameter dual-cycle slice | [Vassena, §2–3](https://arxiv.org/html/2608.27081v2) | New subsidiary question; exact classification advanced |
| Content-dependent stochastic compartments, Open Problem 3.8 | [Primary article](https://link.springer.com/article/10.1007/s11538-026-01659-7) | Explicit open statement; proof cost exceeds the present short pass |
| Critical fragmentation boundary, Conjecture 3.14 | [Primary article](https://link.springer.com/article/10.1007/s11538-026-01659-7) | Recent claimed resolution requires checking; not presented as safely open |
| ACR / zero-divisor converse | [Original preprint](https://arxiv.org/abs/2401.00078) | Recent directly overlapping claim; priority uncertain, parked |
| COBRApy infeasible-solution values #1498 | [Issue](https://github.com/opencobra/cobrapy/issues/1498), [historical PR #1121](https://github.com/opencobra/cobrapy/pull/1121) | Reproduced with two exact LP backends, then parked: exposed infeasible iterates are permitted by existing semantics; changing to NaN needs an API decision |
| RoadRunner delayed initial events #857 | [Issue](https://github.com/sys-bio/roadrunner/issues/857) | Nonzero-start semantics must be fixed before interpreting a trace; not admitted as a new defect |
| RoadRunner conserved mass #1146 | [Issue](https://github.com/sys-bio/roadrunner/issues/1146), [PR #1354](https://github.com/sys-bio/roadrunner/pull/1354) | Closely overlapping merged repair; residual issue unestablished |
| COBRApy multivalued SBO #1454 | [Issue](https://github.com/opencobra/cobrapy/issues/1454) | Scalar SBML term and RDF list semantics differ; no speculative API repair |
| MEMOTE pseudo-reaction consistency #750 | [Issue](https://github.com/opencobra/memote/issues/750) | Existing methods and partition assumptions need comparison; no new algorithm claimed |

Additional duplicate checks excluded MEMOTE #777 (existing PR #780), COBRApy #1463 (PR #1481), COBRApy #1456 / #673 (PR #1264), and RoadRunner #958 / #1193 (merged repairs). We did not create competing patches for those reports.

## Separate earlier engineering result

The [MEMOTE direction audit](../projects/memote-direction-auditor/README.md) predates this screening. Its existing two-line negative-only direction fix is now prepared for upstream review with helper controls and complete ATP-detector regression tests. It is a software correction, not a newly resolved mathematical conjecture.

## What the evidence supports

The endpoint correction is a concise mathematical correction candidate. The symmetry classification is a partial theorem whose priority needs the author's assessment. Upstream acceptance, external peer review, and biological experiments are not supplied by passing exact scripts or continuous integration. Neither result justifies a promise of a major journal paper or a news story.
