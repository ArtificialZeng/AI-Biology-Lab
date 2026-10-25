# Research roadmap

Updated: 26 October 2026.

AI Biology Lab develops small, inspectable research results in mathematical biology. Its first project studies a published phosphorylation model using exact algebra, numerical computation and explicit verification. A result about this model is not an experimental discovery about a living cell.

## Current project and evidence

The starting point is Conradi, Mincheva and Shiu, *Emergence of Oscillations in a Mixed-Mechanism Phosphorylation System* (2019), [arXiv:1809.02886v2](https://arxiv.org/html/1809.02886v2), especially Table 1, equations (2)–(3), Figure 5(c) and Conjecture 6.2. The published model and its previously observed oscillations are prior work.

With the published Table 1 rates and fixed totals

\[
P_{\mathrm{tot}}=8,\qquad S_{\mathrm{tot}}=40,
\]

the present result gives the sufficient bound

\[
K_{\mathrm{tot}}\geq \frac{4660166}{209569}<23.
\]

For every kinase total in this range, the unique positive equilibrium has exactly two eigenvalues in the open right half-plane and four in the open left half-plane on its six-dimensional compatibility class. Thus every \(K_{\mathrm{tot}}\geq23\) is covered. The bound is not asserted to be sharp. Equilibrium uniqueness and the general Hurwitz approach come from prior work; this project's contribution is the explicit uniform bound and its exact certificate.

This supplies a particular witness for Conjecture 6.2(1). It does not establish periodic-orbit existence, stability or invariance in part (2). Equilibrium instability alone does not imply a stable oscillator.

| Evidence | Current status | Inspect |
|---|---|---|
| Finite-kinase equilibrium instability | Exact algebraic argument with a reproducible rational certificate | [Proof](../projects/mixed-phosphorylation/evidence/finite_kinase_proof.md), [certificate checker](../projects/mixed-phosphorylation/manuscript/verify_certificate.py), [short paper](../projects/mixed-phosphorylation/manuscript/main.pdf) |
| Periodic trajectories, periods, extrema and Floquet multipliers | Floating-point evidence; not a validated periodic-orbit certificate | [Independent implementation report](../projects/mixed-phosphorylation/verification.md), [response data](../projects/mixed-phosphorylation/evidence/period_response.json) |
| High-kinase reduction and finite-concentration response | A derivation and numerical research lead, with scope recorded separately from the proved theorem | [Theory notes](../projects/mixed-phosphorylation/theory.md) |
| Novelty and priority | Not confirmed; directly related 2023/2025 conference work has been identified | [Literature audit](../projects/mixed-phosphorylation/literature/novelty.md), [focused theorem comparison](../projects/mixed-phosphorylation/literature/finite_threshold_check.md) |
| External assessment | No external peer review is claimed | The checks in this repository are reproducibility and AI-assisted review evidence, not journal acceptance or expert endorsement |

## Milestone 1: certify a stable periodic orbit

**Question.** Does a computed trajectory correspond to an actual attracting periodic orbit of the published ODE, rather than a long transient or a numerical artifact?

Start with the existing finite-kinase examples at \(K=100,1000,10000\), and the separately derived limiting system. Select the example with the least expensive decisive validation; do not require all four cases to succeed before assessing a useful result.

**Work.** Construct a transverse section or an equivalent phase-fixed boundary-value problem. Use validated integration or another rigorous enclosure method to prove a return-map fixed point and to bound the transverse linearized dynamics. Preserve the conservation laws and positivity conditions of the published model. A small shooting residual or a numerical multiplier below one is a starting guess, not the certificate.

**Acceptance conditions.**

- An explicit enclosed periodic orbit, period interval and positive compatibility class.
- A transverse stability certificate. The neutral time-shift direction must be handled explicitly; roundoff-scale multipliers must not be treated as exact zeros.
- A verifier and a second implementation of the model equations that reproduce the decisive checks.
- A statement limited to the certified orbit and parameters. Uniqueness and global attraction require separate proofs.

If validation fails, record whether the obstacle is an implementation limitation, insufficient precision or a mathematical obstruction. A timeout does not establish nonexistence.

## Milestone 2: persistence for every sufficiently large kinase total

**Question.** Can stable oscillations be proved to persist on an unbounded concentration range?

Use \(\varepsilon=1/K_{\mathrm{tot}}\) to examine the high-kinase reduction. This is not the usual small-enzyme Michaelis–Menten limit. Check the actual singular scaling, the relevant compact positive region and the fast-direction stability before applying a persistence theorem.

**Acceptance conditions.**

- A rigorous limiting system and a justified connection to the full mass-action equations.
- A certified hyperbolic periodic orbit of the limiting system, or an alternative complete route to existence.
- A proof with explicit quantifiers: for the stated fixed rates and totals, a periodic orbit exists for every \(K>K_0\). Give a usable \(K_0\) if the estimates permit it; otherwise disclose that the threshold is nonconstructive.
- A stability conclusion only if its hypotheses and transverse estimates are verified.

The existing equilibrium threshold \(4660166/209569\) is not a periodic-orbit threshold. It must not be reused as one without proof. A result conditional on an unverified limiting cycle must remain labelled conditional.

## Milestone 3: finite-concentration corrections and error bounds

**Question.** How close to concentration-independent is the model's oscillation, and how much kinase is enough for a specified tolerance?

First distinguish the quantities being compared:

- The scalar set projection onto \(S_2\), which is its interval of attained concentrations.
- The period and the minimum/maximum concentrations.
- The output waveform in the model's physical time, compared after phase alignment.

These are not equivalent. The wording of the original conjecture leaves a projection-versus-waveform ambiguity. A changing period alone does not disprove equality of scalar set projections, and two unequal examples do not refute an existential assertion over other totals.

**Target.** Where justified, derive expansions such as

\[
T(K)=T_\infty+\frac{a_T}{K}+R_T(K),
\qquad |R_T(K)|\leq\frac{C_T}{K^2}
\quad(K\geq K_0),
\]

and corresponding formulas for extrema. This displayed form is a research target, not a current theorem. A fitted inverse-power curve does not prove the expansion or its remainder.

**Acceptance conditions.** Certify the coefficients or bound them sufficiently to decide whether they vanish; justify the phase condition and the extrema used; provide explicit remainder constants on a stated range. Translate the bounds into a model-level design question, such as the kinase range that guarantees a chosen period tolerance. Report time in the published model's units, without inventing a biological clock calibration.

Generic sensitivity formulas are established methods. The prospective contribution is a model-specific certified law, useful bounds or a precisely scoped obstruction.

## Milestone 4: interpretable parameter regions

**Question.** Which concentrations or rate combinations preserve oscillation, and which changes destroy it?

Expand beyond \((P,S)=(8,40)\) only after the single-parameter case is understood. Choose a small set of biologically interpretable parameters, explain the choice and freeze the domain before drawing conclusions. Compare the resulting boundary with the published stability surface and known Hopf results.

**Acceptance conditions.**

- An explicitly specified region with certified conclusions, including unresolved boundary portions.
- A distinction between equilibrium instability, periodic-orbit existence and orbital stability.
- A quantitative tolerance or failure condition that can be explained using enzyme binding, sequestration or catalytic steps, with the proposed explanation checked against the equations.
- Verification that any new boundary, sufficient condition or classification is not merely a replot of prior work.

A heat map is exploratory evidence. A complete classification requires coverage and boundary arguments, not a dense grid alone. A sharp theorem on a useful restricted domain may be a better result than an unsupported claim about all parameters.

## Reserve candidates

The following are leads, not accepted projects or confirmed open problems:

- The structured Jacobian stability question in [Vassena, arXiv:2608.27081v2, Conjecture 3.1](https://arxiv.org/html/2608.27081v2). Its exact scope, subsequent literature and current status must be checked before admission.
- Concentration robustness in the reverse mixed mechanism. Grasse and Conradi's [2025 conference abstract](https://biomath.math.bas.bg/biomath/index.php/bmcs/article/download/1787/1594/7828) is directly relevant; the unavailable detailed results must not be ignored.
- Reusable improvements to symbolic Jacobian, stability-certificate or SBML workflows, assessed against real existing tools and identical correctness requirements.

For the first two leads, current open status remains pending fuller documentation and, if later authorized, clarification from the relevant authors. No author contact has occurred. No project will be advertised as resolving an open problem solely because a search did not find a solution.

## Execution and release rules

Use the current computer for exact algebra and narrow numerical experiments before spending on additional hardware or cloud compute. Keep research passes bounded by a named claim and a decisive next test; record actual resource use instead of promising a completion date or cost.

Freeze each candidate statement and its evidence before a fresh mathematical review. Preserve separate labels for theorem, conditional theorem, numerical observation and conjecture. Update the literature comparison around the final statement. Release manuscripts only for the scope that survives review, and keep the full original conjecture unresolved while any required assertion remains unproved.

External specialist review would strengthen confidence and help clarify overlap. It is a future validation goal, not an endorsement already obtained. Media interest is a possible consequence of a result that withstands these checks, not an acceptance criterion for its mathematics.

## 中文简述

当前已得到固定速率、\(P=8,S=40\) 下的精确不稳定性界：\(K\geq4660166/209569\)，因此所有 \(K\geq23\) 均被覆盖。这只是原猜想6.2(1)的一个明确实例，不是稳定振荡或完整猜想的证明；该界也未宣称最优。

后续按四步推进：**认证实际稳定周期轨道 → 证明所有充分大K下周期轨道持续存在 → 给出有限K的周期/振幅校正及误差界 → 扩展到可解释的参数区域**。现有振荡曲线、周期和Floquet数值仍是浮点证据。

2023/2025已有直接相关会议研究记录，优先权尚未确认。项目目前没有外部同行评审或专家背书。其他模型与公开编号猜想仅作备选，最新开放状态仍需补充文献及必要的作者澄清；未进行任何作者联系。
