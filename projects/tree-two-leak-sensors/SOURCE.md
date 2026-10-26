# Immutable source record — 2026-10-26

This is a newly formulated restricted question under a source-authored broad
open research direction. It is not a source-authored conjecture numbered in a
paper.

Primary direction: Meshkat and Shiu, *Structural identifiability of compartmental
models: Recent progress and future directions*, arXiv:2507.04496v1, §5,
Questions 5.2 and 5.4. Q5.4 asks:

> Given a graph structure, what are the minimal configurations of inputs and outputs that make the model generically locally identifiable

https://arxiv.org/html/2507.04496v1

Q5.2 specifically asks for extensions of the single-input, single-output
bidirected-tree classification to multiple inputs or outputs. The source does
not assert that every restricted placement subproblem is open.

Nearest dependency: Bortner, Gross, Meshkat, Shiu and Sullivant,
*Identifiability of linear compartmental tree models and a general formula for
input-output equations*, Advances in Applied Mathematics 146 (2023), 102490,
DOI 10.1016/j.aam.2023.102490; arXiv:2106.08487v2, Theorem 5.2, Lemma 5.5,
Proposition 5.7 and Proposition 2.15.

https://arxiv.org/html/2106.08487v2

These results give generic local identifiability for a bidirected tree with
collocated input and output and at most one leak at any vertex. Proposition
5.7 additionally gives expected coefficient-map dimension for arbitrary leak
sets when input and output are equal or adjacent. No restriction to paths,
stars, equal edge rates or root leaks is present.

Frozen subsidiary question:

For an arbitrary finite labelled undirected tree T, with strictly positive,
independent rates in both directions on every edge, put the sole input and
one output at r, add a second output at s != r, and place exactly two strictly
positive leaks at distinct vertices a,b. Is generic local structural
identifiability equivalent to P(r,s) and P(a,b) sharing an undirected edge?

Scope exclusions: global uniqueness, singular-locus equations, unobserved-root
models, unknown sensor gains, more than two leaks, additional inputs, and
statistical estimation with noisy observations are not part of this question.

Retrieved and inspected 2026-10-26. Current primary-source comparison and
limitations are recorded separately in triage.md. This file is frozen before
the proof handoff.
