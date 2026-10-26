# Source statement and contract

Daniele Cappelletti, Aidan S. Howells and Chuang Xu, *Stability of randomly
switching stochastic reaction networks with asymptotically linear transition
rates*, [arXiv:2507.09421v2](https://arxiv.org/abs/2507.09421v2), revised 12 August
2025. The latest version was refreshed on 26 October 2026; its PDF bytes matched
the earlier primary copy. Conjecture 43 appears in Section 6, printed page 29.

In the source's SAL setting, each environment has a mass-action reaction network
with at-most-monomolecular sources and an asymptotically linear drift matrix Mi.
The finite, independent environment chain Q is irreducible. Conjecture 43 assumes
nonnegative row vectors vi such that (vi Mi)m>0 for every coordinate in supp(vi),
with a nonempty intersection of all these supports. It predicts evanescence for
sufficiently small positive switching rates from initial populations whose
support meets that intersection.

Section 1.1, printed pages 3–4, defines evanescence as positive probability that
the process norm tends to infinity. Section 1.2, pages 4–5, imposes irreducibility
on Q, not on the full molecular chain. Neither the literal conjecture nor SAL
requires inflows, regeneration of every species, product occurrence of every
species, or a common strongly connected component. Product multiplicity two is
allowed because only the source-complex size is restricted.

For the exact model in model.json the vectors are (2,1,0) and (0,1,2), with shared
support B. The initial state is (0,1,0), in either environment. The proof verifies
all hypotheses and establishes extinction for every kappa>0. Thus it contradicts
the source's conclusion for arbitrarily small positive switching rates.

Example 32 and Proposition 33, printed pages 20–21, are the source's
already-known disjoint-support motivation. Their displayed mean matrices are
related to the fixed-type birth/death part of this witness. The present proof
is independent of the source's Lyapunov calculations. Its additional shared
feeder changes the hypothesis that Conjecture 43 proposes to add.

A separate source-contract audit visually checked pages 3,4,5 and29 and found no
extra population-irreducibility premise. This is a literal-statement result;
it does not supply a new general repaired criterion. See literature.md for the
bounded current-source comparison and its limitations.
