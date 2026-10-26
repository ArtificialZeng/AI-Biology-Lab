# Frozen problem and exact contribution scope

Original source question (Tkadlec et al.,2023,Discussion): for every connected simple undirected unweighted graph G on N vertices, G!=K_N, is rho(G,K_N)<1/N for the neutral type-specific birth–death process with a uniform singleton mutant? This original question remains UNRESOLVED by this project.

Admitted natural restricted question: let G(N,m)=K_N minus a matching of size m, N>=3,1<=m<=floor(N/2). Is rho(G(N,m),K_N)<1/N? The matching family with arbitrary m remains UNRESOLVED here. The m=1 formula is already known (source TheoremA.1). The research pass found a sufficient certificate for m=2 only and stopped further scope expansion.

Proposed result-note: for N>=4 and exactly two vertex-disjoint deleted mutant edges, every singleton starting vertex has fixation probability at most

    1/N - 4/[N^2(N-1)(N(N-2)+2)] <1/N.

Uniform singleton averaging inherits this bound. All rates are neutral; residents use K_N, mutants use the graph with the two deleted edges; parents are chosen uniformly and their neighbors uniformly in their own type-specific graph. Both graphs are connected including N=4. No claim concerns death–birth updating, changing fitness, weighted/directed graphs, arbitrary missing-edge sets, or an exact general formula for fixation probability.

Nearest prior X: source TheoremA.1 covers exactly one missing edge with exact hitting function; source Fig4 numerically includes all112 connected graphs at N=6. The predecessor Melissourgos et al.,JCSS2022/1710.07365v3, studies type-specific graphs and bounds for stars/circulants, but the inspected Theorem16 and Corollary18 do not provide the neutral two-edge strict bound. Proposed delta Y is a uniform all-N two-disjoint-edge subclass theorem and a pointwise quantitative bound. Verifier Z is a complete six-case nonnegative-polynomial generator certificate plus a bounded-supermartingale comparison; finite rational state solves are supplemental.

Triage source.md has been copied verbatim and remains immutable. The citation/novelty audit is bounded; absence of a matching result in inspected sources is not a universal priority claim.
