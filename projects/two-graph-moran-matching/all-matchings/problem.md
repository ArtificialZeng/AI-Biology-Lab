# Inherited problem and exact research scope

The immutable source.md records the original 2023 model and arbitrary-G question. Source original: is rho(G,K_N)<1/N for every connected undirected unweighted noncomplete mutant graph G under the neutral type-specific birth–death process? That broader source problem remains unresolved here.

Frozen restricted target carried from R14: for every N>=3 and 1<=m<=floor(N/2), with mutants on K_N minus any matching of m edges and residents on K_N, prove pointwise strict singleton fixation h_v<1/N. This is the same `matching-complement-neutral-two-graph-r14` registry record and `two-graph-neutral-moran` research family, PASS 2 of at most 3, following completed R14 PASS 1. The question is not reset at R15.

Nearest X: source Theorem A.1 gives the m=1 exact formula (N-2)/(N-1)^2; its published N=6 numerics cover matching complements at that size. R14/PR18 gives the all-N m=2 pointwise gap 4/[N²(N-1)(N(N-2)+2)].

New Y proposed for review: for all the restricted target's parameters and every singleton vertex, h_v<=1/N-2m/[N²(N-1)(N(N-2)+2m-2)]. Uniform averaging inherits this inequality. The new mathematical range is unbounded m>=3; m=1 remains credited and m=2 is preserved as the identical specialization of this bound. A related accepted theorem must extend PR18, https://github.com/ArtificialZeng/AI-Biology-Lab/pull/18, rather than create a new small PR for the same family.

Verifier Z: proof.md derives the exact six-rate quotient, supplies a nonnegative falling-factorial certificate for every integer configuration and parameter, and gives the bounded-supermartingale and strict-singleton comparison. certificate.py recomputes the exact algebra, coefficient identities and 12 exhaustive shifted certificates. labelled_check.py separately constructs the full labelled parent-target kernel at varied finite parameters; it imports neither the quotient rates nor the certificate generator. Fresh mathematical review remains required.

The source's connected-graph model restricts this claim to N>=3; removing an edge from K_2 is outside it. Every allowed matching complement H is connected. No all-arbitrary-G result follows merely from this matching theorem, and no edge-monotonicity principle is used.
