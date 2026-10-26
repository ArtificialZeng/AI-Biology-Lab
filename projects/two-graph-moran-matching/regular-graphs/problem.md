# Precise target and admitted scope

Original source question: for every connected noncomplete simple undirected graph G on N vertices, does uniformly initialized neutral type-specific birth–death fixation satisfy rho(G,K_N)<1/N and rho(K_N,G)>1/N? Both broader arbitrary irregular-G questions remain unresolved here.

Proposed new restricted theorem: let G be any connected, simple, undirected, noncomplete d-regular graph on N vertices, and put q=N-1-d>0. For every vertexv, prove

    h_v(G,K_N)<=1/N-q/(dN^3),
    h_v(K_N,G)>=1/N+2q/[dN^4(N-1)].

Every allowed graph has N>=4,2<=d<=N-2; no vertex-transitivity, circulant structure, density limit, or asymptotic population-size limit is imposed. The first argument denotes the mutant graph in this dossier, following the 2023 source. Uniformly averaged bounds are immediate. The bounds are not claimed sharp.

Nearest X: 2023 source single-edge theorem and particular regular-graph families; 2022 predecessor nonneutral asymptotic bounds; labR15 all-matching theorem. DeltaY: an all-regular, all-finite-N, pointwise pair of strict inequalities, certified by a single labelled-state sum-of-squares potential. VerifierZ: derive the general generator identity exactly, establish absorption and stopping, then raw labelled parent-target checks in both color orientations. The explicit reverse inequality is a corollary of the same potential and color complementation, not another project.

Research family two-graph-neutral-moran already completed two passes before this admitted final pass3/3. Cheap triage cycle-complement calculations led to the potential; the full graph-local identity is not inferred by extrapolation. This proposed contribution must augment existing https://github.com/ArtificialZeng/AI-Biology-Lab/pull/18 while preserving its matching theorem. No fourth pass, separate micro-PR, or full arbitrary-G resolution is justified.
