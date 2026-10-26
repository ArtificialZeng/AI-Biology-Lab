# Complete-graph extremality against every regular graph

Proposed exact result-note, awaiting fresh mathematical scrutiny. The source's arbitrary irregular-graph questions remain unresolved. This is the final research pass (3/3) of the existing two-graph neutral Moran family and is intended to augment existing PR18, preserving its all-matching theorem.

## Model and proposed theorem

Let G be a connected, simple, undirected, d-regular graph on N vertices, with G not complete. Consequently N>=4 and 2<=d<=N-2. In the neutral type-specific birth–death process, each discrete step chooses a parent uniformly among all N individuals; it chooses a neighbor uniformly in its type's graph and overwrites that neighbor. Let h_v(G,K_N) be mutant fixation probability from the singleton v when mutants use G and residents use K_N. Write q=N-1-d>0. Then for every vertex v,

    h_v(G,K_N) <= 1/N - q/(d N^3) < 1/N.                 (1)

Conversely let h_v(K_N,G) denote mutant fixation when the complete graph is used by mutants and G by residents. The same proof gives the direct corollary

    h_v(K_N,G) >= 1/N + 2q/[d N^4(N-1)] > 1/N.         (2)

Averaging over a uniformly chosen initial singleton preserves both bounds. These are quantitative bounds, not sharp formulas or an equality of fixation probabilities across regular graphs. No restriction to vertex-transitive graphs is imposed. Complete G gives the neutral equality, as a separate boundary control. Disconnected graphs, weighted or directed edges, unequal fitness, death–birth updating, and irregular G are not claimed here.

## Generator and local statistics

Fix the mutant set M. Put k=|M|, r=N-k, and let e be the number of G-edges with both endpoints in M. For every resident vertex j let g_j be the number of its G-neighbors in M. Set

    B=sum_{j notin M} g_j = d k-2e,
    T=sum_{j notin M} g_j^2,
    U=B/d,    L=k r/(N-1),    Q=N(P-I).

U and L are the total gain and loss rates for mutant number in Q. A resident j becomes mutant at rate g_j/d, because all mutant parents have degree d. A mutant becomes resident at rate r/(N-1), because resident parents use K_N. These observations imply the exact identities

    Qk=U-L,
    Q(k^2)=2k(U-L)+U+L,
    Qe=T/d-2e r/(N-1).                               (3)

For the last identity, gaining mutant j creates g_j internal edges, at rate g_j/d. Losing mutant i removes its internal mutant degree; summing those degrees over mutants gives 2e. These are labelled-state identities; no assumption of a lumped Markov chain on (k,e) is made. In particular T need not be determined by (k,e).

## Sum-of-squares potential

Define

    V(M)= d k(2N-1+k)/[2(N-1)]-e,
    C=d N^2/(N-1)>0,
    D(M)=sum_{j notin M} [g_j-d k/(N-1)]^2.

Substitution of (3), followed by B=dk-2e, gives the exact identity

    QV(M) = -D(M)/d <= 0.                            (4)

For direct verification, expanding the right side yields

    -T/d + 2kB/(N-1) - d k^2 r/(N-1)^2.

Expanding the left side using (3) gives precisely this expression. No estimate, asymptotic argument, fit, or finite enumeration is used in this sign proof. V(empty)=0 and V(all)=C because G has Nd/2 edges. At any singleton, e=0 and V=C/N. At a singleton v, g_j=1 at d neighbors and 0 at q other residents, so

    D({v})=d q/(N-1),    QV({v})=-q/(N-1).           (5)

The formulas are independent of v but not of the graph's degree.

## Absorption and the singleton deficit

Let tau be the first time all vertices share a type. From any transient state, keep one resident parent and choose it repeatedly to overwrite each remaining mutant. Since its graph is complete, each prescribed step has probability 1/[N(N-1)]. Extinction is achieved in at most N-1 steps; padding by resident-to-resident updates if needed gives a uniform positive probability of absorption in each block of N-1 steps. Therefore tau is almost surely finite with finite expectation by a geometric-tail bound. The state space and V are finite, so the stopped V process is a bounded supermartingale. Dominated convergence yields

    C h(M) <= V(M).

Apply this inequality once after a singleton update, using harmonicity of h and (5):

    h({v})=E h(M_1)
           <= E V(M_1)/C
           =[V({v})+QV({v})/N]/C
           =1/N-q/(d N^3).

This proves (1) for every v.

The telescoping version of (4) also gives the explanatory exact identity

    h(M)=V(M)/C - (N-1)/(d^2 N^3)
                     E_M sum_{t=0}^{tau-1} D(M_t).  (6)

Finite expected tau and bounded D justify expectation and summation. This relates the loss in fixation probability to accumulated squared discrepancies in resident neighborhood occupancy. Equation (6) is not needed for either bound.

## Reversing the two types

Color complementation maps a process with mutants on K_N and residents on G, initially a single mutant at v, onto the process with mutants on G and residents on K_N, initially M_0=V(G)\{v}. Thus, path by path,

    h_v(K_N,G)=1-h_G,K_N(V(G)\{v}).                  (7)

For M_0=V(G)\{v}, k=N-1, e=Nd/2-d, V(M_0)/C=(N-1)/N, and D(M_0)=0. All resident-to-mutant conversions lead to fixation, with D=0, and updates without color change also leave D=0. The only updates with nonzero possible D at time one are those where the sole resident v overwrites a mutant u!=v. Each such parent-target event has probability 1/[N(N-1)]. In the resulting state the two residents are v and u, each having mutant-neighbor count d-1_{uv in E(G)}, and k=N-2. Hence

    D(V(G)\{v,u})=2[d/(N-1)-1_{uv in E(G)}]^2.

There are d adjacent choices of u and q nonadjacent choices. Consequently

    E_{M_0} D(M_1)
      =2{d[1-d/(N-1)]^2+q[d/(N-1)]^2}/[N(N-1)]
      =2dq/[N(N-1)^2]>0.                            (8)

Apply (4) over two updates (absorbing states stay fixed), then the global bound h<=V/C:

    h_G,K_N(M_0) <= E V(M_2)/C
                  = (N-1)/N - E D(M_1)/(d N C)
                  = (N-1)/N -2q/[d N^4(N-1)].

Combining with (7) proves (2). This proof does not assert a general complement relation between two singleton fixation probabilities; the complement initial condition has N-1 mutants, as displayed explicitly.

## Scope, predecessors, and verification

The 2023 source proves the single-deleted-edge case and studies regular graph examples numerically and analytically for particular families. The 2022 predecessor's cited regular-graph bounds are nonneutral asymptotic bounds; they do not give the finite r=1 inequalities above. Current bounded primary-source searches located no equivalent all-regular sum-of-squares theorem, but do not establish universal novelty or historical priority. The original arbitrary irregular-G conjectures remain unresolved.

The lab's earlier all-matching-deletion theorem is retained as a distinct complementary subclass; only perfect matching complements are regular. The new theorem includes sparse cycles, arbitrary connected regular graphs, and dense regular graphs, not merely additional matching parameters. All related work should remain in the same existing PR18. `certificate.py` checks the symbolic identities exactly. `labelled_check.py` builds full parent-target kernels independently from graph adjacency on fixed finite fixtures, including non-vertex-transitive regular graphs, complete-graph controls, and both type orientations. These checks supplement the infinite proof and are not a finite-to-infinite extrapolation. Fresh model review is not human peer review, and no formal trusted-runner certificate or priority claim is made.
