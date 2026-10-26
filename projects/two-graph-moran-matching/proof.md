# A uniform fixation bound after deleting two disjoint mutant edges

This is a proposed subsidiary result, awaiting a fresh mathematical referee. The full complete-resident conjecture and the matching-deletion question for three or more deleted edges remain unresolved here.

## Exact statement

Let N>=4. Let H_N be the simple undirected graph obtained from K_N by deleting two vertex-disjoint edges. Residents disperse offspring on K_N; mutants disperse offspring on H_N. At each discrete step choose one of the N individuals uniformly, then replace a uniformly chosen neighbor in that individual's own type-specific graph. Reproduction rates are equal. Write h_v for the probability of eventual mutant fixation when initially only vertex v is mutant. Set C_N=N(N-2)+2. Then, for every vertex v,

    h_v <= 1/N - 4/[N^2 (N-1) C_N] < 1/N.

Consequently the uniformly initialized fixation probability rho(H_N,K_N) satisfies the same bound.

The proposed advance is the infinite N family with two disjoint deleted edges. The one-deleted-edge result rho=(N-2)/(N-1)^2 is Theorem A.1 of Tkadlec et al. (2023), not a new result here. Their published N=6 numerical enumeration includes this graph at that one population size.

## The reduced generator

The four endpoints of the deleted edges are paired. Define a as the number of pairs containing two mutants, b as the number containing exactly one mutant, and c as the number of mutants among the N-4 other vertices. The six possible (a,b) values are (0,0),(0,1),(0,2),(1,0),(1,1),(2,0). Put

    d=2-a-b,  k=2a+b+c,
    T=(2a+b)/(N-2)+c/(N-1),  L=(N-k)/(N-1).

Here T is the total mutant-parent weight seen by a target that is not missing a mutant neighbor. Let P be the one-step transition matrix and Q=N(P-I). The nonzero off-diagonal rates of Q are

| New state | Rate |
|---|---|
| (a,b+1,c) | 2d T |
| (a+1,b-1,c) | b[T-1/(N-2)] |
| (a-1,b+1,c) | 2a L |
| (a,b-1,c) | b L |
| (a,b,c+1) | (N-4-c) T |
| (a,b,c-1) | c L |

Each resident target in an empty pair sees all mutant parents; a resident target in a mixed pair cannot be replaced by its mutant partner. Every mutant has the same resident-parent replacement weight L because the resident graph is complete. The unmatched target rules follow in the same way. Terms with zero rates are omitted at the domain boundary. These observations also establish that the transition rates depend only on (a,b,c), so the reduction is exact.

## A superharmonic potential

Define

    V(a,b,c)=(N-2)k+a+2k(N-k)/[N(N-1)].

At extinction V=0, and at fixation V=C_N. At every singleton state V=C_N/N.

To verify QV<=0, put r=N-4-c. Thus c and r are nonnegative integers and N=c+r+4. With D=N(N-2)(N-1)^2>0, the six identities below give -D QV:

| (a,b) | -D QV after N=c+r+4 |
|---|---|
| (0,0) | 4c(r+4)(c+r+2) |
| (0,1) | 2(c+r+3)(2cr+5c+2r+4) |
| (0,2) | 2(2c^2 r+c^2+2cr^2+12cr+5c+3r^2+11r+2) |
| (1,0) | 4(r+2)(c^2+cr+5c+r+5) |
| (1,1) | 2r(2c^2+2cr+13c+3r+17) |
| (2,0) | 4r(c^2+cr+8c+2r+18) |

Every displayed polynomial has nonnegative coefficients in c,r. This proves QV<=0 throughout the entire state space for every N>=4, without restricting N to the computational baseline range. The identities follow directly by inserting the six transition rates into QV=sum q(x,y)[V(y)-V(x)]. They are checked symbolically in verify_certificate.py and recomputed independently in proof-attack/independent_verify.py.

## Absorption and comparison

Every nonabsorbing population has at least one resident. Since the resident graph is complete, there is a positive-probability path that repeatedly selects a resident parent and replaces a mutant until extinction. It has at most N-1 steps. For fixed N the state space is finite; the minimum over nonabsorbing states of the probability of absorption within the next N-1 steps is positive. Iteration proves that the absorption time tau is almost surely finite. In particular there are no nonabsorbing closed classes.

The process V(X_{t wedge tau}) is a bounded supermartingale. Dominated convergence gives

    C_N Pr_x(mutant fixation) = E_x V(X_tau) <= V(x).

Thus h(x)<=V(x)/C_N at every state x. For a singleton x, apply the Markov property once before this comparison:

    h(x) = E_x h(X_1)
         <= E_x V(X_1)/C_N
         = [V(x)+(P-I)V(x)]/C_N.

For a paired singleton x=(0,1,0), r=N-4, the table gives QV(x)=-4/[N(N-1)]. For an unmatched singleton x=(0,0,1), which exists only when N>=5, r=N-5, the same value follows from the first row. Since (P-I)V=QV/N,

    (P-I)V(x)=-4/[N^2(N-1)]

at every singleton. Substitution of V(x)=C_N/N proves the claimed quantitative strict bound. The N=4 case has only paired singletons and is included.

## Computational role and limits

The proof is the generator identity plus the bounded-supermartingale argument above. Finite computations are supplementary checks: exact rational absorption values were obtained for every N=3..8 and every matching size m<=floor(N/2), with independently constructed full-graph generator residuals checked at every transient labelled configuration. They reproduce the source m=1 formula and give rho(H_4,K_4)=15/68. The general-m baseline is not a uniform theorem.

No all-m theorem is claimed. The naive generalization replacing 2 by m in V fails to be superharmonic: for N=6,m=3,(a,b,c)=(1,2,0), its scaled drift QV is +1/100. This is a counterexample to that potential certificate, not to the matching conjecture or the original conjecture.

No known mathematical gap remains in the stated m=2 argument; correctness and substantive novelty still require the fresh referee. Current bounded primary-source searches found no directly equivalent all-N two-disjoint-edge theorem; they do not establish universal priority.
