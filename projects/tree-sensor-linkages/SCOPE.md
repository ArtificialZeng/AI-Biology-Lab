# Precise reading

T=(V,E) is a finite labelled undirected tree, n>=1. On each directed edge j->i the independent unknown rate is k_ij>0. B is a nonempty subset of V, L=|B|, with independent unknown leak ell_j>0 exactly for j in B. The matrix A has A_ij=k_ij on edges, zero elsewhere off-diagonal, and A_jj=-sum_(i~j) k_ij-1_(j in B) ell_j.

The system is x'=Ax+e_r u. Its data are the labelled rational functions F_o(z)=e_o^T(zI-A)^(-1)e_r, o in O, with known individual sensor gains and r in O. Initial conditions are not additional known constraints. The generic rank means differential rank of the transfer coefficient map on a nonempty Zariski open subset meeting the positive orthant. Global uniqueness, noisy estimation, absent root observation, shared rates, and unknown gains are outside scope.

Let H be the n by L harmonic-extension matrix: H[B,B]=I and, for j outside B, sum_(i~j) k_ij(H_i,*-H_j,*)=0. Thus any harmonic function is h=Hb for its boundary vector b. q is the maximum cardinality of a vertex-disjoint linkage from O to B, including endpoints and permitting a singleton path at an intersection.

The stated result: rank H[O,B]=q for every positive set of edge rates; generic transfer coefficient rank=2n-2+q; generic ambiguity dimension=L-q; full generic local structural identifiability iff q=L. At every positive point, L-q independent harmonic similarities generate an exact local positive ambiguity family. At least L-1 additional sensors are necessary, and L-1 suffice when their positions may be chosen.

Known baseline: root-only generic rank 2n-1 (Bortner et al., Prop. 5.7); diagonal scaling and all-minors matrix-tree arguments are classical. The proposed contribution is their complete placement classification for retained-root trees with arbitrary prescribed leaks, extending the earlier two-leak corollary. Novelty remains bounded-search, not historical-priority certification.
