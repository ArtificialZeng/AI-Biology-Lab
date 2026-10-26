# Theorem and proof: retained-root sensors on a tree

## Statement

Under the hypotheses in SCOPE.md, put m=2n-2+L and let q be the maximum number of pairwise vertex-disjoint paths from distinct observed vertices to distinct designated leak vertices. Endpoints count as vertices; paths of length zero are allowed. Then:

1. The harmonic evaluation matrix H[O,B] has rank q for **every** positive edge-rate assignment.
2. The transfer coefficient map has generic differential rank 2n-2+q and generic fiber dimension L-q. Consequently the model is generically locally structurally identifiable if and only if q=L.
3. If q<L, through every positive parameter point there is an injectively parameterized local (L-q)-dimensional positive family with exactly the same transfer functions at O.
4. Keeping the root sensor, exactly L-1 additional sensors are necessary and sufficient when their positions may be chosen freely.

Only item 1 and the ambiguity family are asserted at every positive point; full differential rank is generic. The theorem includes n=1, a leaking root, all vertices leaking, empty interior, sensors at leaks, and superfluous sensors. It excludes unknown output gains and an unobserved root. It does not assert global uniqueness.

## 1. Harmonic extension

Write the outgoing row Laplacian as

    K_jj = sum_(i~j) k_ij,  K_ji = -k_ij (i~j).

Thus K 1=0, and h is harmonic off B if (Kh)_j=0 for j outside B. Let N=V\B. The Dirichlet matrix K[N,N] is nonsingular: if a harmonic function vanishes on B, a positive maximum at an interior vertex propagates to every neighbor by positivity and eventually to a boundary vertex, a contradiction. Apply this also to the negative minimum. Consequently boundary values determine a unique harmonic h; H is its extension matrix and H[B,B]=I. The maximum principle gives H_jb>=0 and sum_b H_jb=1. In particular the space of harmonic functions has dimension L.

Let M be the n by n matrix whose row j is the j-th row of K for j in N and the j-th identity row for j in B. Then M is invertible and

    H = M^(-1)[:,B].

This also treats N empty; the determinant of an empty matrix is 1.

## 2. Every full harmonic minor is a same-sign forest sum

For any vertex subset W of size L, Jacobi's complementary-minor identity gives

    det H[W,B] = (+/-) det K[V\B,V\W] / det M.       (1)

The sign depends only on the ordered sets B,W. Indeed, the complementary rows of M in this identity lie outside B and therefore equal those of K.

Apply the classical all-minors matrix-tree theorem to the outgoing row Laplacian K (equivalently transpose the column-Laplacian statement). The numerator of (1) is, up to a fixed set-dependent sign, the signed sum over spanning forests with exactly L components, each containing exactly one vertex of B and exactly one vertex of W. Orient each component toward its unique B vertex. Each term has strictly positive weight, the product of its directed edge rates; the only variable sign is the parity of the bijection W→B induced by components. See Zernik, arXiv:1308.2160v1, Definition 1, Theorem 2, Definition 4; this is the classical Chaiken theorem, not a new forest identity.

Every such forest contains L pairwise vertex-disjoint simple W-to-B paths. Conversely, any such linkage extends to a spanning forest of this kind: begin with its disjoint path components and repeatedly attach an unassigned neighboring vertex using one edge to an existing component. Connectedness ensures all vertices can be attached; this procedure never joins two rooted components. Orient toward B.

## 3. A tree has at most one disjoint full-linkage pairing

A vertex in W∩B must match itself in a vertex-disjoint linkage: otherwise its role as a start and an end would place it on two paths. Delete these singleton paths for the following argument.

For any edge e, let U be either component of T-e. Vertex-disjointness allows at most one linkage path to cross e. The number crossing from a start in U toward an end outside U, minus the reverse number, is

    |W∩U| - |B∩U|.

Thus this difference determines whether e is used and, if used, its orientation. (A realizable difference must lie in {-1,0,1}.) Every used edge and its orientation are therefore fixed by W and B. Their directed path components determine the matching uniquely.

It follows that all forests contributing to (1) induce the same pairing and have the same sign. There is no cancellation for positive rates. Combining this with forest existence proves

    det H[W,B] != 0
      iff W can be linked vertex-disjointly to all of B.      (2)

This equivalence holds for every positive rate assignment.

## 4. Arbitrary harmonic rank equals maximum linkage

Consider an s by s minor H[S,J], S⊆O, J⊆B. If it is nonzero then S∩(B\J) is empty, since a row indexed by b in B\J would be identically zero. Put

    W = S ∪ (B\J).

After ordering the rows and columns into blocks, H[W,B] has form

    [ H[S,J]   H[S,B\J] ]
    [   0            I    ].

Hence its determinant equals det H[S,J] up to sign. By (2), W has a full disjoint linkage to B. Every vertex of B\J must match itself, so deleting those singleton paths leaves s disjoint paths S→J. Therefore every nonzero minor has size at most q, proving rank H[O,B]≤q.

For the reverse inequality, choose q disjoint paths from starts S⊆O to distinct leaks. Truncate each path at the first leak it reaches. Their new endpoint set J still has size q, the paths remain disjoint, and they avoid every boundary vertex outside J. Add singleton paths at all b in B\J. The resulting paths link W=S∪(B\J) to all B. By (2), det H[W,B] is nonzero; the block identity then makes det H[S,J] nonzero. Thus rank H[O,B]≥q. This proves item 1.

Notice that neither a separator theorem nor a generic noncancellation assumption is required. The first-leak truncation is essential when considering partial boundary minors.

## 5. Root transfer map and its generic kernel

Let

    D(z)=det(zI-A),  N_r(z)=det((zI-A)[V\{r},V\{r}]).

D is monic degree n and N_r monic degree n-1. Their nonleading coefficient map C_r has 2n-1 coordinates. By Bortner et al., Proposition 5.7, its generic rank is 2n-1 when L>=1. One can also apply their single-leak theorem on the coordinate subspace where all but one prescribed leak equal zero: a nonzero Jacobian minor there is a nonzero polynomial for the full family, so it stays nonzero generically in the open positive orthant. Thus

    dim ker dC_r = m-(2n-1)=L-1                 (generically).  (3)

C_r is dominant onto its (2n-1)-dimensional coefficient space. The universal resultant of the monic pair D,N_r is not the zero polynomial, so its pullback is nonzero. Generically D and N_r are coprime. On this set the root transfer F_r=N_r/D determines its displayed coefficients locally, and the transfer-data differential has the same kernel as dC_r. This prevents an unreduced input-output equation from being mistaken for minimal transfer data.

For n=1 the same conclusions follow directly from F_r=1/(z+ell_r), rank 1, and L=q=1.

## 6. Harmonic diagonal similarities exhaust that kernel

For a harmonic function h with h_r=0, define, for sufficiently small real t,

    d_j(t)=1+t h_j,  Q(t)=diag(d_j(t)),
    A(t)=Q(t) A Q(t)^(-1).

The new directed edge rate is k_ij(t)=k_ij d_i(t)/d_j(t)>0. The diagonal entries are unchanged. At each nonleak j,

    sum_(i~j) k_ij d_i(t)=d_j(t) sum_(i~j) k_ij,

so the new column sum remains zero. At b in B define the new leak as minus its new column sum. It varies continuously from ell_b>0, hence remains positive locally. Thus this is a family in exactly the stated compartment model.

Since d_r=1, diagonal similarity gives

    F_o(z;t)=d_o(t) F_o(z;0),      o in V.             (4)

In particular F_r is fixed. The map from h to its parameter tangent is injective: if all edge derivatives k_ij(h_i-h_j) vanished, connectedness would make h constant, and h_r=0 would make h=0. The harmonic space with h_r=0 has dimension L-1 because H_r,* 1=1. Thus these L-1 tangent directions lie in ker dC_r and, by (3), exhaust it generically.

## 7. Additional sensors give exactly q-1 new directions

Let E be harmonic evaluation h↦h_O. Its rank is q by item 1. Since r∈O and evaluation at r is a nonzero linear functional, restricting E to h_r=0 has rank q-1 and kernel dimension L-q. To see this directly, the image of that restriction is im(E)∩{v_r=0}, a codimension-one subspace of im(E).

For each o, F_o is not the zero rational function. The first nonzero high-frequency Markov term is the positive product of rates along the unique directed r-to-o path. Therefore (4) shows that, on the root-data kernel, the additional transfer differential has exactly the rank of h↦h_O, restricted to h_r=0, namely q-1. Consequently

    rank dC_O = (2n-1)+(q-1)=2n-2+q                (generically).

Here C_O may be taken as the common D and the numerator coefficients of all F_o. Generic coprimality of D,N_r already ensures that these coordinates and rational transfer data have the same local differential rank. Generic local structural identifiability is therefore equivalent to 2n-2+q=m, i.e. q=L.

## 8. Exact ambiguity dimension at every positive point

The space of harmonic h vanishing on O has dimension L-q by item 1. For small b in this space, use Q_b=diag(1+h), A_b=Q_b A Q_b^(-1). The same construction preserves positivity, the prescribed leak support, and every transfer F_o exactly. It is locally injectively parameterized: equality of all directed edge rates for two positive scaling vectors forces their coordinatewise ratio to be constant along the connected tree; both scalings equal 1 at r, so that constant is 1. This supplies the claimed local (L-q)-dimensional ambiguity family through every positive point. Generic regular fibers have exactly that dimension by item 2; no stronger statement about nongeneric fiber dimensions is made.

## 9. Minimal sensor count and the earlier two-leak result

Since q≤|O|, full identifiability needs |O|≥L, hence L-1 additional sensors. Choose a leak b with H_rb>0; one exists because sum_b H_rb=1. Set O={r}∪(B\{b}). If r lies in B, positivity forces b=r and O=B; otherwise these are L distinct vertices. Expanding the resulting harmonic evaluation determinant along its identity boundary rows gives ±H_rb!=0. Hence q=L and L-1 additional sensors suffice. A purely graphical choice is the first leak along a shortest r-to-B path.

For L=2 and O={r,s}, a disjoint linkage between the two endpoint sets exists if and only if P(r,s) and P(a,b) share an edge. Thus the previously audited two-leak criterion is a corollary. This project extends the same structural result rather than claiming an unrelated discovery.

## Dependencies and gap list

Dependencies: Bortner et al., arXiv:2106.08487v2, Proposition 5.7 and the standard Jacobian criterion; Jacobi's complementary-minor identity; the classical all-minors matrix-tree theorem, with its precise sign/pairing statement directly inspected in Zernik arXiv:1308.2160v1, Theorem 2 and Definition 4; finite Dirichlet maximum principle; elementary diagonal similarity.

No unproved mathematical gap is currently known to the proposer. A separate AI mathematical audit is recorded with the verification evidence; external human review is pending. Finite checks supplement it and do not prove the infinite theorem. Historical priority is unestablished beyond the bounded comparison in literature.md.
