# Expected dimension for non-elementary even-order 3-sunlets

## Exact statement and prior result

Let G be a finite abelian group of even order n>=6 and exponent greater than2. For the general group-based model (B={id}) on the source's 3-sunlet network, the complex affine cone has dimension

    dim V_N^G = 5n−4 = 5(|G|−1)+1.

The known odd-order theorem, the parameter-count upper bound, and the tropical dimension bound are credited to Cox–Gross–Martin2025 and their cited predecessors. The contribution here is a uniform lower-rank proof for this even-order class. The exponent-two class is excluded; its remaining dimension problem is not resolved. No new claim of model topology identifiability follows solely from this dimension statement.

Use the source parameterization, indexed by g,h in G and third label −g−h:

    phi_(g,h,−g−h) = w1_g w2_h w3_(−g−h)
                      (w4_(g+h) w5_g + w4_h w6_g).

There are6n independent parameter coordinates. Lemma2.1 of the source gives the upper bound5n−4. Its Section2.2 gives the lower bound dim V>=rank_R A_lambda, where each column of A_lambda is the exponent vector of the uniquely least-weight monomial in the corresponding coordinate. We therefore construct one strict weight vector and determine its rank over characteristic zero.

## One strict tropical chamber

If G is cyclic, fix a generator z and identify G with the integers modulo n using z. Put

    X={1,2,...,floor(n/2)},    Y={floor(n/2)+1,...,n−1}.

If G is not cyclic, let X contain every nonzero element of order2 and exactly one member of each remaining inverse pair {g,−g}. Put Y=G\(X union{0}). Thus, in both cases, every nonzero involution lies in X, and exactly one of g,−g lies in X when g is not an involution.

For compact comparison notation set nu0=4 and nu_g=0 for g!=0; set mu0=−1, mu_g=2 on X and mu_g=−2 on Y. Actual strictly positive weights can be chosen as

    lambda1_g=lambda2_g=lambda3_g=1,
    lambda4_g=nu_g+1,   lambda5_g=3,   lambda6_g=3+mu_g.

Their first-minus-second monomial weight is exactly

    nu_(g+h)−nu_h−mu_g.

For g=0 this equals1, so all g=0 columns choose the second monomial. If g lies in X, the first monomial is chosen for every h except h=−g. If g lies in Y, the first monomial is chosen only at h=0. All comparisons are strict: for g!=0 the nu difference is one of −4,0,4, whereas mu_g is2 or−2; for g=0 the margin is1. In particular, an involution g=−g in X has first choice at h=0, second choice at h=g, and first choice at every other h. No tie-breaking assignment is used.

Denote the resulting6n by n² exponent matrix by A.

## Reduce its entire left kernel to two difference equations

Let (a_g,b_g,c_g,d_g,e_g,f_g) be the coefficients of an arbitrary row functional vanishing on every column of A. These are row-functional coefficients, not original parameter values. The g=0 columns all use the second monomial, so for every h,

    a0+b_h+c_(−h)+d_h+f0=0.

Put

    K=a0+f0,   u_h=c_(−h),   w_h=u_h+d_h,
    r_g=a_g+e_g−K,          s_g=a_g+f_g−K.

Then b_h=−u_h−d_h−K. A first-choice column has zero row-functional value exactly when

    r_g+w_(h+g)−w_h=0.                         (I)

A second-choice column has zero row-functional value exactly when

    s_g+u_(h+g)−u_h=0.                         (II)

We show that w and u must both be constant.

### Noncyclic groups

Every g!=0 generates a proper cyclic subgroup H=<g> of G. Choose a coset t+H different from H; it avoids both0 and−g.

For g in X, equation(I) holds at every h except−g. It therefore holds around the whole coset t+H. Summing these |H| equations gives |H| r_g=0, hence r_g=0 over R (or C). Thus w_(h+g)=w_h at every h except possibly−g. The missing equality also follows by summing differences around H itself: all the other edges in that finite cycle already have zero difference. Hence w is invariant under translation by g.

The set X generates G because every nonzero element belongs to X or is the negative of an element of X. Therefore w is constant on G.

For g in Y, equation(II) holds at every h except0. Summing it on t+H gives |H| s_g=0, so s_g=0. The missing edge at0 is recovered around H exactly as before. Hence u is invariant under translation by every g in Y.

It remains to show Y generates G. Let M={g in G:2g!=0}. Since G has exponent>2, M is nonempty. The set Y has exactly one representative of each inverse pair in M, so <Y>=<M>. If t is an involution, choose v in M. Then v−t is also in M because 2(v−t)=2v!=0, and t=v−(v−t) belongs to <M>. Elements not of order1 or2 already lie in M. Consequently <M>=G, proving u is constant.

### Cyclic groups

Write group coordinates as0,...,n−1. Since n>=6, elements1,2 belong to X and −1,−2 belong to Y. (This argument also works for odd n>=5, whose dimension result is prior work.)

For g=1, equation(I) holds for h=0,...,n−2. Therefore

    w_j=w0−j r_1,  j=0,...,n−1.

For g=2, compare(I) at h=0 and h=n−1, both allowed. The two differences are −2r_1 and (n−2)r_1. They must be equal, so n r_1=0. Therefore w is constant.

For g=−1, equation(II) holds for h=1,...,n−1, giving

    u_j=u0+j s_(−1).

For g=−2, compare(II) at h=2 and h=1. The two differences are −2s_(−1) and (n−2)s_(−1), so n s_(−1)=0. Thus u is constant as well.

## Complete kernel and rank

Write the constants c_g=C and d_g=D. Then b_g=B is constant and K=−B−C−D. Every nonzero g has at least one first-choice and one second-choice column, so (I) and(II) now give

    e_g=f_g=K−a_g,   g!=0.

At g=0, the definition of K gives f0=K−a0, whereas e0 is completely free because no selected column uses w5_0.

Conversely, choose arbitrary a_g for all n elements, arbitrary constants B,C,D, and arbitrary e0. Set K=−B−C−D and define all remaining coordinates by these formulas. Every selected column then has zero row-functional value. These choices are independent: a_g,B,C,D and e0 can all be read off from the resulting functional. Hence the full left kernel has dimension exactly n+4.

It follows that

    rank_R A = 6n−(n+4)=5n−4.

The source's tropical lower bound and parameter-count upper bound now coincide, proving the claimed complex affine-cone dimension for every even-order finite abelian G with n>=6 and exponent>2.

## Checks, boundaries and unresolved part

The executable finite certificate reconstructs all n² coordinate exponent columns from the literal source map, verifies strict weight margins, checks an explicit n+4-vector integer left-kernel basis, and computes exact rational ranks for18 named groups. Positive controls include known odd groups, known finite even examples, and a larger noncyclic example C2×C2×C2×C4. Known elementary-two cases reproduce deficient ranks13,29,61 for orders4,8,16. These finite checks diagnose the construction; they are not the proof of its infinite quantifiers.

The proof uses exponent>2 specifically when showing <Y>=G. For elementary abelian2-groups Y is empty, so it supplies no constraint forcing u constant. Those groups are outside the theorem. We do not claim their variety dimension, a full resolution of the even-order conjecture, a four-state DNA identifiability theorem, a historical-first result, formal proof-assistant certification, or biological experimental evidence.

Proof gap list: none currently known within the stated class. Independent fresh-context review and current priority checking are still required; the researcher does not certify acceptance.
