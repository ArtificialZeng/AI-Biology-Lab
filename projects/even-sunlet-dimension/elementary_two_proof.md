# Elementary-two 3-sunlet dimensions by Hamming order

## New theorem and credited full-scope corollary

Let G=(C2)^r, with integer r≥3 and n=2^r. For the general group-based model B={id} on the Cox–Gross–Martin 3-sunlet network, the complex affine cone has dimension

    dim V_N^G = 5n−4.

The source already established this dimension in the finite cases r=3,4,5 by computations of full-rank chambers. The proposed new contribution is a single explicit strict chamber and complete characteristic-zero rank proof for every r≥3. The finite prefix below is diagnostic, not the proof.

Combined with the published odd-order theorem and the previously submitted R11 theorem for even-order groups of exponent greater than two, this gives

    dim V_N^G = 5|G|−4

for every finite abelian group G with |G|≥5, under B={id}. This is a coherent extension of existing lab PR15 / author PR1, not a separate contribution or split submission. It does not resolve a nontrivial-B interpretation of any broader conjecture, a topology-identifiability question, or a biological experiment. The R11 proof is retained unchanged as `prior-r11-proof.md`; source-notation ambiguity is preserved separately.

## Literal model and dimension bounds

In the elementary-two group identify elements with integers x=0,...,n−1, using XOR for addition. Every element is its own inverse. For g,h∈G the source coordinate is

    phi_(g,h,g+h) = w1_g w2_h w3_(g+h)
                    (w4_(g+h) w5_g + w4_h w6_g).

The 6n parameter coordinates are independent complex variables. The variety is the Zariski closure of this polynomial image, considered as an affine cone. The primary source, Lemma2.1, gives the upper bound 5n−4. Its Section2.2 gives a lower bound by the rank over R of any strict tropical exponent matrix, whose columns are the exponent vectors of the uniquely minimum-weight monomials.

Both uses can be checked directly here. For the upper bound replace w1_g w5_g by alpha_g and w1_g w6_g by beta_g. At g=0 only gamma_0=alpha_0+beta_0 occurs. On the dense parameter set where w2_0,w3_0,w4_0 are nonzero, normalize these three coordinates to1 and absorb their product into alpha,beta,gamma_0. The remaining count is

    (n−1)+(n−1)+(n−1)+2(n−1)+1 = 5n−4.

The image of the dense set has the same closure as the full image. For the lower bound, substitute w_i=t^(lambda_i) in the logarithmic Jacobian. Divide each output column by its uniquely smallest t-power. Its constant column is the selected exponent vector. A nonzero minor of that exponent matrix is therefore the nonzero leading coefficient of the corresponding Jacobian minor. Characteristic-zero generic Jacobian rank equals image dimension. This also explains why no arbitrary tie selection is allowed.

## Strict Hamming-order chamber

Let wt(x) be the number of nonzero bits of x. Put

    nu_x = 4(n wt(x)+x),
    mu_0 = −1,             mu_g = 1 for g≠0.

Actual strictly positive edge-coordinate weights are

    lambda1_x=lambda2_x=lambda3_x=1,
    lambda4_x=nu_x+1,     lambda5_x=2,     lambda6_x=2+mu_x.

The first-minus-second monomial weight is exactly

    nu_(g+h)−nu_h−mu_g.

At g=0 it is1, so all g=0 columns select the second term. At g≠0, the two nu values are distinct multiples of4 and mu_g=1. Hence there is no tie, and the first term is selected precisely when nu_(g+h)<nu_h. Moreover nu induces exactly the total order

    x≺y  iff  (wt(x),x) is lexicographically smaller than (wt(y),y).

Indeed a one-unit Hamming-weight difference contributes n, larger than any possible difference between two integer labels in0,...,n−1. Thus for every g≠0, the paired columns (g,h) and(g,h+g) select different terms: the column whose second label is the larger element selects the first term; the column whose second label is the smaller element selects the second term. The 6n×n² selected matrix is denoted A.

## Scalar Hamming-order lemma

Let F:G→R (or C). Suppose there are numbers t_g for every g≠0 such that whenever x≺y,

    F(y)−F(x) = t_(x+y).                         (L)

Then F is constant for every r≥3.

Proof. Index bits by i=0,...,r−1, with integer label e_i=2^i. If x has bit i equal to0, then x≺x+e_i because the latter has one more nonzero bit. Thus (L) gives

    F(x+e_i)−F(x)=t_(e_i)=a_i.

Adding the nonzero bits of any x successively yields

    F(x)=F(0)+sum_i a_i x_i.                     (A)

For i<j compare the pairs(0,e_i+e_j) and(e_i,e_j). In the first pair0 is smaller by Hamming weight; in the second pair e_i≺e_j by their integer labels. Both pairs have XOR e_i+e_j, so(L) and(A) give

    a_i+a_j = t_(e_i+e_j) = a_j−a_i.

Therefore2a_i=0. Characteristic zero implies a_i=0 for every i<r−1, by taking j=r−1.

Since r≥3, choose two distinct lower indices i,j<r−1 and put k=r−1. The pairs(0,e_i+e_j+e_k) and(e_k,e_i+e_j) have the same XOR. In the second pair e_k≺e_i+e_j because its Hamming weight is1 rather than2, regardless of its larger integer label. Consequently

    a_i+a_j+a_k = t_(e_i+e_j+e_k) = a_i+a_j−a_k.

Thus2a_k=0 and a_k=0. Every a_i vanishes, proving the lemma. No induction on a tested finite prefix is involved.

## Full left kernel

Let(a_g,b_g,c_g,d_g,e_g,f_g) be an arbitrary row functional vanishing on every column of A. These coefficients are independent of the original model parameter values. From the g=0 second-term columns,

    a_0+b_h+c_h+d_h+f_0=0   for every h.

Define

    K=a_0+f_0,  U_h=c_h,  W_h=c_h+d_h,
    r_g=a_g+e_g−K,         s_g=a_g+f_g−K.

Then b_h=−W_h−K. The first-term equation is exactly

    r_g+W_(g+h)−W_h=0,                          (I)

and the second-term equation is exactly

    s_g+U_(g+h)−U_h=0.                          (II)

For any distinct x≺y, let g=x+y. The column(g,y) has g+y=x≺y and selects the first term; equation(I) gives

    W_y−W_x=r_g.

The column(g,x) selects the second term; equation(II) gives

    U_y−U_x=−s_g.

Therefore both W and U satisfy the scalar lemma, with t_g=r_g and t_g=−s_g respectively. Both are constant. It follows that c_g=C and d_g=D are constant, b_g=B is constant, and

    K=−B−C−D.

For every g≠0 there are columns of both types (exactly n/2 of each), so(I) and(II) now force

    e_g=f_g=K−a_g   for g≠0.

At0, f_0=K−a_0; the coordinate e_0 is unused because every g=0 column selects the second term. Thus e_0 remains completely free.

Conversely choose arbitrary a_g for all n group elements, arbitrary B,C,D, and arbitrary e_0. Define K=−B−C−D, set b=B,c=C,d=D, and define the remaining e,f by the displayed formulas. Direct substitution makes every selected column vanish. The n+4 choices are independent, since each a_g,B,C,D,e_0 can be read off individually. This proves that the entire left kernel, not just a displayed subspace, has dimension n+4. Hence

    rank_R A = 6n−(n+4)=5n−4.

Together with the lower and upper bounds, this proves the new theorem.

## Boundary, controls, and full-scope combination

The scalar lemma uses r≥3 exactly when it chooses two lower bits. For r=2, a nonconstant highest-bit function satisfies(L), so this chamber has extra kernel. Its exact matrix rank is14 at n=4; this is not the variety dimension, which is the known exceptional15. For r=1, the selected matrix has rank4. These are diagnostic boundaries, not new dimension claims or counterexamples to the asserted theorem. The integer-order chamber nu_x=4x also has extra kernel and is a failed simpler construction, not an obstruction to every strict chamber.

For the combined corollary, let G be any finite abelian group with |G|≥5. If |G| is odd, apply the published Cox–Gross–Martin theorem. If |G| is even and exponent>2, apply the R11 theorem. If its exponent is2, G is a vector space over F2, hence isomorphic to(C2)^r; |G|≥5 implies r≥3 and the present theorem applies. These cases exhaust G. Isomorphism merely reindexes the coordinates of the parameterization and does not change dimension.

No known gap remains in the stated characteristic-zero B={id} argument. Fresh mathematical/contribution review is requested. No priority certification, external human review, proof-assistant certification, or trusted-runner acceptance is claimed.

## Credited corollary for the source's binary level-1 class

Let N be a binary semi-directed level-1 phylogenetic network in the class used in Cox–Gross–Martin Theorem3.8, with at least three leaves, m edges, and c cycles. Count m after suppressing the root and before any contraction of reticulation outgoing edges; cycles have length at least3 and the source's binary/no-degree-two conventions apply. For every finite abelian G with q=|G|≥5 and B={id},

    dim V_N^G = (q−1)(m−c)+1.

This is the source's existing dimension induction with the odd-order restriction removed by the completed3-sunlet base case. It is recorded as a credited consequence, not a separate new method or publication.

The required published facts are Gross–Krone–Martin2024 Lemma11 (binary-tree dimension), Theorem17 (k-sunlet dimension for k≥4 and |G|>2), and Corollary10 (cut-edge dimension sum). These statements have no odd-order hypothesis in this range. They were checked in the publisher full text at https://link.springer.com/article/10.1007/s11538-024-01314-z .

For completeness, follow Cox–Gross–Martin2025 Theorem3.8's induction on nontrivial cut edges (excluding pendant leaf edges). The indecomposable binary pieces are the3-star tree and k-sunlets. Writing ell=q−1, their dimensions are3ell+1 for the3-star,5ell+1 for the3-sunlet by the combined theorem above, and ell(2k−1)+1 for k≥4 by the published theorem. Each equals ell(m−c)+1 in the stated uncontracted edge convention. When a nontrivial cut edge splits N into N1,N2, the published cut-edge formula is

    dim V_N^G = dim V_N1^G + dim V_N2^G −q.

The pieces satisfy m1+m2=m+1 and c1+c2=c. By induction the right side is

    ell(m1−c1)+1 + ell(m2−c2)+1 −(ell+1)
      = ell(m−c)+1.

This proves the corollary in exactly that source network/model class. Arbitrary nonbinary networks, parallel-edge/2-cycle extensions, nontrivial B, and full topology identifiability are not asserted.
