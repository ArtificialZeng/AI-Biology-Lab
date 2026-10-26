# Precise reading

Let T=(V,E) be a finite tree with n>=2. For every ordered adjacent pair
j->i, k_ij>0 is an independent rate. Let a!=b be designated leak vertices,
with independent leak rates ell_a,ell_b>0. Define A_ij=k_ij on edges, A_ij=0
elsewhere off the diagonal, and

    A_jj = -sum_{i~j} k_ij - 1_{j=a}ell_a - 1_{j=b}ell_b.

The model is x'=Ax+e_r u with measured, individually calibrated states x_r
and x_s, where s!=r. Equivalently, we study the common-denominator transfer
pair H_r(z)=e_r^T(zI-A)^(-1)e_r and
H_s(z)=e_s^T(zI-A)^(-1)e_r. Their generic coefficient-map rank is the
structural identifiability criterion. Initial conditions are not extra known
parameter constraints; this is the usual generic input-output setting in the
cited compartment-model literature.

P(u,v) denotes the unique undirected tree path, including its endpoints. The
criterion concerns an edge, not merely a vertex, in the intersection of two
paths. Vertices r,s may coincide with a,b; only a!=b and r!=s are required.
All edge and designated leak rates are in the open positive orthant.

Nearest prior result X: one output at r and <=1 leak is generically locally
identifiable; with two leaks its coefficient map has generic rank 2n-1.
Proposed delta Y: a complete graphical criterion for the one additional sensor
to recover the missing rank, including a positive exact ambiguity curve when
it cannot. Verification Z: a general harmonic-similarity proof plus exact
finite-field rank checks over all small tree placements and an explicit
rational positive witness on a tree that is neither a path nor a star.

Original review Questions 5.2 and 5.4 remain unresolved in their full scope.
This subsidiary family is a candidate result-note; novelty is not certified.
