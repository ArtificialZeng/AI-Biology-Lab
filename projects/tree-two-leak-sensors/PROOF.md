# Theorem: one additional sensor for a two-leak tree

Under the assumptions in SCOPE.md, the model is generically locally
structurally identifiable if and only if P(r,s) and P(a,b) share at least one
edge. If they share no edge, every positive parameter point lies on a
nonconstant local curve of positive parameters with the same two transfer
functions. The latter is an exact ambiguity, not a numerical rank observation.

## 1. The root-output map has precisely one generic missing direction

There are 2(n-1)+2=2n parameters. Put

    D(z)=det(zI-A),  N_r(z)=det((zI-A) with row/column r deleted).

D is monic of degree n and N_r is monic of degree n-1. Their nonleading
coefficients therefore form a polynomial map C_r with at most 2n-1
coordinates. Its generic rank is exactly 2n-1. One may invoke Bortner et al.,
Proposition 5.7 directly. Alternatively, set ell_b=0: the resulting single-leak
tree with input/output at r is generically locally identifiable by their
Theorem 5.2/Lemma 5.5. A nonzero (2n-1)-minor of that restriction is a nonzero
minor polynomial of the two-leak map. It stays nonzero generically in the
positive orthant. The use of a zero leak here only proves a polynomial is not
identically zero; the theorem's parameter domain has both leaks positive.

The source's Jacobian criterion (Proposition 2.15) applies: the graph is
strongly connected and has an input. The root numerator/denominator are
generically coprime. Indeed, C_r has 2n-1 coordinates and a full-row-rank
Jacobian, so its complex polynomial image is Zariski dense in the space of
monic pairs of these degrees. The universal resultant on that coefficient
space is a nonzero polynomial (for example, z^n+1 and z^(n-1) are coprime).
Its pullback therefore is not identically zero. The complement of that zero
set meets the positive orthant generically. On this set the coefficient map
and root transfer data have the same differential rank, and ker(dC_r) is
one-dimensional.

## 2. An exact positive similarity curve spans that missing direction

For a vector h, call it outgoing-rate harmonic at a nonleak vertex j when

    sum_{i~j} k_ij (h_i-h_j) = 0.                    (1)

There is a unique h satisfying (1) for j outside {a,b} and h_a=0,h_b=1.
Existence and uniqueness follow, for example, by solving the finite absorbing
Markov-chain Dirichlet problem with strictly positive transition weights; or
by the maximum principle, which gives uniqueness and nonsingularity of the
interior Dirichlet matrix. Its coordinates lie in [0,1]. Define

    d_j(t)=1+t(h_j-h_r),  Q(t)=diag(d_j(t)),
    A(t)=Q(t) A Q(t)^(-1).                          (2)

For |t| sufficiently small every d_j(t)>0. The diagonal entries of A(t) are
unchanged; each edge j->i has new positive rate
k_ij(t)=k_ij d_i(t)/d_j(t). For a nonleak column j,

    sum_{i~j} k_ij d_i(t)=d_j(t)sum_{i~j}k_ij

by (1), so its column sum remains zero. At a,b define the new leak as minus
the new column sum. These leaks depend continuously on t and start strictly
positive, hence stay positive on a sufficiently small interval. Consequently
(2) is a curve inside exactly the stated model family.

Because d_r(t)=1,

    H_v(z;t)=d_v(t) H_v(z;0)                       (3)

for every output v, and D(z;t)=D(z;0). In particular H_r is constant. The
curve's parameter derivative is nonzero: if every edge-rate derivative
vanished, h_i=h_j on every edge and connectedness would contradict
h_a=0,h_b=1. Thus at every generic point it spans the one-dimensional kernel
from Step 1. This argument does not require r to be distinct from a or b.

## 3. The voltage difference is exactly the path-intersection test

Let P=P(a,b). Each component hanging off P attaches to exactly one path
vertex. Its harmonic extension is constant at the attachment value: apply
the maximum principle to that component with its one fixed boundary vertex.
Let pi(v) be v's attachment vertex, with pi(v)=v for v on P. Then
h_v=h_pi(v).

Along P, write successive vertices as a=v_0,...,v_m=b. At every interior
v_j, the off-path contributions in (1) vanish, leaving

    k_{v_(j-1),v_j}(h_(v_(j-1))-h_(v_j))
      + k_{v_(j+1),v_j}(h_(v_(j+1))-h_(v_j))=0.

Thus consecutive voltage increments have positive ratios. Their sum is
h_b-h_a=1, so every increment is strictly positive. Therefore

    h_s != h_r  iff  pi(s) != pi(r).

The two attachment vertices differ exactly when the path P(r,s) traverses a
nonempty segment of P(a,b), equivalently when those paths share an edge.
An intersection consisting of one vertex is insufficient.

## 4. The second output removes the missing direction exactly in that case

If h_s=h_r, (3) fixes H_s as well as H_r along the nonconstant positive curve
(2). The model is unidentifiable, including at generic points.

If h_s!=h_r, differentiate (3) at zero:

    d/dt H_s(z;t)|_(t=0) = (h_s-h_r) H_s(z;0) != 0.

H_s is not the zero rational function: the leading nonzero coefficient of
its large-z expansion is the product of the strictly positive rates along
the directed r-to-s path. Hence the second-output data do not annihilate the
one-dimensional root-output kernel. The joint coefficient map has rank 2n
at generic positive points and the model is generically locally identifiable.

The common-denominator description remains legitimate at these generic
points by Step 1; equality of rational functions is equivalent to equality
of their numerator coefficients there. No claim is made about nongeneric
rank-drop loci or global uniqueness. This completes the proof.

## Immediate sensor-design corollary

With the root measurement retained, exactly one additional measurement is
necessary and sufficient to make any fixed two-leak tree generically locally
identifiable, if its location can be chosen freely. Necessity is the known
root-only rank deficit. For sufficiency choose whichever of a,b has a
different projection on P(a,b) from r. More generally, a set of additional
sensors succeeds exactly when at least one has projection different from r.

## Dependency and gap audit

Dependencies: Bortner et al. Theorem 5.2 or Proposition 5.7, their standard
Jacobian/input-output criterion, finite positive Dirichlet maximum principle,
and the resolvent identity for diagonal similarity. All hypotheses used are
listed in SCOPE.md. No unproved mathematical gap is currently known. A separate AI review accepted this mathematical scope. External human review and priority clarification remain pending.
