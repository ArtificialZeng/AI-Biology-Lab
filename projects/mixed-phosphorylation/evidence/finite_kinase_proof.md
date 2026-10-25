# Proposed partial theorem: an explicit finite-kinase instability bound

> Public archive note: this derivation preserves research-stage wording. The current public claim and review status are defined in `README.md` and `claim.json`; internal AI checks are not external human peer review.


Researcher proof for fresh review, 2026-10-26. This file is a mathematical
argument and computational certificate dossier, not an independent referee
report or a publication manuscript.

**Statement.** Consider equations (2) and conservation laws (3) of
Conradi–Mincheva–Shiu, arXiv:1809.02886v2, with the Table 1 rates
(1,1,1,1,100,1,9/10,3,1,100), P_tot=8, and S_tot=40. For every

    K_tot >= B := 4660166/209569 < 23,

the unique positive equilibrium, restricted to its six-dimensional
stoichiometric compatibility class, has exactly two eigenvalues in the open
right half-plane and four in the open left half-plane, counting algebraic
multiplicity. In particular it is unstable for every K_tot>=23.

This establishes the assertion of original Conjecture 6.2(1) with the fixed
witness totals (8,40), with a sufficient finite bound. No optimality of B is
claimed. It does not establish Conjecture 6.2(2), periodic-orbit existence,
periodic-orbit hyperbolicity, or any exact output invariance. The original
two-part problem therefore remains unresolved.

## Equilibrium parameterization and coverage of all kinase totals

Write z=x9 and p=8-(1009/9)z. The steady flux equations give, in the species
order of the source,

    (x1,...,x9) = (r, 200z/r, 100z, 100z, (19/9)z/p,
                   p, (1000/9)z, (101/3)z/p, z),
    G(z) = (2809/9)z + (322/9)z/p,
    r(z) = 40-G(z),
    K(z) = 200z + 200z/r(z).

All nine equilibrium equations and all three conservation identities are
checked as rational-function identities by `src/finite_kinase_exact.py`.
They can also be verified by substitution in the source's nine equations.
The same substitution is separately checked using the directly entered nine
ODEs in `src/verify_finite_kinase.py`.

On 0<z<72/1009, G is strictly increasing, G(0)=0, and G tends to infinity at
the upper endpoint. Thus there is a unique alpha in this interval with
G(alpha)=40. The physical equilibria are precisely 0<z<alpha. On this range
r>0, r'<0, and

    K'(z)=200 + 200/r + 200z*G'(z)/r^2 > 0.

Moreover K(z) tends to zero at z=0 and to infinity as z increases to alpha.
Consequently K maps (0,alpha) bijectively onto (0,infinity), so no finite
positive kinase totals are omitted. The flux equations force every positive
equilibrium to have this parameterization, proving uniqueness among positive
equilibria. Uniqueness of the steady state in the whole compatibility class
is also the source's Proposition 3.1; it is prior knowledge, not a novelty claim.

For the explicit polynomial

    Q(z)=2834281*z^2-568386*z+25920,

one has r(z)=Q(z)/[9(72-1009z)]. With L=7/100 and
U=70122383355/10^12, exact rational evaluation gives

    Q(L)=209569/10000 > 0,
    Q(U)=-6175117557936279/40000000000000000000000 < 0,
    Q'(U)<0,   0<L<U<72/1009.

Since Q' is increasing, Q'<0 throughout [0,U]. The root alpha is therefore
in (L,U), and Q>0 throughout [L,alpha). Direct evaluation gives
K(L)=4660166/209569=B. The remaining task is to count eigenvalues for every
z in [L,alpha), which covers precisely all K>=B.

## Exact characteristic polynomial

Eliminate x2=K-x3-x4, x6=8-x7-x9, and
x8=40-x1-x3-x4-x5-x7-x9, and retain coordinates
(x1,x3,x4,x5,x7,x9). This is an affine coordinate chart on the compatibility
class. Differentiate at fixed K before substituting K=K(z) and the
equilibrium concentrations above. This order is essential: K(z) is a
parameterization of equilibria, not a state-dependent total.

Let J6(z) be this restricted Jacobian and write

    det(lambda I-J6(z)) = lambda^6+c1(z)lambda^5+...+c6(z).

`src/finite_kinase_exact.py` constructs the entire system from the ten reaction
arrows, verifies the conservation matrix, and computes these rational
coefficients. They are serialized completely in
`evidence/finite_kinase_algebra.json`, together with each Hurwitz determinant
Delta_m=det[c_(2i-j+1)]_(i,j=0,...,m-1), where c0=1 and coefficients outside
0,...,6 are zero. Equivalently, the first two Hurwitz rows are (c1,1,0,...) and
(c3,c2,c1,1,...).

The separate checker differentiates the published nine ODEs directly and
verifies the exact rational-function identity

    det(lambda I-J9(z)) = lambda^3 det(lambda I-J6(z)).

Thus the counted roots are precisely the six compatibility-class
eigenvalues; the three conservation zeros are excluded.

## A rational certificate over a continuous interval

Write each determinant as N_j(z)/D_j(z), exactly as stored in the algebra
file. Let A=1009z-72. The denominators, their signs on [L,alpha), and the
degrees of the numerator polynomials are:

| j | D_j | sign(D_j) | degree(N_j) |
|---|---|---|---|
| 1 | 90 A Q | - | 4 |
| 2 | 72900 A^3 Q^2 | - | 10 |
| 3 | 177147000 A^5 Q^3 | - | 17 |
| 4 | 71744535000 A^7 Q^4 | - | 24 |
| 5 | 871696100250 A^8 Q^5 | + | 30 |
| 6 | 94143178827 A^7 Q^6 | - | 34 |

Here A<0 and Q>0 on the physical interval. The certificate proves
-N_j>0 for j=1,...,5 and N_6>0 on the larger closed interval [L,U]. This
uses numerator polynomials only: no determinant is evaluated across the pole
at alpha, and no denominator sign beyond alpha is asserted.

For completeness, the positivity certificate uses the following elementary
identity. If p has degree n and

    p(L+(U-L)t)=sum_(j=0)^n a_j t^j,
    b_k=sum_(j=0)^k a_j*binom(k,j)/binom(n,j),

then

    p(L+(U-L)t)=sum_(k=0)^n b_k*binom(n,k)*t^k*(1-t)^(n-k).

Every basis term is nonnegative on [0,1], and the basis terms sum to one.
Thus all b_k>0 imply p>0 on the entire interval, including both endpoints.
`src/finite_kinase_signs.py` computes these coefficients using only rational
operations. All 125 coefficients (5+11+18+25+31+35), including their full
numerators and denominators, are in
`evidence/finite_kinase_sign_certificate.json`. Every coefficient is strictly
positive; there are no floating-point sign decisions.

The separate checker verifies each claimed expansion in the reverse
direction using Python's standard-library Fraction arithmetic at n+1 distinct
rational points. Both sides have degree at most n, so these exact equalities
establish the polynomial identity. It also checks that each signed numerator
and factored denominator reproduce the stored determinant. No assertion is
inferred merely from a JSON success flag.

It follows that, for every z in [L,alpha),

    (sign Delta_1,...,sign Delta_6) = (+,+,+,+,-,-).

## Root count and conclusion

Since all the Hurwitz determinants are nonzero, the Routh array is regular.
For this monic degree-six polynomial its first column is

    1, Delta_1, Delta_2/Delta_1, Delta_3/Delta_2,
    Delta_4/Delta_3, Delta_5/Delta_4, Delta_6/Delta_5.

Its signs are (+,+,+,+,+,-,+), with exactly two sign changes. The Routh
root-count theorem gives exactly two open-right-half-plane roots. The
regular array with nonzero last determinant excludes imaginary-axis roots;
the other four roots are therefore in the open left half-plane. The
parameterization then proves the stated result for every K>=B. QED.

## Dependencies, reproduction, and limitations

Dependencies: the published reaction equations and rate table; exact rational
linear algebra and differentiation; the displayed elementary Bernstein
identity; the classical Routh root-count theorem. No singular perturbation,
Fenichel theorem, fitted extrapolation, or numerical periodic-orbit data are
needed for this theorem. No known gap in this argument is being deferred;
the fresh referee must still check the formulation, algebra, and root count.

Run, from the project directory:

    python src/finite_kinase_exact.py
    python src/finite_kinase_signs.py
    python src/verify_finite_kinase.py

The last command writes `evidence/finite_kinase_verification.json`. These are
researcher consistency checks, not a separate review context. Current
literature limitations are recorded in `literature/finite_threshold_check.md`.
Neither a manuscript nor a PDF is authorized before independent acceptance.
