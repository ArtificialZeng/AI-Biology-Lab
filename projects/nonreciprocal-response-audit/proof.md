# Finite response jets and complete response curves are different data

## Scope and source

Section8 of Franco, Kepka and Velázquez, *Characterizing the Detailed Balance
Property by Means of Measurements*, ARMA250:25(2026), studies four-state
column Markov generators with complete directed graph. Following its displayed
formulas, the two forward responses here are

    R12(t)=e2^T exp(tA)e1, R13(t)=e3^T exp(tA)e1.       (1)

This avoids any ambiguity with the article's other response-index conventions.
Let F(A) contain their derivatives of orders1,2,3 at0, in that order for each
response. We audit three specific subsidiary assertions in that section:
its full-rank computation at A0=J-4I, the sufficiency of these six coefficients
to characterize the complete responses, and the literal unrestricted-curve
conclusion for nonreciprocal measurements. We also give a limited positive
repair of the local finite-jet argument. The main cut-vertex/pathwise-detailed-
balance theorem is outside this note. The classical real-spectrum consequence
of detailed balance is already implicit in the source's Proposition2.2 and is
not claimed as a new theorem or detection method.

## 1. Correct derivative at the uniform symmetric generator

Parameterize symmetric generators by six positive off-diagonal rates in order
(12,13,14,23,24,34); each diagonal is minus the sum of its column's off-diagonal
entries. Put A0=J-4I, where every entry of J equals1. A symmetric tangent H has
zero row and column sums. Hence JH=HJ=0 and A0H=HA0=-4H. Differentiating a
matrix power gives

    D(A^k)[H] = sum_{l=0}^{k-1} A^l H A^(k-1-l),
    D(A^k)|A0[H] = k*(-4)^(k-1) H.                   (2)

The true six-by-six derivative, in the stated rate order, is therefore

    [ 1  0  0  0  0  0 ]
    [-8  0  0  0  0  0 ]
    [48  0  0  0  0  0 ]
    [ 0  1  0  0  0  0 ]
    [ 0 -8  0  0  0  0 ]
    [ 0 48  0  0  0  0 ],                            (3)

which has rank2. This concerns precisely the symmetric six-rate restriction
used in the printed computation. It is not a claim that the full nine-
dimensional reversible manifold has rank2. Rank failure at one point alone
does not establish that F has no open image.

## 2. A positive finite-jet repair

Take the symmetric generator

    A*=[[-6,  1,  2,  3],
        [ 1,-10,  4,  5],
        [ 2,  4,-12,  6],
        [ 3,  5,  6,-14]].                            (4)

Its rates are(1,2,3,4,5,6). Direct exact differentiation gives

    DF|A* =
    [  1,  0,   0,  0,  0,  0]
    [-18,  3,   4,  1,  2,  0]
    [255,-87,-136,-35,-80, -1]
    [  0,  1,   0,  0,  0,  0]
    [  2,-22,   4, -1,  0,  1]
    [-56,390,-152, 37, -4,-47].                       (5)

Its determinant is800. The ordinary inverse-function theorem therefore makes
F a local diffeomorphism at this strictly positive symmetric point. Its image
contains an open subset of R6. In particular there is no nonzero polynomial
equality in those six coordinates that holds for every positive symmetric
generator: such a polynomial would vanish on that open subset and hence be
identically zero. The determinant polynomial is not identically zero; its
nonzero set is dense in the positive-rate cone, so every neighborhood of the
uniform-rate point contains a point of full rank and maps some smaller open
set onto an open set. This does not prove that F(A0) itself is interior to
every neighborhood image, or that an entire such image is open.

Thus a local finite-jet absence-of-equality-constraints statement has a valid
repair, even though the printed derivative at the chosen point is not correct.
This cannot be promoted to a statement that complete response curves have no
necessary restrictions.

## 3. Same derivatives through order3, different complete curves

Let e1=(1,0,0,0)^T, v=(0,11,-16,7)^T and b=(0,1,-1,0)^T. The matrix A* in(4)
satisfies

    v^T e1=v^T A*e1=v^T A*^2 e1=0,
    v^T A*^3 e1=1112.                                (6)

Define the second generator

    A~=A*+(1/1000)b v^T
      =[[-6, 1,       2,        3],
        [ 1,-9989/1000,498/125,5007/1000],
        [ 2,3989/1000,-1498/125,5993/1000],
        [ 3,5,        6,       -14]].                 (7)

Its column sums are zero because 1^T b=0. All its off-diagonal entries are
strictly positive, so it is an ergodic complete-graph Markov generator.
Equation(6) and induction give

    A~^k e1=A*^k e1 for k=0,1,2,3,
    (A~^4-A*^4)e1=(0,139/125,-139/125,0)^T.            (8)

Consequently even all four state responses have the same first three
derivatives from the same input, while both responses in(1) differ at order4.
Their six coefficients alone do not characterize the full functions when the
underlying generator is unknown.

A* is reversible with uniform stationary distribution. A~ is not reversible:
the product of rates around1->2->3->1 minus the reverse product is1/100.
Detailed balance would make these products equal by cancellation of the
positive stationary weights. Thus these finite data also occur in both a
reversible and a nonreversible four-state realization. This is a source-
specific finite-jet witness, not a claim that every possible finite jet is
ambiguous, or that a complete curve distinguishes every nonreversible model.

## 4. A complete positive network whose forward curves exclude reversibility

Consider

    C=[[-4, 1, 1, 2],
       [ 2,-4, 1, 1],
       [ 1, 2,-4, 1],
       [ 1, 1, 2,-4]].                               (9)

All off-diagonal rates are positive and all column sums are zero. The
characteristic polynomial is s(s+6)(s^2+10s+26). Exact inversion gives

    Laplace(R12)(s)=(2s^2+18s+39)/[s(s+6)(s^2+10s+26)],
    Laplace(R13)(s)=( s^2+13s+39)/[s(s+6)(s^2+10s+26)]. (10)

Each numerator is coprime to the denominator. In particular the remainders
modulo q(s)=s^2+10s+26 are -2s-13 and3s+13, respectively; neither shares a
factor with this irreducible real quadratic. Each measured curve therefore
has genuine nonreal poles -5+i and-5-i. These statements concern the exact
analytic responses, not extraction of poles from noisy finite observations.

Recall the classical obstruction. For a finite reversible column generator B
with positive stationary vector pi and D=diag(pi), D^(-1/2) B D^(1/2) is real
symmetric. It diagonalizes with real eigenvalues; every entry of exp(tB), or
any fixed linear combination of entries, is a finite sum of real exponential
terms. Its rational Laplace transform has only real poles. It cannot equal
either curve(10), by uniqueness of rational functions on their common
right-half-plane domain. Thus one of these exact forward curves already
excludes every finite reversible realization of this measurement, including
ones with more hidden states. Both curves together also exclude one.

This applies the known spectral restriction to a concrete positive source-
domain witness. It refutes the literal assertion that such forward responses
cannot provide any necessary restriction detecting failure of detailed
balance. It does not imply that nonreciprocal measurements are sufficient to
establish detailed balance in every model, and does not affect the source's
reciprocal-measurement theorem.

## Evidence and contribution limits

All matrices, rate inequalities, derivative ranks, determinants, moment
equalities and reduced rational transforms are exactly checkable. The general
steps use elementary differentiation, the inverse-function theorem and the
classical reversible spectral theorem. The contribution is a precise audit
and repair of subsidiary source claims with fully positive witnesses; none
of these general mathematical tools is claimed new. A rank2 observation by
itself would be insufficient; the nonzero determinant, moment witness and
full-curve pole obstruction establish distinct claims with distinct scopes.

The exact mathematical scope received a separate fresh-context AI review.
Substantial AI assistance was used. Trusted runner provenance remains missing;
the actual formal gate failed. No manuscript or PDF is distributed, and no
human specialist confirmation, priority certificate, proof-assistant certificate
or authenticated release-gate PASS is claimed. See provenance/review-receipt.json.
