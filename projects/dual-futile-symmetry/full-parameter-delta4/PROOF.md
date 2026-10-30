# Determinant-conditioned positivity of the fourth Hurwitz determinant

Let all sixteen entries of the parameter vector in `matrix.py` be strictly
positive and independent. Define c0=1,c1,...,c6 by
det(lambda I-G)=sum(k=0,...,6) c_k lambda^(6-k). Let Delta4 be the leading 4x4
Hurwitz determinant. We prove

**Lemma. If c6>=0, then Delta4>0.**

Set

```
A = e4+v5+v6,
B = f10+y11+y12,
L = a+c+u2+u3+x8+x9,
Q = (L+f7)*A^3 + (e1*(v5+v6)+b4*v6)*A^2
  + (L+e1)*B^3 + (f7*(y11+y12)+b10*y12)*B^2
  + A*B*(A+B)^2.
```

Every displayed summand of Q is strictly positive. In particular Q>0; its
expansion has 266 monomials, degree four, and integer coefficient sum 756.

The finite exact algebraic assertion used below is

```
R := Delta4 - Q*c6 has 1,712,869 nonzero coefficients, all positive integers.
[u3*v6^2*x9^3*y12^4] R = 1.
```

Here [m] means the coefficient of monomial m. All sixteen variables remain
formal indeterminates in this assertion. Thus this is coefficient positivity
in the full polynomial ring, not a parameter sample or a positivity claim
after imposing additional relations.

The accompanying solver-independent checker `verify_certificate.py` supplies
a complete reproducible verification of that finite assertion. It constructs
c_k as the sum of the determinants of all k-by-k principal submatrices of -G.
These are the coefficients of det(lambda I-G) by the principal-minor expansion
of the determinant. It constructs the matrix

```
[ c1 c3 c5  0 ]
[  1 c2 c4 c6 ]
[  0 c1 c3 c5 ]
[  0  1 c2 c4 ]
```

and computes its determinant by the 24-term Leibniz formula. It reads the
266-monomial integer Q from `certificate.json`, checks equality to
the compact formula above, computes R, checks every nonzero coefficient, and
extracts and checks the stated strict-positive coefficient witness. FLINT
uses arbitrary-precision integers throughout; this checker imports neither
SciPy nor the discovery cache and invokes no linear-programming solver.
`evidence/exact-verification.json` records its successful run, degrees, term counts and
coefficient sums. A separately implemented discovery path computed the
characteristic polynomial as a 6x6 determinant and used the source's closed
Delta4 formula; it gave the same coefficient counts and identity. This is
computational replication, not external human review.

Because the parameters are strictly positive, every monomial in R is strictly
positive. In particular the exhibited coefficient gives

```
R >= u3*v6^2*x9^3*y12^4 > 0.
```

Consequently Delta4=Q*c6+R>0 whenever c6>=0. This proves the lemma, including
c6=0 within the positive orthant. The proof's computational dependency is the
finite exact polynomial check just specified; no unresolved mathematical gap
is being hidden by numerical tolerances. A separate substantive AI review accepted this partial claim and independently
reconstructed the source model and polynomial identity; see `REVIEW.md`.
External human review and priority confirmation remain pending.

The multiplier also equals the source's negative-part quotient:

```
negative_part(Delta4) = a*c*e4*f10*v6*y12 * Q.
```

The checker verifies this identity too. This identifies Q intrinsically from
the known negative terms rather than only as an opaque fitted LP solution.

Relationship to the original question: the source already proves Delta1,
Delta2,Delta3>0 and c6>0 implies c5>0. This lemma removes Delta4 from the
unproved Hurwitz inequalities on the full determinant-positive region.
Delta5>0 remains unproved, and Conjecture 3.1 remains unresolved. The source
mass-action theorem and prior symmetric classification retain their original
credit. No new conjecture, full stability theorem, or periodic-orbit statement
is asserted here.
