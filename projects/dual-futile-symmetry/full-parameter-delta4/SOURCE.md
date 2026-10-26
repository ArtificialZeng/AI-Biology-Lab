# Primary-source comparison, inspected 26 October 2026

Nicola Vassena, *Sequential and distributive dual futile cycle: Hopf bifurcation
can occur under parameter-rich kinetics but cannot occur under mass action
kinetics*, arXiv:2608.27081v2, revised 14 September 2026.

Primary record: https://arxiv.org/abs/2608.27081
Primary full text: https://arxiv.org/html/2608.27081v2
Primary code: https://github.com/nvassena/dualfutile_hopfnohopf

The arXiv abstract page inspected today lists v2 as its latest revision. Section 3
still states Conjecture 3.1: every strictly positive parameter vector with
det(G_red)>0 gives a Hurwitz matrix. That is the original target of this pass.

The displayed reduced matrix immediately before equation (2.7) is the input.
The row/column order is (Q1,Q2,U,V,Y,X), and the parameter order used here is
(a,b4,b10,c,e1,e4,f7,f10,u2,u3,v5,v6,x8,x9,y11,y12). Direct comparison against
the primary display gives exactly these rows, encoded without relabeling in
`matrix.py`:

```
0,    0,   -u3,                 0,                    y12,                       0
0,    0,     0,                v6,                      0,                     -x9
a,    0, -a-e1-u2-u3,         -e1,                      0,                       0
-b4,-b4,   -e4, -b4-e4-v5-v6,                        -b4,                       0
-b10,-b10,   0,              -b10, -b10-f10-y11-y12,                   -f10
0,    c,     0,                 0,                    -f7,       -c-f7-x8-x9
```

Equation (4.1) defines the monic characteristic coefficients c0=1,c1,...,c6;
the matrix following (4.3) defines the Hurwitz determinants. Equation (4.4)
provides the closed Delta4 formula used only for multiplier discovery.
`verify_certificate.py` instead computes characteristic coefficients by all
principal minors and Delta4 by its defining 4x4 determinant.

Credited prior results: Proposition 4.2 proves c1,...,c4>0 and Delta1,...,Delta3>0
for all sixteen positive parameters, and c6>0 implies c5>0. Its statement does
not include determinant-conditioned positivity of Delta4. Proposition 5.6
proves Delta4>0 under mass-action constraints (MA1),(MA2); its proof reports
the 1,712,913-term expansion, 266 negative terms, and their common factor
a*c*e4*f10*v6*y12. The source's Appendix implementation overview lists a
mass-action-corrected numerator of Delta4. The currently visible primary
`dualfutile_positivity_checks.m` likewise describes that mass-action check.

The present lemma proves Delta4>0 under c6>=0 without mass-action constraints or
parameter identifications. It extends the inspected Delta4 result's parameter
scope; it does not re-claim the mass-action no-Hopf theorem or the lab's prior
eight-parameter equivariant classification. Searches for the exact arXiv ID,
conjecture, and fourth-Hurwitz/futile-cycle keywords found no directly
overlapping extension; this is a bounded comparison, not an absolute priority
guarantee. The source code comparison is to the linked visible primary files; it does
not assert an independently revalidated latest commit hash.

The remaining full-conjecture step is Delta5>0 on c6>0, because Delta6=c6*Delta5.
No conclusion about existence or stability of periodic orbits follows from this
partial algebraic result alone.
