# Exact N=3 ancestral-current classification

## Statement and scope

Fix exactly three Fleming–Viot particles on sites {0,1,2}. Natural motion from any site to each other site has rate1. A particle at site i is killed at arbitrary rate k_i>0 and copies a uniformly chosen other particle, including its ancestral label and its site. Copying at the same site is retained. Use the current orientation and definitions frozen in problem.md:

    Delta_3=F_02-F_20.

For every (k0,k1,k2) in (0,infinity)^3,

    sign Delta_3 = sign[-(k0-k1)(k0-k2)(k1-k2)].

In particular, the complete zero locus is precisely the three planes where at least two killing rates coincide. If k0<k1<k2 then Delta_3>0. The theorem addresses this entire positive three-parameter domain, not a sampled grid or a one-parameter slice. It refines the source's generic-nonzero theorem in the minimal three-site/three-particle case; it does not solve the source's infinite-space nonexplosion question or classify all population sizes.

## Physical probabilities from finite polynomial identities

Let a=k0,b=k1,c=k2. The ten population states and eighteen founder states are enumerated explicitly in evidence/state-ordering.json. The row generator L and founder operator T are the literal source Eq.(3.2) and Eq.(A.4), written in model.py and problem.md. The founder anchor is h(0,(2,0,0))=1/3.

The files population-solution.json and founder-solution.json contain polynomial numerators p_i and h_j and common denominators d_p and d_h. Read-only exact polynomial arithmetic verifies all of the following identities:

    L^T p = 0,          sum p_i = d_p,
    T h = 0,           3h(0,(2,0,0)) = d_h,
    sum_(i:c_i>0) c_i h(i,c-e_i) = d_h    for each of the ten population states.

The 105 coefficients of d_p are strictly positive; all 680 coefficients of d_h are strictly negative. Thus neither denominator vanishes anywhere in the positive domain. Their sign convention is merely a common rational normalization of the archived vector, not a sign claim about the probabilities.

The population chain is finite and irreducible: natural motion connects its count states. Its normalized stationary probability is therefore unique, so p/d_p is pi. The normalized founder system has the unique solution consisting of fixation probabilities (source Proposition3.2). For completeness, uniqueness follows by summing a candidate founder function over the particles of a chosen colour: Eq.(3.7) makes this sum harmonic for the two-colour process. Its boundary values are0 and1, since the all-colour founder sum is a constant population harmonic function fixed at1 by the anchor. Every mixed finite state can reach absorption through positive copying rates, so the finite Dirichlet problem is unique. Hence h/d_h is the actual founder probability, not a fitted surrogate.

No numerical eigensolver, stochastic simulation or floating-point matrix inverse is decisive. Polynomial residual identities supply the exact finite certificate.

## Factorization of the exact current

Form the current numerator directly from the archived p and h vectors according to Eq.(1.2), with the common denominator d_p d_h. The polynomial cross-product identity in verify.py proves

    Delta_3 = -(a-b)(a-c)(b-c) S K Q / (D7 D8),
    S = a+b+c+9,
    K = abc+2(ab+ac+bc)+3(a+b+c).

The explicit expanded numerator, expanded denominator and their factorizations are in evidence/symbolic-current.json. D7 and D8 are the displayed degree7 and degree8 denominator factors; they contain respectively105 and153 monomials, all with strictly positive coefficients. Hence D7D8>0 for all positive rates. Also S>0 and K>0 there.

The residual checker verifies this displayed factorization as a polynomial identity, not by interpolation. It separately verifies that Q is invariant under the swaps a↔b and b↔c; thus Q is symmetric. It remains to show Q>0 away from the full diagonal.

## Positive-coefficient certificate for Q

By symmetry, order the rates and write them as

    a=t, b=t+x, c=t+x+y,    t>0, x>=0, y>=0.

The ordered Q has exactly30 nonzero coefficient polynomials in x,y, listed in the table below. Expanding these coefficient polynomials in t gives exactly94 strictly positive rational (indeed integer) coefficients. This covers the whole polynomial: there is no omitted remainder or truncated series. In particular, the coefficients of x²,xy,y² are all

    C(t)=27(t+2)(t+3)^2(3t²+26t+60)>0.

All the remaining terms are nonnegative. Therefore Q>=C(t)(x²+xy+y²), with equality not needed for the argument. It follows that Q>0 if the ordered triple is not constant, whereas Q(t,t,t)=0 identically. Thus no extra zero of the current can occur when all three killing rates differ.

| x exponent | y exponent | coefficient polynomial in t |
| --- | --- | --- |
| 7 | 0 | `24` |
| 6 | 1 | `84` |
| 6 | 0 | `4*(39*t + 125)` |
| 5 | 2 | `114` |
| 5 | 1 | `12*(39*t + 125)` |
| 5 | 0 | `2*(201*t**2 + 1297*t + 2092)` |
| 4 | 3 | `75` |
| 4 | 2 | `519*t + 1681` |
| 4 | 1 | `5*(201*t**2 + 1297*t + 2092)` |
| 4 | 0 | `3*(171*t**3 + 1670*t**2 + 5432*t + 5848)` |
| 3 | 4 | `24` |
| 3 | 3 | `2*(129*t + 431)` |
| 3 | 2 | `4*(222*t**2 + 1460*t + 2381)` |
| 3 | 1 | `6*(171*t**3 + 1670*t**2 + 5432*t + 5848)` |
| 3 | 0 | `6*(t + 3)*(54*t**3 + 549*t**2 + 1853*t + 2022)` |
| 2 | 5 | `3` |
| 2 | 4 | `3*(18*t + 65)` |
| 2 | 3 | `327*t**2 + 2275*t + 3826` |
| 2 | 2 | `3*(243*t**3 + 2445*t**2 + 8107*t + 8799)` |
| 2 | 1 | `9*(t + 3)*(54*t**3 + 549*t**2 + 1853*t + 2022)` |
| 2 | 0 | `27*(t + 2)*(t + 3)**2*(3*t**2 + 26*t + 60)` |
| 1 | 5 | `3*t + 14` |
| 1 | 4 | `2*(24*t**2 + 191*t + 338)` |
| 1 | 3 | `3*(72*t**3 + 775*t**2 + 2675*t + 2951)` |
| 1 | 2 | `9*(t + 3)*(36*t**3 + 381*t**2 + 1315*t + 1434)` |
| 1 | 1 | `27*(t + 2)*(t + 3)**2*(3*t**2 + 26*t + 60)` |
| 0 | 5 | `(t + 7)*(3*t + 7)` |
| 0 | 4 | `3*(9*t**3 + 104*t**2 + 369*t + 406)` |
| 0 | 3 | `3*(t + 3)*(27*t**3 + 297*t**2 + 1046*t + 1140)` |
| 0 | 2 | `27*(t + 2)*(t + 3)**2*(3*t**2 + 26*t + 60)` |

Combining the positive factors with the alternating Vandermonde factor proves the asserted sign and complete zero locus. If two rates coincide the Vandermonde is zero; all denominators remain nonzero. This handles the equality cases and completes the all-positive-domain classification.

## Relation to the source's perturbation

As a check on the source's exceptional N=3 cubic cancellation, direct substitution of (a,b,c)=(1-e,1,1+e) into the exact rational current gives

    Delta_3 = -24 e^5 (3e²-16)(479e²+38448) /
      [(5e^6+279e^4+4704e²+790272)
       (143e^6+840e^4-137808e²+11757312)].

In particular Delta_3=(89/56010528)e^5+O(e^7)>0 for sufficiently small positive e; the global sign theorem actually gives positivity for every0<e<1. This fifth-order leading term explains why the cubic response in the source vanishes at N=3. It is a consequence of this same exact formula, not a separate inferred discovery.

If one combines this result with the source's explicitly oriented N=4 expansion, whose cubic coefficient is−1/157300 for the same perturbation and Delta=F02-F20, then three and four particles have opposite circulation directions for sufficiently small positive e. This is an asymptotic corollary using the source's N=4 theorem, not a new all-rate N=4 classification.

## Dependencies and remaining limits

The construction uses the source's finite founder representation, a standard finite irreducible-chain stationary/Dirichlet uniqueness argument, and the archived exact polynomial identities/positive-coefficient certificate. A fresh AI referee and independent two-colour implementation are requested; the researcher does not certify their acceptance. No external human confirmation, proof-assistant formalization, trusted-runner gate or historical-first claim is made. The source's supplied supplementary-code claim was not independently located during this algebra pass; bounded source/overlap review is recorded separately. No mathematical gap is currently known in the claimed finite classification, but review can still reject it or identify overlap.
