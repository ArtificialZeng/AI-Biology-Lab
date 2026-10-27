# Ordered mortality does not restore universal fraction stability

This exact addendum received a separate fresh-context AI mathematical/source review on 28 October 2026 local time (27 October UTC), which accepted the scoped witnesses below. This is not external human peer review or formal runner certification. It addresses a newly formulated natural mortality subclass, not a source-declared open conjecture or an additional asserted author error. The common-mortality repair remains sufficient. The main constant-recruitment theorem is not challenged.

## Numbers and fraction fields

Let N=S+x1+x2>0. Use precisely

    S'=nu*N-muS*S-S*(beta1*x1+beta2*x2)/N,
    x'=S*(beta1*x1+beta2*x2)/N*P-M*x.

All newborns enter S with the same per-capita recruitment; neither infected births nor return to S occurs. All postinfection classes remain in N. Off-diagonal entries of M are negative transfer rates, and column sums are mortalities. Put y=x/N, s=1-y1-y2, delta=mu-muS*1 and g=nu-muS-delta^T*y. The quotient rule gives

    F(y)=s*(beta^T*y)*P-M*y-g(y)*y,
    DF=P*(s*beta-(beta^T*y)*1)^T-(M+g*I)+y*delta^T.

The physical triangle is invariant: at yi=0 the infection and transfer inflows are nonnegative; at s=0 the susceptible proportion derivative is nu>0. The total per-capita growth g is bounded on this compact triangle, so positive N persists for every finite time. Fraction equilibria are balanced-growth number solutions, not constant-number equilibria.

## A simple exact model with two sinks and a positive saddle

Take

    nu=1, muS=1/10, mu=(1/10,91/10),
    M=diag(1/10,91/10), P=(1/10,9/10), beta=(135/31,0).

Every mortality is at least muS. Both postinfection nodes are directly reachable because both entries of P are positive. Zero transfers and a zero transmission coefficient are permitted by the source assumptions. With u=y1,v=y2,

    62u'=u*(-35-27u+531v),
    v'=(243/62)*u*(1-u-v)-10v+9v^2.

The disease-free fraction equilibrium has Jacobian

    [[-35/62,0],[243/62,-10]],

with trace -655/62 and determinant175/31, hence is a hyperbolic sink. The two strictly positive equilibria are directly verified by substitution:

* E_saddle: (s,u,v)=(31/54,37/108,1/12), g=3/20. Its Jacobian is [[-37/248,2183/744],[225/248,-2441/248]], trace=-1239/124 and determinant=-37/31. It is a hyperbolic saddle.
* E_sink: (s,u,v)=(31/66,175/396,35/396), g=23/220. Its Jacobian is [[-525/2728,10325/2728],[27/248,-2515/248]], trace=-14095/1364 and determinant525/341. It is a hyperbolic sink.

Thus the simple model has a disease-free sink and a positive sink with distinct local basins, each meeting the strictly positive fraction interior, as well as a positive saddle. This is enough to disprove the proposed universal fraction-GAS conclusion under nonnegative excess mortality. No full equilibrium-count or global basin classification is claimed. Both positive equilibria have positive population growth; this is not a collapse-only example. Positive growth was also present in the old public witness and is not itself a new contribution.

## One final fully strict saddle witness

The following single rational perturbation avoids mortality equality, zero infectivity and zero transfer objections while retaining nu=1 and the same positive equilibrium E_saddle:

    muS=1/10, mu=(11/100,911/100),
    M=[[111/1000,-1/1000],[-1/1000,9111/1000]],
    P=(256241/2504630,2248389/2504630),
    beta=(5008981/1147000,1/1000),
    (s,y1,y2)=(31/54,37/108,1/12).

All postinfection mortalities are strictly greater than muS. Both transfers, both transmission coefficients and both infection-entry probabilities are strictly positive. P sums to1; M has the required sign pattern and mortality column sums. No source assumption is dropped. Exact quotient calculation gives F=0 and g=787/5400>0. The trace and determinant of the reduced fraction Jacobian are

    trace=-387768093698339/38782943235000 <0,
    determinant=-1155754208768357/969573580875000 <0.

The negative determinant is an exact two-dimensional instability certificate: one real eigenvalue is positive and the other negative. Consequently the positive fraction equilibrium is a saddle. Its constant fraction trajectory cannot approach another equilibrium, while this equilibrium is not stable, so the universal GAS conclusion fails even within these strict hypotheses.

For this strict model, direct disease-free linearization gives trace=-1898627083801/179550663125<0 and determinant=817187986178/144798921875>0, a sink. We do NOT claim a second positive sink has been exactly isolated for this strict model. The two-sink claim belongs only to the simple model above. No periodic-orbit or oscillation claim follows from either saddle.

## Reproducibility and contribution limits

`verify_natural.py` uses `fractions.Fraction` and rational dual numbers to reconstruct the original population-number equations and differentiate the quotient. It explicitly checks the source hypotheses, vanishing fraction residuals, balanced-growth number derivatives, exact trace/determinant signs, and agreement with the displayed Jacobian formula. Every gate uses an explicit exception rather than an assertion. Ordinary `verify` performs no writes and compares the certificate with `evidence/natural-mortality-certificate.json`; `archive_tools.py` controls explicit external regeneration. Normal and optimized execution must agree. See the [review summary](review-summary.md) for the review boundary.

The ambient SAIS/two-infective family and the disease-free sink/endemic sink/saddle pattern were already studied by Razvan–Yasaman (2012). This addendum supplies concrete exact witnesses under the present joint source restrictions; it does not establish historical priority for those restrictions or discover a new biological mechanism. The broader prior existence results are not silently imported as a proof at zero-rate boundaries. See [sources.md](sources.md) for the equation-level comparison.

The increment is a coherent strengthening of the existing paper-specific correction: nonnegative, or even strictly positive, infection-associated excess mortality does not by itself repair the proportional-recruitment extension. The common-mortality sufficient repair and the source main constant-recruitment theorem are untouched. There is no claim of oscillations, a complete parameter classification, strict-model bistability, empirical disease validation, public author endorsement, external human review or novelty certification.
