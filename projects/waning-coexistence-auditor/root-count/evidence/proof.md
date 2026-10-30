# Full root classification for the eight-compartment model

This is the proof candidate for research pass2, not a fresh referee decision. The exact quantified statement and explicit parameter tests are in ../problem.md.

## The source polynomial and the feasible curve

Let `J1=i1+eta1*i21`, `J2=i2+eta2*i12`. The actual source equations, in the order `(s,i1,i21,i2,i12,r1,r2,r12)`, are

```text
s'   =Lambda-mu*s-beta1*J1*s-beta2*J2*s+theta1*r1+theta2*r2+theta12*r12,
i1'  =beta1*J1*s-v1*i1,
i21' =beta1*sigma1*J1*r2-v1*i21,
i2'  =beta2*J2*s-v2*i2,
i12' =beta2*sigma2*J2*r1-v2*i12,
r1'  =gamma1*i1-(beta2*sigma2*J2+d1)*r1,
r2'  =gamma2*i2-(beta1*sigma1*J1+d2)*r2,
r12'=gamma1*i21+gamma2*i12-d12*r12.
```

Their unweighted sum is `N'=Lambda-mu*N`. All parameters, p,q,d,Delta,A1,A2,W1,W2 have the definitions in the statement, with `0<d<p<=q`.

At positive coexistence, infected-block singularity forces `r2=(p-s)/e1`, `r1=(q-s)/e2`, and

`i1=s*A1/Delta`, `i21=sigma1*r2*A1/Delta`,
`i2=s*A2/Delta`, `i12=sigma2*r1*A2/Delta`,
`r12=(gamma1*i21+gamma2*i12)/d12`.

These are the source's Cramer formulas. Their numerators are positive on `(0,p)`. A singular Cramer determinant cannot hide a positive solution: eliminating the second amplitude gives `Delta*c1=A1>0`. Conversely, on `(d,p)` every reconstructed coordinate is positive and all seven non-susceptible balances vanish. Delta has its unique zero in `(0,p)` at the displayed d, and is positive on `(d,p)`. Its possible second positive zero lies outside this interval. The d formula is valid also when the quadratic coefficient is zero (`kappa=1`) and when `p=q`.

By the sum identity, the actual susceptible balance after this substitution is exactly

`FG(s)=Lambda-mu*N(s)`.

Since

`N(s)=s+r1+r2+(A1*W1+A2*W2)/Delta`,

the source definition `PG=-Delta*FG` gives precisely

`P(s)=mu*((s+r1+r2-s0)*Delta+A1*W1+A2*W2)`.

There is no s-dependent change of normalization. The linearity of r1,r2,W1,W2 and quadratic degrees of Delta,A1,A2 prove `deg(P)<=3`. Thus the source's nearby quartic description is not the degree of its own defining expression. The checker verifies the original eight balances and this exact identity. The positivity of the Cramer numerators gives `P(d)=mu*(A1(d)*W1(d)+A2(d)*W2(d))>0`, so P is not identically zero. Its zeros in `(d,p)` are exactly the positive coexistence equilibria, bijectively via their susceptible coordinate. Roots of Delta are not included.

## Strict convexity of population

Put

`K(s)=sigma1*sigma2*v1*v2*r1*r2/(gamma1*gamma2*s^2)`
`     =[v1*v2/(gamma1*gamma2*eta1*eta2)]*(p/s-1)*(q/s-1)`

and `H(s)=1/(1-K(s))`. On `(d,p)`, `0<K<1`, `K'<0`, and `K''>0`. One can differentiate the last expression: with its positive constant C,

`K'=C*((p+q)*s-2*p*q)/s^3<0`,
`K''=2*C*(3*p*q-(p+q)*s)/s^4>0`.

Both signs hold for `s<p<=q`, including equality of p and q. Therefore `H>0`, `H'<0`, and

`H''=K''/(1-K)^2+2*(K')^2/(1-K)^3>0`.

The elementary product fact used below is explicit: if f,g are positive, nonincreasing and convex C2 functions, then `(fg)'<=0` and `(fg)''=f''g+2f'g'+fg''>=0`; if one second derivative is strictly positive, the product is strictly convex.

The functions r1,r2 are positive affine decreasing. Also

`T(s)=r1*r2/s=(p*q/s-p-q+s)/(e1*e2)`

is positive, strictly decreasing and strictly convex on `(0,p)`, because

`T'=(1-p*q/s^2)/(e1*e2)<0`,
`T''=2*p*q/(e1*e2*s^3)>0`.

Rewrite the primary infections as

`i1=H*((d1/gamma1)*r1 + (sigma2*v2*d2/(gamma1*gamma2))*T)`,
`i2=H*((d2/gamma2)*r2 + (sigma1*v1*d1/(gamma1*gamma2))*T)`.

These identities follow directly from the Cramer formulas. Each bracket is positive, strictly decreasing, and strictly convex, since all displayed coefficients are positive. The product fact proves `i1,i2` strictly decreasing and strictly convex. The factors `sigma1*r2/s` and `sigma2*r1/s` are positive, strictly decreasing and strictly convex, so

`i21=(sigma1*r2/s)*i1`, `i12=(sigma2*r1/s)*i2`

have the same properties. Their positive linear combination r12 is also strictly decreasing and strictly convex. Finally `s+r1+r2` is affine. Consequently

`N''(s)>0`, for every `d<s<p`.

There is no approximation, symmetry assumption, or restriction on positive ADE/cross-susceptibility. Zero waning is allowed because d1,d2,d12 remain strictly positive. As `s` decreases to d, positive Cramer amplitudes have a pole and `N(s)->infinity`. At p the curve has a finite smooth extension because `Delta(p)=gamma1*gamma2*p^2>0`.

## Exact classification including endpoint and degree degeneracies

Let `ell=p-d`, and set `s=(d+p*t)/(1+t)` for `t>0`. This is a smooth increasing bijection from the positive t-axis to `(d,p)`, with derivative `ell/(1+t)^2>0`. Multiplication by `(1+t)^3>0` gives

`Q(t)=(1+t)^3*P((d+p*t)/(1+t))=a*t^3+b*t^2+c*t+z`,

where the coefficient formulas in the statement follow by expansion at t=0 and at infinity. Here `z=P(d)>0`. The transformation preserves the algebraic multiplicity of every finite feasible root. If P has degree below three, extra roots introduced at `t=-1` are outside the counting domain. In particular the tests also cover all leading-degree cancellations of P.

If `a=P(p)<0`, then `N(p)<s0`. A strictly convex continuous N with left limit infinity and right endpoint below the level s0 crosses that level exactly once. For completeness, if two ordered roots x<y existed, the increasing secant slopes of a strictly convex function would make the slope from y to p strictly positive, contradicting `N(p)<N(y)`. The root is simple: a stationary root would be the strict minimum and would force `N(p)>N(root)`, again a contradiction. This proves the first row.

If `a=0`, the endpoint level is `N(p)=s0`. If `N'(p)<=0`, strict convexity gives `N'(s)<N'(p)<=0` for every s<p, so N(s)>N(p) and there is no interior root. If `N'(p)>0`, N is below its endpoint level immediately to the left of p, while its left limit is infinity. There is one simple interior root by strict convexity and the same slope argument. At a=0,

`b=-ell*P'(p)=-ell*mu*Delta(p)*N'(p)`.

All omitted factors are positive. This proves exactly the two a=0 rows, including the endpoint-tangent case b=0. A root at p is never counted.

It remains to classify `a>0`. Q is then a genuine real cubic with positive leading and constant coefficients. Its discriminant is the displayed Dq. If Dq<0, it has just one real root, which is negative (Q(0)>0 and its sign at minus infinity is negative), hence no positive root. If Dq>0, all three roots are distinct real numbers. Their product is `-z/a<0`; they are either all negative or two positive and one negative. All-negative roots force b,c>0. Conversely b,c>=0 would make every coefficient of Q nonnegative with positive leading and constant terms, excluding a positive root. Thus exactly two roots are positive precisely when `b<0 or c<0`.

If Dq=0, all roots are real and at least one is repeated. A repeated negative root must be accompanied by another negative root because the product is negative; then b,c>0. A repeated positive root is accompanied by a negative root and necessarily has b<0 or c<0, since otherwise Q has positive coefficients on t>0. A triple positive root is impossible because its product would be positive. A triple negative root has b,c>0. Therefore `b<0 or c<0` gives exactly one positive root, and it is double, precisely in this discriminant-zero case. These arguments exhaust all a>0 cases.

The conditions involve only explicit arithmetic, differentiation of a displayed cubic, and the specified positive square root. Hence they are semialgebraic functions of the original parameters. For a literal polynomial presentation one clears the positive parameter denominators, introduces `R>0` with `R^2=(p-q)^2+4*kappa*p*q`, and uses `d=2*p*q/(p+q+R)`; the specified R is unique. The conditions have no equilibrium-solving or unknown critical-point clause.

A double feasible root also satisfies `P''=mu*Delta*N''>0` at that point, so its multiplicity is exactly two. No triple feasible root is possible. No claim about the stability or dynamical bifurcation type of the eight-dimensional ODE follows merely from this scalar fact.

## Consequences and coincident thresholds

When `p<q`, the source invasion formula and endpoint identity give

`P(p)=-gamma2*p^2*mu*(v1+theta1)*q/e2*(R2(E1)-1)`.

Thus ordered invasion implies a unique positive coexistence equilibrium, for the full finite-parameter model. Together with pass1's strict-order converse, if `p<q<=s0`, there is exactly one positive coexistence equilibrium when ordered invasion holds and none otherwise. All subthreshold multiplicity therefore lies outside this both-threshold region. This is uniqueness of equilibria, not a stability or persistence statement.

Coincident thresholds are not obtained by silently extending that strict-order converse. If `p=q`, then `N(p)=p` and

`N'(p)=1-(v1+theta1)/(gamma1*e2)-(v2+theta2)/(gamma2*e1)`.

The full table still applies. In particular at `s0=p=q`, the algebraic invasion expression equals one; there is one simple positive equilibrium exactly when the displayed endpoint derivative is positive, and none otherwise. There need not be a positive one-strain resident at that endpoint. This is a distinct critical-DFE corner.

## Exact examples and boundaries

The following use `mu=gamma1=gamma2=1`, `beta1=2,beta2=1`, `eta1=20,eta2=10`, `sigma1=sigma2=1`, and `theta2=theta12=1`; then `p=1,q=2` and `Delta=(49*s^2+3*s-2)/50`.

1. Set `theta1=535/216`, `Lambda=2011/1440`. The exact source polynomial is

   `P=(3*s-2)*(4*s-3)*(854*s+251)/21600`.

   Its only feasible roots are `2/3` and `3/4`, both simple. Every waning rate is strictly positive. The third root is negative. This illustrates the two-equilibrium branch, not dynamical bistability.

2. Set `theta1=535/262`, `Lambda=872/655`. Then

   `P=(3*s-2)^2*(595*s+169)/10480`.

   The unique feasible root `s=2/3` has multiplicity two, and `N''(2/3)=229095/102704>0`. Its population level is the strict global minimum. Varying Lambda to `872/655-1/100` gives no feasible roots, while `872/655+1/100` gives two simple feasible roots: the endpoint gap at the double-root parameter is `N(p)-Lambda=191/2620>1/100`. This is a transverse scalar double-root example, with `partial P/partial Lambda=-Delta!=0` at the root. We do not claim an ODE saddle-node without the additional dynamical argument.

3. For a coincident-threshold boundary test, set `beta1=beta2=2`, `eta1=eta2=10`, `sigma1=sigma2=1`, all three waning rates one, and `Lambda=mu=gamma1=gamma2=1`. Then `p=q=s0=1`, `d=1/6`, `N'(p)=2/5`, and

   `P=(s-1)*(4*s+1)*(31*s-11)/250`.

   There is one simple feasible root `11/31`; the endpoint root at s=1 is excluded. This prevents a false extension of the strict-order ordered-invasion equivalence to a critical coincident DFE.

The certificate checks the full reconstructed eight-equation balances at the displayed simple and double roots. Independently selected heldout rational cases include zero and positive waning, both threshold orders, coincident thresholds, a=0 endpoints, and the linear-Delta case kappa=1. Exact rational root isolation/Sturm counts are finite verification evidence supporting implementation; the proof above supplies completeness for all parameters.

## Status and limits

This proposes a complete answer to formal actual-model Question2 at its positive-parameter/nonnegative-waning scope, including repeated roots. The separate abstract Conjecture1 remains unresolved. Scalar multiplicity does not establish ODE stability, persistence or a dynamical saddle-node. Cubic reductions and uniqueness in permanent/no-ADE, one-way-zero, and highly transmissible limits are prior work; this contribution is the full finite-parameter classification and strict-convex population argument. Current literature coverage is bounded and the final JTB revision's full text was inaccessible. Fresh mathematical review of this new pass is still required; no formal gate or human peer review is claimed.
