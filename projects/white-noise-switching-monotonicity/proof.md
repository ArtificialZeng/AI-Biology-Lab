# An explicit finite switching-rate decrease

This is the finite example in a technical clarification, within the fixed fast-switching-optimal target branch. A fresh-context AI review found the frozen mathematics valid and required attribution/framing revisions, which are incorporated here. This local review does not establish human peer review, priority or formal runner certification. It is not a claim about the source's ambiguous optimized-target wording or finite-population survival. All populations below are concentrations in the original source Itô model.

## Source model and stationary identities

Let R1 and R2 be positive concentrations in the proliferative and dormant phenotypes, N=R1+R2, and let W1,W2 be independent standard Wiener processes. The original two-population Itô model is

    dR1=[mu1*R1−omega*q*R1+omega*(1−q)*R2]dt+sigma1*R1*dW1,
    dR2=[mu2*R2+omega*q*R1−omega*(1−q)*R2]dt
        +sigma2*R2*[rho*dW1+sqrt(1−rho²)*dW2].

Here mu_i are fixed net-growth means, sigma_i>=0, rho is in[-1,1], omega>0, V>0, and the target is strictly between0 and1. The target dormant proportion q is fixed while the switching speed omega changes. Phenotypes are inherited on division; switching passively redraws the phenotype from the fixed target distribution. This concentration model has environmental noise and no integer demographic extinction state. The independent symbolic derivation in verify_model.py applies their full covariance generator to X=R2/(R1+R2) and log(R1+R2). Write V=sigma1²−2rho*sigma1*sigma2+sigma2² and q=(mu2−mu1+sigma1²−rho*sigma1*sigma2)/V. With the target fixed to this q, it gives

    dX=(q−X)[omega+V*X(1−X)]dt+sqrt(V)*X(1−X)dW,
    h(x)=h(q)−V*(x−q)²/2.

Positive switching points the drift strictly inward at both endpoints. The author solver's density satisfies the zero-current stationary equation exactly. It is positive and normalizable, with exponential barriers at both endpoints; all reciprocal polynomial moments needed below are finite. This is the invariant law of the original ratio diffusion on (0,1). These facts can be justified by the usual one-dimensional speed-density construction, or directly by the normalizable zero-current density and inaccessible endpoints of the positive number SDE. The initial concentrations may be any strictly positive pair; long-run fitness is the stationary expectation of h.

Let lambda=omega/V, R=X/(1−X), nu=(1−2q)(1+2lambda), A=2lambda(1−q), and B=2lambda*q. Changing variables in the stationary density gives

    p_R(r) = C*r^(−nu−1)*exp(−A*r−B/r), r>0.

The endpoint exponential factors justify integrating the generator against log(1−X), by cutoff and dominated convergence. Alternatively, use log R1=log N+log(1−X). The population-number equation gives the stationary identity

    F(omega)=phi1+omega[(1−q)E(R)−q],
    phi1=mu1−sigma1²/2.

Its expectation is finite. Thus no moment of an incorrectly averaged noise coefficient is substituted for the actual Itô fitness.

## Exact source-admissible parameters

Take

    mu1=1, mu2=101/1000, sigma1=1, sigma2=1/10, rho=0.

Then V=101/100, q=1/10, phi1=1/2 and phi2=12/125. Both means and noise amplitudes are positive, the dormant phenotype has smaller mean growth and variance, phi1>phi2, and 0<q<1/4<1/2. The target is exactly the interior fast-switching optimum, not a separately tuned target. Choose

    lambda1=1/1000, lambda2=1/500,
    omega1=101/100000, omega2=101/50000.

All hypotheses of the frozen restricted question are satisfied.

## A rational upper bound for the mean odds at lambda2

At lambda2, nu=502/625, A=9/2500, and B=1/2500. Define

    Nint = integral_0^infinity r^(−nu)*exp(−A*r−B/r) dr,
    Dint = integral_0^infinity r^(−nu−1)*exp(−A*r−B/r) dr.

Then E(R)=Nint/Dint. Both integrals are positive and finite.

For the numerator, substitute v=A*r and drop exp(−AB/v)<=1. Splitting the remaining integral at1 gives

    Nint < A^(nu−1) [1/(1−nu)+exp(−1)].

Here 0<1−nu=123/625<1/5 and A<1. The exact rational fifth-power inequality (25/8)^5>2500/9 gives A^(nu−1)<25/8. Also exp(−1)<3/8, since exp(1)>1+1+1/2+1/6=8/3. Consequently

    Nint < (25/8)(625/123+3/8)=134225/7872 = Nupper.

For the denominator, substitute u=B/r:

    Dint=B^(−nu) integral_0^infinity u^(nu−1)*exp(−u−AB/u) du.

Restrict the integral to delta<=u<=1, with delta=1/100. Since nu<1, u^(nu−1)>=1 there. Since AB/u<=AB/delta=9/62500 and exp(−z)>=1−z for z>=0,

    integral_delta^1 u^(nu−1)*exp(−u−AB/u) du
      > (1−9/62500)[exp(−delta)−exp(−1)]
      > (1−9/62500)(1−1/100−3/8).

Finally nu>4/5, and the integer inequality2500^4>500^5 implies B^(−nu)>500. Therefore

    Dint >500(1−9/62500)(1−1/100−3/8)
         =7686393/25000 = Dlower.

The exact rational inequality Dlower>18*Nupper proves

    E_lambda2(R) < Nupper/Dlower < 1/18.

No special-function evaluation, floating eigenvalue or numerical quadrature is used for this bound.

## Strict finite comparison

At lambda1 the mean odds are strictly positive, so

    F(omega1)>phi1−omega1*q=499899/1000000.

At lambda2 the bound above gives

    F(omega2)<phi1+omega2[(9/10)(1/18)−1/10]
             =499899/1000000.

Thus F(omega2)<F(omega1) for the displayed ordered, positive rational switching rates. The proposed fixed-target nondecreasing assertion is false. verify_bounds.py checks all rational parameter identities, power comparisons and separating bounds, with explicit checks active under normal and optimized Python. It does not formally verify the analytic lemmas; those are given above for independent reading.

No numerical illustration is needed for the separating bound above. This proof concerns concentration fitness only. It does not establish extinction probabilities, new biological mechanisms, an optimized-over-q conclusion, historical priority or the complete resolution of every interpretation of the source prose conjecture.


The general slow-rate mechanism and fast-rate limit are proved in [asymptotics.md](asymptotics.md). See [sources.md](sources.md) for model credit and the precise relation to the source conjecture, and [review-boundary.md](review-boundary.md) for review limitations.

The subsidiary fixed-target assertion is disproved. The complete published prose conjecture remains ambiguous/unresolved; the optimized-target envelope is not addressed. Equivalent stochastic dispersal theory and prior nonmonotonicity are credited in sources.md.
