# Slow-switching obstruction on the fixed optimal-target branch

This is the companion calculation for the precise fixed-target restricted question, not a claimed resolution of every interpretation of Browning–Hamis2025's prose conjecture. Let phi_i=mu_i-sigma_i^2/2, V=sigma1^2-2rho*sigma1*sigma2+sigma2^2>0, and q=(mu2-mu1+sigma1^2-rho*sigma1*sigma2)/V. Assume phi1>phi2 and0<q<1/4. The actual concentration Itô equations, fixed passive redistribution and actual stationary instantaneous log-growth convention are those in the model stated in [proof.md](proof.md).

## Claim

For the target q held equal to this fast-switching optimum, let F(omega)=E_stationary[h(X)], where X=R2/(R1+R2) and

h(x)=mu1+(mu2-mu1)x-[sigma1^2(1-x)^2+2rho*sigma1*sigma2*x(1-x)+sigma2^2*x^2]/2.

Then, as omega decreases to0,

F(omega)=phi1-q*omega+o(omega).

In particular F is not nondecreasing on positive omega. More explicitly F(2omega)-F(omega)=-q*omega+o(omega)<0 for all sufficiently small positive omega. Its fast-switching limit is h(q)=phi1+V*q^2/2>phi1, so the response is nonmonotone overall. This concerns long-run logarithmic concentration growth, not finite-integer-population extinction.

## 1. Original-number Itô reduction and stationary law

The original number equations give

dX=[omega(q-X)+V*X(1-X)(q-X)]dt+sqrt(V)*X(1-X)dB.

Set lambda=omega/V and rescale time byV. The reduced generator is

L f=(q-x)[lambda+x(1-x)] f' + x^2(1-x)^2 f''/2.

The quadratic fitness satisfies h(x)=h(q)-V*(x-q)^2/2. With Y=log(X/(1-X)), another Itô calculation gives unit diffusion and drift

(q-1/2)+lambda*((2q-1)+q*exp(-Y)-(1-q)*exp(Y)).

This drift is locally smooth, points strongly inward at both infinities, and has a normalizable positive zero-current density on the real line. For example the generator applied to1+Y^2 is bounded above by a constant outside a compact set and strictly negative far out, ensuring no explosion; the invariant law is positive recurrent and unique. Thus the original fraction process stays in(0,1), and stationary expectations are legitimate. The stationary Y density is proportional to

exp[-nu*Y-2lambda*(q*exp(-Y)+(1-q)*exp(Y))],

where nu=(1-2q)*(1+2lambda). In odds R=exp(Y), this becomes

p_R(r)=Z_lambda^(-1)*r^(-nu-1)*exp[-lambda*(A*r+B/r)], r>0,

A=2(1-q), B=2q.

Both endpoints have exponential cutoffs for everylambda>0. Every needed odds/log-fraction moment is finite. This density is independently checked against the zero-current law from the original SDE; it matches the actual author stationary implementation rather than relying on the differently displayed supplement antiderivative.

## 2. Exact fitness–odds identity

The original first-population equation has log drift

phi1-omega*q+omega*(1-q)*R.

Since log R1=log(R1+R2)+log(1-X), the stationary generator identity for psi(x)=log(1-x) gives

F(omega)=phi1+omega*((1-q)*E[R]-q).

This identity needs no unsupported pathwise claim about log fractions divided by time. Indeed, zero current gives

integral L psi*p dx = [x^2(1-x)^2*p*psi'/2] at0,1 =0.

The endpoint exponentials kill the boundary terms, and all terms are integrable. The same identity follows by cutoff approximation if desired.

## 3. Moving-order moment asymptotic, without special-function assumptions

Put nu0=1-2q, so1/2<nu0<1 andnu=nu0(1+2lambda). Define

D_lambda=integral_0^infinity r^(-nu-1)*exp[-lambda*(A*r+B/r)]dr,
N_lambda=integral_0^infinity r^(-nu)*exp[-lambda*(A*r+B/r)]dr.

Then E[R]=N_lambda/D_lambda. In the denominator setr=lambda*u, and in the numerator setr=v/lambda:

D_lambda=lambda^(-nu)*I_D(lambda),
I_D(lambda)=integral u^(-nu-1)*exp[-A*lambda^2*u-B/u]du,
N_lambda=lambda^(nu-1)*I_N(lambda),
I_N(lambda)=integral v^(-nu)*exp[-A*v-B*lambda^2/v]dv.

For smalllambda, nu lies in a fixed compact subinterval of(1/2,1). In I_D, exp(-B/u) controls0 and u^(-nu_min-1) controls infinity. In I_N, v^(-nu_max) controls0 and exp(-A*v) controls infinity. These are integrable dominators independent oflambda. Dominated convergence therefore gives

I_D -> B^(-nu0)*Gamma(nu0),
I_N -> A^(nu0-1)*Gamma(1-nu0).

The varying exponent causes no hidden singularity: lambda^(2*(nu-nu0))=exp(4nu0*lambda*loglambda)->1. Consequently

E[R] ~ C(q)*lambda^(2nu0-1) = C(q)*lambda^(1-4q) ->0,
C(q)=A^(nu0-1)*B^nu0*Gamma(1-nu0)/Gamma(nu0)>0.

The exact identity in §2 yields (F(omega)-phi1)/omega -> -q. This proves the stated slow-rate expansion and the separated pair comparison for omega and2omega. A finite numerical quadrature is not needed for the disproof.

## 4. Fast limit

The fraction density can also be written

p_lambda(x) proportional x^(2q-2)*(1-x)^(-2q)*exp[-2lambda*D_q(x)],
D_q(x)=q/x+(1-q)/(1-x)+(1-2q)*log(x/(1-x)).

Its derivative is D_q'(x)=(x-q)/(x^2*(1-x)^2), so D_q has unique strict global minimum atq and diverges at both endpoints. To justify concentration despite the nonintegrable bare prefactor, fix anylambda0>0 and use its already normalizable stationary law as reference. Exponential tilting byexp[-2(lambda-lambda0)D_q] makes the probability outside any neighborhood ofq decay relative to a smaller neighborhood with strictly lowerD_q. Hence X converges in probability toq aslambda->infinity. Bounded continuity ofh givesF->h(q)=phi1+Vq^2/2.

## 5. Fully source-admissible parameter example

mu1=1, mu2=101/1000, sigma1=1, sigma2=1/10, rho=0.

All means and noise amplitudes are positive, correlation is interior, V=101/100 andq=1/10. phi1=1/2>phi2=12/125=.096. The theorem applies. The high-rate fitness limit is10101/20000=.50505, while sufficiently small positiveomega givesF<1/2 and increasingomega to2omega lowersF further. The separate [finite counterexample proof](proof.md) supplies a certified rational pair; the all-smallomega asymptotic above explains the general mechanism.

## Scope

This is a counterexample to nondecreasing fitness on the source's natural fixed-q_inf positive-hedging branch. It does not prove what happens to the reoptimized envelope sup_q F(q,omega), does not address colored environments or survival probability, and does not imply faster switching is always harmful. The original paper's fast-switching optimum and broader bet-hedging conclusions are retained. Published2025 text contains the motivating conjecture, but the unresolved target-quantifier ambiguity must remain explicit. Prior2021 model/optimal-control results must be credited and screened separately before priority claims. No public author confirmation or empirical validation is asserted.

This companion calculation is included coherently with the finite example. Equivalent stochastic dispersal models and general nonmonotone responses are prior work; no general asymptotic priority or standalone journal novelty is established. See [sources.md](sources.md).
