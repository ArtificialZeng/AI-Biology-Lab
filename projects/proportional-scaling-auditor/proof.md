# Proposed scoped correction: proportional recruitment can destroy global stability

This is a research dossier, not a manuscript. It concerns only the unqualified proportional-scaling extension in Earn–McCluskey (PNAS2025), Discussion §Proportional scaling, current author PDF p5. The main constant-recruitment theorem is untouched.

## Exact proportional variant and quotient rule

Let N=S+sum(x), and replace the source recruitment ν by νN and its force of infection by βᵀx/N. Put s=S/N and y=x/N. With δ=μ−μS*1, direct summation and the quotient rule give

N'/N=g(y)=ν−μS−δᵀy,

y'=s(βᵀy)P−My−g(y)y,  s=1−1ᵀy.

Thus arbitrary unequal mortalities produce a state-dependent growth/dilution term. This identity alone would not refute a convergence statement. The following exact witness does.

## Rational witness within the literal source hypotheses

Take two postinfection classes, no transfers between them, and

ν=34/5, μS=5, μ=(1,17), P=(1/10,9/10)ᵀ, β=(24,1)ᵀ, M=diag(1,17).

All vital parameters, infection probabilities, and transmission coefficients are strictly positive. Both postinfection nodes are reachable directly from S because both entries of P are positive. The source explicitly permits zero transfer coefficients and a reducible PIT matrix. There is no return to S. In particular, every source graph/rate condition is satisfied. The witness has μ1<μS; this is allowed by the stated mathematical class but is not presented as a calibrated disease mechanism.

For u=y1,v=y2 and s=1−u−v, the reduced polynomial system is

10u'=−4u+v−64u²+95uv−v²,

10v'=216u−179v−216u²−265uv+111v².

The physical triangle u>=0,v>=0,u+v<=1 is forward invariant: coordinate facets have nonnegative inward flows, and at s=0 the normalized susceptible derivative is ν>0. Since the original vector field has positive finite vital parameters and bounded incidence per capita, positive population numbers exist globally; a balanced-growth solution is legitimate.

At the strictly positive proportion vector (s,u,v)=(4/5,1/10,1/10), both reduced derivatives vanish, while g=1. Its reduced Jacobian is

J = [[−73,103],[1463,−1833]]/100,

with det(J)=−211/125<0 and trace(J)=−953/50. Hence it is a hyperbolic saddle, with one positive and one negative real eigenvalue. The corresponding absolute numbers satisfy (S,x1,x2)=N0*exp(t)*(4/5,1/10,1/10), not a constant-number equilibrium. All stability statements here concern proportions.

The saddle already disproves a universal GAS-equilibrium conclusion for every infected initial fraction: its constant-in-proportions trajectory cannot converge to any other equilibrium, and it is not stable itself. There are also two distinct local attractors, as the following exact calculation shows.

## Bistability in proportions

At the disease-free fraction vector (1,0,0), the reduced Jacobian is

J0=[[-2/5,1/10],[108/5,−179/10]],

with trace=−183/10 and det=5. Both eigenvalues have negative real part. This is strict local asymptotic stability in the physical triangle; its neighborhood contains positive infected fractions. This direct test does not reuse or redefine the source constant-recruitment R0.

Let alpha be the unique root in (13/100,7/50) of

q(v)=41464v²−72863v+9180.

Existence follows because q(13/100)>0>q(7/50), and uniqueness in that interval because q'<0 there. Define

v*=alpha, u*=125/754−37alpha/377, s*=1−u*−alpha.

All three fractions are positive throughout the isolating interval. Substitution in each reduced derivative gives a polynomial multiple of q, hence both derivatives vanish exactly. Reducing trace and determinant of the reduced Jacobian modulo q after that substitution gives

trace = 13405alpha/377−187107/7540,

det = 66266alpha/377−16165/754.

The isolating interval implies the strict rational bounds

trace < −149573/7540 <0,

det > 13302/9425 >0.

Therefore this second positive proportion equilibrium is a hyperbolic sink. It and the disease-free equilibrium have different local basins, both containing positive infected fractions. Thus the proportional variant can have two coexisting stable outcomes in fractions even though all susceptibles are indistinguishable and no infected individual returns to S.

No claim is made about bifurcation direction, periodic orbits, a full parameter classification, or biological prevalence in a real population.

## A sufficient repair already implied by the source theorem

If every mortality equals a common μ>0, then g=ν−μ is constant, and the proportions satisfy

s'=ν−νs−sβᵀy,

y'=s(βᵀy)P−[M+(ν−μ)I]y.

The new PIT matrix has off-diagonal entries −mij and column sums ν>0. This is exactly the source constant-recruitment class with susceptible mortality ν, so its main theorem applies to proportions. This is a sufficient repaired assumption, not a claimed necessary condition. The witness does not decide whether the weaker inequalities μj>=μS suffice.

## Verification and logical scope

`check_candidate.py` independently reconstructs the numbers equations, their total derivative and quotient-rule reduction; verifies every rational witness condition; and checks the second equilibrium using exact polynomial remainders and rational interval bounds. Normal and optimized Python executions produce byte-identical certificates. Discovery script `search_normalized.py` used floating point only to find a lead; no claim depends on its numerical eigenvalues.

The mathematical statement proposed for fresh review is: the source's unqualified §Proportional scaling extension fails for its literal arbitrary-positive-mortality domain. The source main theorem remains valid as far as this work addresses; no independent audit of that main proof is claimed. The paper-specific correction and its priority remain subject to fresh source/scope review. Classical epidemic bistability is not claimed as a new phenomenon.
