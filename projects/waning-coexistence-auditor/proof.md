# Full-model population proof candidate

Research candidate awaiting independent mathematical review. The exact scope is in problem.md.

## Model and population identity

Write J1=i1+eta1*i21, J2=i2+eta2*i12, v_i=gamma_i+mu. The source Eq(6) is

s'=Lambda-mu*s-beta1*J1*s-beta2*J2*s+theta1*r1+theta2*r2+theta12*r12,
i1'=beta1*J1*s-v1*i1, i21'=beta1*sigma1*J1*r2-v1*i21,
i2'=beta2*J2*s-v2*i2, i12'=beta2*sigma2*J2*r1-v2*i12,
r1'=gamma1*i1-(beta2*sigma2*J2+mu+theta1)*r1,
r2'=gamma2*i2-(beta1*sigma1*J1+mu+theta2)*r2,
r12'=gamma1*i21+gamma2*i12-(mu+theta12)*r12.

The sum N of all eight compartments obeys N'=Lambda-mu*N, by direct cancellation of all transfer and infection terms. Therefore every equilibrium has N=s0=Lambda/mu. The accompanying checker independently reconstructs the ODE from 24 elementary mass-action reactions and verifies this identity exactly.

## A positive reconstruction interval

Set p=v1/beta1<q=v2/beta2 and e_i=eta_i*sigma_i. At every strictly positive coexistence equilibrium, the rank-one infected blocks give

r2(s)=(p-s)/e1, r1(s)=(q-s)/e2, 0<s<p,
(i1,i21)=c1*(s,sigma1*r2), (i2,i12)=c2*(s,sigma2*r1).

With d1=mu+theta1,d2=mu+theta2,d12=mu+theta12, the recovery balances become

gamma1*s*c1-sigma2*v2*r1*c2=d1*r1,
-sigma1*v1*r2*c1+gamma2*s*c2=d2*r2.

Define

Delta=gamma1*gamma2*s^2-sigma1*sigma2*v1*v2*r1*r2,
A1=r1*(gamma2*s*d1+sigma2*v2*d2*r2),
A2=r2*(gamma1*s*d2+sigma1*v1*d1*r1).

Cramer's rule gives c1=A1/Delta, c2=A2/Delta. Set r12=(gamma1*sigma1*r2*c1+gamma2*sigma2*r1*c2)/d12. Every reconstructed coordinate is positive exactly when Delta>0 on 0<s<p.

There is exactly one s_Delta in (0,p) with Delta(s_Delta)=0. Indeed Delta=0 is equivalent to

(p/s-1)*(q/s-1)=gamma1*gamma2*eta1*eta2/(v1*v2),

whose left side is strictly decreasing from infinity to zero on (0,p). Thus Delta>0 exactly on (s_Delta,p). A1,A2 are continuous and strictly positive at s_Delta, and Delta tends to zero through positive values from the right.

## Ordered invasion suffices, with arbitrary waning

Let N(s) be the sum of the eight reconstructed coordinates. All seven non-susceptible balance equations vanish identically on (s_Delta,p); this follows from the rank-one kernel, the Cramer equations, and the definition of r12. The population identity therefore implies the exact scalar identity

F_G(s)=s'=Lambda-mu*N(s).

This is the decisive step: the reduced susceptible equation is a population constraint. N(s) is continuous on (s_Delta,p], positive on its open interval, and tends to +infinity as s decreases to s_Delta, because c1,c2 tend to +infinity and their compartment weights are strictly positive. Hence F_G(s) tends to -infinity at that endpoint.

At the right endpoint s=p, r2=c2=r12=0, r1=(q-p)/e2, and i1=d1*r1/gamma1. Therefore

N(p)=p+(v1+theta1)*(q-p)/(gamma1*e2).

The source's invasion number is

R2(E1)=[p+e2*gamma1*(s0-p)/(v1+theta1)]/q.

All denominators are positive. Thus R2(E1)>1 is exactly s0>N(p), equivalently F_G(p)>0. By continuity there exists s* in (s_Delta,p) with F_G(s*)=0. Every reconstructed coordinate is strictly positive there; all eight balance equations vanish. Its population equals s0, so it lies on the model's feasible population simplex. This proves sufficiency for every finite positive beta,gamma,eta,sigma,mu,Lambda and nonnegative waning rates. No extra feasible-root assumption or asymptotic approximation remains. The other ordering follows by exchanging strain labels.

For algebraic checking, put W1=s+sigma1*r2*(1+gamma1/d12), W2=s+sigma2*r1*(1+gamma2/d12). With P_G=-Delta*F_G as in the source,

P_G=mu*[(s+r1+r2-s0)*Delta+A1*W1+A2*W2].

Consequently P_G(s_Delta)=mu*(A1*W1+A2*W2)>0, and P_G(p)=-gamma2*p^2*mu*(v1+theta1)*q/e2*(R2(E1)-1)<0. The feasible root is therefore certified on the correct feasible interval, rather than inferred from an irrelevant root below s_Delta.

## Converse when both resident thresholds are met

Assume a strictly positive coexistence equilibrium and p<q<=s0. Since J2>0 and r1>0, its r1 balance gives i1>d1*r1/gamma1. Summing the other positive compartments yields the strict bound

s0=N>s+(v1+theta1)*r1/gamma1.

Let h=(v1+theta1)/(gamma1*e2)>0. Substituting r1=(q-s)/e2 gives

s0>h*q+(1-h)*s.

If h>=1, use s<p to obtain s0>h*q+(1-h)*p. If h<1, then h*q+(1-h)*p<q<=s0, yielding the same strict inequality. This is algebraically equivalent to R2(E1)>1. Together with sufficiency it proves the claimed equivalence on p<q<=s0, including all waning and ADE parameters.

## Why the resident-existence condition is substantive

Take mu=gamma1=gamma2=1, beta1=2,beta2=1, eta1=30,eta2=10, sigma1=sigma2=1/2, Lambda=23/18 and all waning rates zero. The state

(s,i1,i21,i2,i12,r1,r2,r12)=(1/2,1/3,1/90,1/18,1/60,3/10,1/30,1/36)

is strictly positive and satisfies all eight ODEs exactly. Its population is 23/18. The full state Jacobian determinant is 416/243, so this is a nondegenerate equilibrium. Here p=1,q=2 and s0=23/18. E1 exists with r1(E1)=5/36, but R2(E1)=61/72<1. The lower-fitness resident E2 is not positive because its DFE threshold is 23/36<1. Thus unconditional actual-model coexistence=>ordered-invasion is false if that resident condition is omitted. This is not advertised as a disproof under all abstract A1–A5 assumptions, since the positive E2 required in that setting is unavailable.

## Scope and prior knowledge

The general Conjecture1 for arbitrary rank-one models and Question2's one/many/multiple root parameter classification remain unresolved. We prove no uniqueness or stability. Permanent immunity without ADE is covered by earlier existence/uniqueness results; no priority is claimed there. The source-specific new target is the complete finite-parameter model with arbitrary ADE and waning. The proof uses conservation of transfers and demographic population balance, a classical identity; novelty is the resulting full-model conclusion and closure of the identified feasibility gap, not the identity itself.
