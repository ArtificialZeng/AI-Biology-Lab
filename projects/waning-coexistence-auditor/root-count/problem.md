# Candidate: full actual-model root classification

All `Lambda,mu,beta_i,gamma_i,eta_i,sigma_i` are finite and strictly positive; `theta1,theta2,theta12` are finite and nonnegative. Use the actual eight-compartment source Eq6. Write `v_i=gamma_i+mu`, `d_i=mu+theta_i`, `d12=mu+theta12`, `e_i=eta_i*sigma_i`, and `s0=Lambda/mu`. Relabel strains so `p=v1/beta1<=q=v2/beta2`. The equality case is included, with its distinct endpoint interpretation.

Let `r1=(q-s)/e2`, `r2=(p-s)/e1`,

`Delta=gamma1*gamma2*s^2-sigma1*sigma2*v1*v2*r1*r2`,
`A1=r1*(gamma2*s*d1+sigma2*v2*d2*r2)`,
`A2=r2*(gamma1*s*d2+sigma1*v1*d1*r1)`,
`W1=s+sigma1*r2*(1+gamma1/d12)`,
`W2=s+sigma2*r1*(1+gamma2/d12)`.

Set `kappa=gamma1*gamma2*eta1*eta2/(v1*v2)` and

`d=2*p*q/(p+q+sqrt((p-q)^2+4*kappa*p*q))`, `ell=p-d>0`.

This d is the unique zero of Delta in `(0,p)`; it need not be its only positive zero globally. Define the polynomial

`P(s)=mu*((s+r1+r2-s0)*Delta+A1*W1+A2*W2)`.

It is exactly the source's PG, has degree at most three, and `P(d)>0`. Define

`a=P(p)`, `b=3*P(p)-ell*P'(p)`,
`c=3*P(d)+ell*P'(d)`, `z=P(d)>0`,

`Dq=b^2*c^2-4*a*c^3-4*b^3*z-27*a^2*z^2+18*a*b*c*z`.

These are explicit functions of the original parameters and one specified positive square root. Their sign conditions are semialgebraic after clearing positive denominators, or introducing the unique positive R with `R^2=(p-q)^2+4*kappa*p*q`. No unknown equilibrium or critical point appears in the tests.

Candidate theorem: the distinct feasible roots and multiplicities are exactly:

| Parameter condition | Feasible roots |
| --- | --- |
| `a<0` | one simple root |
| `a=0` and `b<0` | one simple root |
| `a=0` and `b>=0` | none |
| `a>0`, `(b<0 or c<0)`, and `Dq>0` | two simple roots |
| `a>0`, `(b<0 or c<0)`, and `Dq=0` | one root of multiplicity two |
| all other `a>0` cases | none |

There are no triple feasible roots. The positive reconstruction is one-to-one in s, so distinct feasible roots are exactly the distinct strictly positive coexistence equilibria. A scalar double root is not, by itself, a claim about an eight-dimensional dynamical saddle-node or stability.

Structural result: the reconstructed population satisfies `N''(s)>0` throughout `(d,p)`. Thus ordered invasion gives a unique coexistence equilibrium at every finite parameter vector, strengthening pass1's existence. For `p<q<=s0`, coexistence is unique exactly when ordered invasion holds; otherwise none. This strict-order corollary must not be extended to `p=q=s0`: that critical DFE corner may have a positive equilibrium even though the algebraic invasion expression equals one.

The current proposed outcome is a resolution candidate for actual Question2, subject to fresh mathematical review; broad Conjecture1 remains unresolved. No uniqueness of general rank-one systems, stability, persistence, or new ADE mechanism is claimed. Prior no-ADE/permanent and large-transmission uniqueness is credited; the delta is all-finite-parameter full waning classification.
