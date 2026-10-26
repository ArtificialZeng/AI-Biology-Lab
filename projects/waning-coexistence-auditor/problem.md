# Exact candidate scope

Use the actual eight-state mass-action ODE (6), with ordered coordinates (s,i1,i21,i2,i12,r1,r2,r12), demographic inflow Lambda>0, mortality mu>0, transmission beta_i>0, recovery gamma_i>0, secondary infectivity eta_i>0, cross susceptibility sigma_i>0, and arbitrary waning rates theta1,theta2,theta12>=0. Put v_i=gamma_i+mu, d_i=mu+theta_i, d12=mu+theta12, p=v1/beta1, q=v2/beta2, s0=Lambda/mu. The force of infection is J1=i1+eta1*i21, J2=i2+eta2*i12. Define R2(E1)=[p+eta2*sigma2*gamma1*(s0-p)/(v1+theta1)]/q.

Candidate theorem:

1. If p<q and R2(E1)>1, a strictly positive coexistence equilibrium exists, for every finite parameter vector above, without an additional feasible-root hypothesis. The symmetric statement holds when q<p.
2. If p<q<=s0, existence of a strictly positive coexistence equilibrium is equivalent to R2(E1)>1. In particular this gives a complete existence criterion in the natural region where both one-strain endemic boundary equilibria exist (q<s0). It permits all positive eta_i,sigma_i and all nonnegative waning rates.
3. The extra resident-existence condition cannot simply be omitted: the rational source-model example in the certificate has p=1, q=2, s0=23/18, a strictly positive nondegenerate coexistence equilibrium, and R2(E1)=61/72<1. E1 exists but a positive E2 does not. This qualifies the actual-model converse; it does not disprove the general A1–A5 Conjecture1.

This is a result-note candidate. The general original Conjecture1 and full Question2 root-count classification remain unresolved. No uniqueness, local/global stability, periodic orbit or sharp multiplicity claim is made.

Prior delta: Gavish arXiv2405.09813 Theorem2.2 treats a sufficiently large transmission asymptotic regime; Gavish arXiv2412.07536 Theorem2 treats the no-ADE eta1=eta2=1, permanent-immunity case and already gives existence/uniqueness there. We do not claim either prior subcase as new. The proposed delta is the all-finite-parameter sufficiency with arbitrary ADE and waning, and the corresponding both-resident existence equivalence. Literature novelty remains bounded and must be independently reviewed.
