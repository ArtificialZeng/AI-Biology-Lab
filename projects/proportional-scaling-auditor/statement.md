# Exact scope

Consider the proportional-recruitment/standard-incidence variant of the uniform-susceptibility class in Earn–McCluskey (2025): recruitment is nu*N, infection force is beta^T*x/N, susceptible mortality is muS, and postinfection transitions form the source PIT matrix M. All postinfection nodes must be reachable from infection entry; transfer rates may be zero. Every compartment mortality is positive. Let N=S+sum(x), s=S/N and y=x/N.

With nu=34/5, muS=5, M=diag(1,17), P=(1/10,9/10), beta=(24,1), the normalized system has a positive saddle at fractions (4/5,1/10,1/10), a disease-free sink, and another positive sink defined by the isolated quadratic root in the proof. The two sinks have distinct local basins containing positive infected fractions. Consequently this variant does not have a single globally asymptotically stable equilibrium in fractions.

This is a counterexample to the unqualified Discussion section titled “Proportional scaling.” It does not refute the main theorem for constant recruitment. It uses mu1<muS, which is permitted by the source's literal arbitrary-positive-mortality class; it does not decide the subclass mu_j>=muS. Proportions equilibria are balanced-growth states of absolute numbers. No real-world disease claim is made.

Equal mortality in every compartment is one sufficient repair: normalization then produces the original constant-recruitment class with recruitment and every effective mortality equal to nu. Necessity of this repair, a complete parameter classification and oscillations are not asserted.
