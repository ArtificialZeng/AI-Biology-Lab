# Attracting oscillations for all sufficiently large finite kinase totals

This is a mathematical consequence of the validated limiting orbit and standard local singular-perturbation theory. It is not an explicit finite-K threshold, a global-attraction theorem, the full original Conjecture 6.2, or a claim of priority. The source model's oscillations and directly related 2023/2025 research are prior work requiring comparison.

## Precisely stated result

Fix the published Table 1 rates, total phosphatase P=8, and total substrate S=40. There is a finite K0 such that for every finite total kinase K>K0 the six-dimensional reduced mass-action system has a strictly positive periodic orbit that is orbitally asymptotically stable within its compatibility class with the fixed conserved totals (K,8,40). These orbits approach the certified limiting orbit as K tends to infinity. The limiting orbit has a return period in [4.6920959801419082, 4.6920960227201061] in the model's time units. No biological time calibration is asserted.

The numerical input to this theorem is the interval return-map certificate in this folder; the argument below explains how its limited, finite computation yields an unbounded parameter conclusion. The earlier equilibrium-instability bound 4660166/209569 is not a bound for K0.

## Exact singular coordinates

Let epsilon=1/K, x0=epsilon*y and q=(u,v,w,a,b). Set p=8-a-b and s(q)=40-u-v-w-a-b. The original equations become

    epsilon*y' = -y + epsilon*(u+v)*y + u + 100*b,
    u' = y*(1-epsilon*(u+v))-2*u,
    v' = u-v,
    w' = v-100*w*p+a,
    a' = 100*w*p-(19/10)*a,
    b' = 3*(s(q)-epsilon*y)*p-101*b.

All coefficients and coordinate changes are exact. At epsilon=0 the critical graph is y=h0(q)=u+100*b. Substitution into q' gives exactly the five-dimensional vector field in verify_cycle.cpp. The layer vector field has normal derivative -1 at every point of this graph.

## Local persistence argument

1. The interval proof supplies a nonconstant periodic orbit Gamma0 of the limiting slow vector field, whose four transverse return multipliers have modulus strictly below one. The phase direction is removed by the section u=13/4. The initial point is strictly physical. The invariant-domain argument in candidate_lemma.md then proves the entire orbit is physical and strictly positive.
2. Gamma0 is compact. Thus u,v,w,a,b,p and s have strictly positive minima on it, as does h0=u+100*b. Choose a sufficiently small smooth tubular slow neighborhood with compact closure in this strict interior. The critical graph over this neighborhood is uniformly normally attracting because its layer normal derivative is exactly -1. Only a local neighborhood of the orbit is needed; the critical graph need not be compact globally. One can extend the vector field away from this neighborhood to apply a compact local Fenichel theorem without changing the orbit neighborhood.
3. Local Fenichel persistence yields a locally invariant attracting graph y=h_epsilon(q), C1-close to h0, for all sufficiently small positive epsilon. The reduced slow vector field on it is C1-close to the certified limiting vector field. Hyperbolic periodic-orbit persistence (equivalently, the implicit-function argument for a phase-fixed return map) therefore supplies an attracting nearby periodic orbit for every such epsilon. The assertion is uniform on an interval (0,epsilon0), rather than a numerical sequence of K values.
4. Full physical coordinates also remain positive. On these nearby cycles, h_epsilon>0 and s(q)-epsilon*h_epsilon(q)>0 by the positive limiting minima and compactness. The free kinase epsilon^-1-u-v is positive after decreasing epsilon0 further, since u+v stays bounded. Thus x0=epsilon*h_epsilon>0, and all nine reconstructed concentrations satisfy the original conservation laws.
5. Tangential orbital attraction persists from the four slow transverse multipliers. In a coordinate eta=y-h_epsilon(q), the additional normal variational equation has coefficient -1/epsilon+O(1), uniformly along the compact nearby cycles. Its one-period multiplier is strictly contracting for epsilon sufficiently small. Hence the finite-K orbit is orbitally asymptotically stable within the six-dimensional compatibility class. This is not attraction in the ambient nine-dimensional concentration space when conserved totals are allowed to change.

Take K0=1/epsilon0. This proves the quantified assertion, without calculating epsilon0 or K0. The nearby periods converge to the limiting return period by return-map continuity; no nonzero first correction or explicit error constant has been certified.

## Credit and boundaries

The exact reduction and source model are documented in ../theory.md and the 2019 Conradi-Mincheva-Shiu paper (arXiv:1809.02886v2). Fenichel persistence and hyperbolic-cycle continuation are established methods. CAPD provides the interval-integration/return-map machinery. Related Grasse 2023 and Grasse-Conradi 2025 conference materials may already cover the qualitative robustness result; an inspectable model-specific certificate is the concrete deliverable here. No external human peer review, new cellular experiment, exact kinase-independent output waveform, or resolved full source conjecture is claimed.
