# Analytic period-link subaudit

Scope: analytic link only. The numerical sign certificate, original-model reduction, and base-cycle interval proof are inputs to this subaudit. No novelty or source-conjecture resolution is inferred.

## Verdict

The local analytic link is valid. For the displayed polynomial slow-fast system, parameter-dependent Fenichel persistence can be applied with any prescribed finite differentiability order, including r=3. Consequently a jointly C3 slow graph and a C2 (indeed C3) physical-time period branch exist through epsilon=0 on a fixed neighborhood of the base cycle. The artificial family F0+eF1 has the same value and first parameter derivative there and therefore the same T'(0). A rigorously positive T'(0) implies T(K)=T_infinity+T1/K+O(K^-2), approach from above, and eventual dT/dK<0 on the selected local branch. No numerical threshold or remainder constant follows.

The current proof_link.md should cite a parameter-smooth version of Fenichel, explicitly choose r>=3 if its vector-field O(epsilon^2) remainder is meant in C1, and state the short argument below instead of retaining its conditional wording as the final proof. This is a proof-exposition requirement, not an obstruction to the claimed theorem.

## Joint smooth graph through zero

Write F(q,y,e) for the exact physical slow field and g(q,y,e)=-y+e A(q)y+h0(q). The fast-time system is dq/dsigma=e F and dy/dsigma=g. After eta=y-h0(q), it is

    dq/dsigma=e F(q,eta+h0(q),e),
    deta/dsigma=-eta+e[A(q)(eta+h0(q))-Dh0 F(q,eta+h0(q),e)].

At e=0 its critical graph eta=0 consists of equilibria with tangent eigenvalues zero and one normal eigenvalue -1. Take nested slow neighborhoods of the compact positive cycle, with compact closure in the physical interior. Smoothly cut off/extend only the e-dependent terms outside the larger neighborhood onto a compact slow base (for example a sufficiently large torus). At e=0 the compactified graph is the product of that base with eta=0. It is normally attracting and r-normally hyperbolic for every finite r. The vector-field family is C-infinity for signed e near zero. The parameter version of compact NHIM persistence therefore supplies a graph eta=k(q,e) jointly C-r in q and e near zero; its restriction to the smaller original neighborhood gives the desired local slow graph h=h0+k. Taking r=3 is sufficient. Only e>0 is used physically; allowing signed e here is a mathematical extension of a polynomial family, not a physical claim.

This is the parameter smoothness in Fenichel persistence, not the insufficient assertion that a separate graph exists for each fixed positive e. Jones (1995), Geometric singular perturbation theory, Theorem 1 and its graph formulation, explicitly includes finite C-r regularity in e; hypotheses (H1)-(H3) are smoothness on an open parameter interval about zero, compactness, and normal hyperbolicity. Primary chapter DOI: https://doi.org/10.1007/BFb0095239 . The chapter text was inspectable via the Jones text indexed at https://www.scribd.com/document/351009523/Jones-Gspt ; direct publisher/PDF fetches were unavailable in this session. Fenichel original DOI: https://doi.org/10.1016/0022-0396(79)90152-9 .

## First jet and branch dependence

Invariance is e Dh(q,e) F(q,h(q,e),e)=g(q,h(q,e),e). Evaluating at zero gives h(q,0)=h0(q). Differentiating at zero gives

    h_e(q,0)=A(q)h0(q)-Dh0(q)F0(q).

Thus every jointly smooth invariant graph has this same first derivative. Chain rule for F(q,h(q,e),e) gives the asserted F1. With r>=3, Taylor remainder in the reduced field is O(e^2) in C1 on a smaller compact neighborhood. With only r=2, the reduced field is C2 and its first derivative in e exists in C1; scalar period C2 and the scalar O(e^2) conclusion still follow, but one should not silently claim an O(e^2) remainder in C1 from C2 alone.

The phase-fixed return map and return time are C2 jointly in (z,e), by smooth ODE dependence and the transverse-section implicit-function theorem. Invertibility of D_zP-I gives the unique local C2 fixed-point branch z(e), and T(e)=tau(z(e),e). Its derivative is determined entirely by the value and first parameter derivative of the reduced vector field. Hence the artificial polynomial family and the true reduced family have identical T'(0). Scalar Taylor theorem on a compact one-sided parameter interval gives a finite (unspecified) O(e^2) bound.

Slow-manifold nonuniqueness does not change this first jet. More strongly, nearby full periodic orbits lying in a common graph neighborhood lie on every such attracting graph. For eta=y-h(q,e), affine dependence of F and g on y gives exactly

    eta_t=[-1/e+A(q)-Dh(q,e) F_y(q,e)] eta.

The bounded O(1) part is dominated by -1/e for small positive e. A nonzero eta therefore contracts strictly over a full period, contradicting periodicity. The local periodic orbit and its return time are consequently independent of the chosen local graph.

Fenichel is invoked in fast time only to construct the graph. Substitution into q_t=F(q,h(q,e),e) uses the original physical time, so there is no missing factor e in the computed period sensitivity.

## Consequences and limits

For T1>0, continuity of T_e implies T_e(e)>0 for sufficiently small positive e. Therefore dT(K)/dK=-K^-2 T_e(1/K)<0 eventually, and T(K)-T_infinity>0 eventually. The period here is the continued return period selected by the certified phase-fixed map; no minimal-period claim is introduced. No explicit K0, remainder bound, global uniqueness, or statement about arbitrary scalar projections is certified by this argument.
