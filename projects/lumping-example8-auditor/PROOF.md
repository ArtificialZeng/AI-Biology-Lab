# Proof of the Example 8 correction and stability statement

## 1. Derivative correction and exact counterexample

Differentiate the p.38 y5 equation in SOURCE.md with respect to y4 at y1=y4=y5=y6=0. Only the first product contributes, giving −k3*c2. The printed entry −k2*c2 is therefore a transcription error unless k2=k3. The other eight entries agree with direct differentiation of the stated three equations.

At c1=c2=1, k1=k2=km2=7 and k3=k4=km1=km3=km4=1, substitution into the printed matrix gives

    Jp = [-15 -7  0; -7 -3 0; 7 1 -8].

Its characteristic polynomial is (lambda+8)(lambda^2+18 lambda−4). Since 85>81, −9+sqrt(85)>0. Every fast rate and both c1,c2 in this witness are strictly positive. This disproves the footnote's universal extension **when read as referring literally to the printed matrix**. It does not disprove the open-region statement or stability of the actual ODE.

The actual Jacobian at this same witness is

    J = [-15 -7 0; -1 -3 0; 7 1 -8],

with characteristic polynomial (lambda+8)(lambda^2+18 lambda+38). Its roots −8, −9±sqrt(43) are strictly negative.

## 2. Corrected Jacobian is Hurwitz on the full stated domain

Let

    a=k1*c1, b=k2*c2, c=k3*c2, d=k4*c1,
    e=km1, f=km2, g=km3, h=km4.

Here a,b,c,d≥0; e,f,g,h>0. We need not assume a,b,c,d independent, so the argument covers all source parameters. Form the column generator

    Q = [-a-c, e,    g,    0;
          a,  -e-b, 0,    f;
          c,   0,  -g-d,  h;
          0,   b,   d,   -f-h].

Its off-diagonal entries are nonnegative and its column sums vanish. It represents linearized transitions among (P,Ia,Ib,I). In the corresponding directed graph every state can reach P: Ia→P has rate e, Ib→P has rate g, and I→Ia→P has rates f,e. This remains true when c1 or c2 vanishes.

For completeness, the standard spectral facts can be established without any stochastic modeling assumption. Let R=Q^T and choose r larger than every exit rate. Then P0=I+R/r is a row-stochastic matrix with strictly positive diagonal. Uniformization gives

    exp(tR)=exp(-rt) sum_{n≥0} (rt)^n P0^n/n!,

so exp(tR) is uniformly bounded. Gershgorin's theorem puts all eigenvalues of R in Re(lambda)≤0. If Rv=lambda v with Re(lambda)=0, choose an index i with |v_i| maximal and rotate/normalize so v_i=1. Taking real parts in

    lambda v_i=sum_{j≠i} R_ij(v_j-v_i)

forces v_j=1 for every outgoing neighbor j of i. Therefore the right side itself is zero and lambda=0. (If i has no outgoing neighbors, the equation gives this immediately.) Hence no nonzero eigenvalue has zero real part.

The kernel at zero has dimension one. Indeed, for a real vector Rv=0, any maximizing and any minimizing state propagate their value along outgoing edges; since each state can reach P, maximum=minimum=v_P. Thus v is constant; apply this to real and imaginary parts for complex vectors. Boundedness of exp(tR) rules out any Jordan block at zero of size greater than one. Consequently zero is algebraically simple and the other three eigenvalues of Q have strictly negative real parts.

Now change coordinates by H(P,Ia,Ib,I)^T=(P+Ia+Ib+I,Ia,Ib,I)^T, where

    H=[1 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1].

Direct multiplication gives a block matrix whose first row is zero and whose lower-right block is

    J=[-a-b-e, -a,      -a+f;
       -c,     -c-d-g,  -c+h;
        b,      d,      -f-h].

This is exactly the corrected derivative. Thus det(lambda I−Q)=lambda det(lambda I−J), proving that J is Hurwitz for all eight fast rates>0 and c1,c2≥0. This proves the intended full-domain local stability statement, including the boundary cases omitted by an irreducible-only argument.

## 3. Optional source-specific slow-drift obstruction

This is secondary evidence for the original all-time lead, not the core correction. In Equation (32), set all eight fast rates to 1 and k5=epsilon>0, km5=0. Set

    x(0)=(1/2,1/2,13/6,13/12,13/12,13/24).

At epsilon=0 all six derivatives vanish exactly. The source's integrals are mu1=0 and mu2=7, while y1=x3+x4+x5+x6=39/8. Thus the critical reduced solution remains y1=39/8.

For epsilon>0, nonnegativity and the two conserved integrals imply a compact forward-invariant set: mu1+2 mu2=x1+x2+2x3+3x4+3x5+4x6=14. The symmetry x1=x2=u and x4=x5=v is preserved. Put p=x3,w=x6. The equations reduce exactly to

    u'=-u*p+v-u*v+w,
    p'=-2*u*p+2*v+2*epsilon*w,
    v'=u*p-v-u*v+w,
    w'=2*u*v-(2+epsilon)*w.

The conserved mass is u+p+3v+2w=7. The observable h=p+2v+w satisfies h'=epsilon*w≥0 and 7−h=u+v+w≥0. It therefore has a limit. Every point in the nonempty omega-limit set has w=0, because h is constant along trajectories within that invariant set. Invariance then gives uv=0. If u=0,v>0, u'=v>0 and v remains positive locally, contradicting uv=0. If v=0,u>0,p>0, v'=up>0 gives the same contradiction. Thus the only invariant points with w=0 are (u,v,w,p)=(0,0,0,7) and (7,0,0,0). The latter has h=0 and is excluded by h≥39/8. The omega-limit set is therefore the pure-product equilibrium, and h(t)→7.

It follows, for **every** epsilon>0,

    sup_{t≥0}|h_epsilon(t)−h_0(t)|=7−39/8=17/8.

For each fixed finite t, continuous dependence gives h_epsilon(t)→39/8 as epsilon→0. The limits do not commute. This is a classical loss-of-conservation slow-drift mechanism made exact in the source's own model; it is not claimed as a new general singular perturbation principle.

## Scope and verification limits

The model, coordinates and local stability question come from Eilertsen, Romanovski, Schnell and Walcher, following Gijima and Peacock-López. The derivative discrepancy, exact witness and full-domain verification are the specific contribution here; the generator argument and slow drift are standard tools. The source's broad question about useful uniform-in-time criteria remains unresolved. A separate AI review independently reconstructed the reaction system and checked the arguments; this is not external human peer review or a proof-assistant formalization. See SCOPE.md and VERIFICATION.md.
