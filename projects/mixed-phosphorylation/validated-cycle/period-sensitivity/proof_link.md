# Period-sensitivity candidate: formula, enclosure and analytic link

Status: Mac interval development checks; Linux x86/FILIB reproduction and a fresh audit are pending. This is a proposed addition to the previously certified limiting periodic orbit, not a replacement for that certificate.

## Objects

Let q=(u,v,w,a,b), p=8-a-b, s=40-u-v-w-a-b, h0=u+100b and A=u+v. The exact six-dimensional slow-fast system in physical time has

ε y'=-y+ε A y+h0,
q'=(y(1-εA)-2u, u-v, v-100wp+a, 100wp-(19/10)a, 3p(s-εy)-101b).

On a locally invariant slow graph y=h(q,ε)=h0(q)+εh1(q)+O(ε²), the invariance equation ε Dh F=-h+εAh+h0 gives h1=Ah0-Dh0 F0. Thus the slow vector field is F0+εF1+O(ε²), where

F1=(u+10000b-300ps,0,0,0,-3p(u+100b)).

`verify_first_jet.py` verifies the zero and first-order invariance coefficients and exact equality between the implemented C++ family and F0+eF1 using rational polynomial identities. This formal identity does not by itself prove the existence of a smooth invariant graph.

## Poincare sensitivity formula

Use the fixed section u=13/4, positive crossing, and z=(v,w,a,b). Let P(z,e) be the next positive-crossing return point of q'=F0(q)+eF1(q), and τ(z,e) its return time. The existing certified orbit is a fixed point z* at e=0. For its local continued fixed point z(e),

(D_z P-I) z_e=-D_e P,
T1 = D_e τ + D_z τ z_e.

These are derivatives evaluated at (z*,0). The four-dimensional fixed-point derivative is invertible because the previously certified orbit is hyperbolic; the new checker also establishes invertibility with a preconditioned norm bound on its full box.

For the augmented ODE (q,e)'=(F0(q)+eF1(q),0), let M be the derivative of its flow at the return time before including the crossing-time derivative. The section equation gives

D_j τ = -M[0,j]/f_u(return),
DP[i,j] = M[i,j] + f_i(return) D_j τ.

The augmented e column of M and DP supplies the parameter derivative. CAPD computes the flow derivative and `computeDP` includes the crossing correction. No finite difference is used.

## Interval linear solve

A_box encloses D_zP-I and b_box encloses -D_eP at every initial point of X. R is the same explicit approximate inverse as the first certificate. E=I-R A_box and η=||E||∞<1; detR excludes zero. Therefore A at the actual fixed point is invertible and the true derivative obeys z_e=R b+E z_e. The initial enclosure ||z_e||∞ <= ||Rb||∞/(1-η) follows from the Neumann series. Repeated intersections with Rb+E Z preserve a rigorous enclosure of the true z_e. Substitution into the time derivative formula encloses T1.

The checker revalidates a Krawczyk strict inclusion in its own smaller X at e=0. Since X is contained in the earlier radius-1e-8 box, earlier local uniqueness identifies it with the previously certified orbit; no smaller-box location is assumed. The neutral e direction is a constant parameter and is excluded from the four-dimensional fixed-point solve.

## Link to the real finite-K period

The additional analytic assumptions to audit are: a C² slow graph exists jointly in q and ε on a fixed local compact neighborhood of the strictly positive certified limit cycle; the induced slow vector field has the displayed first jet in C¹; and the hyperbolic periodic branch and transverse return time inherit C² parameter dependence. These are the local smooth Fenichel and implicit-function steps, with fast normal derivative -1.

If these hypotheses are accepted, the exact finite-K orbit and the artificial family F0+eF1 have identical value and first derivative at ε=e=0. The periodic-branch implicit derivative is determined by those data, so the checked T1 is the actual coefficient in

T(K)=T∞+T1/K+O(K^-2), K→∞.

The variable t is the original model's time, not fast time. The nonzero T1 establishes eventual period dependence for this fixed parameter choice and this local periodic branch. It does not give an explicit remainder constant, numerical K0, global uniqueness, or a disproof of the source's existential statement about scalar orbit projections. If T1>0, the local branch approaches T∞ from above and, with the stated smoothness, its period eventually decreases as K grows.
