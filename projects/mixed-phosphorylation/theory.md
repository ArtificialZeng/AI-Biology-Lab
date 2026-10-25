# Independent mathematical derivation — research pass 1

> Public archive note: this derivation preserves research-stage wording. The current public claim and review status are defined in `README.md` and `claim.json`; internal AI checks are not external human peer review.


Status: researcher derivation and exact symbolic certificate; **not independently
refereed**. Fixed Table-1 rates throughout. No claim of novelty is made here:
the literature worker has found a close 2025 meeting abstract which must be
compared before a contribution is described as new. The full original
Conjecture 6.2 remains unresolved. In particular, instability does not prove
the existence of an attracting periodic orbit.

## 1. Exact coordinates and singular limit

Write K=K_tot, P=P_tot, S=S_tot and use the five slow coordinates

    q=(a,b,w,u,z)=(x3,x4,x5,x7,x9).

Let p=P-u-z, v=S-a-b-w-u-z, A=a+b, epsilon=1/K and
y=x1/epsilon=K*x1. The three conservation laws give

    x2=K-a-b, x6=p, x8=v-epsilon*y.

Consequently the exact six-dimensional reduced system is

    epsilon*y' = -y + epsilon*A*y + a + 100*z,
    a' = y*(1-epsilon*A)-2*a,
    b' = a-b,
    w' = b-100*w*p+u,
    u' = 100*w*p-(19/10)*u,
    z' = 3*p*(v-epsilon*y)-101*z.                       (1)

These are the original nine mass-action equations after exact conservation
elimination and a coordinate rescaling, not a Michaelis–Menten approximation.

The epsilon=0 critical manifold is y=h0(q)=a+100*z. Its fast normal
derivative is -1. The reduced five-dimensional slow vector field is

    F0(q)=(-a+100*z,
            a-b,
            b-100*w*p+u,
            100*w*p-(19/10)*u,
            3*p*v-101*z).                              (2)

It agrees coordinate-for-coordinate with `src/model.py:limit_rhs`. On any
compact neighborhood in the strict physical interior p>0,v>0,a,b,w,u,z>0,
the critical manifold is uniformly normally attracting. Standard local
Fenichel theory applies for epsilon>0 sufficiently small. An appropriate
compact neighborhood, rather than a global claim at the boundary, is intended.

## 2. First correction to the slow manifold and vector field

Define Lh0=Dh0*F0. Direct differentiation gives

    Lh0=-a-10000*z+300*p*v.

The invariance equation for y=h0+epsilon*h1+O(epsilon^2) gives

    h1=A*h0-Lh0
      =(a+b)*(a+100*z)+a+10000*z-300*p*v.               (3)

Thus the vector field on the invariant slow manifold, in these q coordinates,
is F0+epsilon*F1+O(epsilon^2), where

    F1(q)=(a+10000*z-300*p*v,
            0,
            0,
            0,
           -3*p*(a+100*z)).                            (4)

The direct correction in w'=x5' vanishes in these coordinates; this does NOT
say that w's orbit or period has zero first-order correction. The perturbed
trajectory changes through a and z. Expressions (3)-(4) are checked as exact
polynomial identities through first order by `src/theory_exact.py`.

The O(epsilon^2) expansion means a smooth graph and vector-field expansion
on a specified compact neighborhood. It is not an assertion that arbitrary
initial data are uniformly close to the slow manifold during the initial
fast layer.

## 3. Exact positive equilibrium and large-K expansion

At any positive equilibrium the steady flux equations imply

    a=b=100*z,
    u=(1000/9)*z,
    p=P-(1009/9)*z,
    w=(19/9)*z/p,
    x8=(101/3)*z/p,
    x1=200*z/(K-200*z).                                (5)

Define

    G(z)=(2809/9)*z+(322/9)*z/[P-(1009/9)*z].

The equilibrium is determined by the single equation

    S=G(z)+200*z/(K-200*z),
    0<z<min(K/200,9*P/1009).                            (6)

The right-hand side is strictly increasing in z, vanishes at zero, and
diverges at the right endpoint. Therefore (6) has exactly one positive
solution z_K for every K,P,S>0. This is an independent elementary verification
of an equilibrium-uniqueness fact already known for the original network.

Let alpha be the unique root G(alpha)=S in (0,9*P/1009). Implicit
differentiation of (6) shows z_K is strictly increasing in K. It converges to
alpha as K tends to infinity. Since G'(alpha)>0, the implicit-function theorem
at epsilon=0 gives

    z_K=alpha-[200*alpha/G'(alpha)]/K+O(K^-2),
    x1_K=200*alpha/K+O(K^-2),
    w_K=w_inf-[200*alpha*w'(alpha)/G'(alpha)]/K+O(K^-2), (7)

where w(z)=(19/9)*z/[P-(1009/9)*z] and w'(alpha)>0. Thus the equilibrium
output w_K depends strictly on K, with a nonzero negative first correction
from below. This says nothing by itself about the dependence of a periodic
output: the original conjecture concerns periodic orbits, not equilibria.

## 4. Exact certificate: the P=8,S=40 limiting equilibrium is unstable

For P=8,S=40, G(alpha)=40 is equivalent to

    Q(alpha)=2834281*alpha^2-568386*alpha+25920=0.        (8)

The physical root is the unique root in

    70122383354/10^12 < alpha < 70122383355/10^12.        (9)

Indeed Q has respectively positive and negative endpoint values and Q'<0
throughout that interval; the entire interval lies in (0,72/1009). Equations
(5) with x1=0 produce the limiting equilibrium q_infinity in the strict
five-dimensional slow physical interior.

Let det(lambda*I-DF0(q_infinity))=lambda^5+c1*lambda^4+
c2*lambda^3+c3*lambda^2+c4*lambda+c5. Define Hurwitz determinants Delta_j
with the usual first rows (c1,1,0,...), (c3,c2,c1,...), etc. The exact script
recomputes this Jacobian directly from (2), computes its characteristic
polynomial, and reduces all coefficients and Delta_1,...,Delta_4 in the
quadratic field Q[z]/(Q). Each reduced expression is affine in z. Rational
endpoint evaluation on (9), with no floating-point sign decision, certifies

    c1,c2,c3,c4,c5>0;
    Delta_1>0, Delta_2>0, Delta_3>0, Delta_4<0.          (10)

Display-only decimal values are

    (c1,c2,c3,c4,c5) ~=
      (277.18522755,13847.78576449,33905.06952326,
       27888.45604098,78896.52762001),
    Delta_1 ~= 277.1852,
    Delta_2 ~= 3804496.579,
    Delta_3 ~= 126870873850.6,
    Delta_4 ~= -617736652928122.9.

For a monic fifth-degree polynomial with these nonzero determinants, the
Routh first column is

    1, c1, Delta_2/c1, Delta_3/Delta_2,
       Delta_4/Delta_3, c5,

whose signs are +,+,+,+,-,+. The Routh root-count theorem therefore gives
exactly TWO roots in the open right half-plane. All determinants through
Delta_5=c5*Delta_4 are nonzero, so no roots lie on the imaginary axis.

The complete rational expressions and their rational endpoint values are in
`evidence/theory_exact_certificate.json`. Re-running `src/theory_exact.py`
reconstructs the evidence; the JSON flags are not assumed by that script.

### Proposition: existence of a finite K0 with instability for all K>K0

For the fixed rates and P=8,S=40 there exists K0>0 such that, for every K>K0,
the unique positive equilibrium of the original model, restricted to its
six-dimensional compatibility class, has exactly two eigenvalues with positive
real part and four with negative real part.

Proof. Equation (6) and the implicit-function theorem show that the equilibrium
in the scaled (y,q) coordinates depends smoothly on epsilon near zero and tends
to (h0(q_infinity),q_infinity). Write (1) as epsilon*y'=g and q'=f.
Multiplying its characteristic polynomial by epsilon gives the block
determinant

    det[[epsilon*lambda-g_y, -g_q],
        [             -f_y, lambda*I-f_q]].

At epsilon=0 the first diagonal entry is 1, g_q=Dh0 and the Schur complement
is lambda*I-f_q-f_y*Dh0=lambda*I-DF0. Therefore five characteristic roots tend,
with multiplicity, to the five roots certified in (10). None has zero real
part. The remaining root is real and equals -1/epsilon+O(1): there is only one
root left, and this asymptotic follows also from the trace of the full
Jacobian. It is negative for epsilon sufficiently small. Root continuity now
preserves the two right-half-plane roots and three left-half-plane slow roots.
This proves the proposition. No explicit numerical value of K0 is certified.

Scope: this establishes the EXISTENTIAL total-amount assertion of original
Conjecture 6.2(1), for the particular witness (P,S)=(8,40), subject to a fresh
review of the derivation and certificate. It does not establish part (2),
periodic-orbit existence, or the global conjecture about unstable equilibria.

## 5. Conditional periodic-orbit theorem and period correction

Assume F0 at given fixed P,S has a periodic orbit Gamma0 in the strict physical
interior, of least period T0, with simple unit multiplier and every other
Floquet multiplier strictly inside the unit circle. Then for all K sufficiently
large, the original finite-K system has a nearby attracting periodic orbit.
Reason: the normally attracting slow manifold persists locally; its vector
field is a smooth perturbation (2)-(4); a hyperbolic periodic orbit persists.
The extra fast multiplier is strongly contracting. This is conditional until
the orbit and its hyperbolicity are certified for the chosen P,S. Numerical
shooting alone does not discharge the assumption.

For the first period correction, let psi(t) solve the periodic adjoint equation

    psi'=-DF0(q0(t))^T*psi,
    psi(t)^T*F0(q0(t))=1.

Then the smooth period expansion is

    T(K)=T0+T1/K+O(K^-2),
    T1=-integral_0^T0 psi(t)^T*F1(q0(t)) dt.             (11)

For completeness, put the perturbed orbit on the fixed interval [0,T0]. Its
first variation eta satisfies

    eta'=DF0*eta+F1+(T1/T0)*F0,
    eta(T0)=eta(0),

with one phase condition. Multiplying by psi and integrating gives (11).
If the maximum and minimum of w0(t) are nondegenerate and unique modulo period,
the amplitude derivative is eta_w(t_max)-eta_w(t_min); the shifts of the extrema
times do not contribute because w0' vanishes there. No nonzero value of T1 or
amplitude correction is proved in this file.

## 6. Exact output identifiability: a partial obstruction, and its gap

The following is exact but does not settle the source's x5-only statement.
If two positive periodic solutions at fixed rates,P,S have the same functions
(w(t),p(t))=(x5(t),x6(t)) up to a common time shift, then their kinase totals K
are equal. To see this, u solves u'+(19/10)u=100*w*p. This linear forced ODE
has a unique periodic solution. Next z=P-p-u,
x8=(z'+101*z)/(3*p), b=w'+100*w*p-u, a=b'+b, and

    x1=S-a-b-w-u-x8-z,
    K=a+b+(a'+2*a)/x1.                                 (12)

All denominators are positive. Thus the pair of periodic outputs reconstructs
K uniquely.

For w alone, p may change as K changes. A sufficient missing lemma would be
uniqueness of the positive periodic response of the forced phosphatase
subsystem for every periodic w. I have not proved that lemma. The subsystem
contains nonlinear sequestration, so ordinary uniqueness of solutions from an
initial value does not imply unique periodic response. Consequently this pass
does NOT prove that no continuous family of finite K can share the same
nonconstant w(t), and differences seen at P=8,S=40 alone would not disprove
the source's existential quantifiers.

## 7. Dependencies and remaining decisive work

* Exact algebra: SymPy rational operations and the checked affine sign ranges;
  all formulas should be independently reconstructed for review.
* Analytic tools: implicit-function theorem, polynomial-root continuity,
  Routh root count; Fenichel and hyperbolic-cycle persistence only for Section 5.
* Proven-scope candidate: Section 4's eventual equilibrium instability, with
  explicit algebraic witness totals, plus first-order reduction formulas.
* Unproved: a certified limit cycle, any certified nonzero period correction,
  explicit K0, global periodic-response uniqueness, and original part (2).
* Novelty: compare the 2025 Grasse–Conradi abstract before release. The likely
  new contribution, if any, may need to be the finite-K first correction and
  certified nonzero period dependence, rather than the singular limit itself.
