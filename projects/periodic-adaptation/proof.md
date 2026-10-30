# Periodic forcing of a reversible adaptation motif

## Precise claim and scope

Let a,b,c,d,T be strictly positive real numbers, and let u:R->(0,infinity)
be continuous and T-periodic. Consider, for nonnegative initial concentrations,

    x'=(a+c)u(t)-bx-dxy,
    y'=cu(t)-dxy.                                      (1)

This is the mass-action signalling system S<->X, S<->X+Y with externally
fixed S=u. The arrows are S->X at rate constant a, X->S at b, S->X+Y at c,
and X+Y->S at d. It is the reaction motif in Example3.8 of Franco and
Velázquez (JMB2026), with S=source species1, X=3, Y=2. Its constant-input
adaptation and reversible completion are prior facts, not new contributions.
The published PDF prose labels species3 as adaptive, while its displayed
arrows have conservation vector (1,0,1): species2 is the invariant output.
The claim here follows the explicit equations (1) and uses that output y;
the indexing correction alone is not claimed as a substantive new result.

Write

    h=(a+c)/c, r=ad/(a+c), q=ad/c, y0=bc/(ad).

Claim: there is a unique strictly positive T-periodic solution (xp,yp).
Every nonnegative initial state has a global solution converging exponentially
to it (the bound may depend on the initial state). Along that periodic response,

    z=xp-h(yp-y0)>0,
    y0-average(yp)=h*average((yp-y0)^2/z).                (2)

Here average means integration over one period divided by T. Consequently the
average output equals y0 if and only if the input u is constant; every
nonconstant positive continuous periodic input gives a strictly lower average.
This is a result about one explicit motif, not a resolution of the source's
general periodic-input question. Failure of mean adaptation in other circuits,
and contraction-based entrainment, are established literature and credited.

## 1. Global existence and a compact invariant polygon

Set umin=min u>0 and umax=max u. Given a nonnegative initial state (xinit,yinit),
choose

    M >= max(xinit,(a+c)umax/b),
    L >= max(y0,yinit-xinit/h).

The compact convex polygon

    K(M,L)={(x,y): 0<=x<=M, 0<=y, y-x/h<=L}            (3)

contains the initial state. On x=0 the x derivative is positive; on y=0 the
y derivative is positive; on x=M we have x'<=0. For w=y-x/h, exact cancellation
of the input gives

    w'=r*x*(y0-y).                                     (4)

On the fourth face w=L, y=L+x/h>=y0, so w'<=0. The vector field therefore
points into or along every bounding face. Standard first-exit/uniqueness for a
locally Lipschitz vector field on this convex polyhedron proves positive
invariance. Solutions stay in K, with 0<=y<=Y:=L+M/h, and hence extend for all
t>=0. No finite-time blow-up is possible.

For positivity, scalar comparison gives

    x'>=(a+c)umin-(b+dY)x,
    y'>=c umin-dM y.

Thus both concentrations are strictly positive at every t>0, including when
one or both initial values are zero. These comparisons use M>0.

## 2. Existence of a positive periodic response

Use the fixed canonical polygon K(M0,L0), with M0=(a+c)umax/b and L0=y0.
The time-T solution map is continuous and maps this nonempty compact convex
polygon into itself. Brouwer's fixed-point theorem supplies a fixed point.
Periodicity of the vector field and uniqueness of the initial-value problem
extend its orbit to a T-periodic solution on all of R. Positivity for positive
times, together with periodicity, makes both coordinates strictly positive at
all times. In particular, m=min_{[0,T]}xp>0.

## 3. Incremental contraction, uniqueness and global attraction

Compare any two solutions with the same input; put p=x1-x2 and v=-(y1-y2).
The exact product difference x1y1-x2y2=y1*p-x2*v gives

    p'=-(b+d y1)p+d x2 v,
    v'=d y1 p-d x2 v.                                 (5)

Choose solution2 as a positive periodic solution. Solution1 lies in a polygon
(3), so y1<=Y, while x2>=m>0. Put epsilon=b/(2dY)>0 and

    D=|p|+(1+epsilon)|v|,
    lambda=min(b/2,epsilon*d*m/(1+epsilon))>0.

The absolute-value differential inequality, valid almost everywhere (or for the
upper Dini derivative at zeros), is

    D' <=(-b+epsilon*d*y1)|p|-epsilon*d*x2|v|
       <=-b|p|/2-epsilon*d*m|v| <=-lambda D.            (6)

Integration gives D(t)<=exp(-lambda*t)D(0). Thus every nonnegative initial
state converges exponentially to that periodic solution. A second positive
T-periodic solution would also satisfy (6); periodicity forces its distance
from the first to be zero. This proves uniqueness. This weighted-norm argument
is an application of classical incremental-contraction methods, not a new
general contraction theorem.

## 4. Exact mean deficit and its strict sign

All quantities below refer to the positive periodic solution. Put eta=y-y0
and z=x-h eta. Direct substitution into (1), using ad*y0/c=b, gives

    z'=q*x*eta=q*z*eta+q*h*eta^2.                       (7)

Since z is differentiable and periodic it has a minimum. If that minimum were
nonpositive at t*, eta(t*)=(x(t*)-z(t*))/h>0, so (7) would give z'(t*)>0,
contradicting the derivative zero at a periodic minimum. Hence z>0 everywhere.
Divide (7) by z and integrate over one period. The integral of (log z)' is
zero. Since q>0, this yields (2) exactly.

The right side of (2) is nonnegative. Equality implies eta identically zero
by continuity, hence y=y0. Equation (1) then gives x=a u/b and x'=0, so u is
constant. Conversely constant u has the positive equilibrium x=a u/b,y=y0,
which is the unique periodic response. This proves the claimed strict sign.

## 5. Optional quantified small-sinusoid corollary

Fix u0>0, omega>0 and u(t)=u0+delta*cos(omega*t), with |delta|<u0.
Let x0=a u0/b, K=b*d*x0 and B=b+d*(y0+x0). Then, as delta->0,

    y0-average(y_delta)
      = delta^2 * h*c^2*omega^2 /
        (2*x0*((K-omega^2)^2+B^2*omega^2)) + O(delta^3).  (8)

To justify the expansion, the constant-input equilibrium has Jacobian

    J=[[-b-d*y0,-d*x0],[-d*y0,-d*x0]],

with trace=-B<0 and determinant=K>0. Its eigenvalues have strictly negative
real parts, so exp(J*T)-I is invertible. Smooth parameter dependence of the
finite-time ODE map and the implicit-function theorem provide a twice
continuously differentiable periodic branch near delta=0, uniformly over one
period. Uniqueness from Section3 identifies it with the global periodic orbit.
The first variation of y is Re(H(i omega) exp(i omega*t)), where

    H(s)=c*s/(s^2+B*s+K).

Thus its squared period average is |H(i omega)|^2/2. Substitution into the
exact identity (2), with z=x0+O(delta), proves (8); smoothness of all orders
near the positive equilibrium gives the stated Taylor remainder. The leading
coefficient is maximized over omega>0 at omega^2=K. Indeed, setting v=omega^2,
the derivative of v/((K-v)^2+B^2*v) has numerator K^2-v^2. This corollary is
an infinitesimal-amplitude statement; no optimum at finite amplitude is claimed.

## Verification and remaining review

The invariant-region faces, incremental identity, weighted column bounds,
constant-input Jacobian, mean identity and sinusoidal coefficient can be checked
by exact algebra. Such finite checks do not formalize Brouwer, ODE existence,
the contraction argument or the mean theorem. No floating-point simulation is
needed for the conclusions. The exact scope received a fresh-context AI mathematical review. The primary
comparison is bounded; no human specialist validation, proof-assistant
certificate, authenticated runner release gate or priority is claimed.
Substantial AI assistance was used. See provenance/review-receipt.json.
