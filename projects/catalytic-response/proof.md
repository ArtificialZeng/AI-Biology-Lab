# Exact six-rate classification for the catalytic Fig. 14B model

## Statement

Fix a,b,c,d,e,f>0. For each x>0 the model in source.md has exactly one
nonnegative steady state with A+B+C=1; it is interior and smooth in x.
Its C response is strictly decreasing and parametrizes an interval (0,C0).
Let

    p=a+c, q=b+c, T=a d-p q.

If T<=0, B'(x)>0 for every x>0. If T>0, B has precisely one critical
point x*, with B'(x)>0 for x<x* and B'(x)<0 for x>x*. Thus x* is a
strict global maximum. These assertions hold equally for F_a=2 log x.
In particular, the maximum number of finite-input critical points of B is
exactly one; finite-input degenerate stationary points do not occur.

For T>0 define r=sqrt(p^2+T)>p. The maximum occurs at

    C*=(r-p)/d,
    A*=(r+b-a)/d,
    B*=c/p+(r-p)^2/(d p),
    x*=[c A*-b C*+sqrt((b C*-c A*)^2+4 e f B* C*)]/(2 e C*).

All these coordinates are physical, and the formulas hold at arbitrary positive
e,f. As x→0, B→0. As x→infinity, B→c/(a+c). Consequently, the
overshoot height above the high-input limit is exactly (r-p)^2/(dp).

## Proof

At any nonnegative steady state a zero coordinate is impossible: A=0 forces
aB+bC=0; C=0 forces cA+fB/x=0; B=0 forces C=0 from its own equation,
then contradicts the C equation. Thus all steady states are interior.

Use B=1-A-C in the A equation. With

    D(C)=p+dC,
    M(C)=a+(b-a)C,
    N(C)=c+(d-b-c)C-dC^2,

that equation is equivalent to A=M/D and B=N/D. The polynomial N has
negative leading coefficient and positive constant, hence exactly one positive
root C0. Since N(1)=-b<0, 0<C0<1. For 0<=C<=C0,
M=a(1-C)+bC>0. Hence A>0 and B>0 exactly on 0<C<C0;
at C0, B=0 and A=1-C0>0.

Substituting these A,B into the C equation and multiplying its negative by
xD gives

    H(C,x)=e C D x^2+(b C D-c M)x-f N
          =d(bx+e x^2+f) C^2
           +[a(b+c)x+(a+c)e x^2+f(b+c-d)] C
           -c(ax+f).

For each x>0, H is a quadratic in C with positive leading coefficient and
negative constant. Its two roots have opposite signs; it has exactly one
positive root, at which H_C>0. At C=0, H<0. At C=C0, the A equation
with B=0 gives cA=bC0-dAC0, so

    b C0 D-c M=d A C0 D>0.

Therefore H(C0,x)>0. The unique positive root lies in (0,C0) and gives
an interior steady state. The remaining B equation follows from conservation.
This proves existence and uniqueness without admitting any elimination branch
outside the simplex. Since H_C>0, the state depends smoothly on x.

At the physical root, H considered as a quadratic in x has leading coefficient
eCD>0 and constant -fN<0. Its unique positive root therefore satisfies H_x>0.
Implicit differentiation gives

    C'(x)=-H_x/H_C<0.

For completeness, C maps (0,infinity) onto (0,C0). Monotonicity gives
endpoint limits. At x→0, H(C,x)→-fN(C), so the limit is C0.
At x→infinity, division by x^2 gives e C D→0, so the limit is 0.
It follows that B→0 at the first endpoint, whereas A→a/p, B→c/p at
the second endpoint.

Direct differentiation of B(C)=N(C)/D(C) gives

    dB/dC = S(C)/D(C)^2,
    S(C)=T-2 d p C-d^2 C^2.

Here S'(C)=-2d(p+dC)<0. If T<=0, S(C)<0 for every C>0; with
C'(x)<0 this proves B'(x)>0, including equality T=0.

Suppose T>0. Then S(0)>0. Since N'(C0)<0 and N(C0)=0,
S(C0)=N'(C0)D(C0)<0. Hence S has exactly one zero C* in (0,C0).
Solving S(C*)=0 gives C*=(sqrt(p^2+T)-p)/d. Because C(x) decreases
strictly from C0 to 0, it crosses C* once, at x*. For x<x*,
C>C* and S<0, hence B'(x)>0. For x>x* the two signs reverse, giving
B'(x)<0. At the crossing, B''(x*)=S'(C*)[C'(x*)]^2/D(C*)^2<0;
the critical point is nondegenerate and the global maximum.

At C*, D=r and

    a d+p(a-b)=p^2+T=r^2.

The identity A=(b-a)/d+r^2/(dD) gives the asserted A*. Substituting
into 1-A-C yields B*=c/p+(r-p)^2/(dp). Their positivity follows
from C* in (0,C0), already proved. Solving the C equation as a quadratic
in x yields the displayed positive x*. This completes the proof.

## Scope and gap list

There is no known mathematical gap in this argument, pending fresh independent
review. It proves the stated full-six-rate theorem only for the exact conserved
three-species catalytic Fig. 14B model. It does not classify arbitrary nonlinear
networks, the consumed-catalyst Fig. 14C system, dynamical attractors, or the
parent linear structural-compatibility problem. The novelty comparison remains
a separate source-verification task.
