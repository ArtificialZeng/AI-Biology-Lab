# Arbitrarily many exclusive peaks with one vacant site


## Statement and model

For each integer p >= 2 there are finite, positive rational parameters
K_i, alpha_i, i=1,...,p+1, for which the exclusive concentration

    y(X) = X prod_{i=1}^p (U_i(X)-1) / prod_{i=1}^{p+1} U_i(X), X>0,

has at least p distinct strict local maxima. Here U_i is the unique root >1 of

    U_i^2 + (K_i X-alpha_i-1) U_i - K_i X = 0.

Only site p+1 is vacant. The assertion is a lower bound, not an exact count
or an upper bound for arbitrary numbers of sites. It concerns this equilibrium
model and does not constitute biological validation.

## Source identities, derived on the positive branch

Write q=KX, A=1+alpha and d=sqrt((q-A)^2+4q)>0. The quadratic is negative
at U=1 and has root product -q<0, so exactly one root is positive and that
root is >1. Explicitly U=(A-q+d)/2. All these functions are real analytic
for X>0. Define

    T = q/(U^2+q) = (1+(q-A)/d)/2,
    mu = (U-1)T = (A+q)/(2d)-1/2.

The identities follow by substituting the quadratic. Differentiating it gives
DU=-U(U-1)T, with D=X d/dX. Consequently

    D log((U-1)/U)=-T,    D log(1/U)=mu,
    F(X):=X y'(X)/y(X)=1-sum_{i=1}^p T_i(X)+mu_{p+1}(X).

Since y(X)>0 for X>0, F and y' have the same sign. These identities follow
from the displayed quadratic and do not require a limiting concentration.

## Rational centers and cutpoints

Put c_1=1/10 and c_j=(2j-3)/(2j-1) for j=2,...,p. For j=2,...,p put

    a_j=(4j-7)/(4j-3),    b_j=(4j-5)/(4j-1),
    a_*=p/(p+1),          b_*=3.

All these numbers are positive rational. For j>=2, write them as

    c_j=1-4/(4j-2), a_j=1-4/(4j-3), b_j=1-4/(4j-1).

The denominators are positive. Thus a_j<c_j<b_j, and for j<p,
b_j<a_{j+1}. Also c_1=1/10<a_2=1/5. For j>=3 one has
c_{j-1}<b_{j-1}<a_j, and for j<p one has b_j<a_{j+1}<c_{j+1}.
It follows that exactly j-1 centers lie below a_j, and exactly j centers
lie below b_j. Finally

    b_p=1-4/(4p-1) < 1-1/(p+1)=a_* <1<3,

because 4/(4p-1)>1/(p+1). Every center is below a_* and below 3.
No cutpoint is equal to a center or to 1. The ordered 2p cutpoints are
a_2,b_2,...,a_p,b_p,a_*,3; their p paired open intervals are disjoint.

## The finite rational construction

For a positive rational cutpoint x define the rational number

    E_p(x) = sum_{i=1}^p x*c_i/(x-c_i)^2
             + x*(1+x)/|x-1|^3.

Let Emax be its maximum over the 2p cutpoints, and set

    M=1+ceil(8 Emax),
    alpha_i=M-1  (all i),
    K_i=M/c_i    (1<=i<=p),  K_{p+1}=M.

The cutpoints avoid all denominator zeros, so Emax is finite positive
rational, M is a finite integer >1, and all parameters are positive rational.
This specifies one parameter vector explicitly for every p, with no
existential choice left unresolved.

Under this substitution, for r=X/c_i,

    T_i^M(X) = (1+(r-1)/sqrt((r-1)^2+4r/M))/2,
    mu^M(X) = ((1+X)/sqrt((X-1)^2+4X/M)-1)/2.

At X different from the centers and 1 define the auxiliary function

    F_infty(X) = 1 - #{i: c_i<X} + mu_infty(X),
    mu_infty(X)=X/(1-X) if X<1, and 1/(X-1) if X>1.

These are pointwise limits as M tends to infinity, used only to bound the
finite-M logarithmic derivative. For X<1,
F_infty(X)=1/(1-X)-#{i:c_i<X}. The exact center counts above yield

    F_infty(a_j)=1/4, F_infty(b_j)=-1/4 (j=2,...,p),
    F_infty(a_*)=1,   F_infty(3)=3/2-p <= -1/2.

## Explicit error estimates, valid before taking a limit

For e>=0,

    0 <= 1-(1+e)^(-1/2) <= e/2.

For example this follows by integrating
(1/2)(1+t)^(-3/2)<=1/2 from 0 to e. At X!=c_i take
e=4Xc_i/[M(X-c_i)^2]. The finite T_i formula gives exactly

    |T_i^M(X)-1{X>c_i}| = (1-(1+e)^(-1/2))/2
                          <= Xc_i/[M(X-c_i)^2].

For the vacancy take e=4X/[M(X-1)^2]. Its finite and auxiliary formulas give

    |mu^M(X)-mu_infty(X)|
       = (1+X)/(2|X-1|) * (1-(1+e)^(-1/2))
       <= X(1+X)/[M|X-1|^3].

The triangle inequality therefore proves, at every cutpoint x,

    |F_M(x)-F_infty(x)| <= E_p(x)/M <= Emax/M < 1/8,

where the strict last inequality uses M=1+ceil(8Emax)>8Emax.
In particular F_M(a_j)>1/8 and F_M(b_j)<-1/8; also
F_M(a_*)>7/8 and F_M(3)<3/2-p+1/8<=-3/8.
These are signs for a single finite positive rational parameter vector,
not merely for a singular or vanishing limiting curve.

## From derivative signs to strict local maxima

Fix any one of the p pairs (a,b). The function y is positive and real
analytic on an open neighborhood of [a,b], so it attains a maximum there.
The endpoint signs y'(a)>0 and y'(b)<0 exclude both endpoints as maximizers.
Hence a maximizer z lies in (a,b). The same derivative signs show that y is
not constant. A nonconstant real analytic function on the connected domain
(0,infinity) cannot have y-y(z) vanish on a sequence converging to z.
Thus z is isolated in its local level set, and because it is a maximizer,
y(X)<y(z) for all sufficiently nearby X!=z. It is a strict local maximum.
The p intervals are disjoint, so their maxima are distinct. This proves
the assertion for every integer p>=2.

## Evidence and boundaries

The all-p proof above is decisive. The accompanying p=2,...,6 exact cases
test its implementation and directly isolate the radicals of the original
U-quadratic using rational arithmetic; they are not a proof by extrapolation.
The exact checker validates interval signs from U rather than relying on
the simplified finite-M formulas or the asymptotic estimate.
