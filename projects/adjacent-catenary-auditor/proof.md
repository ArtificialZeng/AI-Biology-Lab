# Candidate theorem: adjacent endpoint measurement on a one-leak path

## Scope and notation

Let n>=2. The n-vertex bidirected path has positive independent forward rates p_i=A_(i+1,i), reverse rates q_i=A_(i,i+1), i=1,...,n-1, and one positive leak e at fixed vertex ell. Diagonal entries are minus total outflow. Sole input is u at vertex1, calibrated sole output is x2. Define

D(s)=det(sI-A), B(s)=det((sI-A)_[2:n,2:n]), R(s)=det((sI-A)_[3:n,3:n]),

with determinant of the empty matrix equal to1. Set a=p1+e if ell=1 and a=p1 otherwise. Cofactor expansion gives the transfer function and recurrence

H(s)=p1 R(s)/D(s),    D(s)=(s+a)B(s)-p1 q1 R(s).                 (1)

The polynomial coefficient map c has all n coefficients of the monic D below its leading term, and all n-1 coefficients of p1R, including its leading coefficient p1. Thus c maps 2n-1 parameters to 2n-1 coordinates. Fixed zero coefficients are omitted. These coefficients are identifiable from the rational transfer function on the open coprime locus. We retain the same polynomial c at special cancellation points, as is standard for the coefficient-Jacobian singular locus; pointwise cancellation is addressed separately below.

**Theorem.** For every fixed ell, at every strictly positive parameter point:

* If ell>=2, Dc loses rank exactly when R(-p1)=0.
* If ell=1, define the monic degree n-2 polynomial T by B(s)=sT(s)+q1R(s). Then Dc loses rank exactly when K=q1R(-(p1+e))-eT(-(p1+e))=0.

At each singular point the rank is exactly2n-2. No boundary/zero-edge assertions are made. For n=2, R=T=1, so ell=1 gives q1=e and ell=2 gives the empty locus.

## Complete proof

### A positive endpoint reconstruction lemma

For an m-vertex irreducible tridiagonal path, write its diagonal decay constants as d_j and the positive neighboring products as b_j. Its full characteristic polynomial F_j for the suffix beginning at j satisfies

F_j(s)=(s+d_j)F_(j+1)(s)-b_j F_(j+2)(s).

The pair of monic polynomials (F_j,F_(j+1)) uniquely and regularly determines all d's and b's: polynomial division gives quotient s+d_j and remainder -b_j F_(j+2). Its nonzero leading coefficient determines b_j, and recursion continues. All divisions in this reconstruction use positive b_j. The forward polynomial map and this rational inverse are local diffeomorphisms throughout the positive domain. Consecutive F_j,F_(j+1) are coprime, because a common root would propagate to the terminal polynomial1. Positive tridiagonal matrices are similar to real symmetric irreducible Jacobi matrices, so each F_j has simple roots (the usual elementary irreducible Jacobi fact; alternatively strict interlacing follows from the same recurrence).

### The physical rate coordinates used below are regular

For ell>=2, keep q1 separately, regarding it as a known boundary loss from the tail2,...,n. The remaining tail physical rates and leak e are locally regularly equivalent to the tail decay constants d2,...,dn and products b2,...,b_(n-1). To invert, start from the left endpoint when j<ell: p2=d2-q1 (unless ell=2), q2=b2/p2, then p3=d3-q2, and continue up to ell. From the right endpoint when j>ell: q_(n-1)=dn (unless ell=n), p_(n-1)=b_(n-1)/q_(n-1), then q_(n-2)=d_(n-1)-p_(n-1), and continue down to ell. At ell the residual diagonal is e. Every divisor is a positive physical rate. The same statement for ell=2 or ell=n uses an empty sweep on that side. Thus (p1,q1,B-coefficients,R-coefficients) are local regular coordinates.

For ell=1, the tail has no leak except its boundary loss q1. Its decay constants/products determine its rates by sweeping from the right endpoint all the way to node2, and finally q1=d2-p2. These are regular positive coordinates. Consequently (p1,e,B-coefficients,R-coefficients) are local regular coordinates. Also B(0)=q1R(0): when q1=0 the tail is closed and has determinant0, and B is affine in q1 with coefficient R. Moreover R(0)>0: every term of the rooted-tree determinant is positive (on a path it is the product q2...q_(n-1)). Hence T=(B-q1R)/s is polynomial, monic of degree n-2.

### Differential for ell>=2

If Dc vanishes on a tangent direction, the numerator coefficients force delta p1=0 and delta R=0. Equation(1) then gives

0=delta D=(s+p1)delta B-p1 delta q1 R.

The polynomial delta B has degree at most n-2. A nonzero kernel exists iff s+p1 divides R, i.e. R(-p1)=0. If so, choose any delta q1 and set delta B=p1 delta q1 R/(s+p1); these are legitimate infinitesimal coordinate changes. Otherwise delta q1=delta B=0. The kernel is exactly one-dimensional when nonzero, proving the assertion including corank.

### Differential for ell=1

Again delta p1=delta R=0. Since q1=B(0)/R(0), delta q1=delta B(0)/R(0). The constant coefficient D(0)=e B(0) yields, for delta D=0,

delta B(0)=-q1R(0) delta e/e,       delta q1=-q1 delta e/e.       (2)

At s=-a, where a=p1+e>0, the variation of(1) is

0=delta e B(-a)-p1 delta q1 R(-a)
 = (delta e/e)[e B(-a)+p1q1R(-a)]
 = (a delta e/e)[q1R(-a)-eT(-a)] = (a delta e/e)K.             (3)

If delta e=0, equations(2) and(1) force delta B=0. Hence a kernel requires K=0. Conversely if K=0, take delta e=1, delta q1=-q1/e, and

delta B(s)=[-B(s)-(p1q1/e)R(s)]/(s+a).

Equation(3) guarantees that this is a polynomial of degree at most n-2. At s=0 it has exactly the constant coefficient in(2), since a=p1+e; therefore it is a valid tangent direction in the regular coordinates and satisfies delta D=0. Again the kernel is exactly one-dimensional. This completes the theorem.

## Exact n=4 specialization (all fixed leak locations)

Let a=p1+e and

R0(s)=s^2+(p3+q2+q3)s+q2q3,
T0(s)=R0(s)+p2(s+p3+q3).

The positive coefficient-Jacobian singular sets are:

ell=1: q1R0(-a)-eT0(-a)=0;
ell=2: R0(-p1)=0;
ell=3: p1^2-p1(p3+q2+q3+e)+q3(q2+e)=0;
ell=4: p1^2-p1(p3+q2+q3+e)+q2q3+e(p3+q2)=0.

All are nonempty. Exact positive witnesses in parameter order(p1,p2,p3,q1,q2,q3,e):
ell1:(1,1,1,1,1,1,1);
ell2:(1,1,1,1,2,2,1);
ell3:(1,1,1,1,1,2,1);
ell4:(1,1,1,1,1,1,1).

The independently computable square determinants in small-path-jacobians.json are polynomials, not sampled evidence: n4 ell1 has prefactor p1^2 p2 q2^2 q3 times K; ell2 has prefactor p1^3 p2 q2^2 q3 times R0(-p1); ell3 has prefactor p1^3 p2^2 q2 q3 times its displayed factor; ell4 has prefactor p1^3 p2^2 p3 q2 times its displayed factor. Every prefactor is strictly positive on the claimed domain, so this is also a stand-alone exact finite proof of the original n4 classification. The general proof does not infer an infinite statement from this finite computation.

## Common-factor cancellation and the data contract

By endpoint coprimality B,R have no common root. Equation(1) therefore gives gcd(D,R)!=1 iff R(-a)=0, and the common factor is exactly s+a (R has simple roots). For ell>=2, this agrees exactly with the coefficient-Jacobian singular condition. For ell=1, cancellation is disjoint from K=0: if R(-a)=0, then B(-a)=-aT(-a)!=0, and hence K=-eT(-a)!=0. Thus the polynomial coefficient map can be regular where a specially canceled zero-initial-condition transfer function has lower degree. For example n4 ell1 at (1,1,1/2,1,3,3,1) has R0(-2)=0 but K=-3/2!=0.

If the intended data map is H itself, with its differential taken as the formal rational function delta H=(D delta N-N delta D)/D^2, its positive differential rank is2n-2 precisely on the union of the coefficient singular set and the cancellation set, and2n-1 elsewhere. Proof: off cancellation, zero delta H implies delta N=delta D=0 by coprimality and degrees. At a degree-one common factor g=s+a, the ambient coefficient-to-rational map has a one-dimensional kernel, consisting of delta D=(D/g)t, delta N=(N/g)t. For ell1, Dc is invertible at cancellation, so the physical differential has corank1. For ell>=2, Dc has corank1 there; its image intersects that ambient kernel trivially. Indeed delta p1=0, and the image relation at s=-p1 is delta D(-p1)=-p1q1 delta R(-p1). An ambient common-factor variation t would instead require B(-p1)t=0 after differentiating(1) in s; B(-p1)!=0, so t=0. Thus its rank remains2n-2. This is a differential-rank statement, not an unsupported assertion that Jacobian rank loss alone proves a global or pointwise positive-fiber classification.

## Scope, novelty and gaps

Mathematical gap list: none identified by the proposer; fresh audit required. Original n4 frozen question is completely classified; the n>=2 result follows by the explicit regular coordinate changes above. Bounded primary comparison found no equivalent adjacent-port exceptional-locus theorem; endpoint continued-fraction inversion itself is standard and is not claimed new. Historical priority is not established, and some older catenary inverse literature was inaccessible. No manuscript or PDF has been created. All review remains AI-assisted until external human review occurs.
