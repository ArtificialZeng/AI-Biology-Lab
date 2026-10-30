# Compatibility under shared scalar site-rate variation

We use the strict homogeneous JC triangle images of Currie et al. (2026),
arXiv:2606.26673v1, Eq. (6) and Theorem 3.1. The result below concerns exact
distributions in this three-taxon model, with no coalescence.

## Theorem

Fix the true triangle orientation `N1`, six positive finite edge lengths and
inheritance probability `delta` in `(0,1)`, with uniform JC root law. Let `R`
have any probability law on `(0,infinity)`, finite almost surely. At each site
multiply **all** six edge lengths by `R`, keeping the inheritance probability
fixed. Let

```
(X,Y,Z,T) = E[(x(R),y(R),z(R),t(R))]
```

be the mixed source Fourier coordinates, in the order
`q_ACC,q_CAC,q_CCA,q_CGT`. Let `M1,M2,M3` be the three strict homogeneous
triangle images, allowing their fitted numerical parameters, including fitted
inheritance, to vary over the source's whole strict domain. Then

```
(X,Y,Z,T) is in M1  if and only if  X*Y*Z-T*T > 0.
```

If `X*Y*Z-T*T <= 0`, the distribution belongs to none of `M1,M2,M3`.
In particular, it cannot belong to a wrong homogeneous orientation while
excluding the true one. The analogous statement holds after leaf relabeling.

## Source coordinates and image conditions

Absorb the JC factor `4/3` into positive lengths `L_i`, and put
`a_i(r)=exp(-r*L_i)`. Source Eq. (6) gives

```
x = a2*a3*a5,
y = delta*a1*a3*a4*a5 + (1-delta)*a1*a3*a6,
z = delta*a1*a2*a4 + (1-delta)*a1*a2*a5*a6,
t = delta*a1*a2*a3*a4*a5 + (1-delta)*a1*a2*a3*a5*a6.
```

For every positive finite `r`, this is a point of `M1`. Theorem 3.1 states that
the full iff conditions for `M1` are `0<x,y,z,t<1` and

```
x-t > 0,
t-x*y > 0,
x*y*z-t*t > 0,
Phi_x = t+x*(t-y-z) > 0.
```

Its additional derived inequalities include `t-x*z>0`, `y-t>0` and `z-t>0`.
The condition `x*y*z-t*t>0` is common to all three images under permutation of
`x,y,z`. These model descriptions concern the full distributions: the JC
probability/Fourier linear symmetries leave precisely the four indicated
nonconstant coordinates. Rate mixing preserves those linear relations.

## Monotonicity needed for the proof

The bounded positive functions `x(r),y(r),z(r)` decrease with `r`: each is a
positive sum of products of decreasing edge factors. The additional function

```
h(r)=y(r)+z(r)-t(r)
```

is positive and decreasing as well. To see this directly from the original
network, delete reticulation edge 6 or edge 4 to obtain the two displayed
trees, with respective fixed weights `delta` and `1-delta`. Their star-edge
factors, after composing degree-two paths, are

```
b^(4) = (a1*a4, a2, a3*a5),
b^(6) = (a1*a6, a2*a5, a3).
```

For either star, `y+z-t = b1*(b2+b3-b2*b3)`. Every `b_i` decreases and lies in
`(0,1)`. The bracket is positive and increasing in each of `b2,b3`, since its
partial derivatives are `1-b3` and `1-b2`. It therefore decreases with `r`.
Multiplication by positive decreasing `b1`, and then mixing with the fixed
inheritance weights, preserves this property.

The path composition uses the standard JC identity
`M(a)=a*I+(1-a)*J/4`, where `J` is the all-ones matrix. Since `(J/4)^2=J/4`,
`M(a)M(b)=M(ab)`. Uniform detailed balance permits rerooting the second tree
at its star center. Thus the path argument keeps the original probability
model; it introduces neither a clock nor new inheritance weights.

## Preservation of the orientation-specific inequalities

For bounded functions `f,g` monotone in the same direction,

```
Cov(f(R),g(R))
  = E[(f(R)-f(R'))*(g(R)-g(R'))]/2 >= 0,
```

where `R'` is an independent copy. Boundedness makes the expectations finite.
This is the standard comonotone covariance identity.

Strict positivity, upper bounds 1, and `X-T>0` survive expectation. Further,

```
T-X*Y = E[t-x*y] + Cov(x(R),y(R)) > 0.
```

The pointwise margin is strictly positive by the source image theorem, and
its expectation is strictly positive. The covariance is nonnegative. The same
argument also gives `T-X*Z>0`.

Since `Phi_x=t-x*h`, the model-specific monotonicity above gives

```
T-X*(Y+Z-T) = E[t-x*h] + Cov(x(R),h(R)) > 0.
```

Again the first term is the expectation of a strictly positive source margin,
and the covariance is nonnegative. A positive almost surely bounded random
variable has strictly positive expectation; no lower bound uniform in rate or
moment condition on `R` is required.

All `M1` iff conditions except `D=X*Y*Z-T*T>0` therefore hold automatically.
If `D>0`, the source iff theorem puts the mixed full distribution in `M1`.
If `D<=0`, its common necessary condition excludes all three strict images.
This proves the dichotomy.

## Nonclosure example

Take `a1=...=a6=1/2`, `delta=1/2`, and rates 1 and 2 with equal probability.
The two coordinate vectors are

```
(1/8,3/32,3/32,1/32),
(1/64,5/512,5/512,1/1024).
```

Their average is `(9/128,53/1024,53/1024,33/2048)`, and

```
D = (9*53*53 - 32*33*33)/2^27 = -9567/134217728 < 0.
```

It lies outside every strict homogeneous triangle image. Direct original
ancestral-state summation gives the following probabilities for each individual
pattern in its equality class:

| Class | Number of patterns | Probability of each |
| --- | ---: | ---: |
| AAA | 4 | 1657/65536 |
| ACC | 12 | 1101/65536 |
| CAC | 12 | 1025/65536 |
| CCA | 12 | 1025/65536 |
| CGT | 24 | 879/65536 |

They are strictly positive and their weighted sum is 1. The accompanying exact
checker reconstructs all 64 entries without using contracted stars, then checks
both Fourier directions and the optional path contractions.

## Interpretation and limits

The theorem gives an exact compatibility criterion, not orientation uniqueness.
It permits loss of homogeneous fit and ambiguity among orientations. It makes
no claim about maximum-likelihood projections when every strict image is
incompatible, finite-data error, parameter accuracy or boundary models.
The fitted homogeneous inheritance need not equal the generating inheritance.
Lineage-specific rate variation, rate-dependent inheritance, coalescence,
larger networks and other substitution processes are not covered.

The homogeneous image inequalities, covariance method and tree-topology scalar-
rate precedents are credited in [sources.md](sources.md). The bounded new
consequence is the triangle-orientation compatibility dichotomy. Independent
human confirmation and exhaustive priority assessment remain pending.
