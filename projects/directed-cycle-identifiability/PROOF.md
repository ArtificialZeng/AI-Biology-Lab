# Exact candidate: endpoint counterexample and an interior-leak repair

Status: exact mathematical contribution checked in a separate AI review. No external researcher has confirmed it. No manuscript or PDF has been created.

## Original assertion fails at the output-site leak

Use a three-compartment cycle with input at 1, output at 2, and a sole leak at 2. Set k1=k21, k2=k32, k3=k13, e2=k02. Then

    x1' = -k1*x1 + k3*x3 + u
    x2' =  k1*x1 - (k2+e2)*x2
    x3' =  k2*x2 - k3*x3
    y   =  x2.

The monic-denominator transfer function is

    N(s) = k1*(s+k3),
    D(s) = (s+k1)*(s+k2+e2)*(s+k3) - k1*k2*k3.

Its nonconstant coefficient map, ordered (C2,C1,C0,B1,B0), is

    C2 = k1+k2+k3+e2,
    C1 = k1*k2+k1*e2+k1*k3+k2*k3+e2*k3,
    C0 = e2*k1*k3,
    B1 = k1,
    B0 = k1*k3.

The four rows (C2,C0,B1,B0) of its Jacobian have determinant

    -k1^2*k3.

Consequently this Jacobian is full column rank at every positive parameter point. In particular at

    (k1,k2,k3,e2)=(2,1,3,1)

all parameters are positive, d2=k2+e2=k1=2, but the selected minor is -12 and the rank is 4. The model has n>=3, a=1 in [p-1] minus Leak, ell=p=2 in Leak intersect[p], and at most one leak; every hypothesis of Conjecture4.6 is met. Its predicted hyperplane therefore is not contained in the singular locus.

At this point N(s)=2s+6, D(s)=s^3+7s^2+16s+6 and gcd(N,D)=1, so the counterexample does not exploit transfer-function cancellation.

A stronger direct check is the rational inverse

    k1 = B1,
    k3 = B0/B1,
    e2 = C0/B0,
    k2 = C2 - B1 - B0/B1 - C0/B0.

It is defined whenever B1*B0 !=0, including every positive k1,k3. This inverse is a property of this specified three-node model; it is not a claim that all cycle models are globally identifiable.

## Proposed corrected interior-leak inclusion

Claim: for a generically locally identifiable cycle with input1/outputp, every a in [p-1] minus Leak and every ell in Leak intersect[p-1] satisfies

    {d_ell=k_a} subset singular locus.

Proof: set K=product_{i=1}^n k_i, kappa=product_{i=1}^{p-1} k_i, and Q(s)=product_{i=p+1}^n(s+d_i). The transfer numerator and denominator are

    N(s)=kappa*Q(s),
    D(s)=product_{i=1}^n(s+d_i)-K.

Consider the parameter tangent vector v with only these potentially nonzero entries:

    delta k_a   =  k_a,
    delta k_ell = -k_ell,
    delta e_ell =  k_ell-k_a.

These coordinates are distinct because a is not a leak and ell is a leak. Since a,ell<p, both k_a and k_ell occur once in kappa and K. The product rule gives delta kappa=delta K=0 identically, as polynomial identities; no division by a parameter is required. Neither a nor ell contributes to Q, hence delta N=0. Also delta d_a=k_a, delta d_ell=-k_a, and all other delta d_i=0. Therefore

    delta D(s)=k_a*(d_ell-d_a)*product_{i notin{a,ell}}(s+d_i),

which vanishes on d_ell=d_a=k_a. Every coefficient has zero directional derivative, so J*v=0 there. On the positive parameter region v is nonzero, since k_a>0, and J fails to have full column rank. The same conclusion extends to the complete real hyperplane: all maximal minors are polynomials vanishing on its dense relatively open subset k_a!=0, and hence vanish identically on the hyperplane. This proves the stated inclusion. Generic local identifiability is used only to give the source's singular-locus interpretation.

The endpoint ell=p is excluded precisely because k_ell then lies outside kappa. The cancellation delta kappa=0 no longer holds, explaining why the source's endpoint extension fails.

## Verification scope and limits

`python3 verify.py` uses exact SymPy integers/polynomials. It derives the three-node transfer function by determinant/adjugate and by direct ODE elimination, computes all maximal minors, composes the rational inverse, and checks the interior-kernel identity in all 28 eligible cases with n=3,4 and one or two leaks. These finite checks support the written general proof; they do not replace it.

The work does not challenge the source's main generic-identifiability classification. It gives a counterexample to one separately identified conjecture and a proof of its restricted interior-leak version. A separate AI review accepted this mathematical scope; priority and original-author confirmation remain pending. There is no basis here to promise journal acceptance or AI media coverage.
