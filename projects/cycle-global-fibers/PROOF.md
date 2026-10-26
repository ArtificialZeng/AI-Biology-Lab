# Theorem: exact fibers of two-leak directed cycles

This exact mathematical scope passed a separate AI mathematical review on 26 October 2026; external human review remains pending. Its
novelty is unverified beyond the bounded primary-source comparison in SOURCE.md.

## Statement

Fix n>=3, 2<=p<=n-1 and 1<=r<p<=t<=n. In the directed cycle
1->2->...->n->1 let k_i denote the rate on the edge leaving i. Only r and t
have leaks e_r,e_t. There is one input at 1 and one output at p. Put
d_i=k_i+e_i at the two leaking vertices and d_i=k_i otherwise. Set

    kappa = product_{i<p} k_i,
    C = product_{i=1}^n k_i,
    T(s) = product_{i>p}(s+d_i),
    U(s) = product_{i<=p}(s+d_i),
    D(s) = U(s)T(s)-C,  N(s) = kappa T(s).

These are precisely the denominator and numerator of the input-output equation
of Ahmed et al., Lemma 3.10. The labels of compartments and parameters are fixed.
The fiber below refers to the coefficient map of that equation over the affine complex extension C^(n+2), where either leak coordinate may vanish. On the nonzero-parameter torus a special fiber can omit zero-leak branches, while its generic degree is the same.

Suppose the observed data belong to the image, kappa*C != 0, and U and T have
only nonzero simple roots. Write their factorizations as

    U(s)=product_{j=1}^p(s+alpha_j),
    T(s)=product_{j=1}^q(s+beta_j), q=n-p.

Then:

1. The entire complex coefficient fiber consists of exactly p!q! distinct
   parameter vectors. It is given by the inverse branches below. Consequently
   the generic algebraic identifiability degree is p!(n-p)!.
2. If the data arise from strictly positive edge and leak parameters, all these
   complex preimages are real and all their edge rates are positive. Set

       A=U(0), B=T(0), a=C/(kappa B), b=A/kappa,
       m=#{j: a<alpha_j<b}.

   Exactly m(p-1)!q! preimages have both leak rates strictly positive. Thus the
   possible positive fiber cardinalities are precisely

       (p-1)!q!, 2(p-1)!q!, ..., p!q!.

   Every one occurs on a nonempty open subset of the strictly positive
   parameter space satisfying the simple-root assumptions. This is a regional
   real count, not a replacement for the generic complex degree.

## Recovery of invariants from data

The numerator has degree q and leading coefficient kappa != 0, so the data
recover kappa and the monic T=N/kappa. Since q>=1, polynomial division is
available. On image data, dividing D=UT-C by T leaves the constant -C. Thus C,
then U=(D+C)/T, and then A and B are uniquely recovered. No knowledge of a
particular parameter preimage is used in this step.

There is also no transfer-function cancellation on this domain: at any root
z of T, D(z)=-C != 0. Hence gcd(D,N)=1. A monic D and N are uniquely determined
by the rational transfer function. The coefficient-fiber statement therefore
also applies to transfer-function fibers on this domain; it does not make a
claim about special unknown-initial-state degeneracies outside the standard
coefficient-map setting.

## Complete inverse branches

Choose an arbitrary bijective assignment of alpha_1,...,alpha_p to
d_1,...,d_p and independently of beta_1,...,beta_q to d_{p+1},...,d_n.
There are p!q! choices. At nonleak vertices set k_i=d_i. At the two leak
vertices set

    k_r = kappa / product_{i<p,i!=r} d_i,
    k_t = C / (kappa * product_{i>=p,i!=t} d_i),
    e_r = d_r-k_r, e_t = d_t-k_t.

Every division is legitimate because all assigned roots and kappa are nonzero.
These formulas force product_{i<p}k_i=kappa and then product_i k_i=C.
They reproduce the assigned d_i by construction, so they reproduce D and N.
Thus every assignment gives a parameter vector in the fiber.

Conversely, any parameter vector in the same coefficient fiber has the same
kappa,T,C,U by the preceding recovery. Its labeled d_i therefore give exactly
one of the assignments. At every nonleak vertex k_i=d_i, and the two product
identities uniquely force k_r and k_t to be the displayed values. There are
no other preimages, including ones with a zero edge rate, because C != 0.
Distinct assignments yield distinct vectors: d_i is recovered from a vector
as k_i at nonleaks and k_i+e_i at leaks. This proves both completeness and the
p!q! count.

The hypotheses define a nonempty Zariski-open condition on the image. For
example, choose any distinct nonzero d_i and nonzero kappa,C and use the inverse
formulas to construct a preimage. Locally near simple roots, each labeled root
varies analytically with the coefficients, and each inverse branch above is
analytic. Thus these generic branches are distinct and unramified; the count
also equals the generic degree counted with multiplicity.

## Positivity

For positive source parameters, all d_i, kappa and C are positive. Thus all
alpha_j and beta_j are positive real numbers. Every inverse branch has positive
edge rates because its formula uses products and quotients of positive numbers.
Let alpha=d_p be the root assigned to the output. Using A=product_{i<=p}d_i
and B=product_{i>p}d_i gives

    k_r/d_r = kappa*alpha/A = alpha/b,
    k_t/d_t = C/(kappa*alpha*B) = a/alpha.

These formulas hold also when t=p: then k_p=C/(kappa B)=a, so its leak is
positive exactly when alpha>a. Since d_r,d_t are positive, the inequalities
e_r>0 and e_t>0 are equivalent to a<alpha<b. Positivity is therefore independent
of the assignment of the remaining p-1 first-block roots and all q tail roots.
For each eligible output root there are exactly (p-1)!q! such assignments.
The original positive parameter vector supplies one eligible root, so m>=1;
trivially m<=p. Roots equal to an interval endpoint yield a zero leak and are
excluded by the strict inequalities. No assumption that such equalities never
occur is needed for the formula itself.

## Realizing every positive count on open parameter regions

Fix any allowed n,p,r,t and any m in {1,...,p}. Choose first-block diagonal
values 2,3,...,p+1 and tail diagonal values p+2,...,n+1. Assign d_p=2 and assign
the other first-block values in any order. Set

    a=3/2, b=m+3/2,
    A=2*3*...*(p+1), B=(p+2)*...*(n+1),
    kappa=A/b, C=kappa*B*a.

Apply the inverse formulas to that assignment. The two edge/diagonal ratios are
2/b<1 and a/2=3/4<1. All rates and both leaks are therefore strictly positive.
The interval (a,b) contains exactly the first m values among 2,...,p+1. All
diagonal roots are nonzero and simple, and no root equals an endpoint. Root
separation, strict positivity, and the interval inequalities persist under
small perturbations of the parameter vector. Thus the same count occurs on a
nonempty open positive parameter region, proving attainability without merely
testing finitely many sizes.

## Scope and payoff

By Ahmed et al.'s existing local-identifiability classification, rotating the
single input to 1 puts every locally identifiable single-input/single-output
two-leak cycle in this family. If the output is also 1, two leaks cannot
interlace with distinct input/output positions. If p=n, interlacing would force
one leak at n, giving their exceptional unidentifiable family. For 2<=p<n,
interlacing is equivalent to having one leak in [1,p-1] and one in [p,n].
This last reduction is credit to the prior local theorem, not a new local
classification. The proposed new content is the complete inverse and exact
complex/positive fiber count.

For n=4,p=3, every locally identifiable leak placement has complex degree 6;
its positive fiber count can be 2,4,or6. For n=3,p=2, the complex degree is 2
but a positive-data region can have exactly one positive preimage. Thus
positivity can remove an algebraic ambiguity, and the size of that removal
depends on the observed data through an explicit root interval.

Known gap list in the algebraic argument: none identified by the proposer.
A separate AI mathematical review accepted this exact scope. Priority and original-author confirmation remain pending; no journal publication or external endorsement is claimed.
