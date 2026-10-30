# Global dynamics of two competitors under density-triggered serial transfer

This is an analytic argument. Finite script checks support its identities but do not prove its quantified claims. Current checking status and source-search limitations are stated in README.md and source.md.

## Original question and model

Schreiber (2025), DOI 10.1007/s00285-025-02283-0, Sections 2.2 and 5.2, equation (26), asks whether the strict arithmetic-mean inequality gives exclusion from every positive initial composition, rather than merely local attraction of the resident-only two-cycle.

Fix 0<delta<1. For phase j in a finite period m>=1, growth rates a_j,b_j are strictly positive constants. The phase begins immediately after dilution with x_1+x_2=delta. Both species grow exponentially, with rates a_j,b_j, until total density first reaches 1. Both are then multiplied by delta and phase advances. The original question has m=2, a_j=rho_1(j), b_j=rho_2(j). The phase cycle is fixed in advance and is independent of the populations; waiting times depend on composition. Let L=log(1/delta)>0.

For post-transfer frequency u=x_1/delta in (0,1), the unique waiting time tau solves

    delta [u exp(a tau)+(1-u) exp(b tau)]=1.                (1)

Existence, uniqueness and smooth dependence follow because the left side is strictly increasing from delta to infinity and its tau derivative is strictly positive. Bounds L/max(a,b)<=tau<=L/min(a,b) hold, and are strict for interior u if a!=b. Put

    v=delta u exp(a tau), 1-v=delta(1-u)exp(b tau),
    z=log[u/(1-u)], T_j(z)=log[v/(1-v)].

Then

    T_j(z)=z+(a-b)tau(u).                                (2)

## Lemma 1: pointwise invasion bounds

For every interior composition,

    L(1-b/a) <= T_j(z)-z <= L(a/b-1).                    (3)

Both inequalities are strict if a!=b; both are equalities (with value zero) if a=b.

Proof. If a>b, (1) is a strict convex combination of exp(a tau) and exp(b tau), so exp(b tau)<exp(L)<exp(a tau). Hence L/a<tau<L/b, and multiplying by a-b>0 gives (3). If a<b, the inequalities are reversed before multiplication by a-b<0, with the same result. If a=b then tau=L/a and (2) is zero. No invasion approximation or small-composition expansion is used.

In particular, along any prescribed sequence of positive growth-rate pairs, even without periodicity, summing (3) gives

    z_0 + L sum_{k<n}(1-b_k/a_k) <= z_n
         <= z_0 + L sum_{k<n}(a_k/b_k-1).                (4)

This envelope is pathwise; it does not by itself assert a limit for sequences whose cumulative bounds do not diverge with the required sign.

## The original strict two-phase question

Suppose (a_0/b_0+a_1/b_1)/2<1, and set

    gamma=L[2-a_0/b_0-a_1/b_1]>0.

At every complete two-transfer period, (4) gives

    z_{2n}<=z_0-n gamma,
    x_{1,2n}/x_{2,2n}<= (x_{1,0}/x_{2,0})exp(-n gamma). (5)

The bound is valid for either starting phase because the period sum is unchanged. Since x_1+x_2=delta after every transfer, it follows that x_{1,2n}->0 and x_{2,2n}->delta. The intervening phase has a bounded log-odds increment by (3), so the same limits hold along every transfer, not only a subsequence. This proves exclusion from every positive point of the source post-transfer state space. It supplies the affirmative mathematical implication in the frozen strict two-phase question. Source comparison and historical-priority limits are separate from this implication.

Between transfers, tau is uniformly bounded above and below by positive constants because the phase set is finite. At elapsed s in a phase, odds multiply by exp[(a_j-b_j)s], a uniformly bounded factor. Thus the species-1 frequency, and hence its absolute density (total density<=1), tend to zero throughout continuous time. The survivor traverses a growth-dilution cycle; no fixed continuous-time survivor-density limit is asserted. Distance to the resident boundary flow-kick cycle tends to zero. No finite-time accumulation of transfer events is possible, by the positive lower waiting-time bound.

## Lemma 2: order and strict contraction in log odds

Each T_j is a strictly increasing C^1 diffeomorphism R->R. If a=b it is the identity. Otherwise 0<T_j'(z)<1 everywhere.

Proof. Write A=exp(a tau), B=exp(b tau), D=uA+(1-u)B=1/delta, and G=a u A+b(1-u)B>0. Implicit differentiation of (1) gives

    tau_z=-u(1-u)(A-B)/G,
    T_j'=1-(a-b)u(1-u)(A-B)/G
         =[a u+b(1-u)]D/G
         =[a u+b(1-u)]/[a v+b(1-v)]>0.                 (6)

If a!=b, (a-b)(A-B)>0 since tau>0, proving T_j'<1. Boundedness of T_j(z)-z from (3) implies the limits at +/-infinity are +/-infinity, hence onto; strict positive derivative supplies the smooth inverse. As z->-infinity, u->0 and tau->L/b; as z->+infinity, u->1 and tau->L/a. Therefore the two displacement limits are exactly the upper and lower expressions in (3).

## Complete finite-period classification (one coherent corollary)

Define

    U=L sum_{j=0}^{m-1}(a_j/b_j-1),
    V=L sum_{j=0}^{m-1}(1-b_j/a_j),
    F=T_{m-1} o ... o T_0, H(z)=F(z)-z.

If all a_j=b_j, every T_j and F are the identity and all compositions are preserved. Otherwise 0<F'<1, so H is continuous and strictly decreasing. Because each finite prefix differs from z by a bounded amount, Lemma 2 gives

    lim_{z->-infinity} H(z)=U,
    lim_{z->+infinity} H(z)=V.                          (7)

The nonneutral condition ensures V<U. Indeed each summand difference is

    (a/b-1)-(1-b/a)=a/b+b/a-2=(a-b)^2/(ab)>=0,

and at least one is positive. A useful stronger fact is that V>=0 and U<=0 cannot both occur unless all ratios equal one: writing r_j=a_j/b_j, Cauchy gives (sum r_j)(sum 1/r_j)>=m^2, with equality only for constant ratios, and the two sum inequalities force equality and ratio one.

All nonneutral cases fall into these three disjoint alternatives:

1. U<=0: H(z)<U<=0 for every finite z. Iterates of F strictly decrease. If they had a finite lower limit ell, continuity would give F(ell)=ell, impossible. Thus z_{mn}->-infinity and species 2 excludes species 1 from every interior composition. When U<0, z_{mn}<=z_0+nU gives the explicit exponential odds envelope. The threshold U=0 is included without claiming a uniform exponential rate.
2. V>=0: H(z)>V>=0, and the symmetric argument gives z_{mn}->+infinity and global exclusion of species 2. If V>0, z_{mn}>=z_0+nV. No uniform exponential claim is made at V=0.
3. V<0<U: by (7) and strict decrease there is exactly one finite fixed point z*. If z<z*, strict increase of F and H(z)>0 imply z<F(z)<z*. If z>z*, z*>F(z)>z. Thus every orbit converges monotonically at period times to z*. Applying continuous T_j gives a unique globally attracting interior m-phase cycle for the fixed starting phase. With m distinct phase labels the augmented state has period m; the population-composition sequence can repeat a shorter block if the phase schedule has repeated equivalent labels. No uniform contraction constant below one on all of R is asserted or needed.

Changing the starting phase cyclically conjugates the period map by a diffeomorphism and leaves U,V unchanged. The conclusions hold for every starting phase. The continuous-time boundary argument above also applies to the exclusion cases. For coexistence, composition and waiting times converge phasewise to the unique resident cycle; elapsed-time convergence to the cycle as a set follows from bounded phase lengths and continuous exponential flows, not from a claim that all solutions synchronize to the same absolute clock.

## Limits and proof status

This uses common dilution, fixed common threshold, exactly two species, strictly positive constant within-phase growth, and a prescribed periodic environment sequence. It does not prove multispecies coexistence, density-dependent growth, stochastic demographic extinction, unequal survival fractions or arbitrary feedback-driven switching. Formula (4) is only a stated envelope for nonperiodic sequences; no general nonperiodic classification is claimed.

The original source's broad invasion framework and its other ecological examples are not changed. Prior invasion/coexistence criteria and the Yi–Dean model are credited. A current primary-source comparison is recorded in source.md; no 2013 supplementary document was established or represented as read. No historical priority, biological experiment, author endorsement, human peer review or proof-assistant certification is asserted.

Separate checking passes and their limitations are recorded in provenance/review-receipt.json. They are AI checks, not external human peer review or proof-assistant certification. No trusted-runner manuscript release is asserted.
