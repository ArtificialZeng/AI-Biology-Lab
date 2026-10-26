# Proof candidate: affine-rate mass-action isola

A fresh substantive AI review accepted the mathematics at the scope below. This is not external human peer review or formal release-gate certification.

Take lambda in (0,infinity) and species X,Y. Deterministic mass action uses:

| Reaction | Constant | Propensity |
|---|---:|---|
| X -> emptyset | 7 | 7x |
| 2X -> 3X | 4 | 4x^2 |
| 3X -> 2X | 1 | x^3 |
| X+Y -> 2X+Y | 4 | 4xy |
| X+2Y -> emptyset | 1 | xy^2 |
| X+2Y -> X+4Y | 1 | xy^2 |
| emptyset -> Y | lambda | lambda |
| Y -> emptyset | 1 | y |

The paired X+2Y reactions cancel their net Y contribution and produce the net X term -xy^2. Thus f=(x[-7+4x-x^2+4y-y^2],lambda-y). All constants are positive, and the only varying constant enters affinely.

**Chemicality.** Every propensity is a smooth nonnegative monomial, monotone in its reactants, with the required reactant factor for every consumed species. Even the strict literal net-consumption support reading of source Definition 2 is satisfied: in every column with a negative net component, all rate dependencies are among the net-consumed species. The X+2Y loss consumes both X,Y; its compensating column has no negative components. At x=0, x'=0; at y=0, y'=lambda>0. Full stoichiometric rank follows from the X and Y outflows, so there are no compatibility-class constraints.

**Siphons and boundary.** The Y inflow prevents every set containing Y from being a siphon. The set {X} is a siphon because every reaction producing X has X as a reactant. Therefore the full nonempty siphon list is exactly {{X}}; adjoining the empty set gives a finite two-element lattice. No boundary equilibrium has y=0, because y'=lambda>0 there. The sole boundary equilibrium branch is E0(lambda)=(0,lambda).

**Complete equilibrium classification.** Equilibrium requires y=lambda. For x>0 the remaining condition is

1-(x-2)^2-(lambda-2)^2=0.

Conversely every such point is an equilibrium. Consequently the entire strictly positive equilibrium set is the image of the circle under theta -> (2+cos(theta),2+sin(theta),2+sin(theta)), with coordinates (x,y,lambda). It is a connected compact circle, all coordinates in [1,3]. Its closure is itself and cannot meet a state boundary. This is the entire positive set, not merely a selected arc or numerically traced segment.

**Boundary stability.** At E0(lambda), J=diag(-7+4lambda-lambda^2,-1)=diag(-3-(lambda-2)^2,-1). Both eigenvalues are strictly negative for every lambda>0. Its scalar transversal block is Metzler. Thus there are no boundary losses of hyperbolicity; the universally conditional hypothesis on any such loss is satisfied. A regular splitting is F=4lambda, V=7+lambda^2. Its invasion number is R=4lambda/(7+lambda^2)<1 because V-F=3+(lambda-2)^2. No R=1 boundary relay edge exists.

**Isolation and ordinary folds.** At fixed lambda there are exactly two positive equilibria x=2±sqrt(1-(lambda-2)^2) when 1<lambda<3, one double equilibrium at lambda=1 and 3, and none otherwise. All are isolated. At positive equilibria the eigenvalues are 2x(2-x) and -1. The only nonhyperbolic positive points are (2,1,1),(2,3,3), with Jacobians [[0,4],[0,-1]], [[0,-4],[0,-1]]. For both use right nullvector p=(1,0); normalized left nullvectors are q=(1,4),(1,-4). The parameter coefficients q^T f_lambda are 4,-4, and q^T D_z^2f[p,p]=-4 in both cases. Both are nonzero, so these are ordinary saddle-nodes. The upper-x equilibria are locally asymptotically stable and the lower-x equilibria are saddles.

**Boundedness.** At x=4, x'=-12-4(y-2)^2<0; y(t)=lambda+(y(0)-lambda)e^{-t}. Every trajectory remains bounded and exists for all forward times. For lambda in [1/2,7/2], the common compact rectangle [0,4]^2 is positively invariant and contains the full isola in its interior. This supplies compactness near the isola even though the source conjecture does not demand a common compact domain for all unbounded parameter values.

**Contradiction.** The positive equilibrium set has exactly one connected component C; its closure is disjoint from every boundary equilibrium. There are no boundary relay edges at any admitted parameter. Therefore neither the full chain conclusion nor the weaker connectivity conclusion holds.

Scope: this is a single-parameter family. No assertion is made about connectedness after independently varying all rate constants. Isolas are classical; the contribution is the exact application refuting this particular universal 2026 conjecture. The local relay theorem and Feng–Gavish worked family are unaffected.

The accompanying exact script reconstructs the field from integer reactant/product vectors, exhaustively enumerates all siphons on two species, and checks the polynomial identities and fold coefficients. It is computational confirmation of the displayed argument, not external peer review or priority evidence.

## Explicit path in the full kinetic-rate space

Let the X-loss constant a vary while the other fixed rates remain unchanged and lambda=2. Set

x=t, y=2, a=4+4t-t^2, 0<=t<=3.

The first equilibrium equation becomes exactly -a+4t-t^2+4=0, while the second is zero. All kinetic constants remain strictly positive: a=4+t(4-t) lies in [4,8]. At t=3, a=7 and the path meets the original circle. As t decreases to zero it reaches the boundary equilibrium (0,2) at a=4, which is an admitted positive kinetic value. There the transversal eigenvalue is 4-a=0. The scalar expansion x[(4-a)+4x-x^2] has nonzero x-parameter derivative -1 and quadratic coefficient 8, giving the ordinary transcritical crossing. Thus this exact network fails as a counterexample if all reaction constants vary independently. This is an explicit limitation and possible stronger reformulation of the source conjecture, not a result about the truth of that reformulation for arbitrary networks.
