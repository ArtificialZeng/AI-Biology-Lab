# Pointwise fixation deficit for every matching complement

Proposed complete proof of a subsidiary theorem; awaiting fresh mathematical review. The arbitrary connected noncomplete mutant-graph question from the 2023 source remains unresolved.

## Statement and credit

Let N>=3 and 1<=m<=floor(N/2) be integers. Let M be any matching of m edges on N labelled vertices and put H=K_N\M. Residents use K_N and mutants use H in the neutral, discrete-time, type-specific Moran birth–death process: choose a parent uniformly, choose one neighbor uniformly in that parent's type-specific graph, and replace that neighbor by an offspring. If h_v is eventual mutant fixation probability starting with only vertex v mutant, then, writing C=N(N-2)+2m-2,

    h_v <= 1/N - 2m/[N²(N-1)C] < 1/N.

The same inequality holds for the fixation probability averaged uniformly over the initial singleton. The m=1 case is already covered by the exact formula of Tkadlec et al., Theorem A.1 (2023); their N=6 numerical enumeration covers all matching complements at that size. The m=2 case and its identical stated gap were proved in R14 and published in lab PR18. The new proposed contribution is the entire unbounded matching family m>=3, including perfect matchings, with a pointwise quantitative bound. This is an extension of PR18, not a separate micro-result. The bound is not claimed sharp.

## Exact state reduction

For a configuration let a be the number of missing-edge pairs occupied by two mutants, b the number occupied by one mutant, and d=m-a-b the number occupied by no mutants. Let c and r be the numbers of mutants and residents at vertices outside the matching. Then

    a,b,d,c,r >=0 are integers,
    m=a+b+d,  N=2m+c+r,  k=2a+b+c.

The pair states and the unmatched vertices can be permuted by automorphisms of H and K_N. Define

    T=(2a+b)/(N-2)+c/(N-1),  L=(N-k)/(N-1),  Q=N(P-I).

The nonzero off-diagonal rates in Q, expressed in (a,b,c), are

| target state | rate |
|---|---|
| (a,b+1,c) | 2d T |
| (a+1,b-1,c) | b[T-1/(N-2)] |
| (a-1,b+1,c) | 2a L |
| (a,b-1,c) | b L |
| (a,b,c+1) | r T |
| (a,b,c-1) | c L |

Every mutant at a matched vertex has degree N-2; every unmatched mutant has degree N-1. Each resident target in an empty pair sees all mutant parent weights, totaling T. A resident target in a mixed pair is inaccessible to exactly its mutant partner, removing 1/(N-2). Unmatched resident targets see weight T. Every mutant target sees resident-parent weight L because residents use K_N. Target multiplicities are respectively 2d,b,2a,b,r,c. These facts prove the rate table and exact lumpability for every parameter in the statement. Terms with zero multiplicity at the boundary contribute zero; when b>0, T>=1/(N-2), so the mixed-pair gain rate is nonnegative. Denominators are positive for N>=3.

## Uniform superharmonic certificate

Set

    W=[N-2+(m-2)/N]k+a+m k(N-k)/[N(N-1)],
    D=N(N-2)(N-1)² > 0.

We prove QW<=0. For a nonnegative integer x, write x_j=x(x-1)...(x-j+1) for its falling factorial, with j=2,3,4. Each x_j is nonnegative: it vanishes if x<j and is a product of nonnegative integers otherwise.

Substitute the six transition rates into QW=sum_y q(x,y)[W(y)-W(x)], and then substitute N=2(a+b+d)+c+r and m=a+b+d. Direct expansion gives the exact polynomial identity

    -D QW = F00 + c F10 + r F01 + c² F20 + cr F11 + r² F02
             + 2(a+b+d)(c²r+cr²),

where the following expressions use the falling factorial notation just defined:

    F20 = 4ad + b_2 + 6bd + 4d²
    F02 = 4a + b² + 2bd + b
    F11 = 2[4a²+7ab+8ad+2b²+2b_2+9bd+2d²+2d_2]

    F10 = 16a²d+4a b_2+36abd+24ad²
          +4b_3+5b_2+24b_2 d+28b d_2+38bd+8d_3+16d_2

    F01 = 8a_3+28a_2+8a+12a_2 b+8a b_2+24ab
          +8a_2 d+20ad+12abd+4b_3+9b_2+2b
          +12b_2 d+8b d_2+18bd

    F00 = d(16a_3+72a_2+40a)
          +bd(32a_2+84a+24)+d_2(16a_2+24a)
          +b_2(2+48d)+b_3(6+2a+14d)+2b_4
          +28a b_2 d+b d_2(36+24a)+20b_2 d_2+8b d_3.

Every displayed expression is a sum of products of nonnegative numbers for all admissible integer a,b,d. Since c,r>=0, the identity proves -D QW>=0 and hence QW<=0 for the entire infinite parameter domain. This is an exact algebraic certificate, not an extrapolation from finite enumeration.

The identity is reproducible by expansion alone. `certificate.py` recomputes the generator identity and all eight c,r coefficient identities using exact symbolic arithmetic. It also verifies a second sign representation: all coefficients are nonnegative in all 12 exhaustive substitutions a=0 or 1+A, b=0 or 1 or 2+B, and d=0 or 1+E, where A,B,E,c,r>=0. `evidence/uniform-certificate.json` stores the expanded polynomial and each finite coefficient certificate; these checks supplement the displayed self-contained identity.

## Normalization, strictness, and absorption

At extinction W=0. At fixation a=m,b=0,c=N-2m,k=N, so W=C=N(N-2)+2m-2>0. For either singleton orbit (a,b,c)=(0,1,0) or (0,0,1), W=C/N. The second orbit exists only if N>2m. Direct substitution in the same generator formula gives, for every existing singleton,

    QW=-2m/[N(N-1)],
    (P-I)W=-2m/[N²(N-1)].

These exact identities are checked symbolically. N=2m has only matched singleton vertices and is included.

The absorption time tau is almost surely finite. Indeed, from any nonabsorbing state choose one resident and repeatedly choose that same resident parent to overwrite any specified remaining mutant. Each prescribed step has probability 1/[N(N-1)] because the resident graph is complete. Extinction is reached within N-1 steps. Thus, from every transient state, the probability of absorption within N-1 steps is at least epsilon=[N(N-1)]^{-(N-1)}>0. Iterating blocks gives a geometric upper bound on the survival probability and proves almost-sure absorption.

For fixed N,m the state space is finite. The stopped process W(X_(t wedge tau)) is a bounded supermartingale by QW<=0. Dominated convergence therefore gives

    C h(x)=E_x W(X_tau)<=W(x)

for every initial configuration x. At a singleton x, apply the Markov property for h and the preceding comparison once:

    h(x)=E_x h(X_1)
        <=E_x W(X_1)/C
        =[W(x)+(P-I)W(x)]/C
        =1/N-2m/[N²(N-1)C].

All denominator factors are positive and m>=1, proving strict inequality. Averaging over singleton locations gives the corresponding rho bound.

## Failure diagnosis and computational limits

The R14 naive all-m potential omits the new linear correction (m-2)k/N. It has QV=1/100 at N=6,m=3,(a,b,c)=(1,2,0), and at (a,b,c)=(m-1,1,N-2m) its drift equals (m-2)(N-2m)/[N(N-2)(N-1)], positive whenever m>=3,N>2m. That invalidates only that potential. The new correction cancels this entire latter obstruction and retains exact singleton/fixation calibration; at the former N6 state its corrected drift is -1/150.

`labelled_check.py` separately constructs every parent-target event in the full labelled process, without importing quotient rates or the certificate generator. It exactly compares 6152 configurations and 575480 parent-target events at N=3..11, including m=1,2,3,4,5 and perfect matchings, to the symbolic polynomial. It confirms every singleton drift and nonpositive drift throughout those finite checks. This is supplemental implementation verification; uniform validity comes from the general rate derivation and displayed sign identity.

The all-matching subsidiary question is resolved by the proposed argument, subject to fresh review. The source's arbitrary connected noncomplete G question is still unresolved. No exact fixation formula, optimal gap, unequal fitness, death–birth update, weighted graph, arbitrary missing-edge set, or current universal priority claim is asserted. No external human review, public action, or formal-proof-assistant certification has occurred in this pass.
