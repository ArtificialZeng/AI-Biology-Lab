# Four particles permit reversal with a fixed order of killing rates

## Exact scope

Use the complete three-site Fleming–Viot motion graph with q(i,j)=1 for i≠j, uniform copying among the other N−1 slots, positive killing rates and the source orientation Delta_N=F02−F20. Same-site copying is included. Among N≥2, the smallest population for which two killing vectors in the same strict ordering chamber have opposite nonzero current is N=4.

More precisely, for exactly four particles on the entire path k=(1,t,t+1), t>1, there is a unique threshold rho such that Delta4 is negative for1<t<rho, zero att=rho, and positive fort>rho. It is the unique positive root of the degree19 polynomial P below, and

    4769/446 < rho < 3133/293.

This is approximately10.69283. The rational interval is the certificate; the decimal is explanatory only. Ratios change along this path: only the order k0<k1<k2 is held fixed. No fixed-ratio scale reversal or complete all-positive N4 classification is asserted.

## Finite exact identities identify the physical current

The literal general-N count and founder operators are source Equations3.2 and3.6, with founder copying coefficient2/(N−1), donor exclusion, and anchorh(0,(N−1,0,0))=1/N. For N4 the two systems have15 and30 states. four_particle_model.py constructs them exactly; common-denominator vectors p/dp and h/dh are stored with their state ordering.

The primary read-only verifier verify_four_particle.py checks all15 stationary residuals, all30 founder residuals, the anchor, sum(p)=dp and all15 founder normalizations. After t=1+u, every coefficient of dp is positive and every coefficient of dh is negative. Neither denominator vanishes foru>0. The source's finite-chain stationary/Dirichlet uniqueness argument therefore identifies these vectors as the actual stationary and founder probabilities; they are not interpolated fits. The separate independent checker verifies all 96 mixed two-colour Dirichlet residuals directly from victim/donor events, including same-site cross-colour copies.

Forming the source current from these vectors and checking the polynomial cross-product gives

    Delta4(1,t,t+1) = 2t(t−1)(3t²+13t+8) P(t) / D(t).

D is the degree27 polynomial in evidence/four-particle/symbolic-current.json; all its coefficients are strictly positive, so D(t)>0 fort>1. The factorization is exact for the whole path. The remaining polynomial P is:

| degree | coefficient |
| --- | --- |
| 19 | 17798400 |
| 18 | 1536534720 |
| 17 | 60567213600 |
| 16 | 1420684660560 |
| 15 | 21344899045872 |
| 14 | 195208427960028 |
| 13 | 584253437592398 |
| 12 | -12565959488966207 |
| 11 | -246248681426819416 |
| 10 | -2533044309774322902 |
| 9 | -18238231413999297101 |
| 8 | -98691048418285812541 |
| 7 | -410901931250671582674 |
| 6 | -1320965946402702265528 |
| 5 | -3247192332959739375159 |
| 4 | -5965979589385014298578 |
| 3 | -7868434215092961460880 |
| 2 | -6944980494592538947760 |
| 1 | -3588761435808030789504 |
| 0 | -782299069753995118080 |

P has seven positive leading coefficients followed by thirteen negative coefficients, with no omitted coefficients. Descartes' rule of signs gives at most one positive root, counting multiplicity. Its constant term is negative and its leading coefficient positive, so the intermediate value theorem gives a positive root. It is thus unique and simple. Exact rational evaluation gives P(1)<0 and

    P(4769/446)<0<P(3133/293).

Hence the unique root rho lies in the displayed interval, above1. The other numerator factors and D are positive fort>1. This proves the complete sign classification along the unbounded path.

For concrete independent controls, exact substitution gives

    Delta4(1,2,3) = -7553807516755111581622 / 1096028058625091839860416349,

    Delta4(1,100,101) = 458226106942025952304289497342294745709790415205
                      / 4105486436463376678337301093835196669201815807466961979.

Both vectors have the same strict ordering. They supply a reversal even without the stronger path classification. Finite controls alone would not prove the path theorem.

## Why four is the minimum in this fixed motion model

The source's two-particle balance theorem makes every current zero for N2. The companion N3 classification proves sign(Delta3)=sign[−(k0−k1)(k0−k2)(k1−k2)] over the whole positive domain. It has a fixed nonzero sign on each strict ordering chamber, so reversal inside one chamber is impossible for N3. The two four-particle witnesses above prove existence at N4. These three facts establish the claimed minimum, which is restricted to this fixed unit-rate triangle.

## A local zero surface away from equal-rate planes

The source proves the current is analytic in positive rates. At t=rho the displayed path current has nonzero derivative: its other factors are positive and P has a simple crossing from negative to positive. Thus at k=(1,rho,rho+1), the full gradient of Delta4 has nonzero dot product with(0,1,1). The implicit function theorem gives a local smooth zero surface inside the strict ordering chamber, with both current signs in adjacent open regions. This is a local existence statement, not a classification of the entire N4 zero set.

## Prior work and limits

[Lin–Lu2610.12589v1](https://arxiv.org/html/2610.12589v1) supplies the model, founder representation, rationality, N2 balance, generic circulation and local N4 cubic; these are prior results. The N3 global classification is the companion lab result, not their generic theorem. The nearby [Brownian paper](https://arxiv.org/html/2610.18802v1), AppendixA.3, already shows that neutral copying rates can affect ancestry direction. No first general direction-reversal mechanism is claimed here. The added content is the specified discrete within-chamber minimum and exact one-path threshold.

The source's advertised supplementary Python/Lean files were not found in the inspected archive; historical priority and external human verification remain unconfirmed. The full N4 sign/zero classification, other motion graphs and all larger populations remain outside this result. No infinite-space nonexplosion result or biological experiment is claimed.
