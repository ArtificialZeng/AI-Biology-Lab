# Candidate extension by stopped growth stages: beta > 7/8

Research pass 2 for the existing beta-mixture contribution. This is a new proof candidate, not a referee acceptance or a new paper. The existing beta>10/11 result and the public endpoint results remain credited and unchanged until review.

## Exact theorem proposed

In the exact equation-(7) compartment model and source mixture

psi_beta((e,s),.)=(1-q(s))psi_b((e,s),.)+q(s)kappa_p((e,s),.),
q(s)=c/(1+s)^beta,

fix beta in (7/8,1), c,p in (0,1). From every finite state with positive total substrate, the minimal CTMC has strictly positive probability of finite-time explosion. There is no almost-sure or sharp-threshold assertion. The interval 0<beta<=7/8 is not classified here.

The proof uses empty budding only as a baseline, then couples the exact mixture until its first marked alternative allocation. It never assumes the post-alternative process has just one substrate-bearing compartment.

## 1. Exact baseline rates and credited occupancy lemma

In the empty-budding process, start with one substrate-bearing compartment (the tag). Write e,s for its enzyme and substrate counts, C for total compartments, and B for the number of outside compartments containing enzyme. The tag remains unique. Its enzyme and substrate counts never decrease. As in Yoo v1.1.0, the exact predictable rates are

s -> s+1 at rate es; C -> C+1 at rate s;
C -> C-1 at rate C(C-1);
B -> B+1 at rate C-1-B;
B -> B-1 at rate B(B+1).

Tag-carrier mergers occur at rate 2B and each imports at least one enzyme. Direct tag enzyme immigration has rate 1. These are functions of the complete process; no Markov-projection or independence assertion is needed.

We use the following lower-occupancy estimate, credited to the public preprint's Lemma 3. Suppose an adapted nonnegative nearest-neighbor count Y has upward rate at least lambda while below 2R, and downward rate at most mu=lambda/4 while between 1 and 2R. For any initial Y>=0 and T>=t_fill=8R/lambda, its probability of failing to satisfy Y>=R throughout [t_fill,T], before any actual-process explosion, is at most

5/(2R)+mu T 4^(-R).                                    (1)

This holds under an adapted guard up to guard failure, whose probability is added separately. A short proof is included to make the application inspectable. Couple Y from below with a birth-death chain on {0,...,2R}, initially 0, whose upward/downward attempt clocks have rates lambda,mu. At t_fill their Poisson means are 8R,2R. Chebyshev bounds the probability that the upward count is below 6R or downward count above 4R by 2/R+1/(2R)=5/(2R). Otherwise the chain has reached 2R. From 2R-1 its embedded walk, with probabilities 4/5 and 1/5, reaches R-1 before 2R with probability 3/(4^(R+1)-1)<=4^(-R). The expected number of departures from 2R by T is at most mu T. Strong Markov and a union bound prove (1). The adapted comparison can be built at current-state equality and stopped at guard failure, so this proof does not condition on future guard success. All actual-process couplings may first be made on finite complete-state localizations and continued only until explosion.

## 2. One growth stage with substrate exponent eight

Let x be a power of two, x>=1024. Suppose initially

s=x^8,       e>=x/1024,                                (2)

with any finite configuration of outside compartments. Put S=x^8 and let theta be the first time s reaches (2x)^8=256S. Define

t_*=1/x,        T_x=2^20/x,        e_next=x/512.

**Bounded-substrate nonexplosion.** The baseline cannot explode before reaching a fixed finite substrate level. Indeed, let E be its total enzyme count and V=1+C+E+s. Before s exceeds a fixed M, its complete-state generator satisfies

LV=C+s+es-C(C-1)<=C+M+ME<=(M+1)V.

The finite sublevel sets of V contain only finitely many complete compartment configurations. Stopped Dynkin, Gronwall and Markov's inequality bound the probability of reaching V>=n by time t and before the substrate target by V(0)exp((M+1)t)/n. Letting n increase rules out explosion before that target. This also follows by dominating compartment births by a rate-M Poisson process and then total enzyme immigration by its finite compartment count. Consequently explosion before theta has probability zero in each stage. The countable union of these exceptional events over all stages is null; the stage-time and occupation bounds below are asserted outside that null set.

If the baseline explodes during the stage, call the stage successful. Otherwise we prove that, except with probability at most

p_x=128/x,                                             (3)

theta<=T_x and e(theta)>=e_next. The bounds are uniform over the finite outside configuration and over every initial e satisfying (2).

### Filling total compartments

Set R_C=x^4/4 and lambda_C=x^8. While C<2R_C, its birth intensity s is at least x^8. While 1<=C<=2R_C, its death intensity C(C-1) is at most x^8/4=lambda_C/4. Applying (1), the probability of failure of

C>=R_C throughout [2/x^4,1/x]

before explosion is at most

10/x^4 + (x^7/4)4^(-x^4/4).                            (4)

### Filling outside enzyme carriers

From time 2/x^4 apply the guarded estimate to B, with guard C>=R_C. Put R_B=x^2/16 and lambda_B=x^4/8. When B<2R_B and the guard holds,

C-1-B >= x^4/4-1-x^2/8 >= x^4/8.

For B<=2R_B,

B(B+1) <= (x^2/8)(x^2/8+1) <= x^4/32=lambda_B/4.

These inequalities hold for x>=4. Its filling time is 4/x^2, and

2/x^4+4/x^2 <=1/(2x)                                   (5)

for x>=16. Apart from the failure in (4), the probability of failure of B>=R_B throughout [1/(2x),1/x] is at most

40/x^2 + (x^3/32)4^(-x^2/16).                          (6)

### Delivering enzymes by time t_*

While these occupancy guards hold on [1/(2x),1/x], tag-carrier merger intensity is at least 2R_B=x^2/8. Couple its counting process from below by constant-rate Poisson attempts until guard failure. Over the interval those attempts have mean x/16. Chebyshev gives probability at most 64/x of fewer than x/32 attempts. Since each actual merger imports at least one enzyme and x/32>=x/1024, the tag reaches e_next by t_*, apart from (4), (6), the import failure 64/x, or prior explosion. This argument counts actual mergers; it does not assume independent per-enzyme import delays or bounded enzyme cargo.

### Substrate cannot escape the next level too early while enzymes lag

Let eta be the first time e reaches e_next. If e initially exceeds this threshold, eta=0 and the following bad event is empty. Until eta, substrate birth intensity es is at most (x/512)s. Thus s is dominated up to eta by a Yule process Y with initial population S and per-capita rate x/512. At time t_*=1/x,

E Y <=2S,       Var(Y)<=2S,

because exp(1/512)<2. Hence

P(theta<eta, theta<=t_*, no prior explosion)
 <= P(Y(t_*)>=256S)
 <= 2S/(254S)^2 <1/S.                                  (7)

No conditioning on future nonexplosion is used; the comparison stops at eta or an actual explosion. If theta happens after t_* but e(theta)<e_next, then eta>t_*, already bounded by the enzyme-delivery estimate. Combining these observations bounds failure of e(theta)>=e_next by the sum of (4), (6), 64/x and 1/S.

### The upper time bound

Until theta, substrate birth intensity is at least (x/1024)S. Over T_x the corresponding lower Poisson attempt process has mean 1024S. Reaching the target requires 255S births. Unless actual explosion occurs first,

P(theta>T_x) <= 1024S/(769S)^2 <1/S.                   (8)

This is a lower count comparison up to the target, not an explosive upper comparison process.

### Summable error

For every integer R>=4, 4^R>=R^4: equality holds at 4 and the induction ratio is at least 4(4/5)^4>1. Consequently for dyadic x>=16,

(x^7/4)4^(-x^4/4)<=1/x^2,
(x^3/32)4^(-x^2/16)<=1/x^2.

For the second inequality it suffices that x^3>=2048, and for the first that x^7>=64. Summing (4), (6), (7), (8) and the import failure gives

64/x+42/x^2+10/x^4+2/x^8 <=128/x

for x>=16. This proves the stage statement (3), with all threshold counts integral for x>=1024 dyadic.

## 3. Positive-probability explosion with finite diversion occupation

Choose any dyadic x_0>=1024 and start one tag with s=x_0^8 and e>=x_0/1024, with any finite outside configuration. Set x_n=2^n x_0. Apply the stage estimate at successive substrate hitting times. Prior explosion is treated as permanent success. The strong Markov property and a union bound for the first failed stage give a success event G with

P(G)>=1-sum_n 128/x_n=1-256/x_0>=3/4.                  (9)

On G, either explosion occurs during a finite stage or all substrate levels are reached by total time

sum_n T_{x_n}=2^21/x_0<=2048.

In the latter case infinitely many substrate births occur by this finite accumulation time, so it too is explosion.

More importantly, until the next level is reached within stage n, the substrate is bounded above by (2x_n)^8. Thus for every delta in (0,1/8), on G,

integral_0^zeta s(t)^delta dt
 <= sum_n (2^20/x_n)(2x_n)^(8delta)
 = [2^(20+8delta) x_0^(8delta-1)]/[1-2^(8delta-1)]
 =: K_delta < infinity.                               (10)

If explosion occurs inside a stage, only a prefix of that stage is included and later stages contribute zero. Formula (10) is precisely the additional estimate needed beyond a bare explosion stage proof. It relies on stopping each stage at its substrate target; a lower growth bound alone would not control the integral above.

## 4. Finite compensator and marked-mixture transfer

Set delta=1-beta in (0,1/8). For the baseline fragmentation times t_j, let

W=sum_{t_j<zeta} q(s(t_j-)),
A=integral_0^zeta s(t)q(s(t))dt.

Since s>=1, s q(s)<=c s^delta. Equation (10) gives A<=cK_delta on G.

We need the pathwise finite-compensator implication, rather than incorrectly conditioning a compensator expectation on the future event G. For any finite K, stop the baseline when its accumulated nonnegative intensity integral first reaches K, and also at deterministic time and complete finite-state exits. The weighted stopped fragmentation count has mean equal to that stopped integral, at most K. Passing to increasing localizations by Fatou shows that the weighted sum before that stop is finite almost surely. On {A<K} the intensity stop is never reached before explosion; hence W<infinity almost surely there. Taking K=cK_delta+1 shows W<infinity almost surely on G.

Independently attach uniforms to each baseline fragmentation, marking it with probability q(s(t_j-)). Up to the first marked event, use the same full states, clocks and empty-budding transitions for the mixture. At its first mark allocate using kappa_p and thereafter use the exact mixture generator. Before that mark the common-state fragmentation law is exactly (1-q)psi_b+q kappa_p. Thus this constructs the requested CTMC until divergence, without a claim about its dynamics after a diversion.

Conditionally on a complete baseline path, the probability of no marks before its explosion is

product_{t_j<zeta} [1-q(s(t_j-))].

On G this is strictly positive almost surely: q<=c<1 and

-log(1-q)<=q/(1-c),

so finiteness of W makes the logarithmic sum finite. Since P(G)>0, the expectation of this positive product times 1_G is strictly positive. On that event the mixture follows the same infinite baseline path to the same finite explosion time. Therefore the exact mixture explodes with positive probability.

This transfer works for any admissible replacement allocation in place of kappa_p. That observation is a proof remark only; the theorem above retains the exact source family.

## 5. Every finite positive-substrate initial state

Given an arbitrary finite state with positive total substrate, choose a dyadic x_0>=1024 with x_0^8 at least that total substrate. There is a finite positive-probability reaction path that coagulates all compartments into one tag, immigrates enzymes until e>=x_0/1024, and produces substrate until s=x_0^8. Every desired event has positive rate at its finite current state, and a finite sequence has positive probability even when enormously long. Competing fragmentations can be absent along that specified finite path; no infinite suppression event is required. At the target the tag is unique and the outside configuration is finite (indeed empty), so the construction in Sections 2–4 applies. The strong Markov property proves the stated positive-probability conclusion from the original state.

## Scope and credit

Yoo v1.1.0 already proves empty-budding explosion using lower-occupancy and actual-process stages. This argument credits that method and the exact carrier rates. The proposed delta is the exponent-eight staged cone, stopped substrate upper bound and finite occupation estimate (10), which transfer explosion to beta>7/8 within the formerly source-unclassified range. It enlarges the same contribution previously proved for beta>10/11. It does not newly solve original Open Problem 3.8, prove an almost-sure mixture result, or classify beta<=7/8. Finite code checks below support constants and rate algebra; they do not establish the infinite-state theorem.
