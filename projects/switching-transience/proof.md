# A shared-species branching skeleton resolves source Conjecture 46

## Exact statement and scope

Let Z(t)=(X(t),E(t)) be exactly the mass-action process of Example44 in
Cappelletti, Howells and Xu, arXiv:2507.09421v2, printed pp.29–30.
The complete reactions, unit inflows and independent symmetric environment
switching are recorded in model.json and source.md. For every alpha>0 and
every 0<kappa<=1, the chain is nonexplosive, irreducible and transient.
Indeed, from every finite initial state, X1(t)+X2(t)+X3(t) tends to infinity
almost surely as t tends to infinity.

This implies source Conjecture46, which asks only for some alpha>0 and all
sufficiently small kappa. It does not settle the general common-support
Conjecture43. The bound kappa<=1 is sufficient and is not asserted sharp.
The proof uses no simulations, stationary truncation or annealed finite-mean
branching assumption. The known source mean matrices and branching-process
interpretation are credited; the new step is the pruned shared-species skeleton.

## 1. Construction and nonexplosion of the original process

All non-inflow reactions have one reactant. Conditional on the environment
path, give every molecule a label and independent rate clocks: reactions
replace that molecule by its product molecules with fresh labels. Inflow
molecules start new independent families. This is the standard individual
construction of a linear mass-action branching process with immigration.
Families do not interact; they share only the independent environment path.
It realizes exactly the propensities in model.json.

Write N=X1+X2+X3. In environment1,

    L N = 4 X1 - alpha X3 + 1,

and in environment2,

    L N = 4 X3 - alpha X1 + 1.

The only positive population jumps are +4 at rate at most N and +1 at
rate1. Thus N is dominated by a nonexplosive pure-birth process having
these two jumps/rates. For example, dominate both by +4 jumps at rate N+1;
its successive holding means are 1/(N0+4j+1), whose sum diverges.
The environmental clock has the fixed finite rate kappa. Degradations and
conversions have total intensity bounded by a constant times N, so on each
bounded time interval the finite upper population bound also bounds their
integrated intensity. Consequently, no reaction or environmental explosion
can occur. This also justifies the individual-family construction for all time.

## 2. Kill-at-switch skeleton on S2

Start with one labelled S2 molecule in either environment and ignore all
immigrants and all other initial molecules. Let tau0=0, tau1,tau2,... be the
successive environment switching times, with IID holding times
Ti=tau(i+1)-taui distributed Exp(kappa), independent of all family clocks.
At every switching time delete all tracked molecules except those of type S2.
Deleted labels and their future descendants are never tracked again.

In environment1 the tracked active pair (A,B)=(S1,S2) evolves by

    A -> 4 A + B at rate1 per A,
    B -> A at rate1 per B.

In environment2 exactly the same active pair process is obtained by putting
(A,B)=(S3,S2). The other tracked type is absent because of the preceding
deletion. In particular, alpha acts only on untracked inactive molecules.
Within each holding interval the total number A+B of tracked molecules
never decreases: its jumps are +4 and 0. At the switch itself A is pruned.

Let Yn be the number of tracked S2 molecules just after taun, with Y0=1.
Conditional on the whole sequence T=(T0,T1,...), different S2 ancestors in
generation n have independent descendant families. Each contributes an IID
number of S2 molecules at the next switch with the same law as B(Tn)
from initial (A,B)=(0,1). The two environments have identical offspring
laws under this relabeling. Thus Yn is a one-type branching process in an
IID environment. This is a conditional assertion: offspring of different
ancestors are generally dependent after averaging over their shared holding
time. The proof below never treats them as unconditionally independent.
The pruned process is pathwise a labelled subpopulation of the full process.

For deterministic t>0, write m(t)=E[B(t)] and v(t)=Var(B(t)). The active
pair is nonexplosive and all its moments used below are finite. Its mean
matrix, with column population vectors, is

    A0 = [[3,1],[1,-1]].

Put r=sqrt(5), lambda_plus=1+r, lambda_minus=1-r,
c_plus=(1-2/r)/2 and c_minus=(1+2/r)/2. Directly solving the two-dimensional
mean equation from (0,1) gives

    m(t)=c_plus exp(lambda_plus t)+c_minus exp(lambda_minus t).

Both coefficients are strictly positive. The event of no conversion of the
initial B has probability exp(-t), so

    exp(-t) <= m(t) <= exp(4t).

For the upper bound, the active family total is dominated by the pure-birth
process U with U(0)=1 and jumps +4 at rate U. Its first two moments are

    E U(t) = exp(4t),
    E U(t)^2 = 5 exp(8t)-4 exp(4t).

Indeed their generator equations are a'=4a and b'=8b+16a. The same
pathwise domination implies

    v(t)/m(t)^2 <= E B(t)^2/m(t)^2 <= 5 exp(10t).             (1)

These rough bounds are enough; no infinite-state numerical approximation
or finite unconditional expectation of m(T) is required.

## 3. Strict positive logarithmic mean for kappa<=1

The exact exponential formula yields

    log m(t) >= log(c_plus)+lambda_plus t.

We have lambda_plus>3 and c_plus>1/20. The latter follows from
sqrt(5)>20/9, equivalently 405>400. Also log20<3, since the first nine
terms of the positive exponential series give

    exp(3) > sum(j=0..8) 3^j/j! = 89641/4480 > 20.

Consequently, for T~Exp(kappa) and 0<kappa<=1,

    mu := E log m(T)
       >= log(c_plus)+lambda_plus/kappa
       > -log20+3
       > 0.                                                   (2)

The bounds exp(-t)<=m(t)<=exp(4t) give |log m(t)|<=4t, hence
E|log m(T)|<infinity. Every conditional offspring mean is finite, although
the unconditional mean E[m(T)] is in fact infinite on this parameter range.
The latter fact is harmless because only quenched finite moments and an
integrable logarithmic mean will be used.

## 4. Self-contained quenched survival lemma

For completeness we prove the needed survival implication directly, rather
than invoking a branching-in-random-environment theorem as a black box.
The classical Smith–Wilkinson1969 criterion provides historical context,
but is not a dependency of this argument. Standard SLLN, Borel–Cantelli
and Cauchy–Schwarz are the only probability tools in this step.

Fix a holding-time realization. Write

    mn=m(Tn), vn=v(Tn), un=vn/mn^2,
    M0=1, Mn=product(i=0..n-1) mi.

Conditional branching independence gives

    E[Yn | T] = Mn,
    E[Y(n+1)^2 | T] = mn^2 E[Yn^2 | T]+vn E[Yn | T].

Since Y0=1, division and iteration yield the exact formula

    E[Yn^2 | T]/Mn^2 = 1+sum(i=0..n-1) ui/Mi.                 (3)

By the SLLN applied to the IID integrable log m(Ti),
log Mn/n tends to mu>0 almost surely. Also Tn/n tends to0 almost surely:
for every epsilon>0,

    sum(n>=1) P(Tn>epsilon n)
      =sum(n>=1) exp(-kappa epsilon n)<infinity,

and Borel–Cantelli, first on rational epsilon, gives the assertion.
By(1), un<=5 exp(10Tn). Therefore the series

    C(T)=1+sum(i>=0) ui/Mi

is finite almost surely. Explicitly, eventually log Mi>mu i/2 and
10Ti<mu i/4, so ui/Mi<=5 exp(-mu i/4).

For every n, conditional Cauchy–Schwarz now gives

    P(Yn>0 | T) >= (E[Yn | T])^2/E[Yn^2 | T] >= 1/C(T)>0.

Zero is absorbing for Yn, so the events {Yn>0} decrease with n. Their
intersection is the survival event, whose conditional probability is at
least 1/C(T)>0 for almost every holding-time realization. Averaging proves
positive unconditional survival from Y0=1.

## 5. Survival forces Yn to tend to infinity

Let p0(t)=P(B(t)=0 | (A(0),B(0))=(0,1)). The path where the initial B
converts to A and that A has no birth before t has probability t exp(-t).
Hence

    p0(t)>=t exp(-t)>0 for every t>0.

For an integer K>=1 and 1<=Yn<=K, the conditional-on-past probability
of extinction at the next generation is

    E[p0(T)^Yn] >= E[p0(T)^K]
                  >= E[T^K exp(-KT)]
                  = kappa K!/(kappa+K)^(K+1) =: delta_K>0.     (4)

Future holding times and fresh family clocks are independent of the past.
At successive visits to {1,...,K}, the probability of avoiding the offered
next-step extinction is at most (1-delta_K) per visit. By induction over
these stopping times, survival together with infinitely many such visits
has probability0. Thus, on survival, each finite positive set is visited
only finitely often. Intersecting over integer K shows Yn tends to infinity
almost surely on survival. This elementary step avoids any additional
branching growth theorem or martingale-limit integrability assumption.

## 6. Continuous-time escape and transience of the source chain

For every t in [taun,tau(n+1)), the tracked active family total is at least
Yn, because its only population changes before the next pruning are +4
or0. Consequently the full original population satisfies

    X1(t)+X2(t)+X3(t) >= Yn throughout [taun,tau(n+1)).          (5)

The switching times tend to infinity almost surely. By positive survival
and the preceding section, the original process started from a single S2
has positive probability of tending to infinity in continuous time.
Equation(5), rather than a bound only at switching times, excludes
arbitrarily frequent finite-population dips between switches on that event.

Here is an explicit check of the source's irreducibility assertion. From
any state, enter environment1, prescribe the finitely many conversions
needed to turn all S2 into S1 and the degradations needed to remove all S3,
with no births or inflows during this finite prescribed sequence. Then
switch to environment2 and degrade all remaining S1 with no other reaction
before reaching zero. Each prescribed finite sequence has strictly positive
probability because alpha,kappa>0 and each required rate is positive when
its reactant is present.

Conversely, to reach any target (x1,x2,x3) from zero, first enter
environment2 and immigrate one S1. In environment1 prescribe exactly x2
S1-birth reactions, producing x2 molecules of S2 and 1+3x2 molecules of
S1, and prescribe x3 S3-immigration events, with no unwanted events.
In environment2 degrade all 1+3x2 S1 molecules, keeping S2 and S3
unchanged during this finite sequence, and immigrate x1 S1 molecules.
The target is reached. An additional switch gives either target environment.
This proves irreducibility of the full chain on Z_nonnegative^3 x {1,2}.

A recurrent state cannot have a positive probability of escaping all finite
sets, because recurrence entails infinitely many visits to that state.
Therefore the single-S2 state is transient. By irreducibility all states
are transient. For a nonexplosive countable-state CTMC with strictly
positive finite holding rates, transience implies that each state, and
hence every finite set, is visited only finitely often almost surely.
Here the total holding rate is always at least 1+kappa, so there are no
absorbing states or infinite last holding times. Finite-population sets
are finite after including the two environments. Taking their countable
union proves the asserted almost-sure divergence of the total population
from every initial state. This completes the source Conjecture46 proof.

## Dependencies, exact checks and limitations

The stochastic argument above is a proof candidate pending fresh independent
mathematical review. It is self-contained apart from standard SLLN,
Borel–Cantelli, conditional Cauchy–Schwarz and irreducible CTMC classification.
The exact checker verifies the source reaction drifts, active-pair mean
identities, comparison-moment identities, rational exponential-series bound,
conditional second-moment algebra and concrete versions of the generic
irreducibility paths. Such finite algebra/path checks are not a formalization
of the probability proof and cannot establish its universal statements alone.

No mathematical gap is currently known within the exact stated scope.
General Conjecture43, the optimal switching-rate threshold, other source
examples, experimental realization and biological data are not settled.
Current-priority screening remains a separate source-audit task.
