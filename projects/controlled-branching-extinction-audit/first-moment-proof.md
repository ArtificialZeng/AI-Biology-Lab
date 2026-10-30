# Exponential extinction moments with unbounded offspring

This strengthens the corrected bounded-offspring theorem in `proof.md`. The two
literal-hypothesis counterexamples in that file are unchanged. The result below
does not repair the original hypotheses silently: strict negativity of the
limiting spectral abscissa and accessibility of zero are explicit assumptions.

## Statement

Consider the minimal continuous-time Markov process on nonnegative integer
vectors of dimension d, with a type-i event at state x replacing that individual
by a finite nonnegative integer offspring vector ell. Its rate and offspring law
are gamma_i(n) x_i and pi_(i,n), where n=x_1. Self-replacements may be retained as
null events; they do not change the generator. Zero is absorbing. Assume:

1. There are constants 0<g_-<=g_+<infinity with
   g_-<=gamma_i(n)<=g_+ for every type and every nonnegative integer n.
2. K=sup_(i,n) E_(i,n)|ell|_1 is finite. Offspring support may be unbounded,
   and no second moment is required.
3. For some epsilon>0, P_(i,n)(ell_1>=1)>=epsilon for every i>=2 and n.
   This condition is vacuous if d=1.
4. Writing M_ij(n)=E_(i,n)ell_j and
   A(n)=diag(gamma_i(n))(M(n)-I), the matrices A(n) converge to a Hurwitz
   Metzler matrix H as n tends to infinity.
5. From every finite state a finite path of positive-rate state changes leads
   to zero.

Then the process is nonexplosive. There are theta>0 and C<infinity, independent
of the finite initial state, such that

    E_x exp(theta T_0) <= C(1+|x|_1),

where T_0 is the first hitting time of zero. In particular absorption occurs
almost surely. No exponential moment of the offspring distribution is assumed.

These are explicit sufficient hypotheses, not a claim that the unrepaired
Assumptions A/B of Mehra--Taylor imply assumptions 4 and 5. A pointwise
carrying-capacity spectral sign is not needed once these conditions hold.

## 1. Nonexplosion and finite-time first moments

Let U(x)=1+|x|_1. Conditional on an event, the absolute change of U is at most
|ell|_1+1, so its conditional expectation is at most K+1. The total event rate
is at most g_+|x|_1, and

    QU(x) = sum_i gamma_i(n)x_i (E_(i,n)|ell|_1 - 1)
          <= b U(x),                  b=g_+K.

This includes K=0, when the displayed upper bound is zero. Construct the
minimal event process and let tau_m be its first entrance into |x|_1>=m,
with m>|x|_1 initially. Before tau_m the event rate is bounded by g_+m.
There cannot be infinitely many events before t and before tau_m. More
explicitly, a rate-g_+m Poisson clock, with state-dependent thinning before the
stop, constructs this stopped process; the exit event is included and then
the state is held fixed.

The expected total absolute change of U in these stopped events up to t is at
most g_+m t (K+1). To see this without an integrability assumption, first sum
only the first J clock events, use their conditional offspring expectation
bound K+1, and then use monotone convergence on the nonnegative absolute
increments. The expected number of clock events is g_+m t. This argument also
shows that the unbounded exit overshoot has finite expectation. Compensation
or conditioning on successive clock events therefore gives stopped Dynkin's
identity for U, with no unbounded-function interchange left unjustified.

Consequently

    E U(X_(t wedge tau_m))
      <= U(x) + b integral_0^t E[U(X_s) 1_(s<tau_m)] ds
      <= U(x) + b integral_0^t E U(X_(s wedge tau_m)) ds
      <= U(x) exp(bt).

On {tau_m<=t}, the stopped population is at least m. Hence

    P_x(tau_m<=t) <= U(x) exp(bt)/m.

Explosion while the population remains bounded is impossible by the bounded
event rate there. Any explosion by time t must therefore cross every finite
population threshold by t; the last probability bound rules it out. Now the
process is nonexplosive, and for fixed t it has only finitely many events,
each of finite offspring size. Thus X_(t wedge tau_m) converges to X_t as m
tends to infinity. Fatou gives

    E_x U(X_t) <= U(x) exp(bt).                         (1)

This proof uses first moments only, including at the unbounded overshoot.

## 2. A cutoff weight gives drift on the entire state space

Since H is Hurwitz Metzler, choose a strictly positive vector xi with H xi<0.
For example xi=integral_0^infinity exp(Ht) 1 dt has this property. Matrix
convergence supplies an integer Y>=1 and a>0 with

    A(n) xi <= -a xi            for every n>=Y.

Write L(x)=xi dot x, m=min_i xi_i>0, M=max_i xi_i, D=MK, and B=g_+D/m.
Directly from the mean generator,

    QL <= BL everywhere,             QL <= -aL when x_1>=Y.       (2)

Choose any eta satisfying

    0<eta<= a/[4g_+Y(1+D/m)],

and put R=1+eta, w(n)=R^(-min(n,Y)), and V(x)=w(x_1)L(x). At every finite
state the generator sums for L and V are absolutely convergent: w is bounded
between R^(-Y) and 1, L(after)<=L(before)+M|ell|_1, and the offspring first
moment is finite. The product identity is therefore legitimate:

    QV(x)=w(n)QL(x)
          + sum_i gamma_i(n)x_i E[(w(n')-w(n)) L(x')].           (3)

The weight can increase only in a type-1 event with ell_1=0 and 1<=n<=Y.
There w(n')-w(n)=eta w(n). Although L(x') has no pathwise upper bound, the
needed conditional integral has the uniform bound

    E[1_(ell_1=0) L(x')] <= L(x)+M E|ell|_1 <= L(x)+D.

Its rate multiplier is at most g_+Y. Since L>=m at nonzero states, all positive
terms in the correction in (3), taken together, are bounded above by

    g_+Y eta w(n)(L+D) <= (a/4)V.                              (4)

Every other correction is nonpositive. In particular x_1>=Y implies
QV<=-3aV/4, including x_1=Y where a controller loss can increase w.

Next suppose n<Y and let N=sum_(i>=2)x_i. At a non-controller event with
ell_1>=1 we have w(n')<=w(n)/R. Also, regardless of the offspring size,

    L(x') >= L(x)-xi_i >= L(x)-M.

If L>=2M, such an event contributes at most
-eta w(n)L/(2R) to the weight correction. Their total rate is at least
g_- epsilon N, by assumption 3. Using (2) and (4) in (3) yields

    QV <= (B+a/4)V-cNV,             c=g_- epsilon eta/(2R)>0.  (5)

For d=1 this step is unnecessary, since n<Y already describes a finite set;
one may take any positive epsilon in the definitions below and N=0.
Choose an integer N_0>=max{2M/m,(B+a)/c}. Outside

    C_0={x: x_1<Y, sum_(i>=2)x_i<=N_0},

which is finite and contains zero, (2)--(5) imply QV<=-aV/2. This covers the
whole state space, including arbitrarily large non-controller populations
when the controller count is below Y.

Let v_*=R^(-Y)m and W=1+V/v_*. Then W>=1, W is comparable to 1+|x|_1, and
W>=2 at nonzero states. Outside C_0,

    QW <= -(a/2)(W-1) <= -delta W,             delta=a/4.       (6)

## 3. Exponential return, with explicit localization

Let tau_C be the first hitting time of C_0. If x lies in C_0 the following
bound is immediate. Otherwise stop at

    rho=t wedge tau_C wedge tau_m.

Before rho, population size is less than m and time at most t. The expected
sum of absolute increments of W at the stopped event times is finite: W is
bounded by a constant times U, and for a pre-stop state the expected absolute
increment is bounded by a constant times m+K+1. The same Poisson-clock
argument from Section 1 applies, also after multiplying by exp(delta t).
Thus the time-dependent stopped Dynkin identity is valid, including the
unbounded exit overshoot. Equation (6) gives

    E_x[exp(delta rho) W(X_rho)] <= W(x).

Nonexplosion allows m to tend to infinity. Fatou gives

    E_x[exp(delta(t wedge tau_C)) W(X_(t wedge tau_C))] <= W(x).

Since W>=1, monotone convergence for exp(delta(t wedge tau_C)) now yields

    E_x exp(delta tau_C) <= W(x).                              (7)

In particular tau_C is finite almost surely; this is proved before using any
return/absorption cycle. Equation (1) and W<=constant*U also give

    sup_(z in C_0) E_z W(X_1) < infinity.                      (8)

Neither (7) nor (8) assumes a finite offspring second moment.

## 4. From finite-set return to exponential absorption

For a process beginning at a point of C_0 set

    sigma=1+tau_C after time 1.

The strong Markov property at deterministic time 1, (7), and (8) give

    sup_(z in C_0) E_z exp(delta sigma)
       <= exp(delta) sup_(z in C_0) E_z W(X_1) < infinity.     (9)

For each nonzero z in C_0, choose a finite positive-rate path to zero.
All its finitely many waiting rates are finite, and each selected transition
has positive conditional probability. There is therefore positive probability
that this path is completed within one time unit. Finiteness of C_0 supplies
a common p>0 such that P_z(T_0<=1)>=p for every nonzero z in C_0.
If C_0 contains only zero, (7) already proves the theorem.

Otherwise, for nonzero z,z' in C_0 define the finite matrix

    K_theta(z,z') = E_z[exp(theta sigma);
                         X_sigma=z', T_0>sigma].

Its row sums at theta=0 are at most 1-p. Domination by exp(delta sigma) in
(9), and finiteness of the matrix, imply that some 0<theta<delta makes every
row sum at most rho_0<1. The one-cycle absorption contribution

    b_theta(z)=E_z[exp(theta T_0); T_0<=sigma]

is uniformly finite by (9). Repeating the cycles, their survival probabilities
are at most (1-p)^j after j cycles. Every cycle is finite almost surely and
each lasts at least one unit, so this proves eventual absorption. Decomposing
according to the first cycle in which absorption occurs, nonnegative
expectations give exactly the renewal sum

    E_z exp(theta T_0) = sum_(j>=0) (K_theta^j b_theta)(z)
                       <= max b_theta/(1-rho_0) < infinity.

This series proof does not assume exponential integrability in advance.
Finally use the strong Markov property at tau_C and (7), with theta<delta.
The common finite bound on C_0 gives E_x exp(theta T_0)<=C'W(x), and
W<=constant*(1+|x|_1) proves the stated estimate.

## 5. Consequence for the source's unbounded-offspring scope

The source's Assumption A includes uniformly bounded first and second
offspring moments and positive bounded individual rates. If its B is read
as the prose probability bound, assumption 3 holds directly. If B is read as
its displayed mean bound E ell_1>=epsilon_m>0, let
K_2=sup_(i>=2,n) E ell_1^2<infinity. Cauchy--Schwarz gives

    epsilon_m^2 <= (E ell_1)^2
                 <= E ell_1^2 P(ell_1>=1)
                 <= K_2 P(ell_1>=1).

Thus the probability bound holds with epsilon=epsilon_m^2/K_2. For d>=2,
epsilon_m>0 ensures K_2>0. For d=1 B is vacuous.

Consequently A and either reading of B, **together with an explicitly Hurwitz
limiting mean matrix and zero accessibility**, imply the exponential-moment
conclusion without the bounded-support Assumption C. This establishes the
exponential-moment assertion for that repaired unbounded-offspring scope.
The two extra conditions are not consequences of literal A/B, as the two
unchanged examples in `proof.md` demonstrate. The literal original conjecture
must therefore not be reported as proved unconditionally.

## 6. An unbounded-support, infinite-second-moment illustration

This example shows that the first-moment theorem covers laws beyond both the
earlier bounded-offspring theorem and the source's second-moment hypothesis.
It is an illustration, not a novel branching mechanism or biological fit.

Take two types and individual rates one. At controller level n a type-1 event
produces (2,0) with probability 1/(n+1), and (0,0) otherwise. A type-2 event
produces (1,Z), where

    P(Z=0)=3/4,
    P(Z=k)=1/[k(k+1)(k+2)] for k>=1.

Telescoping gives sum_(k>=1)P(Z=k)=1/4 and EZ=1/2. Moreover
k^2 P(Z=k)=k/[(k+1)(k+2)]>=1/(6k), so EZ^2=infinity.
The total offspring first moments are at most max{2,3/2}=2. The mean matrix
and its limiting generator are

    M(n)=[[2/(n+1),0],[1,1/2]],
    H=[[-1,0],[1,-1/2]],

so H is strictly Hurwitz and every non-controller event creates a controller.
Zero is accessible: choose successive type-2 events with Z=0 until there are
no type-2 individuals, then choose zero-offspring type-1 events (each has
positive probability at every positive controller count). This is a finite
positive-rate path from any finite state.

The theorem therefore applies, despite the infinite offspring second moment.
It is not claimed that this example meets source A's second-moment requirement.
Its mean spectral sign is sign(1-n), but that does not restore the missing
second moment. The existing finite checker in this packet concerns the two original
witnesses. This heavy-tailed illustration and the infinite-state theorem use
the analytic arguments above, not a finite prefix of computations.

## Status and methods

This is an analytic theorem with explicit hypotheses. The cutoff function is
the same model-specific ingredient as in the earlier bounded-offspring note;
the strengthening replaces upper jump bounds by expectation bounds and gives
the unbounded-jump localization argument. General Lyapunov, stopped Dynkin,
Metzler-matrix and finite-set renewal methods are classical and are not claimed
as new. No finite simulation, truncation, empirical biological validation,
sharp extinction exponent, human peer review or priority certification is
asserted. Substantial AI assistance was used. This note updates the same
research packet; it is not a separate result count or a new paper.
