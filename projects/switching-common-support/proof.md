# Extinction despite intersecting positive-drift supports

This exposition concerns the literal statement of Conjecture 43 in Cappelletti,
Howells and Xu, *Stability of randomly switching stochastic reaction networks
with asymptotically linear transition rates*,
arXiv:2507.09421v2, Section 6, p. 29. The model and argument below are the fixed
candidate considered in the separate source-contract and mathematical reviews.
Those reviews support the stated counterexample. They are AI-assisted checks,
not external human peer review or a certification of priority. No claim is made
about a strengthened conjecture that adds population irreducibility or a
species-regeneration condition.

## Model and conclusion

Let the molecular population be \(X(t)=(X_A(t),X_B(t),X_C(t))\in\mathbb N_0^3\).
Independently of the reaction clocks, let the environment \(J(t)\in\{1,2\}\)
switch from either state to the other at rate \(\kappa>0\). The generator of
the environment is therefore

\[
Q=\kappa\begin{pmatrix}-1&1\\1&-1\end{pmatrix}.
\]

In environment 1 the reactions and rate constants are

\[
A\xrightarrow{1}2A,\qquad B\xrightarrow{1}A,\qquad
C\xrightarrow{2}\varnothing.
\]

In environment 2 they are

\[
A\xrightarrow{2}\varnothing,\qquad B\xrightarrow{1}C,\qquad
C\xrightarrow{1}2C.
\]

The constants over the arrows are mass-action rate constants: a reaction with
source species \(S\) has propensity equal to that constant times \(X_S\).
There are no inflows. Start from \(X(0)=(0,1,0)\), with either initial
environment. For every \(\kappa>0\), the molecular population reaches zero in
finite time almost surely. Nevertheless, the two networks have nonnegative
left vectors with strictly positive drift on their respective supports, and
the initial species belongs to the intersection of those supports.

## Drift, support, and nonexplosion

Using column population vectors, sum the rate constant times the reaction
increment in the column corresponding to each source species. This gives

\[
M_1=\begin{pmatrix}1&1&0\\0&-1&0\\0&0&-2\end{pmatrix},
\qquad
M_2=\begin{pmatrix}-2&0&0\\0&-1&0\\0&1&1\end{pmatrix}.
\]

For the nonnegative row vectors
\(v_1=(2,1,0)\) and \(v_2=(0,1,2)\), direct multiplication yields

\[
v_1M_1=v_1,\qquad v_2M_2=v_2.
\]

Thus \((v_iM_i)_m>0\) whenever \(m\in\operatorname{supp}(v_i)\), while

\[
\operatorname{supp}(v_1)=\{A,B\},\qquad
\operatorname{supp}(v_2)=\{B,C\},\qquad
\operatorname{supp}(v_1)\cap\operatorname{supp}(v_2)=\{B\}.
\]

The initial support meets this intersection. All six rate constants are
strictly positive and all sources are monomolecular. The drifts are exactly
linear, so the asymptotically linear drift requirement is satisfied. The
finite environment chain is irreducible for every \(\kappa>0\).

Write \(N=X_A+X_B+X_C\). The total molecular jump rate is \(X_A+X_B+2X_C\)
in environment 1 and \(2X_A+X_B+X_C\) in environment 2, hence is at most
\(2N\). Every upward population jump has size one, and its total rate is
\(X_A\) in environment 1 and \(X_C\) in environment 2, hence at most \(N\).
Consequently the population can be dominated by a linear pure-birth process
with transition rate \(n\mapsto n+1\) equal to \(n\). That process is
nonexplosive. On a finite time interval the resulting population bound is
finite, and the total rate of molecular and environment events is bounded
by \(2N+\kappa\). This also excludes accumulation of events with bounded
population. The joint process is therefore nonexplosive.

## Reduction to a single lineage

No reaction produces \(B\). In either environment the only reaction available
to a \(B\) particle consumes it at rate one. Until the initial \(B\) converts,
there are no \(A\) or \(C\) particles and hence no other molecular reactions.
Its conversion time \(\tau\) can thus be realized by a single rate-one clock
independent of the environment. In particular, \(\tau<\infty\) almost surely.

At time \(\tau\), the population consists of a single \(A\) if the conversion
occurs in environment 1, or a single \(C\) if it occurs in environment 2.
The two independent clocks do not have simultaneous jumps almost surely,
so this description is unambiguous. Thereafter every descendant retains that
same species: \(A\) reactions produce only \(A\), and \(C\) reactions produce
only \(C\). For either possible species, call its birth environment the good
environment and the other the bad environment. Per particle, births occur
at rate one in the good environment and deaths at rate two in the bad
environment. There are no births in the bad environment and no deaths in
the good environment.

## Conditional first moment and extinction

Condition on the full environment path and on \(\tau\). These data also
determine which species starts the lineage. Its initial size is one.
Conditional on these data, the postconversion population is a one-type
linear birth-death process with the prescribed time-dependent per-particle
birth and death rates. Independent individual reaction clocks give its
finite-time conditional mean

\[
m(t):=\mathbb E\!\left[N(\tau+t)\mid J(\cdot),\tau\right]
=\exp\!\left(\int_0^t
\bigl(b(\tau+s)-d(\tau+s)\bigr)\,ds\right).
\]

Here the integrand is \(+1\) in the good environment and \(-2\) in the bad
environment. Equivalently, if \(G(t)\) and \(D(t)\) are the amounts of good
and bad environment time in \([\tau,\tau+t]\), then
\(m(t)=\exp(G(t)-2D(t))\) and \(G(t)+D(t)=t\).

The symmetric irreducible two-state environment has stationary distribution
\((1/2,1/2)\). Its occupation fractions converge almost surely to those
values from either initial environment. Shifting the start by the finite
time \(\tau\) does not change those limits: the portion removed or added
has length at most \(\tau\), and its ratio to \(t\) tends to zero. This
applies to both possible assignments of the good environment. Therefore

\[
\lim_{t\to\infty}\frac{\log m(t)}{t}
=1\cdot\frac12-2\cdot\frac12=-\frac12
\quad\text{almost surely}.
\]

In particular, \(m(t)\to0\) for almost every conditioned environment path
and conversion time. For integer offsets \(n\), conditional Markov
inequality gives

\[
\mathbb P\!\left(N(\tau+n)>0\mid J(\cdot),\tau\right)
\leq \min\{1,m(n)\}\longrightarrow0.
\]

Zero molecular population is absorbing because there are no inflows.
Consequently the events \(\{N(\tau+n)>0\}\) decrease with \(n\).
Continuity of conditional probability from above shows that their
intersection has conditional probability zero for almost every environment
path and conversion time. Averaging gives unconditional probability zero
as well. Thus, almost surely, some finite integer \(n\) satisfies
\(N(\tau+n)=0\). Since \(\tau<\infty\), molecular extinction occurs in
finite time and the population remains zero afterwards. The environment
continues switching.

## Relation to the literal conjecture and to the checker

The source defines evanescence using a strictly positive probability that
the process norm tends to infinity. Since the environment takes only two
values, divergence of the mixed-process norm is equivalent here to divergence
of the molecular population norm. Almost-sure finite-time extinction gives
probability zero to that event. Hence the initial state
above is not evanescent for any \(\kappa>0\), including arbitrarily small
positive switching rates. Together with the verified drift and support
conditions, this contradicts the literal conclusion of Conjecture 43.

The source-contract review checked that the relevant irreducibility
requirement concerns the environment; the literal conjecture does not impose
irreducibility of the full population chain. The present population chain
is reducible. A transience observation for the initial \(B\) state would
not establish the stronger evanescence conclusion in that setting. If a
different source premise were imposed that excluded this reducible model,
the counterexample label for that modified statement would not be ready;
the model and extinction calculation would not resolve it.

The exact checker derives the matrices and rate identities from the six
arrows in `model.json`, verifies the support and lineage structure, and
checks exact coefficientwise rate bounds. Its small-tuple tests and altered
input controls verify computation behavior. These finite checks do not
prove the nonexplosion argument, the conditional mean formula, the ergodic
theorem, or the probability limit leading to extinction. Those steps are
the mathematical argument above, subject to its separate review. This
result does not modify the separate Example 44/Conjecture 46 process, whose
active reactions regenerate the common species.
