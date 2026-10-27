# Extinction hypotheses and an exponential-moment theorem for bounded offspring

## Scope

The source process on \(\mathbb Z_{\ge0}^d\) has jumps
\[
x\longmapsto x+\ell-e_i
\quad\hbox{at rate }\gamma_i(x_1)x_i\pi_{i,\ell}(x_1).
\]
Zero is absorbing. Offspring are nonnegative integer vectors.

Mehra and Taylor, arXiv:2609.12567v1, Section 2, impose a pointwise carrying-capacity spectral sign in Assumption A, uniform boundedness and convergence of rates and first/second offspring moments, and a positive lower bound on individual rates. Assumption B describes a uniformly positive probability of producing a type-1 offspring at a non-type-1 death; its display instead bounds that offspring mean. Section 2.1 then asserts that the limiting mean matrix is Hurwitz as a consequence. Theorem 3.1 asserts almost-sure extinction and all polynomial moments under A and B.

Two examples below test these literal hypotheses separately. They do not refute a theorem which additionally assumes the missing conditions. The positive result explicitly assumes a Hurwitz limit and accessibility of zero and, conservatively, the source's bounded-offspring Assumption C. It is not an unconditional assertion about the unrepaired source model or an experimentally calibrated infection model.

## 1. Zero need not be accessible, even with a Hurwitz limit

Take \(d=2\), \(\gamma_1=\gamma_2=1\), and type-2 offspring \((1,0)\) deterministically. For type 1, use offspring \((2,0)\) with probability \(m(z)/2\) and \((0,0)\) otherwise, where
\[
m(z)=
\begin{cases}
2,&0\le z\le1,\\
3-z,&1<z<3,\\
0,&z\ge3.
\end{cases}
\]
This is a bounded continuous probability specification; no differentiability is needed in A/B or Theorem 3.1. The mean matrix and its limiting generator matrix are
\[
M(z)=\begin{pmatrix}m(z)&0\\1&0\end{pmatrix},
\qquad H=\lim_{z\to\infty}(M(z)-I)
=\begin{pmatrix}-1&0\\1&-1\end{pmatrix}.
\]
The spectral abscissa is \(m(z)-1\), with the sign of \(2-z\). Thus A holds with integer carrying capacity two, including a strictly Hurwitz limit. B holds with probability one. All offspring have total size at most two and fit C with \(\Xi=4\). First and second moments are bounded and have the required limits.

Starting in \((1,0)\), the process has the closed positive class
\[
(1,0)\ \xrightarrow{1}\ (2,0),\qquad
(2,0)\ \xrightarrow{1}\ (1,0),\qquad
(2,0)\ \xrightarrow{1}\ (3,0),\qquad
(3,0)\ \xrightarrow{3}\ (2,0).
\]
Its generator on these three states is \(\begin{pmatrix}-1&1&0\\1&-2&1\\0&3&-3\end{pmatrix}\). It never reaches zero. In fact, once the first coordinate is positive it can never become zero: at one type-1 individual the only type-1 offspring is two, and all non-type-1 events add a type-1 individual. The same non-extinction conclusion therefore holds from \((1,1)\). Positive individual event rates do not make every nonzero state transient.

## 2. A pointwise subcritical tail need not have a strictly subcritical limit

Keep \(d=2\), both individual rates one and the same type-2 offspring. Replace the type-1 mean by
\[
m(z)=1+\frac{1-z}{(1+z)^3}.
\]
Again choose offspring two or zero with probabilities \(m(z)/2\) and \(1-m(z)/2\). For \(0\le z\le1\), \(1\le m(z)\le2\); for \(z\ge1\), \(3/4\le m(z)\le1\). All these probabilities are valid. Literal A has carrying capacity one, B holds with probability one and C holds as above. But \(m(\infty)=1\), so
\[
H=\begin{pmatrix}0&0\\1&-1\end{pmatrix}
\]
is singular. Pointwise negative spectral abscissae above the carrying capacity have not supplied a negative limiting gap.

On the type-1 axis, set
\[
a_n=\frac{n-1}{(n+1)^3},\qquad
\lambda_n=\frac n2(1-a_n),\qquad
\mu_n=\frac n2(1+a_n),\quad n\ge1.
\]
The holding rate is \(n\), jumps are \(\pm1\), and the embedded upward probability is at most \(1/2\). The process is nonexplosive, for example by domination of population increases by a linear Yule process. The finite-interval birth-death hitting formula, obtained by solving
\(\lambda_i(h_{i+1}-h_i)=\mu_i(h_i-h_{i-1})\), gives
\[
\mathbb P_1(T_N<T_0)=
\left(\sum_{j=0}^{N-1}\prod_{k=1}^j\frac{\mu_k}{\lambda_k}\right)^{-1}
\le\frac1N.
\]
A process avoiding zero cannot remain forever in a finite interval whose boundaries are accessible. Letting \(N\) tend to infinity proves \(T_0<\infty\) almost surely.

Let \(u_n\) be the expected occupation time at state \(n\) before absorption, starting at one. The expected number of \(1\to0\) jumps is one, so compensation gives \(\mu_1u_1=1\). Upcrossings and downcrossings of each edge agree pathwise before absorption. Induction, using compensation and nonnegative integrals, gives
\[
u_n=\frac1{\mu_n}\prod_{k=1}^{n-1}\frac{\lambda_k}{\mu_k}.
\]
Every individual \(u_n\) is finite; this argument does not assume that their sum is finite.

For \(n\ge1\), \(0\le a_n\le1/4\), and
\[
\sum_{n\ge1}a_n
\le\frac1{27}+\sum_{j\ge4}\frac1{j^2}
\le\frac1{27}+\int_3^\infty t^{-2}\,dt
=\frac{10}{27}.
\]
For every finite \(N\), the elementary product inequality gives
\[
\prod_{k=1}^N\frac{1-a_k}{1+a_k}
\ge\prod_{k=1}^N(1-2a_k)
\ge1-2\sum_{k=1}^Na_k\ge\frac7{27}.
\]
As \(\mu_n\le n\), \(u_n\ge7/(27n)\). Tonelli's theorem now yields
\[
\mathbb E_1T_0=\sum_{n\ge1}u_n=\infty.
\]
Zero is accessible from every finite two-type state: first choose successive type-2 deaths, then successive zero-offspring type-1 deaths. In fact, absorption is almost sure from every finite state: type 2 has no births and its finite pure-death population vanishes in finite time almost surely; nonexplosion leaves a finite type-1 population, to which the axis argument applies. From state one on the axis the first-step equation implies \(\mathbb E_2T_0=\infty\); from \((1,1)\), a first type-2 event has probability \(1/2\) and reaches \((2,0)\), so its mean extinction time is also infinite. This example therefore separates the limiting-gap issue from inaccessibility of zero.

## 3. Corrected bounded-offspring theorem

Assume the following explicit conditions for the displayed process:

1. \(0<g_-\le\gamma_i(n)\le g_+<\infty\).
2. Every offspring vector has total size at most a fixed \(\Xi<\infty\).
3. For every \(i\ge2,n\), \(\mathbb P_{i,n}(\ell_1\ge1)\ge\varepsilon>0\).
4. \(A(n)=\operatorname{diag}(\gamma_i(n))(M(n)-I)\) converges to a Hurwitz Metzler matrix \(H\).
5. Zero is reachable by a finite positive-rate path from every finite state.

Then there are constants \(\theta>0,C<\infty\), independent of the finite initial state, such that
\[
\mathbb E_x e^{\theta T_0}\le C(1+\lVert x\rVert_1).
\]
The exact carrying-capacity sign change is not needed for this corrected theorem. Conditions 4 and 5 are not deduced from the source's A/B wording.

For the source's mean version of B under C, \(\ell_1\le\Xi\) and
\(\mathbb E\ell_1\ge\varepsilon_m\) imply condition 3 with
\(\varepsilon=\varepsilon_m/\Xi\). Thus the probability/mean discrepancy can be accounted for without identifying the two constants.

### A global drift inequality

Since \(H\) is Hurwitz Metzler, choose a strictly positive vector \(\xi\) with \(H\xi<0\). Matrix convergence gives an integer \(Y\ge1\) and \(a>0\) such that
\[
A(n)\xi\le-a\xi\quad(n\ge Y).
\]
Put \(L=\xi\cdot x\), \(m=\min_i\xi_i\), \(M=\max_i\xi_i\), \(D=M\Xi\), and \(B=g_+D/m\). The generator satisfies
\[
QL\le BL\quad\hbox{everywhere},\qquad QL\le-aL\quad(x_1\ge Y).
\]
Choose
\[
0<\eta\le\frac{a}{4g_+Y(1+D/m)},\quad R=1+\eta,\quad
w(n)=R^{-\min(n,Y)},\quad V(x)=w(x_1)L(x).
\]
For a jump from \(x\) to \(x'\),
\[
Q V=w(x_1)QL+
\sum_{x'}q(x,x')\,[w(x'_1)-w(x_1)]L(x').
\]
An increase of \(w\) is possible only at a type-1 event with zero type-1 offspring and \(1\le x_1\le Y\). Its relative increase is \(\eta\), its event rate is at most \(g_+Y\), and \(L(x')\le L+D\). Because \(L\ge m\) at nonzero states, all positive corrections are bounded by
\[
g_+Y\eta\,w(L+D)\le(a/4)V.
\]
Every other correction is nonpositive. Hence for \(x_1\ge Y\),
\(QV\le-3aV/4\).

Now let \(n=x_1<Y\) and \(N=\sum_{i\ge2}x_i\). At each non-type-1 event with \(\ell_1\ge1\), \(w(x'_1)\le w(n)/R\). If \(L\ge2M\), \(L(x')\ge L-M\ge L/2\). Such events have total rate at least \(g_-\varepsilon N\), giving
\[
QV\le(B+a/4)V-cNV,\qquad
c=\frac{g_-\varepsilon\eta}{2R}>0.
\]
Choose an integer \(N_0\ge\max\{2M/m,(B+a)/c\}\). Outside the finite set
\[
\mathcal C=\{x:x_1<Y,\ \sum_{i\ge2}x_i\le N_0\},
\]
which contains zero, we have \(QV\le-aV/2\). This covers the entire infinite state space, including the high-parasite/low-controller strip.

Let \(v_*=R^{-Y}m\) and \(W=1+V/v_*\). Then \(W\ge1\), \(W\) is comparable to \(1+\lVert x\rVert_1\), and \(W\ge2\) at every nonzero state. Therefore
\[
QW\le-(a/4)W\quad\hbox{outside }\mathcal C.
\]

### From geometric return to exponential absorption

Write \(\delta=a/4\) and \(\tau_{\mathcal C}\) for the first entrance time in \(\mathcal C\). Stopped Dynkin's formula and Fatou's lemma for
\(e^{\delta(t\wedge\tau_{\mathcal C})}W(X_{t\wedge\tau_{\mathcal C}})\)
give
\[
\mathbb E_x e^{\delta\tau_{\mathcal C}}\le W(x).
\]
Localization is legitimate: bounded offspring and linear total event rates give nonexplosion and finite population moments on every finite time interval, by comparison with a bounded-batch linear birth process. The same bounds give
\(\sup_{c\in\mathcal C}\mathbb E_c W(X_1)<\infty\).

For a process starting in \(\mathcal C\), put
\(\sigma=1+\tau_{\mathcal C}\circ\vartheta_1\).
The preceding return bound implies
\[
\sup_{c\in\mathcal C}\mathbb E_c e^{\delta\sigma}
\le e^\delta\sup_{c\in\mathcal C}\mathbb E_c W(X_1)<\infty.
\]
By zero accessibility and finiteness of \(\mathcal C\), a finite positive-rate absorption path can be completed within one time unit with probability at least a common \(p>0\), uniformly over its nonzero states. Zero remains absorbing thereafter.

For \(c,c'\in\mathcal C\setminus\{0\}\), define
\[
K_\theta(c,c')=
\mathbb E_c[e^{\theta\sigma};
X_\sigma=c',\ T_0>\sigma].
\]
At \(\theta=0\), each row sum is at most \(1-p\). The uniform exponential bound and dominated convergence, on this finite matrix, give some \(0<\theta<\delta\) for which every row sum is at most \(\rho<1\). The contribution when absorption occurs during a cycle is bounded by \(\sup_c\mathbb E_c e^{\theta\sigma}<\infty\). Summing the substochastic renewal series yields a finite bound on \(\mathbb E_c e^{\theta T_0}\), uniformly over \(c\in\mathcal C\). It also proves eventual absorption, without invoking the unrepaired source theorem.

Finally apply the strong Markov property at \(\tau_{\mathcal C}\). Its exponential bound gives
\(\mathbb E_x e^{\theta T_0}\le C W(x)\), proving the stated result.

## Status and limits

The generator, cutoff weight and finite-set argument are analytic. Finite-state truncation, simulation and numerical fitting are not used to infer an infinite-state theorem. The two source-contract examples and the corrected theorem must be read at their distinct stated scopes. No biological parameterization, universal extinction claim under unrepaired A/B, quasi-stationary distribution, sharp exponent or unbounded-offspring theorem is asserted.

The construction used substantial AI assistance. Separate source/value and fresh mathematical AI contexts reviewed the stated scope. These checks are not human peer review or priority certification. Classical birth-death occupation identities, Metzler-matrix Lyapunov theory and regeneration arguments are credited as methods; the proposed contribution is the source-specific hypothesis audit and the explicit bounded-offspring model theorem.
