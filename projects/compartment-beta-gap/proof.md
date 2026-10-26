# Candidate partial classification: explosion for 10/11 < beta < 1

This is a frozen proof candidate for fresh mathematical review, not an acceptance or public-priority claim. All rate conventions and kernel parameters are those in source.md/problem.md. The public endpoint results of Yoo v1.1.0 are credited prior work.

## Exact statement

Fix beta in (10/11,1), c in (0,1), and p in (0,1). Under the allocation kernel

psi_beta((e,s),.) = (1-q(s)) psi_b((e,s),.) + q(s) kappa_p((e,s),.),
q(s)=c/(1+s)^beta,

the model of current equation (7) has strictly positive probability of finite-time explosion from every finite state with positive total substrate. No almost-sure claim is made. Zero total substrate is the already-known nonexplosive class. No conclusion is drawn for 0<beta<=10/11.

The argument in fact permits any admissible replacement kernel in place of kappa_p, because it only uses the actual mixture before its first marked replacement event.

## 1. Empty-budding variables and exact rates

Start an empty-budding process with exactly one substrate-bearing compartment, called the tag. Empty budding preserves this property. Let e,s be the tag's enzyme and substrate counts, C the total compartment count, and B the number of outside compartments containing at least one enzyme. Thus s>=1 and C>=B+1. The four counts need not form a Markov chain; all calculations are for functions of the full state.

The exact predictable rates are:

- e gains one by direct tag immigration at rate 1.
- s gains one at rate es.
- C gains one by fragmentation at rate s.
- C loses one by coagulation at rate C(C-1).
- B gains one by immigration into an outside empty compartment at rate C-1-B.
- B loses one by merging two outside enzyme carriers at rate B(B-1), or by a carrier merging with the tag at rate 2B. Thus its total loss rate is B(B+1).
- At each tag-carrier merger, the tag receives k>=1 enzymes. The total rate of such events is 2B.

Correlated changes in e,C,B are retained by an additive test function; no independence of these counts or molecule import times is assumed.

Put

F = 1/log(e+2) + 1/log(s+2) + 1/(C+1) + 1/(B+1).

This is bounded and strictly positive at every finite state. The inequality log(1-u)<=-u for 0<u<1, the fact that importing k>=1 decreases the first summand by at least the decrease from importing one, and the exact rates give

LF <= 5/2 - T,
T=A+D+H+J,
A=s/(C+2)^2,
D=C/(B+2)^2,
H=2B/[(e+3)log(e+3)log(e+2)],
J=e/[4log(s+3)^2].                                      (1)

For completeness, before the final simplification the bound is

LF <= (C-1)/(C+1) + 1_{B>=1}
      -s/[(C+1)(C+2)]
      -(C-1-B)/[(B+1)(B+2)]
      -2B/[(e+3)log(e+3)log(e+2)]
      -es/[(s+3)log(s+3)log(s+2)].

The direct immigration term in e has been discarded because it is negative. The identity

-(C-1-B)/[(B+1)(B+2)]
= -C/[(B+1)(B+2)] + 1/(B+2)

and s/(s+3)>=1/4 yield (1). In particular, all positive contributions from correlated coagulation events are included.

## 2. A polynomial occupation drift

Let delta=1-beta in (0,1/11), r=1-11delta=11beta-10 in (0,1), and define the explicit constants

L=4096/r^2,     a=1/(L+1).                              (2)

We claim

LF <= -s^delta whenever F<a.                           (3)

Indeed F<a implies e,B,C>=3 and t=log(s+3)>L. For these states,

A D^2 H^4 J^4
= (s/16) [C/(C+2)]^2 [B/(B+2)]^4 [e/(e+3)]^4
  /[log(e+3)^4 log(e+2)^4 log(s+3)^8]
>= 2^(-14) s/[log(e+3)^8 log(s+3)^8].                  (4)

Every term of T is at most T, so A D^2 H^4 J^4 <= T^11.

If e<=s, (4) implies

T^11 >= 2^(-14) s/t^16.                                (5)

For all t>=L we have exp(rt)>=2^38 t^16. To verify this uniformly, the function rt-16log(t) is increasing there because L>=16/r. At L, use log(4096)<12, log(1/r)<=1/r, and log2<1 to obtain

rL-16logL-38log2
> 4096/r - 16(12+2/r) -38
=4064/r-230>0.

Since s>=exp(t)/4 and r<1, it follows that

s^r >= exp(rt)/4 >= 2^36 t^16.

Substituting this into (5) yields T>=4s^delta.

If e>s, then J>=s/(4t^2). For t>=4096, exp(t/2)>=64t^2: the logarithmic difference t/2-2logt-6log2 is positive at 4096 and increasing thereafter. Since 1-delta>1/2,

s^(1-delta) >= exp((1-delta)t)/4 >= exp(t/2)/4 >=16t^2.

Thus T>=J>=4s^delta in this case as well. Equation (1), together with s^delta>=1, gives LF<=5/2-4s^delta<=-s^delta, proving (3).

The exponent 1/11 comes only from the weights (1,2,4,4) in the product (4). No sharpness is claimed.

## 3. Stopping, explosion and the integrated intensity bound

Let z be a finite one-tag initial state with F(z)<a. For the empty-budding process write zeta for its minimal-chain explosion time, and sigma for its first exit from D={F<a} before explosion, setting sigma=infinity if no such exit occurs. Let tau=min(sigma,zeta). Exhaust the countable full state space by finite sets, with exit times tau_m increasing to zeta, and stop at t∧sigma∧tau_m.

Because F is bounded and rates are bounded within each finite set, stopped Dynkin and (3) give

E F(X_{t∧sigma∧tau_m})
 + E integral_0^{t∧sigma∧tau_m} s(u)^delta du <= F(z).

Fatou and monotone convergence give

E integral_0^tau s(u)^delta du <= F(z),   E tau <= F(z).  (6)

Here s^delta>=1. The boundary value is at least a when exit occurs. Letting m,t increase in the same stopped inequality yields

P(sigma<zeta) <= F(z)/a,
P(zeta<infinity, zeta<=sigma) >=1-F(z)/a.                (7)

There is no coercivity assumption on F. Although the projection does not control every outside enzyme count, finite-state localization is performed on the complete CTMC state, not just the four displayed counts. Every state has a finite total jump rate. Thus the limit of the finite-state exits is the actual minimal-chain explosion time; finite rates preclude infinitely many jumps within one finite set.

The substrate counts in time integrals are pre-jump values (the distinction at individual times does not change the integrals). The number of empty-budding fragmentations has predictable intensity s(u). Consequently its weighted stopped compensator identity, first on finite localizations and then by monotone convergence, gives

E sum_{fragmentation times t_j<tau} q(s(t_j-))
 <= E integral_0^tau s(u) q(s(u)) du
 <= c E integral_0^tau s(u)^delta du
 <= c F(z).                                            (8)

The inequality uses beta=1-delta, s>=1, and
s/(1+s)^beta <= s^(1-beta)=s^delta.
On each finite localization the compensated stopped sum including a possible terminal fragmentation has mean equal to the integrated intensity. Discarding a terminal event only lowers the sum. Fatou and monotone convergence give the displayed bound for events strictly before tau. A finite expectation proves that the weighted sum is finite almost surely, even on the finite-time explosion event.

## 4. Coupling the exact mixture

Construct the empty-budding trajectory first, up to its explosion time. Independently attach a uniform mark U_j to each of its fragmentation events. Declare an event marked if U_j<=q(s(t_j-)). Up to the first marked event construct the mixture process with exactly the same immigration, substrate production, coagulation and fragmentation clocks and jumps as the empty-budding process. At the first mark, draw daughters using kappa_p instead; thereafter continue with the ordinary mixture generator.

Before the first mark, the two processes have identical full states. At a fragmentation in that common state the mixture uses psi_b with probability 1-q(s) and kappa_p with probability q(s), exactly as required. No hypothesis about the process after the first mark is needed. Marking an alternative allocation that happens to produce empty budding is allowed; it only weakens the lower bound.

By the conditional union bound and (8),

P(at least one marked fragmentation before tau) <= cF(z).

On the event that the baseline explodes before exit and there is no preceding mark, the mixture follows the identical infinite trajectory up to the identical finite explosion time. Therefore (7) gives

P_mixture(zeta<infinity) >= 1-F(z)/a-cF(z).              (9)

Equivalently, conditionally on a baseline path, the no-mark probability is the infinite product of 1-q(s(t_j-)); (8) makes that product strictly positive. The union-bound form (9) also gives a quantitative lower bound without conditioning on future survival.

## 5. Explicit finite initial state and extension to every positive-substrate state

With a from (2), set

M=ceil(64/a),    N=2^M.

Take one tag containing N enzymes and N substrate molecules and N outside compartments each containing exactly one enzyme and no substrate. Then B=N, C=N+1, and

F_N=2/log(N+2)+1/(N+2)+1/(N+1)
   <4/M+2/N <=6/M <=3a/32 <a/4.

We used log2>1/2 and 2^M>=M. Since a<1 and c<1, (9) gives explosion probability strictly greater than 1/2 from this explicitly specified finite state.

For example beta=11/12 gives r=1/12, L=589824, a=1/589825, M=37748800 and N=2^37748800. These are deliberately loose proof constants; no simulation or construction at that size was performed.

Finally, from any finite state with positive total substrate there is a finite positive-probability path to a state with F<=F_N. Coagulate all initial compartments into one. If necessary immigrate enzymes in that compartment until e>=N, then produce substrate until s>=N. Choose N consecutive empty-budding fragmentation outcomes, creating N outside empty compartments. Give each outside compartment exactly one enzyme by its immigration reaction. Each desired reaction has positive rate; every desired empty-budding outcome has probability at least 1-q(s)>0. A finite sequence of such choices has strictly positive probability despite its enormous length. No competing unwanted event is forced to occur. The resulting state has e,s>=N, B=N and C=N+1, so F<=F_N. Apply the strong Markov property and (9). This proves strictly positive explosion probability from every finite positive-substrate state.

## Scope and remaining gap

The original equation-(7) endpoint issue already has a stronger public solution; it is not newly solved here. The proposed contribution is only the interval 10/11<beta<1 inside the cited preprint's explicitly unclassified range. No absolute priority, almost-sure explosion, sharp threshold, or classification of 0<beta<=10/11 is asserted. Finite exact checks verify the stated event/count algebra and constants; the infinite-state probabilistic claims rest on this analytical argument and require fresh review.
