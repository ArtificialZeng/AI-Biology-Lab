# A sharp count-coupling class for mixed Moran updating

Shawn Zeng (Zijian Zeng), AI Biology Lab, UCSI University.

Proposed subsidiary proof/code note, 28 October 2026. Author validation and priority confirmation have not been obtained.

## Model and claim

Let G be the complete multipartite graph with q>=2 nonempty parts of integer sizes n_1,...,n_q and total N. There are no self-loops, weights, mutation, or type-dependent edges. The native discrete process independently chooses Bd with fixed probability lambda in[0,1], or dB with probability1-lambda. Residents have fitness1 and mutants constant fitness r>0. Poissonize these same steps at rate1; this leaves fixation probabilities and the update mixture unchanged.

Let X_i be the mutant count in part i. On the finite product of integer intervals, use coordinatewise order. Attractivity means that any two ordered initial count states of the same parameter process admit a pathwise order-preserving coupling for the native generator L=P-I in common physical time (equivalently its semigroup preserves stochastic order). The coupled copies need not share one prescribed Poisson attempt schedule.

For fixed G and lambda, the following are equivalent:

(A) This count process is attractive for every r>0.
(B) lambda=0, or q=2, or all part sizes are equal.

Moreover, when q=2 or all parts are equal, for every fixed lambda, all0<r1<=r2, and initial sets S1 subset S2, fixation satisfies
rho(G,lambda,r1,S1)<=rho(G,lambda,r2,S2).
The stronger sufficient condition is dominance of their per-part mutant counts. The lambda=0 comparison is a known pure-dB case and also follows from this argument.

For unequal parts with q>=3 and lambda>0, the proof gives a same-clock count-attractivity obstruction at large finite r, NOT a reversal of fixation probability. The historical general-graph joint conjecture remains unresolved. Alternative time changes or proof methods are not excluded.

## 1. Exact native lumping

Vertices inside each part are exchangeable and every vertex in part j has degree d_j=N-n_j. Therefore the count vector is a Markov lumping of the exact binary-state native process. Put

K=sum_j k_j, W=N+(r-1)K,
t_i=sum_(j!=i) k_j, F_i=d_i+(r-1)t_i.

The rate for k -> k+e_i is

b_i(k;r)=(n_i-k_i)[lambda*r/W*sum_(j!=i) k_j/d_j
                 +(1-lambda)/N*r*t_i/F_i].

The rate for k -> k-e_i is

a_i(k;r)=k_i[lambda/W*sum_(j!=i)(n_j-k_j)/d_j
            +(1-lambda)/N*(d_i-t_i)/F_i].

These are the original one-step change probabilities; the Poisson generator uses them as off-diagonal rates and absorbs all type-preserving steps into its diagonal. The Bd term sums parent probability times the number of eligible targets in part i divided by the parent's degree. The dB term first selects an eligible target uniformly, then a neighboring parent with fitness-proportional probability. All W,F_i,d_j are strictly positive, including at boundaries, for every finite r>0. Empty and full count vectors are absorbing.

## 2. Sufficient order inequalities

If q=2 or all parts have equal sizes, for every target i all outside parts have the same degree, call it c_i. Thus both sums above reduce to t_i/c_i and (d_i-t_i)/c_i.

Fix k_i=u. Write H=n_i+(r-1)u>0, so W=H+d_i+(r-1)t with t=t_i. Increasing other counts increases t. The needed scalar derivatives are

d/dt [r*t/W] = r*(H+d_i)/W^2 >0,
d/dt [(d_i-t)/W] = -(H+r*d_i)/W^2 <0,
d/dt [r*t/F_i] = r*d_i/F_i^2 >0,
d/dt [(d_i-t)/F_i] = -r*d_i/F_i^2 <0.

The factors n_i-u,u,lambda/c_i,(1-lambda)/N are nonnegative. Consequently, whenever x<=y and x_i=y_i,
b_i(x;r)<=b_i(y;r), a_i(x;r)>=a_i(y;r).
At coordinate boundaries, weak inequalities are all that is required.

For lambda=0 the terms involving c_i vanish; the last two derivatives prove the same inequalities for EVERY multipartite graph, without equal part sizes.

At a fixed count vector, births are nondecreasing and deaths nonincreasing in r. Indeed

d/dr [r/W]=(N-K)/W^2 >=0, d/dr [1/W]=-K/W^2 <=0,
d/dr [r/F_i]=(d_i-t_i)/F_i^2 >=0,
d/dr [1/F_i]=-t_i/F_i^2 <=0.

These fitness derivatives apply even before the special graph restriction.

## 3. Constructive coupling and fixation ordering

For two count states x<=y and fitnesses r1<=r2 in one of the sufficient classes, the preceding inequalities imply, whenever x_i=y_i,

b_i(x;r1)<=b_i(y;r2), a_i(x;r1)>=a_i(y;r2).

Construct a continuous-time pair chain. For each coordinate and each jump sign, move both coordinates together at the minimum of the two corresponding marginal rates; assign the remaining rate to a move in the appropriate marginal alone. This gives the exact marginal generators. All rates are bounded because the state space is finite. A simultaneous same-sign move preserves order. At equal coordinates, a lower-only birth and an upper-only death have zero excess rate. At unequal integer coordinates, a single one-step move can meet but cannot cross. Hence the pair stays ordered at all times.

Setting r1=r2 proves count attractivity in the sufficient classes. For joint fixation ordering, both marginal processes absorb almost surely: in every nonabsorbing configuration, connectedness and positive fitness supply a finite positive-probability sequence to either monochromatic configuration. Finiteness gives the absence of another closed class. If the lower count process reaches the full vector, the ordered upper process is already full. Thus the event of lower fixation is contained in the upper fixation event.

Within-part exchangeability identifies these count-process fixation probabilities with those for any original vertex sets having the given counts. In particular S1 subset S2 gives ordered counts. This proves the stated joint fitness/initial-set comparison, for all positive fitnesses and all fixed mixtures in the sufficient graph classes.

This construction uses the exact rate-one Poissonized process, not a state-dependent clock that would alter the fixed update mixture.

## 4. Necessity of the count-attractivity class

Suppose q>=3 and not all part sizes are equal. Choose distinct parts j,l with d_j>d_l, and a third part i. Start at x with one mutant in part l and none elsewhere, and y=x+e_j. Both states have zero mutants in target part i. They are ordered and feasible because every part is nonempty.

The difference in target-i birth rates for the pure Bd component is

D_Bd(r)=n_i*r*[N/d_j+(r-1)(1/d_j-1/d_l)]
        /[(N+r-1)(N+2(r-1))].

Its limit as r -> infinity is
(n_i/2)*(1/d_j-1/d_l)<0.

The pure dB difference is

D_dB(r)=n_i*r*d_i/[N*(d_i+r-1)*(d_i+2(r-1))]>0,

and its limit is0. For every fixed lambda>0,

lim_(r->infinity) [b_i(y;r)-b_i(x;r)]
 =lambda*(n_i/2)*(1/d_j-1/d_l)<0.

Thus the difference is strictly negative at some finite r>1. For the increasing test function f(k)=1 if k_i>=1 and0 otherwise, f(x)=f(y)=0 but Lf(x)=b_i(x;r)>b_i(y;r)=Lf(y). The first-order small-time expectations violate stochastic order. No order-preserving coupling of these same-clock count-process laws for all ordered starts is possible at this r. This proves necessity of(B) for(A).

The obstruction concerns the transition law, NOT the all-mutant absorption event. It gives no fixation-order counterexample on unequal multipartite graphs.

## Scope and prior boundary

The old ITCS manuscript explicitly proposes general joint monotonicity; the later distinct journal manuscript proves a neutral-versus-advantage special comparison and supplies other bounds/formulas. This candidate only resolves the natural subsidiary count-coupling classification and resulting complete-bipartite/equal-part fixation comparison. General graphs, unequal multipartite fixation order, weights, directed edges, mutation, time-varying mixtures and alternative clocks are not resolved. Existing pure-rule results and source star/cycle formulas are not claimed as new.

No finite parameter scan proves this theorem. Exact native kernel examples may corroborate the algebra; a separate fresh AI context reconstructed the argument and checked symbolic identities and finite native examples. These are not external human peer review, formal proof-assistant certification, or priority clearance.
