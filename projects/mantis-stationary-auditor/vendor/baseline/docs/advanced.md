# Advanced analysis (structural robustness & analytic stochastics)

Beyond the classical Deficiency Zero / Deficiency One theorems, `mantis-delta`
provides a pipeline that, for a given reaction network, certifies global stability
and robustness from structure alone and emits the exact stochastic stationary
distribution. This page documents those capabilities with the theory and worked
examples; see the [API Reference](api.md) for full signatures.

All features degrade gracefully when their optional dependency is absent
(`pip install "mantis-delta[all]"` installs them all).

---

## Exhaustive steady states (completeness guarantee)

`CRNetwork.steady_states()` integrates the ODEs and runs multi-start least-squares —
it can **silently miss** unstable fixed points that forward integration flows away
from. At mass-action equilibrium the positive steady states are the positive real
roots of a polynomial system, so they can instead be found *exhaustively*.

`all_steady_states()` rationalises and non-dimensionalises the rate polynomials,
eliminates the conservation laws, triangularises the system with a lexicographic
Gröbner basis, recovers every complex root by back-substitution (`numpy.roots`),
polishes each with Newton's method, and returns the real non-negative ones. For a
zero-dimensional ideal this finds **all** solutions (bounded by the Bézout number).

```python
# Schlögl model: three equilibria, the middle one unstable.
rn = CRNetwork.from_string(
    ["A + 2 X <-> 3 X", "X <-> B"],
    rates={"A + 2X -> 3X": 6.0, "3X -> A + 2X": 1.0, "X -> B": 11.0, "B -> X": 6.0},
    chemostatted={"A": 1.0, "B": 1.0},
)
for ss in rn.all_steady_states({"X": 0.5}):
    print(ss.concentrations["X"], "stable" if ss.is_stable else "UNSTABLE")
# 1.0 stable / 2.0 UNSTABLE / 3.0 stable
```

The multi-start solver typically returns only the two stable branches; the
exhaustive engine recovers the unstable middle branch as well. Pass
`backend="phcpy"` to use polynomial homotopy continuation (PHCpack) for larger or
stiffer systems (`pip install "mantis-delta[homotopy]"`).

---

## Global stability via the Horn–Jackson Lyapunov function

For a **complex-balanced** network with positive equilibrium $c^*$ (the *Birch
point*), the pseudo-Helmholtz function

$$V(c) = \sum_i \left[ c_i\left(\ln\frac{c_i}{c_i^*} - 1\right) + c_i^* \right]$$

is a strict Lyapunov function on every positive compatibility class: $\partial
V/\partial c_i = \ln(c_i/c_i^*)$, hence $dV/dt = \sum_i \ln(c_i/c_i^*)\,\dot c_i \le
0$ with equality only at $c^*$. So $c^*$ is **globally** asymptotically stable — a
strictly stronger statement than local eigenvalue stability. Every weakly reversible
deficiency-zero network is complex-balanced for *all* rate constants (Horn–Jackson).

```python
rn = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1.0, "B -> A": 2.0})
rn.is_complex_balanced()                       # True
rn.complex_balanced_equilibrium({"A": 3.0})    # Birch point {'A': 2.0, 'B': 1.0}
cert = rn.certify_global_stability({"A": 3.0})
print(cert)                                    # readable certificate
cert.globally_stable                           # True
cert.lyapunov_function, cert.lyapunov_derivative   # symbolic V(c) and dV/dt
```

---

## Absolute Concentration Robustness (ACR)

A network has **ACR** in species $S$ when $S$'s steady-state concentration is
identical in *every* positive steady state, regardless of the initial totals — the
robustness property at the heart of biosensor and homeostatic design. Shinar &
Feinberg (*Science* 2010) give a structural sufficient condition: a **deficiency-one**
network with two **non-terminal** complexes that differ in a single species $S$ has
ACR in $S$.

```python
rn = CRNetwork.from_string(["X + Y -> 2 Y", "Y -> X"],
                           rates={"X + Y -> 2Y": 3.0, "Y -> X": 6.0})
res = rn.detect_acr({"X": 5.0, "Y": 5.0})
res.species          # ['X']
res.acr_values["X"]  # 2.0 = k2/k1, independent of the initial totals
```

---

## Injectivity (ruling out multistationarity beyond D1T)

An **injective** network admits at most one positive steady state per compatibility
class, so it cannot be multistationary — and unlike the Deficiency One Theorem this
works for higher-deficiency networks. Following Craciun & Feinberg (*SIAM J. Appl.
Math.* 2005/06), `is_injective()` forms the square steady-state map (rank-many
independent rows of $N\cdot v(c)$ plus the conservation totals) and tests whether its
Jacobian determinant is **sign-definite** in the positive concentrations and rate
constants. Sign-definiteness ⇒ the determinant never vanishes ⇒ injective.

```python
CRNetwork.from_string(["A <-> B"], rates={"A -> B": 1, "B -> A": 1}).is_injective().injective  # True
schlogl.is_injective().injective                                                                # False (multistationary)
```

A `True` verdict is rate-independent; `False` is inconclusive (the test is sufficient,
not necessary).

---

## Pseudo-arclength continuation (fold / Hopf detection)

The log-scan `bifurcation()` recomputes steady states *independently* at each parameter
value, so at a saddle-node **fold** — where the branch turns back on itself — it loses
the branch and never reports the fold. `continuation()` parameterises the branch by
*arclength* $s$ along the solution curve $(c(s),\lambda(s))$ and solves the augmented
system

$$
G(c,\lambda)=0,\qquad
\langle t_c, c-c_{\text{prev}}\rangle + t_\lambda(\lambda-\lambda_{\text{prev}})-\mathrm{d}s = 0,
$$

with a predictor–corrector (Keller 1977). The arclength row keeps the augmented Jacobian
regular *through* the fold, where the steady-state block $G_c$ alone is singular, so
folds are traversed smoothly. They are detected from a sign change of $\det G_c$ (refined
to $\det = 0$); **Hopf** points from a change of the unstable-eigenvalue count with no
fold (a complex pair crossing the imaginary axis).

```python
res = schlogl.continuation("X -> B", (1.0, 20.0), {"X": 0.5})
res.folds()              # both folds bounding the bistable region — (X,k) ≈ (1.35,10.72), (2.53,11.15)
res.species_branch("X")  # X along the full S-curve (arclength order, non-monotonic in λ)
```

On the Brusselator (`A`, `B` chemostatted) continuing the rate of `B + X -> Y + D`
detects the Hopf at the analytic $kB = 1 + A^2$ threshold, with a purely imaginary
critical eigenvalue.

---

## Multistationarity parameter regions (Conradi–Feliu–Mincheva–Wiuf 2017)

Injectivity answers *whether* multistationarity is possible;
`multistationarity_region()` maps *where in parameter space* it occurs (Conradi, Feliu,
Mincheva & Wiuf, *PLoS Comput. Biol.* 2017). It forms the **critical function**
$\varphi = (-1)^s \det J$ ($s=\operatorname{rank} N$) and reads it as a polynomial in the
steady-state concentrations with coefficients in the rate constants:

* every coefficient positive for all positive rates ⇒ $\varphi>0$ everywhere ⇒ the
  network is **monostationary** for *all* parameters (an injectivity certificate);
* otherwise each coefficient that can go negative is a face of the **multistationarity
  region** — the explicit "where".

```python
res = schlogl.multistationarity_region()
res.monostationary           # False — Schlögl is genuinely multistationary
res.region_conditions        # coefficient expressions whose negativity ⇒ multistationarity
res.region_boundary          # exact discriminant boundary (when reducible to one polynomial)
res.multistationary_at_rates # exact ≥2-root verdict at the supplied rates, when available
```

Reading the coefficient signs of $\varphi$ over the *whole* positive orthant is only
**necessary**: $\varphi<0$ might occur at a $c$ that is not a steady state. When the
steady states reduce to a univariate polynomial $p(x)$ — directly for
$\operatorname{rank}N=1$ networks like Schlögl, or after eliminating the conservation
laws — the test is put **on the steady-state variety** and becomes exact:
multistationarity $\iff p$ has $\ge 2$ positive roots, whose boundary is the discriminant
$\operatorname{disc}_x(p)=0$. For Schlögl this is

$$
\Delta(k) = -4k^3 + 36k^2 + 648k - 6156,
$$

zero at $k = 10.72$ and $11.15$ — exactly the folds the continuation reports — so the
bistable window is $\Delta(k) > 0$. `region_boundary` returns this discriminant and
`multistationary_at_rates` the exact positive-root count.

For a fully definitive verdict on a *specific* compatibility class (any network, no
parametrisation needed), `is_multistationary(initial_conditions)` enumerates the positive
steady states and checks for two or more:

```python
schlogl.is_multistationary({"X": 0.5})   # True at k = 11 (three steady states)
```

---

## Exact stochastic stationary distributions (Anderson–Craciun–Kurtz)

For a complex-balanced network the chemical master equation has an *exact* stationary
distribution in closed form (ACK, *Bull. Math. Biol.* 2010): a product of independent
Poissons with means $\lambda_i = N_A V c_i^*$, supported on the reachable state space.
With conservation laws it is that product conditioned on the conserved totals — a
Binomial for a single moiety. You can therefore **write down the stationary
distribution from structure alone, with no simulation**.

```python
rn = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
sd = rn.stationary_distribution({"A": 40, "B": 0}, volume_L=1e-15, initial_as="count")
sd.poisson_means()             # {'A': λ_A, 'B': λ_B}
sd.expected_counts()           # exact stationary means
vals, probs = sd.marginal("A") # exact marginal pmf (matches Binomial to machine precision)
```

---

## Finite State Projection (CME without noise)

FSP (Munsky–Khammash, *J. Chem. Phys.* 2006) truncates the count-state space to the
reachable set, assembles the sparse CME generator $Q$, and solves the master equation
directly — accurate exactly where the SSA is noisiest. It returns the transient
distribution $p(t) = e^{Qt}p_0$ with a rigorous **leaked-mass error bound**, or the
stationary distribution (null space of $Q$).

```python
rn = CRNetwork.from_string(["A <-> B"], rates={"A -> B": 2.0, "B -> A": 3.0})
stat = rn.fsp({"A": 30, "B": 0}, volume_L=1e-15, initial_as="count")        # stationary
pt   = rn.fsp({"A": 30, "B": 0}, volume_L=1e-15, t=0.1, initial_as="count") # transient
pt.expected_counts(); pt.truncation_error
```

For complex-balanced networks the FSP stationary distribution reproduces the ACK
product form to machine precision — an independent cross-check against the closed
form and the SSA.

---

## SBML interoperability

With the `[sbml]` extra, models round-trip through SBML Level 3:

```python
rn.to_sbml("model.xml")
rn2 = CRNetwork.from_sbml("model.xml")
```
