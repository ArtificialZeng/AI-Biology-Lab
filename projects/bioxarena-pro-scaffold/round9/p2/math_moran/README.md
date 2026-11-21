# Exact consensus time of the neutral Moran process on the star graph

This directory contains a complete analytic proof and dual independent computational checks resolving Conjecture (9) of Brewster, Svoboda, Roscow, Chatterjee, Tkadlec, and Nowak (*Maintaining diversity in structured populations*, *PNAS Nexus* 4(8):pgaf252, 2025; `arXiv:2503.09841v2`).

## Main results

Let $S_n$ be the simple undirected star graph on $N = n + 1$ vertices ($n \ge 1$ leaves of degree $1$, $1$ central vertex of degree $n$, no self-loops), initialized with $N$ pairwise distinct types. Under neutral discrete-time updating where each step selects directed edge $(0 \to i)$ with probability $p_{\mathrm{in}}$ and $(i \to 0)$ with probability $p_{\mathrm{out}}$ ($p_{\mathrm{in}} > 0$, $p_{\mathrm{out}} \ge 0$, $n(p_{\mathrm{in}} + p_{\mathrm{out}}) = 1$), the expected number of steps until consensus (homogeneity) is exactly

$$T_N(p_{\mathrm{in}}, p_{\mathrm{out}}) = n H_n + n(n - 1)\frac{p_{\mathrm{out}}}{p_{\mathrm{in}}}, \qquad H_n = \sum_{k=1}^n \frac{1}{k}.$$

Specializing to the standard update rules yields:

* **Birth-death (`bd`, Conjecture Eq. (9)):** $p_{\mathrm{in}} = \frac{1}{Nn}$, $p_{\mathrm{out}} = \frac{1}{N}$, giving
  $$T_N^{\mathrm{bd}} = n^3 - n^2 + n H_n.$$
* **Death-birth (`db`):** $p_{\mathrm{in}} = \frac{1}{N}$, $p_{\mathrm{out}} = \frac{1}{Nn}$, giving
  $$T_N^{\mathrm{db}} = n H_n + n - 1 = (N - 1) H_{N-1} + N - 2.$$
* **$\lambda$-mixed ($\lambda \in [0, 1]$):** $\lambda$ fraction `bd`, $(1-\lambda)$ fraction `db`, giving
  $$T_N^{(\lambda)} = n H_n + n(n - 1)\frac{n\lambda + 1 - \lambda}{\lambda + n(1 - \lambda)}, \qquad T_N^{(1/2)} = n^2 - n + n H_n.$$

See [proof.md](proof.md) for the complete coalescing-walk duality and recurrence solution, [REPORT.md](REPORT.md) for the task-by-task audit report, and [sources.md](sources.md) for primary source definitions and bounded literature checks.

## Reproduce checks

The exact rational Markov solvers (`checks/exact_markov_checker.py`) use only the Python 3 standard library (`fractions`). The symbolic checker (`checks/independent_mc_and_symbolic.py`) uses SymPy (tested with `sympy==1.13.3`). Run from `round9/p2/math_moran/`:

```sh
python3 -B checks/verify.py
python3 -B -O checks/verify.py
```

Both normal and `-O` modes retain all explicit `require(...)` checks and verify live results against [checks/verification_summary.json](checks/verification_summary.json):

* **Checker 1 (`checks/exact_markov_checker.py`):** Exact $\mathbb{Q}$-arithmetic linear solves on (i) the unreduced $N^N$ vertex-labeled forward chain ($n=1,2,3$), (ii) the type-isomorphism-only forward chain ($n=1\dots 5$), (iii) the full $(2^{n+1}-1)$-subset backward coalescing chain with explicit strong-lumpability checks ($n=1\dots 6$), and (iv) the $(2n-1)$-state $(k, c)$ chain ($n=1\dots 20$) across `bd`, `db`, and $\lambda$-mixed rules.
* **Checker 2 (`checks/independent_mc_and_symbolic.py`):** Independent SymPy symbolic identity verification for arbitrary $n, k, r, \lambda$, plus literal forward-time Monte Carlo simulations ($17$ suites across `bd`, `db`, and $\lambda=1/2$).
* **Negative controls (`checks/negative_controls.py`):** $8$ adversarial controls confirming both checkers reject wrong formulas, self-loop graph perturbations, corrupted transition kernels, and swapped update rules.

## Scope and status

The theorem holds for all $n \ge 1$ ($N = n+1 \ge 2$) on simple connected star graphs without self-loops under neutral fitness and discrete single-replacement steps. Adding self-loops alters the exact value (explicitly demonstrated in `checks/negative_controls.py`).

AI assistance was used for literature retrieval, proof derivation, implementation, and verification. Bounded API searches through 2026-11-21 across arXiv, OpenAlex, Crossref, and Europe PMC do not establish universal priority, and no external human peer review or formal proof-assistant certification has been obtained. No third-party paper PDF or author code is redistributed.
