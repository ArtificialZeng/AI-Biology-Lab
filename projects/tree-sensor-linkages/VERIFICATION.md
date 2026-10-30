# Reproduction and evidence limits

From the repository root:

```
python projects/tree-sensor-linkages/verify.py
```

Dependencies: Python 3.10+, NetworkX, SymPy. The recorded run used Python 3.13.5, NetworkX 3.4.2, SymPy 1.13.3. The CLI interpreter is an environment detail, not part of the mathematical claim. For a portable installation, invoke the script with a Python interpreter supplying those packages.

The checker uses a brute-force dynamic program over vertex masks to maximize disjoint paths. This is independent of the harmonic linear solve. It forms each physical A with its true dependent diagonals, and computes sensitivities of the Markov parameters (A^j)[o,r], j=1,...,2n-1, recursively from the individual edge and leak directions. It does not infer physical rank from harmonic rank.

Why these moments suffice for tangent-rank checking: F_o=N_o/D has denominator degree n. Its parameter derivative has denominator D^2, leading coefficient 1, and numerator degree at most 2n-2 (the z^-1 coefficient is parameter-independent). Consequently its first 2n-1 potentially varying Laurent coefficients determine that derivative. The Jacobian of the listed moments and the complete transfer-data differential have the same kernel. This remains a finite exact algebraic calculation.

All linear algebra is modular at prime 1,000,003. Nonzero modular minors certify corresponding nonzero rational minors at the recorded positive integer parameter points, when the Dirichlet denominator is invertible modulo the prime. The script fails explicitly on a singular Dirichlet solve; none occurred. Sampled deficient rank is not treated as proof of a generic upper bound. Four deterministic seeds are allowed for a physical lower-rank witness; two retries were used in the recorded campaign.

Coverage: all nonempty leak and observation sets for every representative of each unlabelled tree on 1–7 vertices, giving 204,625 harmonic placements. Every root retained in every observation set for all trees on 1–5 vertices, giving 8,497 physical Jacobians. Relabelling does not change the theorem, but the parameter assignment is one deterministic positive integer witness per case, with retries when necessary; this is not exhaustive parameter coverage. A held-out attack selected 100 random labelled trees on 6–16 vertices with independent root/leak/observation placements. The random sequence and all parameter seeds are fixed in verify.py.

The rational seven-vertex branching-tree witness serializes H exactly, includes four leaks, and exhibits ambiguity dimensions 3,2,1,0 as sensors are added. Its physical-rank labels follow the written theorem and are explicitly marked as such. The independent physical Jacobian campaign is the separate lower-bound check. The final all-leak observation row retains root 1; no unobserved-root physical rank is asserted.

The infinite theorem is established by PROOF.md, with a separate AI audit; external human review remains pending. These computations are counterexample searches and finite witnesses, not formal verification, external peer review, or historical novelty certification. No PDF/manuscript is created in this research phase.

The public verifier uses explicit exceptions, so `python -O` retains every check. It emits deterministic JSON to stdout by default and writes evidence only when an output directory is explicitly requested. Ordinary replay never regenerates or alters the archived evidence.
