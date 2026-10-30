# Directed-cycle identifiability: checking an output-leak endpoint

A positive-integer counterexample to the endpoint case of Conjecture 4.6 in Ahmed et al., *Identifiability of Directed-Cycle and Catenary Linear Compartmental Models*, SIAM J. Applied Dynamical Systems 25(1), 304–350 (2026), [DOI](https://doi.org/10.1137/25M1728636), [arXiv v3](https://arxiv.org/html/2412.05283v3#S4.SS1).

For the three-node directed cycle with input 1, output 2, and sole leak 2, take

```
(k21, k32, k13, k02) = (2, 1, 3, 1).
```

The point lies on the conjectured singular hyperplane `k32+k02=k21`, but a 4×4 coefficient-Jacobian minor is `−k21² k13 = −12`. The Jacobian therefore has full column rank. An explicit rational inverse also recovers all four parameters throughout the positive domain of this particular model. This contradicts the stated hyperplane inclusion; the paper's main generic-identifiability theorem is unaffected.

The [proof](PROOF.md) additionally proves the restricted inclusion for a leaking vertex **strictly before the output** (`ell<p`) using a polynomial tangent-kernel identity. It does not characterize the complete singular locus, or assert that all output-leak cases fail.

## Reproduce

From the repository root, with Python 3.12 or 3.13:

```sh
python -m pip install -r projects/directed-cycle-identifiability/requirements.txt
python projects/directed-cycle-identifiability/verify.py --self-test
python -O projects/directed-cycle-identifiability/verify.py --self-test
python projects/directed-cycle-identifiability/independent_reconstruction.py
```

The exact SymPy checker reconstructs the transfer function from the ODE matrix and by elimination, verifies the minor and inverse, and checks the interior kernel in 28 finite cases. Those cases supplement the general written proof. Three deliberately invalid claims must be rejected. The script prints a receipt; `--output PATH` writes a separate receipt without changing archived evidence.

[Source and hypothesis audit](SOURCE.md) · [proof](PROOF.md) · [verification and limits](VERIFICATION.md) · [archived exact receipt](evidence/certificate.json) · [credit](THIRD_PARTY.md)

## Status

Prepared with AI assistance for derivation, coding, source checking, and a separate skeptical review. Computational checks and a fresh AI review are not external human peer review. Original-author confirmation and priority clarification are pending. This is a conjecture-boundary correction, with no experimental biological claim or promise of publication or media coverage.
