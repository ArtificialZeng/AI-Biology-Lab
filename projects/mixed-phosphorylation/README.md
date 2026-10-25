# Mixed-mechanism phosphorylation: exact instability and validated attracting cycles

**Status:** exact equilibrium certificate and computer-assisted periodic-orbit/period-correction proof; novelty and external human review pending.

## Model and results

Use the Conradi–Mincheva–Shiu mixed-mechanism mass-action model with Table 1 rates `(1,1,1,1,100,1,9/10,3,1,100)`, total phosphatase 8, and substrate 40.

1. For every finite total kinase K ≥ 4660166/209569 (hence K≥23), the unique positive equilibrium has exactly two unstable and four stable eigenvalues on the six-dimensional compatibility class. The sufficient bound is not claimed sharp. [Exact proof](evidence/finite_kinase_proof.md), [certificate](evidence/finite_kinase_sign_certificate.json), [earlier equilibrium-only PDF](manuscript/main.pdf).
2. There is a K0 such that every finite K>K0 has a nearby strictly positive attracting periodic orbit within its fixed compatibility class. The limiting orbit is certified by an interval Poincare map; local normal attraction and hyperbolic-cycle persistence prove the unbounded-K assertion. K0 is not numerically computed.
3. On this local branch, in physical model time, T(K)=T∞+T1/K+O(K^-2), with T∞∈[4.69209598,4.69209603] and T1∈[7.65771717,7.68396573]. Thus T1>0. No explicit remainder constant is given. [Full new proof, code, interval output and independent audits](validated-cycle/README.md).

The old instability threshold is not a periodic-orbit threshold. The full original Conjecture 6.2, global attraction, exact scalar-output invariance and priority remain unresolved. The source paper already numerically observed oscillations; directly related 2023/2025 conference work may already establish qualitative robustness. The candidate increment is the inspectable, model-specific nonzero period correction and strict certificate. No external human referee or biological experiment is reported.

## Reproduce

The earlier rational equilibrium checks are unchanged:

```sh
python -m pip install -r requirements-exact.txt
python scripts/reproduce.py verify
```

The new interval integration uses the separate [C++/CAPD reproduction route](validated-cycle/README.md), including Linux x86/FILIB runs and independent exact postprocessing. A successful equilibrium-only Python run does not by itself rerun the periodic-orbit integration.

The earlier manuscript and archived derivation/evidence are retained. Historical `theory.md` wording predates the new cycle certificate; this README, `claim.json` and the new validated-cycle proof define the current scope.
