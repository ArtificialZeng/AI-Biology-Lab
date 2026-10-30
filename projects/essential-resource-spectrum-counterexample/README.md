# An exact counterexample for essential-resource exchange

An explicit three-consumer network in the published Gibbs et al. (2022) model has an effective interaction matrix B with three real, strictly negative eigenvalues, yet its positive coexistence equilibrium is unstable at low resource supply. The same fixed network is locally stable above an exact positive supply threshold. This refutes the paper's proposed real-negative-spectrum protection conjecture in its equal-abundance analytical setting. It does not address the separate sufficient theorem requiring symmetric matrices.

For every n>0 the equilibrium is R=(1,1,1), N=(n,n,n), with mortality 4/39 and strictly positive inflow

```
rho = n * (1043/90, 16697/2970, 4142/1485).
nc = 1639836184859336 / 29522765162582466287
```

B is unchanged as n varies. The full Jacobian has exactly two eigenvalues in the open right half-plane for 0<n<nc and is Hurwitz for n>nc. The critical point n=nc is excluded from these strict claims. [proof.md](proof.md) gives the model, exact calculation, and a two-consumer obstruction limited to the same canonical normalization.

Install Python 3.10 or later and the pinned dependency, then run from this directory:

```sh
python -m pip install -r requirements.txt
python verify.py
python -O verify.py
```

The verifier uses rational arithmetic, differentiates the original component equations on the proven strict limiting-resource branch, verifies positive equilibrium and biomass balances, computes the full six-dimensional characteristic polynomial, cross-checks a quadratic matrix pencil, and certifies the entire two open intervals by Routh determinants. It compares the resulting certificate with [certificate.json](certificate.json) and validates every archived input hash. Ordinary verification writes no files and requires no network, API, GPU, R, or Julia. The `PASS` output means these explicit checks succeeded; it is not formal proof-assistant certification or external peer review. The dependency installation itself may use the network.

The archive manifest covers the eight listed inputs and excludes itself to avoid a circular hash. Changed mathematical inputs require a fresh derivation and review of witness, proof, certificate and manifest together; merely updating hashes does not validate a changed claim. This packet contains no regeneration command that overwrites archived evidence.

This is an independent reconstruction of the published two-pathway equations, not a bug patch to the authors' current simplified R implementation. All included code and prose are original; the MIT license applies to this directory's files. Third-party papers and author code are linked, not redistributed.

The contribution is a bounded model-specific correction, not a new general principle about block spectra. Meacock, Amicone and Mitri (2026, v2) previously analyzed full population/environment block spectra and the role of environmental damping; Liu et al. (2025) studied destabilization from mismatched resource consumption and growth. [sources.md](sources.md) explains their relationship and links the primary records. A bounded primary-source comparison located no equivalent native Gibbs example; priority has not been established exhaustively.

The model is dimensionless and uncalibrated, and all three resources are supplied externally at a fixed ratio. The result concerns local equilibrium stability only: no periodic orbit, extinction, global persistence, experimental outcome, or general minimum-dimension claim follows. Development and checking were AI assisted. Independent human confirmation and external peer review remain pending; no author endorsement is implied.
