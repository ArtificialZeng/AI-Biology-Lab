# Independent numerical verification

Status: floating-point evidence only; not a rigorous periodic-orbit or hyperbolicity certificate.

The implementation in `src/independent_numeric.py` was constructed directly from the ten chemical reaction arrows in Conradi, Mincheva and Shiu, arXiv:1809.02886v2, equations (1)–(3) and Table 1. It imports no primary implementation, and no primary model code was read. Both implementations nevertheless share the published model and use SciPy; this is an independent implementation pass, not peer review.

## Setup

Species: (S0,K,S0K,S1K,S2,P,S2P,S1,S1P). Rates: (1,1,1,1,100,1,0.9,3,1,100). Totals: P=8, S=40. The explicit 9-by-10 stoichiometric matrix and three conservation rows are serialized in `evidence/independent_numerics.json`; their integer-coefficient product is zero.

The primary independent trajectories retain all nine species, rather than eliminating conservation laws. Radau starts with all substrate free (S0=40, Kfree=K, Pfree=8). BDF starts from the separate interior state (2,K-5,3,2,12,5,2,18,1). Both run to t=400. Radau uses rtol=1e-10, atol=1e-12; BDF uses rtol=2e-11, atol=2e-13. Peak events use the S2 derivative and a negative crossing.

The infinite-kinase slow model independently replaces instantaneous S0 rebinding by the flux S1P -> S0K+P; S0K dissociation returns immediately to S0K and cancels. Its six slow species are (S0K,S1K,S2,S2P,S1,S1P), with substrate conservation. Primary slow trajectories use Radau and BDF; periodic shooting and variational integration use five coordinates and DOP853.

## Periodic shooting and Floquet evidence

| K | period | S2 minimum | S2 maximum | peak-to-peak | largest nontrivial Floquet modulus |
|---:|---:|---:|---:|---:|---:|
| 100 | 4.778280750337 | 0.223009957951 | 4.137256868428 | 3.914246910477 | 0.589280116519 |
| 1000 | 4.699851454354 | 0.243469684238 | 3.916313140754 | 3.672843456516 | 0.629796977729 |
| 10000 | 4.692863922764 | 0.245548456781 | 3.894660377022 | 3.649111920241 | 0.633655052935 |
| infinity | 4.692096001431 | 0.245779663343 | 3.892262760496 | 3.646483097153 | 0.634081418609 |

Finite-K shooting uses an independent conservation reconstruction only after the full 9D trajectory has supplied a guess. A transversal phase condition is imposed using the vector field at the initial guess; it is not imposed at an S2 peak. The variational equations use complex-step Jacobians of the independently built vector field. One neutral multiplier near 1 is removed before reporting the largest nontrivial modulus.

Maximum shooting closure residual: 7.2387e-14. Maximum neutral multiplier deviation from 1: 3.9213e-13.

Multipliers at approximately 1e-17 or below are numerical roundoff-scale outputs and must not be interpreted as exact zeros, accurate signs, or rigorous bounds. In particular the serialized monodromy determinants are not reliable evidence of invertibility.

## Cross-checks and limits

Across eight trajectories, the maximum absolute conservation discrepancy was 2.5102e-10; all 20,001-point sampled trajectories were nonnegative (the free initial state contains exact zeros). The independently chosen interior initial conditions stayed positive in the sampled data. The maximum Radau/BDF period difference was 3.9614e-10; maximum peak-to-peak difference was 2.4287e-09.

For K=100,1000,10000, both period and amplitude differ from each other by much more than cross-solver discrepancies. Their convergence toward the separately implemented slow system is consistent with a leading inverse-K correction. Scaled finite differences are serialized in `evidence/independent_crosschecks.json`; three finite-K values alone neither prove an asymptotic expansion nor determine a coefficient exactly.

These data support a stable periodic response for the stated parameter example and two tested initial conditions. They do not prove periodic-orbit existence, uniqueness, global attraction, persistence for every sufficiently large K, or any assertion under different existential quantifiers in Conjecture 6.2. Any mathematical persistence theorem must retain its hyperbolicity hypothesis unless a rigorous certificate is supplied.

Reproduction command: `python src/independent_numeric.py` from the project directory. SciPy 1.15.3; NumPy 2.2.3.
