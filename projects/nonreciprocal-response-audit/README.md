# Finite response jets and full biochemical response curves

This code/proof note audits three subsidiary assertions in Section8 of
Franco, Kepka and Velázquez, ARMA250:25(2026), and supplies a limited local
repair. It concerns four-state column Markov generators with strictly positive
off-diagonal rates and the displayed forward curves e2^T exp(tA)e1 and
e3^T exp(tA)e1.

| Exact finding | Scope |
|---|---|
| The printed symmetric six-rate derivative at J-4I has rank2 | Includes the required diagonal derivatives; does not describe the full nine-dimensional reversible tangent |
| At rates(1,2,3,4,5,6), its determinant is800 | Repairs existence of an open finite-jet image locally at that point |
| A reversible and a nonreversible positive generator share all state derivatives through order3 but differ at order4 | Six low-order coefficients do not determine complete response curves |
| A positive complete four-state network has forward curves with uncancelled poles -5±i | The classical real-pole restriction excludes every finite reversible realization of those exact curves |

The source-specific contribution is the exact scope correction, positive
witnesses and finite-jet repair. The paper's main cut-vertex theorem is outside
this audit. Rank failure at one point alone would not prove an image has no
interior. Complete idealized curves, six derivatives, and finite noisy
measurements are different data; this is not a practical noisy-data classifier.

The reversible spectral theorem and oscillatory evidence for nonequilibrium
are established prior knowledge, including the source's Proposition2.2,
[Qian–Elson2004](https://doi.org/10.1073/pnas.0305962101) and
[Amann–Schmiedl–Seifert2010](https://doi.org/10.1063/1.3294567).
No new oscillation criterion, general identification theorem or biological
discovery is claimed. See [literature.md](literature.md) for observable
distinctions and the bounded correction/priority search.

Read [proof.md](proof.md), [source.md](source.md) and [problem.md](problem.md).
All witnesses have positive rates and obey the source's complete-graph domain.

## Reproduce

Python3.12 or newer, with pinned exact-algebra dependencies:

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
python -B verify.py regenerate --output-dir /absolute/new/external-directory
```

The wrapper verifies complete inventory/hash binding, runs the owner and
separate checker outside the read-only input tree, and compares their full
deterministic certificates. Worker optimization is forced to0. Regeneration
requires an exclusive new external absolute directory and creates a complete
replayable packet. Exact checks cover derivative matrices, determinant800,
positive moment twins, cycle imbalance, coprime resolvents, eighteen
directional identities and six owner mutation controls. General mathematical
steps are justified in the written argument, not certified by finite tests.

The work used substantial AI assistance. A separate context accepted the
mathematical scope for a modest source-specific code/proof note. Its overall
trusted-runner review remains `revise` because genuine job provenance was not
supplied; the actual formal gate failure is preserved. No human author
confirmation, priority certificate, proof-assistant certificate or formal
release PASS is claimed. No manuscript/PDF, third-party paper or private
correspondence is distributed. No GPU, API or paid computation is needed.
