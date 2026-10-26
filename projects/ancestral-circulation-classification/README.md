# Ancestral current: N3 sign classification and N4 reversal

For the complete three-site Fleming–Viot motion graph with every motion rate one, uniform copying among the other N−1 slots, positive killing rates, and current Delta_N=F02−F20, this project proves two related results.

- **N=3, whole positive domain:** sign(Delta3)=sign(−(k0−k1)(k0−k2)(k1−k2)); the zero locus is exactly the three equal-rate planes.
- **N=4, the entire path k=(1,t,t+1), t>1:** the current changes from negative to positive at the unique simple positive root rho of an explicit degree-19 polynomial, with 4769/446<rho<3133/293. Combining this with N2 balance and the N3 theorem proves that four is the minimum population admitting reversal within one strict killing-rate ordering chamber of this fixed motion model. A regular local zero surface lies inside that chamber.

Killing-rate ratios change along the N4 path. This is not an all-positive N4 classification, fixed-ratio scale reversal, or a first general direction-change mechanism. The nearby Brownian result already establishes copying-rate dependence of ancestral direction. The source's infinite-space question, other motion graphs, and a full classification for larger populations remain outside these results.

See [statements](statement.md), [model and current](problem.md), [N3 proof](proof.md), [N4 proof](four-particle-proof.md), [primary sources and limits](sources.md), and [review status](review-summary.md). Fresh AI mathematical checks are not external human peer review, a priority certificate, proof-assistant formalization, or trusted-runner certification. No manuscript or PDF is supplied.

## Reproduce

Tested with CPython 3.13.5, SymPy 1.13.3, and python-flint 0.9.0; CPU only. Keep environments outside the project.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
python -B verify_independent.py verify
python -O -B verify_independent.py verify
python -B verify_four_particle.py verify
python -O -B verify_four_particle.py verify
python -B verify_four_particle_independent.py verify
python -O -B verify_four_particle_independent.py verify
```

The two N3 drivers check the whole three-parameter polynomial identity and positive-coefficient sign proof. The two N4 drivers check the whole unbounded path, not a finite grid. The N4 primary driver uses the 15-state count and 30-state founder systems in `four_particle_model.py`; the independent driver imports no primary mathematical code and reconstructs all 96 mixed two-colour states from victim/donor events. Same-site copies between different colours are retained. The N4 data are isolated under `evidence/four-particle/`; the N3 mathematical model, closed forms, and original certificates are preserved.

Every `verify` command checks the v2 manifest, recomputes its deterministic certificate, and leaves the archive unchanged. Each driver also accepts `compute` to emit its own certificate and `regenerate` to propose **all four** certificates plus a manifest in a new absolute directory outside the project:

```sh
python -B verify_four_particle.py compute
python -B verify.py regenerate --output-dir /absolute/outside/path/new-circulation-proposal
```

The parent directory must exist and the output must not exist. Regeneration checks the supplied closed forms; it does not silently replace them. Changed mathematics needs review as well as new hashes. The manifest covers every public project input, excluding itself and Python bytecode.

A separate standard-library Fraction solver retains fifteen rational-rate N3 regression cases and reproduces the source Appendix A tables:

```sh
python -B colour_reference.py > /absolute/outside/path/colour-reference.json
cmp evidence/colour-reference.json /absolute/outside/path/colour-reference.json
```

Finite regression cases do not establish the global theorem. See `THIRD_PARTY.md` for attribution and rights.
