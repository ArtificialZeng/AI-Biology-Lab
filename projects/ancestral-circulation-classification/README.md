# Three-particle ancestral circulation: complete sign classification

**Fresh AI mathematical review found no gap; external human review and historical priority remain pending.** For exactly three Fleming–Viot particles on a symmetric three-site graph with unit motion rates and arbitrary positive killing rates, the theorem classifies the direction and complete zero locus of the surviving-label current.

With current orientation Delta=F02−F20, its sign is the sign of −(k0−k1)(k0−k2)(k1−k2). It vanishes exactly when two killing rates coincide. The source already proves generic nonzero circulation and gives one exact witness; the added target is the complete three-parameter classification.

See [model and current](problem.md), [full proof](proof.md), [sources](sources.md) and [review status](review-summary.md) for assumptions, prior work and limits. This does not classify all population sizes, solve infinite-space nonexplosion, or establish a biological experiment. Historical priority and external human review remain pending.

For rates (1−e,1,1+e), the new N=3 result gives a positive quintic response. Combined with the source’s negative N=4 cubic response, the two population sizes circulate in opposite directions for sufficiently small positive e. No numerical threshold or all-rate N=4 classification is claimed.

## Reproduce

Tested with CPython 3.13, SymPy 1.13.3 and python-flint 0.9.0; CPU only. Keep environments outside the project. The pinned native dependency provides CPython 3.13 Linux and macOS wheels.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
python -B verify_independent.py verify
python -O -B verify_independent.py verify
```

The primary checker verifies polynomial residuals for the stationary and founder closed forms, the exact current factorization, and its global positivity certificate. The second checker independently constructs the two-colour victim/donor process and verifies all 36 mixed-state Dirichlet identities over the whole polynomial domain. It shares candidate closed-form data, and imports no primary model or checker code. These are exact finite identities, not stochastic simulations or a sampled proof of the positive domain.

The archive manifest binds every public project input; each verify command recomputes its certificate without modifying the archive. `compute` emits a verification certificate to stdout. Explicit regeneration proposes both certificates and a manifest only in a new absolute directory outside this project:

```sh
python -B verify.py regenerate --output-dir /absolute/outside/path/new-circulation-proposal
```

The parent directory must exist and the output must not exist. Regeneration checks the supplied closed forms; it does not silently rediscover or replace them. Changed mathematics needs review as well as new hashes. Third-party rights are separate from the directory-scoped MIT license.

An additional standard-library Fraction reference solver checks fifteen rational-rate cases and reproduces every source Appendix A table:

```sh
python -B colour_reference.py > /absolute/outside/path/colour-reference.json
cmp evidence/colour-reference.json /absolute/outside/path/colour-reference.json
```

These finite cases are regression controls; the polynomial identities and positivity proof establish the whole-domain theorem. No external human review, proof-assistant formalization or trusted-runner certification is claimed. The archived proof records its research-stage status; the review summary records subsequent checks.
