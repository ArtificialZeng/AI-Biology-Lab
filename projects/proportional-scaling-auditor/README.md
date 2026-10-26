# Proportional recruitment: a scoped stability counterexample

**A fresh, substantive AI mathematical and source review accepted this scoped correction on 26 October 2026.** The original authors have not confirmed it; external human peer review and historical priority are not established. This project audits only the proportional-scaling paragraph of Earn–McCluskey (2025). Its main constant-recruitment global-stability theorem is untouched. See the [review summary](review-summary.md) for the exact reviewed scope and evidence bindings.

The two-infective model family and bistability are prior work: Razvan and Yasaman (2012) already studied both. Our contribution is only the exact, paper-specific correction to the unqualified 2025 extension, not a new model or biological phenomenon. The [source comparison](sources.md) identifies that prior family and its relation to this example.

For the literal arbitrary-positive-mortality domain, the exact normalized example has a disease-free sink, a positive saddle and a second positive sink. Thus fractions have two different locally attracting outcomes. The proof uses rational arithmetic, an isolated quadratic root and exact two-dimensional Hurwitz tests. No simulation or numerical eigenvalue is decisive.

The example uses one postinfection mortality below susceptible mortality. It is a mathematical scope test, not a calibrated disease mechanism. Fraction equilibria describe balanced growth of population numbers. It does not decide the subclass in which every postinfection mortality is at least susceptible mortality. Equal mortality is a sufficient repaired assumption already covered by the main source theorem. See [statement](statement.md), [complete proof](proof.md) and [primary sources](sources.md).

## Reproduce

Python 3.12 or 3.13 with pinned SymPy 1.13.3; CPU only. Keep the environment outside this directory.

```sh
python -m pip install -r requirements.txt
python verify.py verify
python -O verify.py verify
python verify_independent.py verify
python -O verify_independent.py verify
```

The first checker invokes the preserved exact research calculation and compares its archived JSON. The second starts separately from the three population-number equations, derives fractions, eliminates one variable by a resultant and verifies both attractors and the saddle. Each ordinary check verifies hashes for every public input and byte-compares its certificate; it writes no files. Historical candidate-status text in the preserved calculation/certificate is provenance, not the final review status. These programs are not proof assistants or external human peer review.

Explicit regeneration proposes both certificates and a manifest in a new, absolute directory outside this complete project:

```sh
python verify.py regenerate --output-dir /absolute/outside/path/new-scaling-proposal
```

The parent must exist and the output must not exist. No archive is overwritten. Changed mathematics needs review; replacing hashes alone does not validate it. `compute` emits a recomputed certificate to stdout. The [MIT license](LICENSE) is directory-scoped; [third-party rights](THIRD_PARTY.md) remain separate. No third-party paper or private research process is included.
