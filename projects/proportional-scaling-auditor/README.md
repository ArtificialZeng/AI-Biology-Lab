# Proportional recruitment: a scoped stability counterexample

This project corrects only the unqualified proportional-scaling paragraph of Earn–McCluskey (2025). Its main constant-recruitment global-stability theorem is untouched. A fresh-context AI mathematical/source review accepted the original witness on 26 October 2026. A separate fresh-context AI mathematical/source review accepted the ordered-mortality extension on 28 October 2026 local time (27 October UTC), within the exact scope below. This is a scoped mathematical review, not a formal runner-certified release pass. No external human peer review, public author endorsement or historical priority is claimed. See the [review summary](review-summary.md).

The two-infective model family and bistability are prior work: Razvan and Yasaman (2012) already studied both. The contribution here is an exact, paper-specific correction and a strengthening of its demonstrated parameter scope, not a new model, epidemic mechanism or general dynamics theorem. The [source comparison](sources.md) explains the overlap and the limits of the prior-work check.

The [original proof](proof.md) gives a disease-free sink, a positive saddle and another positive sink in the literal arbitrary-positive-mortality class. Its witness has one postinfection mortality below susceptible mortality. The [ordered-mortality addendum](natural-mortality.md) removes that limitation: a simple exact model with every postinfection mortality at least susceptible mortality also has two sinks and a positive saddle. A separate fully strict model has every postinfection mortality greater than susceptible mortality, both infectivities and both transfer directions positive, and an exact positive saddle. A second positive sink is **not** claimed for that strict model. These witnesses have positive total-population growth; the original witness already did too.

All stability statements concern fractions; their equilibria are balanced-growth states of population numbers. No calibrated disease prediction or oscillation claim is made. Common mortality remains a sufficient repaired assumption. See [statement](statement.md), [proof](proof.md), [addendum](natural-mortality.md) and [primary sources](sources.md).

## Reproduce

Python 3.12 or 3.13 with pinned SymPy 1.13.3; CPU only. Keep the environment outside this directory.

```sh
python -m pip install -r requirements.txt
python verify.py verify
python -O verify.py verify
python verify_independent.py verify
python -O verify_independent.py verify
python verify_natural.py verify
python -O verify_natural.py verify
```

The first checker invokes the preserved original exact calculation. The second starts separately from the population-number equations and uses exact elimination for the original witness. Their code and certificates are unchanged. The third uses standard-library rational dual numbers to differentiate the original population equations for the ordered-mortality examples, checks balanced growth and compares the resulting Jacobian with the explicit quotient formula. Its gates use explicit exceptions and remain active under `-O`.

Every ordinary verification checks hashes for every public input and byte-compares its recomputed certificate. It writes no files. Historical candidate-status text inside the preserved original calculation/certificate is provenance, not the current review status. Arithmetic checks are neither proof assistants nor external human peer review.

Explicit regeneration proposes **all three** certificates and a manifest in a new absolute directory outside this complete project:

```sh
python verify.py regenerate --output-dir /absolute/outside/path/new-scaling-proposal
```

The parent must exist and the output must not exist. No archive is overwritten. Changed mathematics needs review; replacing hashes alone does not validate it. `compute` emits a recomputed certificate to stdout. The [MIT license](LICENSE) is directory-scoped; [third-party rights](THIRD_PARTY.md) remain separate. No third-party paper or private correspondence is included.
