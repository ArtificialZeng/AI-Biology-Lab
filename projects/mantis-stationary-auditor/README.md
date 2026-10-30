# Finite conditioned stationary-law support in Mantis

This reproducible engineering audit fixes two errors in the finite, fully enumerated initial-state-conditioned stationary-distribution API of [Mantis](https://github.com/EmilioVenegas/mantis-delta), pinned at `0423b4ce556e987a2367fd293ca8db8ae0479927`.

The original implementation can assign positive probability outside its declared finite support, and can normalize over states that satisfy linear conservation laws but belong to different communicating classes. For `2A <-> 2B`, initial counts `(2,0)`, the unreachable state `(1,1)` receives probability `1/2`. The corrected law assigns probability `1/2` to each reachable endpoint. The variance of A changes from `1/2` to `1`; its mean remains `1`.

This is a software contract correction using established theory, not a new biological finding, a new mathematical theorem, or a contradiction of the Anderson–Craciun–Kurtz result. The original law in this example is a valid stationary mixture; exact stationarity alone does not identify the initial communicating class.

## Reproduce offline after installing dependencies

Python 3.12 is recommended because it is an upstream CI target. The original local run used Python 3.13.5. After `python -m pip install -r requirements.txt`, run:

```sh
python -B verify.py verify
python -B -O verify.py verify
```

These commands check every archived input and certificate, import the actual vendored baseline/candidate `mantis` classes, and replay the independent exact CTMC checks. They do not install Mantis, use a GPU, contact a network service, or write inside this project. Runtime caches are placed in temporary directories outside it. Both commands produce identical deterministic JSON. Raw measured floating-point errors are excluded from default certificates; the archived certificate records exact references and explicit comparison tolerances.

The first oracle supplies 352 helper/product-law cases and five public API cases. A separate fresh-context reviewer supplied an independently implemented oracle covering 990 two-species and 87 multispecies/catalyst cases. These finite checks support the bounded patch; they are not an infinite-family theorem or external human peer review.

## Optional upstream tests

From this project directory, on a writable checkout, with cache directories outside the archive:

```sh
PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=/tmp/mantis-numba PYTHONPATH=vendor/candidate python -B -m pytest -p no:cacheprovider vendor/candidate/tests/test_stochastic_stationary.py vendor/candidate/tests/test_fsp.py vendor/candidate/tests/test_stochastic.py vendor/candidate/tests/test_chemostatted.py -q
PYTHONDONTWRITEBYTECODE=1 NUMBA_CACHE_DIR=/tmp/mantis-numba PYTHONPATH=vendor/candidate python -B -m pytest -p no:cacheprovider vendor/candidate/tests -q
```

To expose the four original regression failures, use `PYTHONPATH=vendor/baseline` with the candidate stationary tests and the pytest selection `-k "outside_support or discrete_communicating_class or below_required_count"`. The candidate tests are deliberately exercised against the actual pinned baseline package.

The related suite passed 37 tests in the original local run and in the packaged replay. The complete local candidate suite had **166 passed, 1 skipped, 1 failed**. The unchanged failure `test_cha_unique_positive_state_matches_integration` was reproduced twice on pinned upstream and twice on the candidate under Python 3.13.5 / SciPy 1.15.3. It compares concentrations near `1e-12` with an absolute tolerance of `1e-15`; its call-path implementations and test were unchanged. This is not an all-green full-suite claim. Optional `python-libsbml` / `numba` availability can change skip counts. Upstream CI tests Python 3.10–3.12, so its result must be checked separately.

## Archive contents

- `vendor/baseline/`: pinned original author package, tests and relevant docs.
- `vendor/candidate/`: the same files with the bounded patch.
- `upstream.patch`: five-file upstream patch, including documentation and regressions.
- `scripts/`: two independent Fraction-based finite CTMC oracles, actual-library audits.
- `certificates/`: deterministic replay outputs; `archive_manifest.json` hashes inputs and certificates.
- `SOURCE.md`, `SCOPE.md`, `THIRD_PARTY.md`: provenance, exact scope and attribution.

Hashes detect archive corruption; they are not a trusted runner attestation. Maintainer confirmation remains pending. No private contact or internal research log is included.
