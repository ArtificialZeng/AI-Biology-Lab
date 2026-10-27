# msprime terminal hazard auditor

At a demographic jump placed at the final requested time, pinned msprime main can return the preceding instantaneous coalescence rate. With one population changing from size1000 to10000 at time100, a query ending at100 returns0.0005 instead of0.00005; appending101 changes the shared answer at100. The included two-file patch evaluates the terminal time through the existing interior path, including the existing positive-time lineage movements.

This package contains the exact frozen patch, unchanged neutral regression and analytic-oracle inputs, source hashes, native bootstrap and paired replay. Source contracts are in source.md. It prepares msprime commit `45cfa59d828757352f7ab1a55c19b55102d40d22` and the exact tskit gitlink `f1b139edb376fbffb2978a5da2c9e4b948930b0c`; downloaded archives, build trees and native binaries stay outside the public project.

## Actual local evidence

Separately downloaded and compiled baseline/candidate source trees were replayed using CPython3.13.5 on macOS26.5 arm64 and an existing Python/GSL environment:

| Replay | Baseline | Candidate |
| --- | --- | --- |
| Frozen new regressions | 20 failed, 10 passed | 30 passed |
| Independent analytic cases | 42 failed of68 | 68 passed |
| Existing trajectory cases | 16 passed | 16 passed |

`verify_variant.py` requires these exact counts, zero skips and no missing calls. An expected baseline pytest exit1 is accepted only after its30 cases give exactly20 failures and10 passes. Errors, extra skips, missing cases, changed scientific-input hashes or a different failure count fail the wrapper. Full stdout/stderr/JSON reports are written to the external output directory. The analytic oracle is preserved byte-for-byte. Its environment-specific expected_sources.json is generated beside a work copy after the immutable source hash is verified, binding the actual newly built extension hash without pretending binary hashes are platform independent.

Compact host evidence and retained failed-attempt descriptions are in evidence/. The original two candidate failures with overly strict survival tolerance remain documented: upstream's five-term Taylor approximation gave a1.3915e-11 mixture-survival difference, below the2.0792e-11 analytic remainder bound. Both variants use the same survival rtol2e-12/atol5e-11, while rate checks retain rtol2e-12/atol2e-15. The tenfold terminal rate witness is separate from that allowance.

## Build and replay

Use CPython3.13 with a C compiler, Python headers, Git and GSL development files. Ubuntu native prerequisites include `build-essential` and `libgsl-dev`; the CI template uses setup-python3.13. On macOS an installed GSL prefix can be passed with `--gsl-prefix`. `requirements-observed.txt` pins the observed direct Python dependencies, including SciPy; it is not a universal wheel/hash lock. Pip installation needs network access.

Run from the package directory, choosing two fresh external work directories:

```sh
python3.13 bootstrap.py --work /tmp/msprime-base --variant baseline
python3.13 bootstrap.py --work /tmp/msprime-candidate --variant candidate
```

Each bootstrap writes `bootstrap-result.json` with absolute `python` and `source` paths. Use the corresponding values to run:

```sh
/absolute/venv/bin/python verify_variant.py --variant baseline --source /absolute/base/source --output-dir /tmp/base-evidence
/absolute/venv/bin/python verify_variant.py --variant candidate --source /absolute/candidate/source --output-dir /tmp/candidate-evidence
```

A deliberate `--reuse-python /absolute/venv/bin/python` avoids changing that environment; it is the mode used in the actual local packaging replay. `--prepare-only` downloads, validates and stages sources without compiling or testing. The native build's SCM version1.4.5.dev2 is metadata only; commit and file hashes establish source identity. A candidate build applies the frozen patch only inside its downloaded tree and verifies the resulting demography hash.

`neutral_replay.py` is the outcome-neutral pytest runner. It records actual package/native import paths and hashes, selected node IDs, collected counts and outcomes. Its default is the16 preexisting trajectory cases. `verify_variant.py` additionally executes the unchanged30-case standalone tests and68-case independent analytic oracle in fresh processes. No module substitutes or monkeypatches implement the proposed correction.

Install `ci/native-replay.yml` at repository `.github/workflows/msprime-terminal-hazard.yml` for the public destination `projects/msprime-terminal-hazard-auditor`. It installs native Ubuntu prerequisites, builds both variants, and runs all three evidence groups. This workflow has not been dispatched by this packaging worker.

## Scope and limits

The patch does not repair inherited time-zero duplicate processing, ancient-sampling ties, growth-survival grid approximation or the finite Taylor matrix exponential. For positive-time jumps it reuses existing interior event order. Unsupported empty-query input naturally changes from an internal exception to empty arrays; this incidental behavior is disclosed, not promoted to a new API guarantee. Floating-zero survival retains its established NaN-rate convention.

Local compilation emitted GSL deployment-target warnings. No fresh pip environment, Ubuntu/macOS portability matrix, full upstream/CUnit suite, stochastic Monte Carlo validation, maintainer approval or external peer review is claimed. Root owns the independent technical/source/value review and public-release decision. Source archives and native libraries are not vendored. See THIRD_PARTY.md and the included licenses.

Final test-only revision: a separately reviewing context identified that `(P*C)/P` need not equal `C` bit for bit. The size-rate regression now uses explicit relative tolerance1e-12 (zero absolute allowance). The runtime correction and independent68-case oracle are unchanged; earlier validation evidence is retained in the private dossier. This avoids treating rounding as a biological contract.
