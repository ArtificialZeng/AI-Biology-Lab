# msprime terminal hazard auditor

At a demographic jump placed at the final requested time, the pinned pre-fix msprime baseline can return the preceding instantaneous coalescence rate. With one population changing from size1000 to10000 at time100, a query ending at100 returns0.0005 instead of0.00005; appending101 changes the shared answer at100. The included two-file patch evaluates the terminal time through the existing interior path, including the existing positive-time lineage movements.

This package contains the exact frozen patch, unchanged neutral regression and analytic-oracle inputs, source hashes, native bootstrap and paired replay. Source contracts are in source.md. It prepares msprime commit `45cfa59d828757352f7ab1a55c19b55102d40d22` and the exact tskit gitlink `f1b139edb376fbffb2978a5da2c9e4b948930b0c`; downloaded archives, build trees and native binaries stay outside the public project.

## Upstream acceptance (27 October 2026)

[msprime PR #2547](https://github.com/tskit-dev/msprime/pull/2547) was approved by Jerome Kelleher at 13:59:54 UTC and merged by him at 14:09:01 UTC as `7954a7285c1738e1395a11f9e41446ec8f3ce610`. Aaron Ragsdale also left positive maintainer feedback. All 13 upstream checks succeeded on reviewed head `56427173f6c013dc11903a9e0d641f83e514ea51`. At the recorded 14:48:49 UTC check, upstream main contained the merge and its demography SHA-256 `415c7b43d43043e0966af33817dab59ce073a1a95b18d96019324ec817d9bb37` matched the frozen candidate. See `evidence/upstream-acceptance.json` for public review/check links and the verification scope.

The lab's earlier PR #33 head `94b4e1f6a1d1c69986d964fc65478dc0d7732ff2` passed all 9 recorded hosted checks on 27 October 2026, including both native paired replay jobs, by 10:26:09 UTC ([native hosted run](https://github.com/ArtificialZeng/AI-Biology-Lab/actions/runs/40488186021)). Those checks apply to that earlier head. The acceptance metadata was prepared before this subsequent revision's checks ran; current checks are shown on [lab PR #33](https://github.com/ArtificialZeng/AI-Biology-Lab/pull/33). The original packaging-worker flags remain historical records of that worker's execution.

This is upstream maintainer acceptance of a bounded numerical diagnostic correction. It does not establish experimental biological validation or journal peer review.

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

Install `ci/native-replay.yml` at repository `.github/workflows/msprime-terminal-hazard.yml` for the public destination `projects/msprime-terminal-hazard-auditor`. It installs native Ubuntu prerequisites, builds both variants, and runs all three evidence groups. The original packaging worker did not dispatch this workflow; the later hosted run recorded above succeeded for the dated lab head.

## Scope and limits

The patch does not repair inherited time-zero duplicate processing, ancient-sampling ties, growth-survival grid approximation or the finite Taylor matrix exponential. For positive-time jumps it reuses existing interior event order. Unsupported empty-query input naturally changes from an internal exception to empty arrays; this incidental behavior is disclosed, not promoted to a new API guarantee. Floating-zero survival retains its established NaN-rate convention.

Local compilation emitted GSL deployment-target warnings. The original local packaging replay did not use a fresh pip environment or execute a full upstream/CUnit suite, Ubuntu/macOS portability matrix, or stochastic Monte Carlo validation. Later upstream CI and maintainer acceptance are recorded separately above. No experimental biological validation or journal peer review is claimed. Source archives and native libraries are not vendored. See THIRD_PARTY.md and the included licenses.

Final test-only revision: a separately reviewing context identified that `(P*C)/P` need not equal `C` bit for bit. The size-rate regression now uses explicit relative tolerance1e-12 (zero absolute allowance). The runtime correction and independent68-case oracle are unchanged; earlier validation evidence is retained in the private dossier. This avoids treating rounding as a biological contract.

Upstream formatting synchronization (27 October 2026): the frozen patch and standalone regression now include Ruff's requested assertion line wrap from msprime PR #2547. Both test ASTs, the correction, and the independent analytic oracle are unchanged. The updated patch was applied to the original pinned baseline and matched both upstream files exactly; see `evidence/format-synchronization.json`. This is maintenance of the same result.

Archive integrity note: the earlier packaging file map included an ignored, untracked `frozen/__pycache__/test_terminal_contract.cpython-313-pytest-9.1.1.pyc` entry whose local hash no longer matches. `archive_manifest.json` preserves that historical map and records the discrepancy; its current `files` map covers tracked project files plus the new acceptance evidence, excluding the archive manifest itself. No scientific input or cache file was changed by this metadata update.
