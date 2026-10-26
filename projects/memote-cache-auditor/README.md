# MEMOTE sequential API cache audit

**A local candidate resets model-dependent support caches before each sequential `test_model` run.** It corrects demonstrated stale results after changes between calls while retaining cache reuse within a run. This is a software audit; no biological discovery or default single-run CLI error is claimed.

For reactions A→B and A→2B, positive conserved masses would require both mA=mB and mA=2mB, which is impossible. The original API nevertheless passes the modified model when an earlier call cached A→B alone. A fresh model copy correctly fails. Changing the second reaction's biomass SBO annotation also leaves stale classifications and collection parameters across calls.

| Sequential API scenario | Frozen original | Candidate |
| --- | --- | --- |
| Add inconsistent internal reaction between calls | Incorrect pass | Correct failure |
| Remove inconsistent internal reaction between calls | Incorrect failure | Correct pass |
| Add biomass annotation between calls | Stale failure and collection parameters | Correct pass and biomass parameter |
| Exception after caching, then modify and retry | Incorrect pass | Correct failure |
| Stable queries in two tests during one run | One classification | One classification |

The [production patch](patches/clear-caches-between-api-runs.patch) adds an explicit reset for all ten existing model-key caches and calls it before pytest collection. The API entry documentation states the lifetime policy. Upstream maintainer acceptance is pending; independent AI review is a computational review, not external human peer review.

## Reproduce

Use Python 3.11 and Git on Linux or macOS, install the pinned runtime stack, then run from the repository root:

```sh
python3.11 -m venv .venv-cache
. .venv-cache/bin/activate
python -m pip install -r projects/memote-cache-auditor/requirements.txt
python -B projects/memote-cache-auditor/verify.py
```

The verifier checks seven relevant installed source files against frozen MEMOTE commit `bf76477b4087feefc8ad429ebd492e30ef10eebc`, copies the package to disposable directories, applies only the production patch, and checks actual imports before running the real public API/core tests. Baseline and candidate run in separate processes with an explicit seven-key environment allowlist and plugin autoload disabled. Every pytest entry uses `--tb=short`; no environment or locals are dumped. No network is needed after dependency installation, and installed packages are never modified.

Default output is JSON on stdout. `--output artifacts/memote-cache-report.json` additionally writes a report. Archived evidence is checked by manifest and remains read-only. [Exact oracle](evidence/exact-oracle.json), [recorded red/green evidence](evidence/red-green.json), [source pin](source-manifest.json), and [license notices](THIRD_PARTY.md) are included.

## Costs and limits

A five-trial local diagnostic on iJO1366 (2,583 reactions, 1,805 metabolites) measured about **0.240 s** to call all ten helpers once after reset, including about **0.205 s** for the quadratic duplicate-metabolite classifier. Clearing all populated caches took about 0.121 ms. These are classification diagnostics, not end-to-end report timings or production bounds; individual cold times must not be summed because dependencies share caches. Each run intentionally loses cache reuse from earlier runs, including unchanged models. Within-run reuse remains. Full measurements are in [performance.json](evidence/performance.json).

The automatic change covers sequential, non-overlapping API entry boundaries. Direct support calls after mutation need explicit `clear_support_caches()`. Same-suite changes, temporary-context cache contamination, arbitrary annotation rollback, concurrent/nested API sessions and immediate resource release remain unaddressed. The global reset affects caches for every model instance; new decorated support functions require an inventory update. No full genome-scale test suite or alternate solver was validated; the local API checks use GLPK.

Caching is intentional prior work: [issue #325](https://github.com/opencobra/memote/issues/325), [PR #420](https://github.com/opencobra/memote/pull/420), and [PR #428](https://github.com/opencobra/memote/pull/428). Their performance motivation is preserved. Searches were bounded; no priority, uniqueness, or first-discovery claim is made. The upstream candidate is separate from the reaction-direction repair. Frozen primary source: [public API](https://github.com/opencobra/memote/blob/bf76477b4087feefc8ad429ebd492e30ef10eebc/src/memote/suite/api.py), [internal-reaction helpers](https://github.com/opencobra/memote/blob/bf76477b4087feefc8ad429ebd492e30ef10eebc/src/memote/support/consistency_helpers.py).
