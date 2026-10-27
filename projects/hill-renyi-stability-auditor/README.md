# Hill / Renyi large-finite-order stability auditor

This packet documents a numerical correction to scikit-bio's public Hill and Renyi diversity APIs at commit `2f573b9a9411c3c41465d9abcd41c3ad3a7dc077` (`0.7.5-dev`). For the supported integer abundance input `[1, 1]` at order 2000, the unmodified APIs return positive infinity; the defining values are Hill diversity 2 and natural-log Renyi entropy log(2). A uniform 1024-taxon community at order 128 also reproduces the failure. The supplied patch uses a maximum-probability-scaled power sum for finite orders >=2. This is a scoped software correction based on established mathematics. The frequency and practical impact of affected workflows remain unestablished; no new ecological theorem or biological discovery is claimed.

This is the own-lab reproduction packet. The upstream report is [issue #2651](https://github.com/scikit-bio/scikit-bio/issues/2651), followed by [draft PR #2652](https://github.com/scikit-bio/scikit-bio/pull/2652) while the issue-first discussion is pending. Local technical and value checks were AI-assisted computational reviews; human and upstream maintainer review remain pending.

## Reproduce

Use Python 3.13 (tested: 3.13.5), a C compiler/toolchain, `patch`, network access to the official GitHub codeload endpoint and the Python package index, and enough disk space for a complete source build. The driver accepts Python >=3.10, but other versions and platforms have not been verified with these exact dependencies.

```sh
python3.13 reproduce.py
```

An optional work directory must not already exist:

```sh
python3.13 reproduce.py --work-dir /tmp/hill-renyi-new-replay
```

The default driver downloads and SHA256-verifies the pinned official archive, checks the prepatch source hashes, creates a fresh work-directory-owned virtual environment, installs the exact observed dependency versions and builds the full pinned source editable with runtime dependencies enabled. It invokes genuine public APIs in a fresh baseline subprocess, verifies and applies the frozen patch, checks the postpatch hashes, and invokes the candidate in another fresh subprocess. There is one native build because the patch changes Python only. Imported Python/native paths and hashes are recorded. No stubs, extracted-function substitutes or monkeypatches are used. `-I -B`, source checks and bytecode removal guard import provenance. Reproduction paths are derived from the packet and fresh work directories; existing environments are not installed into. `uv` is used when available, otherwise the fresh environment's pip is used.

`audit_api.py` checks analytic uniform values, 110 direct 100-digit Decimal definition cases and 20 extreme finite-order cases against an explicitly stated limiting bound, plus empty-community controls. The limiting-bound cases are not direct finite-order Decimal oracle calculations. The patched upstream relevant test file is also run. Failure checks are explicit conditional exceptions and subprocess exit checks, rather than removable `assert` statements.

## Evidence and limits

The normal default driver passed after a portable export into an unrelated temporary directory. See [VALIDATION.md](VALIDATION.md) and the compact evidence files. `pin.json` binds the archive, frozen patch and pre/post source files. `MANIFEST.json` binds the public packet files; hashes establish content identity, not mathematical correctness or external review.

Public baseline/candidate records are labeled sanitized derivatives of the original local execution records, with private runtime paths replaced by `<REPLAY_WORKDIR>` and original-record SHA256 bindings retained. `evidence/portable-export.json` records the tested export and the byte identity of executable/reproduction inputs. `evidence/portable-summary.json` preserves the execution method and test outcomes. `evidence/optimized-failure-summary.json` records the separate failed optimized candidate attempt. Internal research summaries and raw local logs are not distributed.

Normal Python execution is the verified route. A fresh `-O` candidate attempt failed during an unrelated, unmodified upstream import before the audit or diversity API calls. No workaround was applied and no optimized-mode support is claimed. `ci/reproduce.yml` is an opt-in, normal-Python-3.13 template; The lab also supplies an active workflow at `.github/workflows/hill-renyi-stability.yml` for Ubuntu/macOS using Python 3.13. Hosted results are pending at publication; consult the current PR checks.

The patch does not address accuracy near order one, count-normalization/sum overflow, or every platform. Negligible scaled terms may still underflow; global NumPy `under=raise` is not guaranteed. Exact dependency versions describe one observed environment, not a supply-chain lock with wheel hashes. Future package-index availability and build compatibility remain external constraints. No population-level impact estimate, exhaustive accuracy guarantee or absolute novelty claim is made.

## Sources, licensing and contribution

Primary sources are the [pinned implementation](https://github.com/scikit-bio/scikit-bio/blob/2f573b9a9411c3c41465d9abcd41c3ad3a7dc077/skbio/diversity/alpha/_base.py), [Hill API](https://scikit.bio/docs/latest/generated/skbio.diversity.alpha.hill.html), [Renyi API](https://scikit.bio/docs/latest/generated/skbio.diversity.alpha.renyi.html) and [original API PR #2074](https://github.com/scikit-bio/scikit-bio/pull/2074). Definitions and upstream implementation credit belong to scikit-bio and its contributors. Original packet scripts are MIT-licensed (`LICENSE`, copyright Zijian Zeng). The scikit-bio-derived patch retains its BSD-3-Clause notice in `THIRD_PARTY_LICENSE.txt`; downloaded upstream source retains its own licenses.

The [contributor guide](https://scikit.bio/contribute.html), inspected 2026-10-27, asks for issue discussion before code contribution, tests, a changelog and maintainer review. The [pinned PR template](https://github.com/scikit-bio/scikit-bio/blob/2f573b9a9411c3c41465d9abcd41c3ad3a7dc077/.github/PULL_REQUEST_TEMPLATE.md) requires disclosure and attribution of external content and compatible licensing. Any later contribution should recheck current guidance and source freshness.

This patch and reproduction packet were developed with substantial AI assistance (OpenAI Codex), including code drafting, test generation and local auditing. These checks are computational evidence, not external human peer review. The human submitter remains responsible for reviewing and validating any contribution.

The upstream submission has byte-identical implementation and tests. Its changelog additionally links issue #2651; `evidence/upstream-submission.json` records this documentation-only difference from the frozen reproduction patch.
