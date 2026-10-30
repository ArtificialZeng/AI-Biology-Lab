# tskit trait-correlation: large-offset native variance audit

A bounded numerical-software defect witness and minimal C correction. On the pinned upstream baseline, valid nonconstant traits with a large common offset can produce incorrect or nonfinite trait correlations because variance is computed by subtracting large raw second moments. Anchoring each column and then summing centered squared residuals removes this failure on the included exact-reference cases.

This packet contains original audit/reproduction code and a four-file upstream patch (C, twelve regression cases, and two changelog entries). It does not redistribute the upstream source archive. See [sources.md](sources.md), [note/statement.md](note/statement.md), and [status.json](status.json).

## Reproduce

Use Python 3.11 or newer, a working C compiler, and the `patch` command. Create a new environment; do not alter a preexisting research environment:

```sh
python3 -m venv /tmp/tskit-trait-audit-env
/tmp/tskit-trait-audit-env/bin/python -m pip install -r requirements.txt
/tmp/tskit-trait-audit-env/bin/python reproduce.py --workdir /tmp/tskit-trait-native-replay
```

Run from this packet directory, or use absolute script and requirements paths. Choose a new or empty work directory. The default mode downloads the public archive of commit `f154d1190e0142c0499d2698d76eed0fe6429876`, checks the baseline C SHA-256 in both independently extracted trees, applies `native.patch` to the candidate only, and invokes the current Python interpreter's `setup.py build --build-base` separately for each tree. Safe extraction preserves `python/lib -> ../c`. Build logs, patch log, two full audit receipts and `summary.json` go to the selected work directory. No msprime or scipy dependency is needed.

Success requires exactly 396 cases in each isolated native replay, zero candidate failures, and at least one baseline failure. The observed reference machine had 84 baseline passes / 312 failures and 396 candidate passes / zero failures; baseline failure counts can differ with platform. All runner gates use explicit exceptions, and child `-E -B` interpreters receive the requested native path through `sys.path`, so inherited optimization settings cannot disable checks or redirect native imports.

For a quick byte-integrity and Python-syntax check with standard-library dependencies only:

```sh
python3 reproduce.py --verify-packet-only
```

This option does not execute scientific checks and is not a correctness certificate.

For a fast replay against two existing, separately built native libraries:

```sh
python reproduce.py --workdir /tmp/tskit-trait-existing-replay \
  --existing-baseline-build /path/to/build-baseline \
  --existing-candidate-build /path/to/build-candidate
```

Pass the build-base directories containing `lib.PLATFORM/tskit` and `_tskit` extensions. This mode checks import provenance and replays the audit, but cannot certify where those preexisting binaries came from; the receipt identifies it explicitly. Fresh source builds are the default reproducibility route.

The test is a standard floating-point stability correction. It makes no claim of mathematical novelty, human external verification, universal correctness for all finite doubles, or a new biological result. Extreme dynamic ranges can still overflow or underflow centered residuals or squared residuals. Upstream's larger optional-dependency regression suite is separate from this small standalone packet.

The frozen v2 packet was also replayed from a fresh public source download under parent `python -O`; both independent native builds completed and the expected baseline/candidate split reproduced. Portable receipts are in [validation.json](validation.json). Generated third-party source archives, source trees and native build directories are not bundled.

Upstream correction: [tskit PR3494](https://github.com/tskit-dev/tskit/pull/3494), commit `369ff4fce5eda3ab4dd25366c8812af9df3a8149`. The default work directory is `/tmp/tskit-trait-native-replay`; choose `--workdir` to retain separate runs.

## Upstream validation and review status

On 27 October 2026, all 17 upstream PR checks passed at candidate commit `369ff4fce5eda3ab4dd25366c8812af9df3a8149`, including the 13-job [Tests workflow](https://github.com/tskit-dev/tskit/actions/runs/40505067326), lint and documentation builds. This hosted validation supplements the bounded packet audit; it does not establish correctness for every floating-point input.

In a [public comment](https://github.com/tskit-dev/tskit/pull/3494#issuecomment-6442764401), maintainer Jerome Kelleher gave an initial positive assessment and invited original code author Peter Ralph to review. This is a preliminary comment, with formal author review and merge still pending. The scientific witness, patch, and reproduction code are unchanged.
