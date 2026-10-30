# Cogent3 missing-site distance audit

The public legacy calculator and `fast_slow_dist` application in Cogent3 2026.10.22 can reuse a zero-distance pair as if its sequences were identical. With a = A-, b = AA and c = AC, pairwise deletion gives distances 0, 0 and 1/2. The baseline app can instead report all zero, depending on row order.

The frozen patch limits reuse to complete encoded-row equality with a positive valid-site count. Other zero-distance pairs retain their own statistics. A fixed-size off-diagonal count guard now avoids a full-row equality scan when an observed mismatch already rules out equality. Equal encoded rows can only contribute diagonal counts, so the guard removes no valid alias. The separate modern `alignment.distance_matrix` path is unchanged. This is an engineering correction, not a new biological theorem or maintainer acceptance.

## Explicit setup and offline verification

Use Python 3.13 in an isolated environment. Installation is the only step that needs package-network access:

```sh
python3.13 -m venv /tmp/cogent-audit-env
/tmp/cogent-audit-env/bin/python -m pip install -r requirements.txt
/tmp/cogent-audit-env/bin/python -B verify.py verify
```

All 20 dependencies are pinned to tested versions. Verification uses the actual installed Cogent3 package, checks its baseline source identity, copies it into external temporary directories, and overlays the frozen baseline or candidate module. It then runs the unchanged public calculator/app replay, exact-column exhaustive script, fresh independent checker, edge controls and 53 regression tests. There is no AST substitute or network request during verification. Scientific workers always run with normal assertions, `-B`, an external Numba cache and no pytest source cache, including when the outer verifier uses `-O`.

All 504 ordered triples of distinct length-two rows over A, C and - give 1,512 pair values per method. The baseline differs from direct-column counts in 648 values each for p-distance and Hamming; the candidate differs in 0. Fresh checking performs 6,204 numeric comparisons and 138 edge checks. The 84 all-encoded-identical `include_duplicates=False` StopIteration cases are inherited and explicitly paired to the baseline. The original targeted source run had 208 passing tests and 2 intentional unsupported-model skips; this portable package reruns 53 passing tests and 2 skips, plus the independent checkers.

Limits: finite checks do not establish universal package correctness. The inherited exact-alias pair length of 0, unsupported nonzero Hamming/p-distance variance and invalid stderr, and all-identical exclusion exception remain. The full upstream suite, all Python versions and performance across workloads were not tested. The package pins Python 3.13; its CI recipe covers macOS and Linux.

## Focused performance refinement

The second source pass responds to an actual inline review of the same upstream PR. NumPy row equality scans data whose size grows with alignment length. The new guard examines a matrix whose size depends only on the fixed alphabet, before the existing row-equality and valid-observation conditions. It preserves all recorded correctness outputs from v1, including the inherited exceptions.

A finite local benchmark used actual source calculators, preconverted DNA rows, one terminal mismatch, 3 warmups and 17 alternating repetitions. At 5,000,000 positions, the median calculator time was 3.79 ms for v1 and 2.07 ms for v2. This is a workload-specific measurement, not a universal speedup or an end-to-end alignment-loading claim. There is no timing-threshold test. Detailed source-pass evidence is identified by the scientific snapshot below.

## Deterministic certificates

`certificate.json` stores categorical results, counts and source identities, without raw floating output, paths or timestamps. Every scientific worker must complete and actual numerical comparisons must agree before certificate comparison. `archive_manifest.json` lists every public input; changed, missing and extra files and symlinks are rejected.

```sh
python -B verify.py regenerate --output /tmp/new-cogent-certificates
```

The output directory must be new and external. Regeneration checks actual calculations and writes a byte-identical certificate and input manifest there. It never overwrites archived inputs; reusing an output directory fails.

## Provenance and license

Pinned upstream baseline: commit 9d1bd6b2441f13b49e8e73ef1bd8ced9448bf464, release 2026.10.22. The focused v2 scientific snapshot is 1005f00f170a60f7f088d2f05cd1de3cd7a669c62677e3eaf657a6c5e99d0138; candidate source SHA-256 is cb6bd15dac34e9ca6074da73ff69c627e1c6ca50cf78302897110ff26dcceb52. `SCIENTIFIC_SNAPSHOT.json` preserves its 30-file inventory. This updates the same PR 2789 correction rather than creating a separate result.

Cogent3-derived source, tests and patch retain the verified BSD 3-clause license in `vendor/LICENSE-Cogent3.txt`. They are an explicit exception to the lab's root MIT license. See THIRD_PARTY.md. Upstream acceptance is pending. AI assistance and a separate AI review were used; neither is human peer review.
