# Cogent3 missing-site distance audit

The public legacy calculator and `fast_slow_dist` application in Cogent3 2026.10.22 can reuse a zero-distance pair as if its sequences were identical. With a=A-, b=AA, c=AC, pairwise deletion gives distances0,0,1/2; the baseline app can instead report allzero, depending on row order.

The frozen patch limits reuse to complete encoded-row equality with a positive valid-site count. Other zero-distance pairs retain their own statistics. The separate modern `alignment.distance_matrix` path is unchanged. This is an engineering correction, not a new biological theorem or maintainer acceptance.

## Explicit setup and offline verification

Use Python3.13 in an isolated environment. Installation is the only step that needs package-network access:

```sh
python3.13 -m venv /tmp/cogent-audit-env
/tmp/cogent-audit-env/bin/python -m pip install -r requirements.txt
/tmp/cogent-audit-env/bin/python -B verify.py verify
```

Dependencies are fully pinned to the tested versions. Verification uses the actual installed Cogent3 package, first checks its baseline source identity, copies it into an external temporary directory, and overlays the frozen baseline or candidate module. It then runs the unchanged public calculator/app replay, exact-column exhaustive script, fresh independent checker, edge controls and53new regression tests. There is no AST substitute and no network request in verification. Scientific workers always run with normal assertions, `-B`, an external Numba cache and no pytest source cache, even if the outer verifier is invoked with `-O`.

All504ordered triples of distinct length-two rows over A,C,- give1512pair values per method. The baseline differs from direct-column counts in648values each for pdist/Hamming; the candidate in0. Fresh checking performs6204numeric comparisons and138edge checks. The84all-encoded-identical `include_duplicates=False` StopIteration cases are inherited and explicitly paired to the baseline; they are not hidden or called new regressions. The original full targeted source run was208passed/2intentional unsupported-model skips; this portable package reruns53passed/2skipped regressions plus both independent checkers.

Limits: finite checks do not establish universal package correctness. The inherited exact-alias pair length0, unsupported nonzero Hamming/pdist variance/invalid stderr, and all-identical exclusion exception remain. Full upstream suite, all Python versions and performance were not tested. The package pins Python3.13; its CI recipe covers macOS and Linux.

## Deterministic certificates

`certificate.json` stores categorical results/counts/source identities, not raw floating output, paths or timestamps. Every actual scientific worker must complete successfully and numerical comparisons must agree before certificate comparison. `archive_manifest.json` lists every public input; changed/missing/extra files and symlinks are rejected.

```sh
python -B verify.py regenerate --output /tmp/new-cogent-certificates
```

The output directory must be a new external path. Regeneration checks actual calculations, writes a byte-identical certificate and input manifest there, and never overwrites archived files. Reusing an output directory fails. It does not modify source, certificates or manifests in this package.

## Provenance and license

Pinned upstream develop commit9d1bd6b2441f13b49e8e73ef1bd8ced9448bf464, release2026.10.22. Scientific snapshotd79e93edca83cfa650e944e1de5545197537d5d510d93678185258e33ced7675. `SCIENTIFIC_SNAPSHOT.json` preserves that upstream research inventory. Vendored source/module/test bytes retain Cogent3's verified BSD3-clause license in `vendor/LICENSE-Cogent3.txt`; these files are an explicit exception to the lab's root MIT license. See THIRD_PARTY.md. Upstream has not reviewed this patch. AI assistance and a separate AI review were used; neither is human peer review.
