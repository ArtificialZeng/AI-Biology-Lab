# Sources and provenance

- Upstream project: [tskit-dev/tskit](https://github.com/tskit-dev/tskit), MIT license retained in fetched source.
- Fixed baseline: [f154d1190e0142c0499d2698d76eed0fe6429876](https://github.com/tskit-dev/tskit/commit/f154d1190e0142c0499d2698d76eed0fe6429876).
- Public source archive: `https://codeload.github.com/tskit-dev/tskit/tar.gz/f154d1190e0142c0499d2698d76eed0fe6429876`.
- Baseline `c/tskit/trees.c` SHA-256: `2c5e9f77ad5901d6c05801abf7d46d2e8782187fe68551b7556138d27f21a616`.
- Relevant API: [`TreeSequence.trait_correlation`](https://github.com/tskit-dev/tskit/blob/f154d1190e0142c0499d2698d76eed0fe6429876/python/tskit/trees.py) and [native C implementation](https://github.com/tskit-dev/tskit/blob/f154d1190e0142c0499d2698d76eed0fe6429876/c/tskit/trees.c).
- `native.patch` is the frozen diff from that baseline to candidate commit `369ff4fce5eda3ab4dd25366c8812af9df3a8149`, affecting `c/tskit/trees.c`, `python/tests/test_tree_stats.py`, `c/CHANGELOG.rst`, and `python/CHANGELOG.rst`. It includes upstream context lines for application; full third-party files are not bundled.
- `audit.py` is original independent reference code, formatted from the R40 fresh-review oracle. Expected values come from `fractions.Fraction`, fixed descendant sets, explicit edge lengths, and allele partitions; tskit traversal/statistic/genotype routines and NumPy variance/correlation routines do not derive them. NumPy handles only API input/output arrays and tolerance comparison.

The MIT license in this packet covers original audit/reproduction material and original changes. It does not replace the upstream license, whose notices remain in downloaded source. No third-party paper or private correspondence is included.
