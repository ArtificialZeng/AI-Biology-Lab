# Frozen upstream contract

Base: msprime commit `45cfa59d828757352f7ab1a55c19b55102d40d22`, with exact tskit gitlink `f1b139edb376fbffb2978a5da2c9e4b948930b0c`. Source archives and file hashes are pinned in manifest.json.

The [current-source API equation](https://github.com/tskit-dev/msprime/blob/45cfa59d828757352f7ab1a55c19b55102d40d22/msprime/demography.py#L4239) defines the conditional mean pairwise coalescence rate at each requested time. `_pop_size_and_migration_at_t` selects the epoch to the right of a demographic boundary. The preexisting `TestDemographyTrajectories.test_one_pop` checks the post-change rate at the boundary; `test_two_pop` checks positive-time lineage movement; `test_nan_coalrate` requires NaN when computed survival is zero.

The base implementation separately handles the last augmented time using the previous rate matrix and without its ordinary point-movement operation. The patch uses the normal evaluation path for the last time and propagates only when a subsequent time exists. Its duplicate-zero guard preserves the previously executed zero-time point operations.

This is a bounded scientific-software correction. It is not a new coalescent theorem, a simulator-wide accuracy guarantee, or a claim that a published scientific result is affected. No maintainer or external human confirmation is claimed here.
