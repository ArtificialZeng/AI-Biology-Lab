# Validation history and scope

The mechanism was already reported by Fabian Fröhlich in pyPESTO issue1794.
The contribution here is a bounded proposed correction and reproducible tests.

Early local setup attempts failed because the private build-tool directory was
absent from PATH, a noise placeholder was omitted, AMICI reserves parameter `k`,
different initial states shared a compiled-model cache, an empty PEtab condition
was invalid, and a pre-created model-output directory suppressed compilation.
These were audit-fixture/setup errors, corrected before the decisive runs;
they are not defects attributed to the upstream project.

An early nominal baseline pytest command collected a test inside candidate,
causing pytest to prepend candidate's package path. Its passing baseline result
was rejected. The authoritative paired replay pre-imports and verifies the
selected source before starting pytest, using a neutral test copy.

An auxiliary initial assertion expected an absent inner-parameter return to be
present as None; return filtering removes it. That assertion was removed from
the final contract. Consequently four baseline Hessian-containing requests that
already return inf/NaN are reported as preserved behavior, not newly discovered
false-finite objective failures. Only three direct public requests and three
NaN utility cases fail the final baseline regression.

The original native run produced6baseline failures/33passes and39candidate
passes, with no skips; nine selected existing solver tests also passed.
Existing solver cases emitted zero-sigma warnings. The original neutral pytest
configuration produced cache-write warnings under /dev; the portable replay
uses an external work-directory config and disables pytest caching instead.

Archived supplemental controls recorded19passing candidate assertions,
including actual integration failure, empty sensitivity requests and supported
adjoint requests; the baseline had two additional false-finite failures.
These supplement the native matrix but are not silently counted as39-case CI.
No complete upstream test suite or external biological experiment was run.

Fresh independent review identified a low-level edge outside the measured-NaN
repair: manually supplying sigma=0 at a perfect fit changes `compute_nllh` from
−inf to NaN. The actual profiled perfect-fit path raises `AssertionError` in
both versions while estimating sigma. This boundary is retained explicitly;
no claim of identical zero-sigma edge behavior or a solved perfect-fit policy
is made.
