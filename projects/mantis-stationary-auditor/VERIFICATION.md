# Verification record and reproducibility limits

The package was replayed using `/opt/anaconda3/bin/python`, CPython 3.13.5, NumPy 2.2.3, SciPy 1.15.3, SymPy 1.13.3 and networkx 3.4.2. Python 3.12 remains the recommended upstream CI target; the archive does not install, download, or assume an editable Mantis distribution.

Checks performed on this package:

- Deterministic certificate verification succeeds normally and with `-O`; JSON outputs are identical.
- The archive was copied to a different directory containing spaces, then made read-only; normal and optimized replay both succeeded without project writes.
- A copy with a modified candidate implementation was rejected by the manifest before executing an audit.
- Related actual-library upstream/regression tests: 37 passed.
- Four selected new regressions reproduce failures against the pinned baseline package.
- The full candidate suite retains the disclosed local baseline CHA precision failure: 166 passed, 1 skipped, 1 failed. No unrelated precision tolerance or implementation was changed.

Certificates represent exact finite CTMC reference laws and successful floating comparisons at explicit tolerances, not exactness of floating equilibrium computations. The 352-case owner family uses a separate standard-library Fraction oracle; five owner cases exercise the public API end to end. The fresh AI reviewer implemented an independent oracle and covered 990 two-species plus 87 multispecies/catalyst cases, with additional zero-rate, positive-underflow, completeness and constructor controls. These are computational checks and an AI review, not external human peer review, maintainers' acceptance, or a formal proof-assistant result.

The raw private process logs and measured floating errors are not part of this public-ready archive. The full vendored baseline and candidate tests allow the disclosed outcomes to be checked directly. The immutable input/certificate manifest is a corruption check, not an external attestation.
