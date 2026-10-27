# Validation state — 2026-10-27

The default portable replay passed using Python 3.13.5 on macOS ARM64. A tar.gz export was extracted into an unrelated temporary directory without Git metadata, the research environment or local research source. The driver downloaded and verified the official pinned archive, created its own fresh environment, installed the observed dependency versions, built the complete pinned source editable once, and passed `uv pip check`. It reproduced the two baseline positive-infinity failures using genuine public APIs, applied the byte-identical frozen pure-Python patch, checked all postpatch source hashes and started the candidate in a new subprocess.

The candidate passed 110 direct 100-digit Decimal definition cases, 20 extreme finite-order limiting-bound cases, analytic uniform cases and empty-community controls. The relevant upstream `test_base.py` passed 36 tests and 32 subtests. These are test counts, not a count of independent theorems. Baseline and candidate native-extension hashes match, consistent with one native build. A second independent native build was not performed or claimed.

The original research stage also recorded broader alpha/driver and lint checks. Those are inherited research evidence, not checks rerun by the portable packet. The public execution claim above is restricted to the portable replay and its relevant upstream test file.

## Compact evidence

- `evidence/portable-baseline.json` and `evidence/portable-candidate.json`: labeled sanitized derivatives of the actual local records. Private runtime prefixes are replaced by `<REPLAY_WORKDIR>`; source/native hashes, observed values, warnings and case counts are retained. Each derivative binds its original raw record by SHA256. Sanitization changes record bytes and is not presented as an unchanged raw log.
- `evidence/portable-summary.json`: compact default replay method, archive/source/patch bindings, one-build baseline/candidate process sequence and relevant test outcomes, with the original record bindings.
- `evidence/portable-export.json`: recorded export method and byte identity of the tested and final executable/reproduction inputs. README and licensing clarification occurred after that original export. Documentation edits do not imply another native execution.
- `evidence/optimized-failure-summary.json`: separate unsuccessful `-O` candidate attempt and its original evidence bindings.

`pin.json` fixes the upstream commit, archive hash, patch hash and the three changed source files' pre/post hashes. The public manifest lists only distributed public files and excludes itself and runtime outputs. The core reproduction inputs remain byte-identical to the tested export. Raw runtime logs, private local paths, internal research decisions and prompts are excluded from this distribution.

## Optimized-mode failure

A genuine fresh `python -I -B -O` candidate subprocess failed before any audit or diversity API execution. The unmodified upstream `skbio.sequence._substitution` import raises `TypeError: 'numpy.float64' object cannot be interpreted as an integer` in `_vector_to_matrix` (line 494). No unrelated fix or monkeypatch was applied. This blocks optimized-mode verification for the exact pinned upstream source. The normal default replay passed. The auditor's explicit failure gates do not establish that the entire upstream package supports `-O`; no successful optimized audit or full optimized driver replay is claimed.

## Remaining review and scope

The supplied template and active lab workflow use normal Python 3.13. Hosted execution is pending at publication; consult the current PR checks. Cross-platform and other-Python testing, a second independent native build, external maintainer review and human peer review are pending. Local technical/value review was AI-assisted computational review. A public release, upstream issue or pull request is not evidence of maintainer acceptance.

The correction targets finite orders >=2 under the existing supported normalization contract. Accuracy near one, normalization/sum overflow and global NumPy `under=raise` are outside its scope. Individual negligible scaled terms may underflow. The independent extreme-order bound is stated in the auditor and is not a direct Decimal oracle at those orders. The frequency of affected scientific workflows remains unknown. Established diversity formulas and numerical scaling do not constitute a new ecological theorem or biological discovery.
