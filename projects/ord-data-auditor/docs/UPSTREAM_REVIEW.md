# Review packet: CAS check-digit warning

Submitted second with the lab owner's authorization on 2026-10-26: [upstream PR #1165](https://github.com/open-reaction-database/ord-schema/pull/1165). Head commit: `4e524972a7d4b321d65369f76b7f049655b5079d`. Open at submission; acceptance remains pending.

Target: `open-reaction-database/ord-schema` commit `a9c1603a4fb040f251c1f174548d1747838fb247`.

## Problem and proposal

The current CAS-specific branch warns on malformed syntax but accepts a well-formatted value even when its check digit disagrees with the prefix. A frozen 42-file sample contains one such case in reaction `ord-04b7893fc4e441518687949508e32cdc`. The arithmetic witness and exact record path are in the project README and `evidence/sample_audit.json`.

The candidate adds a `context.warn` in the existing format-valid branch. It retains existing format behavior and severity, does not mutate messages, does not upgrade this finding to an error, and does not guess the identity-correct replacement. Targeted upstream tests cover known valid examples and wrong check digits; the test helper already checks that validation leaves the message unchanged.

## Validation supplied

- Offline source hashes and patch-application check.
- Actual extracted CAS-branch execution: 9 warnings before and 10 after on the ten finding fixtures.
- Positive/negative/boundary controls and 1,200 artificial arithmetic cases, checked against a separate arithmetic formulation.
- Full frozen convenience-sample rerun: 42 files, 20,720 reactions, 11,949 CAS occurrences; each input hash verified.

The offline checker isolates the CAS-specific body. Before submission, a separate development environment was synced from the upstream lockfile: all 322 tests in `ord_schema/validations_test.py` passed, including the five added cases and the helper's unchanged-message checks. Ruff check/format and `ty check ord_schema` passed. The full suite was attempted with live resolvers disabled and stopped on a missing Git LFS dependency in a history fixture (72 passed, 1 setup error). Full-suite success and complete semantic validation of the frozen corpus are not claimed. See `evidence/upstream_validation.json`. The candidate patch was formatted by the upstream pinned Ruff version before submission; its arithmetic is unchanged.

## Known overlap and unresolved identity

Issue https://github.com/open-reaction-database/ord-data/issues/300 already discusses OCR/name defects in the same dataset. Its inspected text does not mention this checksum, but that does not establish priority. Existing dataset-review tools also overlap the broad auditing task. The proposed contribution is a small standard validation improvement, not a new chemistry discovery.

The correct identifier for the compound remains unverified. In particular, changing only the final digit would produce an arithmetic-valid string without proving it denotes the intended compound.
