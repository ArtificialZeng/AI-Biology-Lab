# Source credit and licenses

Our new auditor and regression code are covered by the lab's MIT license.

Proposed modifications to the upstream `ord-schema` source and tests are supplied under Apache-2.0, matching that project's license.

The files in `fixtures/upstream/ord_schema/` are preserved copies from [Open Reaction Database's ord-schema](https://github.com/open-reaction-database/ord-schema), commit `a9c1603a4fb040f251c1f174548d1747838fb247`. Their original copyright headers are retained; the Apache-2.0 license is included as `fixtures/upstream/LICENSE-APACHE-2.0.txt`. The patch documents proposed changes rather than modifying those archived originals.

The finding records and retirement metadata in `fixtures/data/` derive from [ord-data](https://github.com/open-reaction-database/ord-data), whose datasets and descriptive metadata are CC-BY-SA-4.0. That license is retained as `fixtures/data/LICENSE-CC-BY-SA-4.0.txt`; these derived data fixtures retain it. The findings fixture projects CAS values and record locations from the pinned official mirror. It does not contain complete reaction records. The sample manifest links the immutable original files rather than redistributing all 42 Parquet datasets.

Dataset-derived metadata and record locations in `evidence/sample_manifest.json` and `evidence/sample_audit.json` likewise retain CC-BY-SA-4.0 attribution and terms.

Dataset mirror: https://huggingface.co/datasets/open-reaction-database/ord-data/tree/b6c89761dacf809d4082acbdbebe72042a08db91 . Source identifiers, hashes, and sample-selection limits are recorded under `evidence/`. The original experimental authors remain the source of the measurements; no new experiments were performed.

The check-digit rule is established CAS documentation: https://www.cas.org/ja/training/documentation/chemical-substances/checkdig . Passing its arithmetic test does not verify a registered identity or license the substitution of an identifier.
