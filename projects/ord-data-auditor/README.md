# ORD Data Auditor: checksum validation without guessing chemical identity

**A standard-library identifier guard, a reproducible frozen-sample audit, and a proposed upstream warning-level patch.**

[中文简报](docs/BRIEF.zh-CN.md) · [Review packet](docs/UPSTREAM_REVIEW.md) · [Patch](patches/ord-schema-cas-checksum.patch) · [Sample results](evidence/sample_audit.json) · [Offline regression](evidence/verification.json)

## Result

We inspected a frozen convenience sample from the official Open Reaction Database mirror:

| Quantity | Count |
| --- | ---: |
| Active dataset files | 42 |
| Reaction records | 20,720 |
| CAS identifier occurrences | 11,949 |
| Format anomalies already covered by ORD's existing warning | 9 |
| Well-formatted value with a mismatching check digit | 1 |

The sample consists of all files below 250,000 bytes in the retrieved mirror tree at commit `b6c89761dacf809d4082acbdbebe72042a08db91`. It is **not** a representative sample; these counts must not be extrapolated to the full database. [The manifest](evidence/sample_manifest.json) pins file URLs, byte sizes, and SHA-256 hashes.

The checksum witness is reaction `ord-04b7893fc4e441518687949508e32cdc`, field `inputs[Ligand].components[0].identifiers[2]`. The recorded value `977739-46-3` has a weighted prefix sum of 251, inconsistent with its stored check digit. This proves an arithmetic inconsistency, **not** which compound or replacement identifier is correct. The tool never changes the data or proposes an identity replacement.

The existing [ORD CAS branch](https://github.com/open-reaction-database/ord-schema/blob/a9c1603a4fb040f251c1f174548d1747838fb247/ord_schema/validations.py#L994-L1000) checks format, but not the check digit. This project integrates a standard, established check; it is not a new chemical algorithm. The same dataset already has an [OCR/name curation issue](https://github.com/open-reaction-database/ord-data/issues/300). No priority or maintainer-acceptance claim is made.

## Offline checks: no dependencies beyond Python, Git

From the lab repository root:

```sh
python projects/ord-data-auditor/verify.py
```

Checks include known positive and negative controls, 1,200 artificial arithmetic cases, the ten real finding records, source integrity, patch applicability, and execution of the **actual extracted CAS-specific branch** before and after patching. Synthetic checksum cases do not represent registered compounds.

The proposed patch retains the upstream format expression and warning severity, adding only check-digit diagnostics plus targeted upstream tests. The standalone guard uses an ASCII-strict format; the narrower upstream patch intentionally leaves legacy Unicode-format behavior unchanged. Neither performs a registry lookup.

Example diagnostic (findings produce exit status 1):

```sh
python projects/ord-data-auditor/cas_guard.py --value 7732-18-5 --value 977739-46-3
```

For identifier JSONL files with `value` and optional `identifier_type`, `dataset_file`, `reaction_id`, and `path` fields:

```sh
python projects/ord-data-auditor/cas_guard.py --input identifiers.jsonl --output report.json
```

## Reproduce all 42 frozen sample files

Use a separate environment for the optional data reader:

```sh
python3 -m venv .venv-ord
source .venv-ord/bin/activate
python -m pip install -r projects/ord-data-auditor/requirements-sample.txt
python -m pip install --no-deps ord-schema==0.9.6
python projects/ord-data-auditor/audit_sample.py --download
```

The download is about 2.4 MB. The reader uses only the schema's generated Protobuf classes; the `--no-deps` installation avoids installing unrelated chemistry tooling. The auditor verifies every file hash, excludes retired dataset IDs using the frozen reference, checks row/Protobuf IDs, visits CAS fields, and compares aggregate counts with the recorded sample. Outputs and downloaded data stay under ignored `artifacts/`. An existing directory can be supplied with `--data-dir`.

Fastparquet is used because the earlier PyArrow 19 reader failure was an already known upstream compatibility issue, not an ORD data defect. Schema 0.9.6 parses the checked fields; unknown newer fields are not semantically validated.

## Proposed upstream patch

The patch targets `ord-schema` commit `a9c1603a4fb040f251c1f174548d1747838fb247` and adds warning-level checks and targeted test cases. Inspect it before applying in a local upstream checkout:

```sh
git apply --check /absolute/path/to/AI-Biology-Lab/projects/ord-data-auditor/patches/ord-schema-cas-checksum.patch
git apply /absolute/path/to/AI-Biology-Lab/projects/ord-data-auditor/patches/ord-schema-cas-checksum.patch
```

In a compatible upstream development environment, the proposed tests can be run with `pytest ord_schema/validations_test.py -k 'cas_number or compound_identifier_format'`. The full upstream suite and complete RDKit/ORD semantic validation have **not** been run here; our offline tests isolate the CAS branch. No reaction balancing, biological inference, missing-product accusation, or model-performance claim follows from this audit.

## Licensing and status

Our code is MIT-licensed. The pinned schema code/tests are Apache-2.0, and the derived data records/metadata are CC-BY-SA-4.0. See [THIRD_PARTY.md](THIRD_PARTY.md) and the retained notices. No upstream issue/PR or external message has been sent, no author endorsement is reported, and no chemical identities were automatically repaired.
