# Paper-Code Auditor: executable evidence and scoped reconstruction

**A reproducibility engineering case study, with offline CPU checks and an applicable patch.**

[中文简报](docs/BRIEF.zh-CN.md) · [Review packet](docs/UPSTREAM_REVIEW.md) · [Patch](patches/reproduction-fixes.patch) · [Recorded checks](evidence/verification.json)

We audited the documented execution route in `ziul-bio/DMS-SiteEffect-PLM`, pinned to commit `e168a2c825be0c3e28c300b8f6277dda1d999c4d` (2026-09-21). The repository accompanies *Intrinsic dataset features drive mutational effect prediction by protein language models*.

## Findings and corrections

| Finding in the pinned public code | Evidence | Proposed correction |
| --- | --- | --- |
| Cellular wrapper loops over an undefined array | Clean Bash control-flow check: zero scheduled calls; intended grid contains 198 | Use its defined cellular array |
| Both wrappers invoke a file absent from the current tree | Complete GitHub tree, source hashes, and public deletion history | Recover the original partial trainer **and its matching model module** from public history; invoke as a module |
| Available LoRA setup supplies unsupported keywords | Extracted actual setup raises `TypeError`; corrected call passes on a real-assay identifier fixture | Pass the existing `test_size` keyword |
| Recovered runtime selects logical GPU 1 while a wrapper exposes only one GPU | Configuration inconsistency in historical source | Logical GPU 0; preserve the user's visibility override and default to a single GPU |

The primary patch touches five files, including the two restored files. It preserves the historical partial-finetuning model implementation; it does **not** route the wrappers to LoRA or to a different current model. Hardware selection is a separately documented adaptation. The smaller [local-fixes patch](patches/local-fixes.patch) contains only the two independently verifiable local corrections; do not apply both patches to the same checkout.

The missing source is recoverable. Commit `8cd64d647a4841b56d4aa9f6d160c05d03b35b09` deleted it and its dependency; its parent `e0c2d6ae06903cb356abb7e1b579058115e5623a` contains both. [Restoration provenance](evidence/restoration_provenance.json) gives immutable URLs and hashes. Recovering public historical code does not establish which exact scripts produced the paper's published results.

## Run the audit

From the lab repository root, using Python 3.12 or 3.13:

```sh
python -m pip install -r projects/paper-code-auditor/requirements.txt
python projects/paper-code-auditor/audit.py
```

The check is offline and requires NumPy, pandas, Git, and Bash. It verifies hashes, reproduces the scoped failures, applies each proposed patch in a temporary copy, checks shell syntax and argument compatibility, and checks three seeded data splits. Output goes to ignored `artifacts/paper-code-auditor/verification.json`.

Only selected, reviewed data-handling AST nodes and reconstructed loop control execute. The original training shell body, model libraries, and training jobs do not run. The real-assay fixture retains only row and mutation identifiers; sequences and phenotype values are not needed for these checks. No GPU, model download, API key, or paid service is required.

## Apply the reconstruction to a local upstream checkout

First inspect the patch. In a local checkout at the pinned upstream commit:

```sh
git apply --check /absolute/path/to/AI-Biology-Lab/projects/paper-code-auditor/patches/reproduction-fixes.patch
git apply /absolute/path/to/AI-Biology-Lab/projects/paper-code-auditor/patches/reproduction-fixes.patch
```

This is a local proposed reconstruction, not an upstream-accepted fix. Compatible training dependencies, memory, and an actual GPU run still need verification. The original alternative custom checkpoint refers to author-specific files and is outside this reconstruction. An output directory is not evidence of successful scientific reproduction.

## Limits and negative control

- The original Lasso and LoRA splitting helpers agree on the three checked seeds. The initial allegation of unfair Lasso/LoRA folds was rejected and is not a finding.
- These checks establish software control flow, call compatibility, source provenance, and selected split invariants. They do not reproduce training, predictions, metrics, or the paper's biological conclusions.
- No author endorsement, upstream acceptance, or external human review is reported. Absence from the inspected issue tracker is not exhaustive priority clearance.

## Source credit and licensing

The pinned code and recovered historical code are third-party MIT-licensed material; the original copyright notice is retained in [fixtures/upstream/LICENSE.txt](fixtures/upstream/LICENSE.txt). [THIRD_PARTY.md](THIRD_PARTY.md) records provenance. Our audit code is covered by the lab's MIT license. This project is a software audit, not a claim of a new biological mechanism.
