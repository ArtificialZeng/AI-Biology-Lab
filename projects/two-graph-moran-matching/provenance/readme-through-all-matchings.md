# A fixation bound for every matching complement

For integers N≥3 and 1≤m≤⌊N/2⌋, residents reproduce on K_N and neutral mutants on K_N minus a matching of m edges under type-specific Moran birth–death. The proposed theorem bounds fixation from **every singleton vertex** v:

    h_v ≤ 1/N − 2m/[N²(N−1)(N(N−2)+2m−2)] < 1/N.

Uniform averaging over the starting vertex inherits the bound. The proof covers unbounded m, including perfect matchings. It combines an exact superharmonic identity with absorption; finite labelled-state audits are supplementary. The arbitrary connected noncomplete mutant-graph conjecture remains unresolved. The bound is not claimed sharp.

[All-matching proof](all-matchings/proof.md), [precise scope](all-matchings/problem.md), [source/model credit](all-matchings/source.md), and [bounded literature search](all-matchings/literature.md) give the scientific argument. The new range is m≥3. The exact m=1 formula is due to Tkadlec, Kaveh, Chatterjee and Nowak (2023), [Evolutionary dynamics of mutants that modify population structure](https://doi.org/10.1098/rsif.2023.0355); their N=6 numerical enumeration already includes all matching complements at that size. Their official [Figshare deposit](https://doi.org/10.6084/m9.figshare.16910170) is credited. No author GitHub submission endpoint was verified. No third-party PDF or source code is included.

The earlier m=2 [proof](proof.md), scientific checkers, and both certificates remain byte-identical. The old proof/problem/literature/provenance describe the historical state in which m≥3 was unresolved. That limitation is superseded by this proposed all-matching argument; those historical bytes are retained as evidence. This is an augmentation of the existing lab PR18, not a separate result split. `all-matchings/proof.md` names the owner scripts `certificate.py` and `labelled_check.py`; their public locations are `science/all_matching_owner.py` and `science/all_matching_labelled.py`.

The derivation and checks used AI. A separate fresh AI context independently derives the drift from event counts and audits held-out labelled configurations. Mathematical/package review status is recorded in accompanying receipts, and neither AI review nor successful reproduction is human peer review or formal proof-assistant certification. The formal workflow gate lacked trusted runner provenance and did not certify this work; no synthetic job IDs are supplied. Current priority is limited to the inspected primary-source set.

Use Python 3.12 with the pinned dependencies:

    python -m pip install -r requirements.txt
    python -B verify.py verify
    python -B -O verify.py verify

The wrapper runs all five scientific workers with PYTHONOPTIMIZE=0 in external disposable directories, retaining their assertions even when the wrapper is optimized. It recomputes the two historical certificates, the new symbolic certificate, the owner's separate labelled audit, and the fresh independent certificate, then compares exact bytes. The complete public inventory is hash-checked before and after computation; missing, unexpected, changed, or symlink inputs are rejected. Verification writes no project files. Direct scientific-script execution in this archive is unsupported because the original checkers write beside themselves.

Regeneration requires a new absolute directory outside this project with an existing parent:

    python -B verify.py regenerate --output-dir /absolute/new-directory

It writes all five certificates and an entire project manifest, and refuses existing or in-project destinations. All six outputs must compare byte-for-byte with the archive. The manifest covers every public input except itself; it is a consistency inventory, not a signature or mathematical proof. No GPU, network, author software, or special arithmetic backend is required. New original code is MIT licensed. The existing Linux/macOS Python3.12 workflow remains applicable.
