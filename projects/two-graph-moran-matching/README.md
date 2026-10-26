# A fixation bound for two deleted mutant edges

For every integer n≥4, residents reproduce on K_n and neutral mutants on K_n minus two vertex-disjoint edges under type-specific birth–death. Every singleton starting vertex satisfies

    rho_v ≤ 1/n − 4/[n²(n−1)(n(n−2)+2)] < 1/n.

The proof is an all-n symbolic superharmonic certificate, followed by an absorption argument; finite state audits are supplementary. The general arbitrary-graph conjecture and matching deletions m≥3 remain unresolved here. Positive drift of a proposed potential at n=6,m=3 is a failure of that method, not a counterexample to the conjecture.

Tkadlec, Kaveh, Chatterjee and Nowak (2023), [Evolutionary dynamics of mutants that modify population structure](https://doi.org/10.1098/rsif.2023.0355), supplies the model, the known one-deleted-edge Theorem A.1 and numerical n=6 evidence. Their official [Figshare deposit](https://doi.org/10.6084/m9.figshare.16910170) is credited; no author GitHub submission endpoint was verified. No third-party paper PDF or source code is redistributed. See literature.md for the bounded priority search; priority is not established.

The original proof and checkers were AI-assisted. A fresh AI context accepted the mathematical scope and independently derived the raw parent-target drift. This is not human peer review. The formal workflow gate could not certify the work because trusted runner provenance was absent. `provenance/research-claim.json` preserves the original intake status and its private evidence inventory; it is provenance, not a public file inventory or release certificate. Publication readiness requires separate mathematical and package acceptance.

Use Python 3.12 with pinned dependencies:

    python -m pip install -r requirements.txt
    python -B verify.py verify
    python -B -O verify.py verify

Both wrapper modes explicitly run the two scientific workers with PYTHONOPTIMIZE=0, retaining their assertions. This does not claim that their assertions survive direct optimized execution. Workers write only in a disposable directory. Verification recomputes both complete certificates, compares exact bytes, and hashes every public file. Unexpected files, missing files, changed bytes and symlinks are rejected. Do not run science/owner_checker.py directly in the archive: the frozen original writes beside itself.

Regeneration requires a new absolute directory outside this project, whose parent exists:

    python -B verify.py regenerate --output-dir /absolute/new-directory

It recomputes both certificates and an entire project manifest; it never overwrites an existing output directory. Compare all three outputs with their archived counterparts. The manifest covers every public input except itself and is a consistency inventory, not a signature or proof of mathematical correctness.

The owner checker is byte-identical to the frozen research input. The independent checker preserves all mathematical functions and assertions; only private snapshot/path/provenance plumbing was removed. Its original byte-identical source and adaptation mapping are in provenance/. New code is MIT licensed under the repository license. No special arithmetic backend, GPU, network access or author software is required for verification.
