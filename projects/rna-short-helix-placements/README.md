# RNA short-helix placement results

Two structural results for a simplified RNA inverse-folding model, with exact finite verification:

- **Three-helix classification.** For motif-free targets with exactly three maximal helices in total, each of length 2, the proof identifies precisely two geometries that prevent a modulo-2 separated coloring. Compressing each positive unpaired run to one dot gives 288 target classes: 4 failures of this coloring certificate and 284 constructive successes.
- **Comparable-short corollary.** If all length-2 helices lie on one ancestor chain, the published F/Q construction gives a separated coloring with no bound on their number. Other helices may branch and have any length at least 3. This is a short structural corollary of Jogalekar's local construction, not a new general folding algorithm.

A successful coloring yields a sequence with the target as its unique maximum-pair compatible fold, by the prior separated-coloring theorem. A failed modulo-2 certificate does **not** establish that a target is undesignable.

The model permits only A-U/C-G pairs, forbids pseudoknots, uses zero minimum base-pair span, and assigns energy minus the number of pairs. These statements make no thermodynamic, kinetic, experimental, or biological-realizability claim. [Exact scope](SCOPE.md) · [Proof](PROOF.md) · [Sources and credit](SOURCE.md) · [Literature screen](LITERATURE.md).

A separate AI mathematical review found no gap in the stated scope and requested attribution and release-framing changes, which are reflected here. Independent human expert review and priority confirmation remain pending. The new statements and Python checks have **not** been Lean formalized. This package is not an accepted runner/release-gate record.

## Reproduce

Python 3.12 or 3.13; standard library only:

```sh
python projects/rna-short-helix-placements/replay.py
```

The command verifies the archive manifest, reruns all three checkers in ordinary and optimized (`-O`) Python, and compares their complete deterministic JSON results with the archived receipts. It prints a compact replay receipt and does not rewrite the archive. All mathematical checks use explicit errors, so optimization cannot remove them. [Verification details](VERIFICATION.md).

The raw dot-bracket checker exhausts 729 colorings for every normalized target and independently computes the maximum-pair optimum and number of optimal folds for all 284 positive sequences. Six supplementary examples have 3, 4, 5, 7, 9, or 12 short helices, with lengths 45–170. The proof, rather than those examples, establishes the unbounded corollary.

Third-party papers are linked and are not redistributed. [Third-party material and licensing](THIRD_PARTY.md).
