# Verification and review limits

On 26 October 2026, a separate AI review accepted the precise endpoint counterexample and strict-interior inclusion in `PROOF.md`. This is not external human peer review, an author response, or a priority certificate.

The review independently derived the transfer function from the state matrix and an input/output cofactor, checked every hypothesis against the source, obtained the minor `-k1^2*k3`, and checked the rational inverse. The witness has all five maximal minors `(6,18,0,-12,-60)` in lexicographic row order and no transfer-function cancellation. The general inclusion proof was checked on the complete real hyperplane, including its zero-rate portion by polynomial density.

`independent_reconstruction.py` constructs the state matrix without importing `verify.py`. It verifies the endpoint exactly and differentiates numerator/denominator determinants by cofactor expansion for 60 held-out evaluations (58 distinct tuples), with `n=3,...,6`, every `p>=3`, swapped first/last interior positions, and one-, two-, and `(n-1)`-leak patterns. Each tuple is tested at one positive integer point satisfying the hyperplane; the finite tests supplement the written all-n proof. Both implementations use SymPy, so this is independence of construction, not independence from a shared CAS.

The main checker also passes under `python -O`, rejects three deliberately corrupted claims, and was tested with stdout and a separate `--output` path. Root command `python scripts/verify_open_biology.py` recomputes receipts without changing the archives.

The published conjecture's endpoint was independently checked through the publisher's search-indexed primary full text. Direct page fetches timed out in that review; it cannot exclude an unindexed same-day amendment. The accessible arXiv v3 agrees. Bounded correction/erratum and counterexample searches found no direct matching correction; unpublished or unindexed work is not excluded. Ask the authors to assess the interpretation and priority before promoting the result.
