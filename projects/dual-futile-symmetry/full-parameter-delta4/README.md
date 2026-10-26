# Fourth Hurwitz determinant for the full sixteen-parameter dual futile cycle

For the reduced Jacobian in Nicola Vassena's [arXiv:2608.27081v2](https://arxiv.org/html/2608.27081v2), all sixteen parameters remain independent and strictly positive. Write `c6=det(G)` and let `Delta4` be the fourth leading Hurwitz determinant. The exact partial result is:

**If `det(G)>=0`, then `Delta4>0`.**

A positive quartic `Q` gives `Delta4=Q*det(G)+R`, where the exact polynomial `R` has 1,712,869 nonzero positive integer coefficients. Its coefficient of `u3*v6^2*x9^3*y12^4` is one, proving strictness throughout the positive orthant, including `det(G)=0`.

This uses neither mass-action relations nor parameter pairing. The source already proves a separate mass-action no-Hopf theorem and `Delta4>0` under MA1/MA2 without a determinant-sign premise. The result here enlarges the kinetic class on the determinant-nonnegative region; it does not globally subsume that mass-action proposition. It complements the [existing eight-parameter classification](../README.md).

The full sixteen-parameter Conjecture3.1 remains unresolved. Positivity of `Delta5` on `det(G)>0` remains to be established. No full stability, periodic-orbit, or nonnegative-parameter-boundary statement is made.

## Reproduce

Python3.12 or Python3.13; no optimizer, NumPy, network access, or large cache is required. From this directory:

```sh
python -m pip install -r requirements.txt
python -B verify_certificate.py
python -B -O verify_certificate.py
python -B independent_reconstruction.py
python -B -O independent_reconstruction.py
python -B -O integrity.py --self-test
```

Each command prints a deterministic JSON receipt and changes no package input or archived evidence. Receipts contain no current manifest digest or elapsed-time field, avoiding a self-hash cycle. The first implementation computes characteristic coefficients by principal minors and evaluates the defining 4x4 Hurwitz determinant. The second reconstructs the full nine-species reaction Jacobian, reduces by conservation laws, obtains characteristic coefficients through exact trace recursion, and derives a generic Hurwitz determinant with SymPy. All checks use explicit exceptions and remain active under `-O`.

The compact `Q` formula is in [the proof](PROOF.md), and its complete 266-term integer serialization is in [certificate.json](certificate.json). [manifest.json](manifest.json) lists and hashes every package file, including archived evidence, except the manifest itself and Python bytecode. Both checker entry points enforce the complete inventory, SHA-256 hashes, and a no-symlink rule before calculation. The checkers construct the remainder afresh; no multi-million-term remainder or discovery cache is distributed. Local receipts are in [evidence](evidence/); the CI matrix covers Ubuntu/macOS and Python3.12/3.13.

[Source and scope](SOURCE.md) · [Review summary](REVIEW.md) · [Third-party credit](THIRD_PARTY.md) · [License](LICENSE)

Prepared with AI assistance and independently reconstructed in a separate substantive AI review. External human review and absolute priority confirmation remain pending.
