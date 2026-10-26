# Dual futile cycle: exact stability on an eight-parameter symmetry slice

For the positive Jacobian family in Nicola Vassena's [arXiv:2608.27081v2](https://arxiv.org/html/2608.27081v2), impose the eight paired-parameter equalities specified in the [proof](PROOF.md). They make the reduced six-dimensional matrix equivariant under the involution `(Q1 Q2)(U X)(V Y)`.

For all eight remaining parameters strictly positive, define

```
Delta = w*(f+v+z) - f*z.
```

The exact spectral classification on this slice is:

| Delta | Eigenvalues |
| --- | --- |
| positive | all six have strictly negative real part |
| zero | one simple zero; five have strictly negative real part |
| negative | one simple positive real eigenvalue; five have strictly negative real part |

Thus positive determinant is equivalent to Hurwitz stability **on this slice**. The proof reduces the matrix to two cubic blocks and gives explicit positive Routh expressions. No member of the slice has a nonzero purely imaginary eigenvalue.

This is a subsidiary result. The unrestricted sixteen-parameter Conjecture 3.1 remains unresolved. The source already proves the separate mass-action no-Hopf theorem. This project neither reproves that as a new discovery nor excludes periodic orbits in nonlinear systems.

## Reproduce

The main checker needs only the Python standard library:

```sh
python projects/dual-futile-symmetry/verify_symmetry_slice.py --self-test
python -O projects/dual-futile-symmetry/verify_symmetry_slice.py --self-test
```

It checks the change of basis, full characteristic-polynomial factorization, coefficient identities, positive Routh expressions, and all three threshold regimes using exact rational polynomial arithmetic. Three corrupted claims must be rejected.

A second algebra implementation is optional:

```sh
python -m pip install -r projects/dual-futile-symmetry/requirements-optional.txt
python projects/dual-futile-symmetry/derive_symmetry_slice_sympy.py
python projects/dual-futile-symmetry/independent_reconstruction.py
```

[Proof](PROOF.md) · [source record](SOURCE.md) · [verification and limits](VERIFICATION.md) · [exact receipts](evidence/) · [credit](THIRD_PARTY.md)

Prepared with AI assistance and a separate skeptical AI review. No external human review or priority confirmation has been obtained. A bounded search did not locate this precise threshold; that does not establish novelty. The natural first upstream question is whether the author already knows this symmetry corollary.
