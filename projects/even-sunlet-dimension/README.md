# General group-based 3-sunlet dimensions

For every finite abelian group G with order q≥5, the general group-based model (B={id}) on the source 3-sunlet has complex affine-cone dimension **5q−4**. The published odd-order theorem is credited to Cox–Gross–Martin. The contribution here supplies strict tropical proofs for the even-order groups: [exponent greater than two](proof.md), and [elementary abelian two groups](elementary_two_proof.md).

The elementary-two proof uses a Hamming-order chamber and its complete characteristic-zero left kernel for every r≥3. The source already computed the cases r=3, 4, 5; those finite cases are credited rather than claimed as new. The same source cut-edge induction consequently gives dimension **(q−1)(m−c)+1** for its binary semi-directed level-1 network class, counting m in the root-suppressed, uncontracted representation. The model with nontrivial B, groups of order≤4, and topology identifiability are outside this uniform theorem.

This is one coherent extension of the existing even-order research contribution. [Statement](statement.md), [sources](sources.md), and [review status](review-summary.md) separate the known results, new proof and limitations. Finite matrices are diagnostics, not proofs of infinite quantifiers.

## Reproduce

Use CPython 3.13 and python-flint 0.9.0 for the full exact-rational matrices; the scalar and separate literal-monomial checkers use only the standard library.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B verify_independent.py verify
python -B verify_elementary.py verify
python -B verify_scalar.py verify
python -B verify_elementary_independent.py verify
```

Each command also works with `python -O -B`. Verification hashes every input and compares recomputed certificates without writing files. `compute` prints a certificate. Explicit regeneration proposes all certificates and the manifest in a new absolute directory outside the project:

```sh
python -B verify.py regenerate --output-dir /absolute/outside/path/new-sunlet-proposal
```

The output must not exist and its parent must exist. Mathematical changes need new review as well as new hashes. `prior-r11-proof.md` is an identical dependency copy of `proof.md`, retained to preserve the frozen elementary-two proof's file reference. Original additions use the directory-scoped MIT license. No author Julia source, paper or private correspondence is redistributed.
