# Even-order 3-sunlet dimensions

For every finite abelian group G of even order n≥6 and exponent greater than two, the general group-based model (B={id}) on the source's 3-sunlet network has complex affine-cone dimension **5n−4**.

The [proof](proof.md) constructs a strict tropical chamber and parameterizes its entire left kernel, uniformly for this infinite class. The known odd-order theorem, weight pattern, upper bound and finite even examples are credited to Cox–Gross–Martin. Elementary abelian 2-groups and order four are excluded. This is a partial result on the even-order problem; dimension alone does not establish topology identifiability or a biological experiment.

Fresh AI mathematical review found no gap. Historical priority, author confirmation and external human review remain pending. No trusted-runner or proof-assistant certificate is claimed. See [statement](statement.md), [sources](sources.md) and [review status](review-summary.md).

## Reproduce

The primary exact-rational checker uses CPython 3.13 and python-flint 0.9.0. The independent checker uses only the Python standard library; it constructs literal monomial exponents separately and performs integer elimination over Q.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
python -B verify_independent.py verify
python -O -B verify_independent.py verify
```

These finite matrices are diagnostics, not the proof of the infinite quantifiers. Every verify command hashes all project inputs and compares its recomputed certificate without writing files. `compute` prints a certificate to stdout. Explicit regeneration proposes both certificates and a manifest in a new absolute directory outside this project:

```sh
python -B verify.py regenerate --output-dir /absolute/outside/path/new-sunlet-proposal
```

The parent must exist and the output must not exist. Changed mathematics requires review as well as new hashes. Original additions have a directory-scoped MIT license; no author Julia source, paper or private correspondence is redistributed.
