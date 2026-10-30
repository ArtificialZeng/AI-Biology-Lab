# Group-based 3-sunlet dimensions

This archive proves exact complex affine-cone dimensions for two classes of group-based models on the source 3-sunlet:

| Model | Range | Dimension |
|---|---|---|
| General model, B={id} | Every finite abelian G, q=|G|≥5 | 5q−4 |
| Cyclic inversion model, G=C_q, B={id,inversion} | Every integer q≥5 | 5 floor(q/2)+1 |

For the general model, the published odd-order theorem is credited to Cox–Gross–Martin. The even-order proof is supplied by [the exponent-greater-than-two argument](proof.md) and [the elementary-two argument](elementary_two_proof.md). The published finite elementary-two cases r=3,4,5 are credited prior work. A consequence for B={id}, using the source's existing cut-edge induction, is dimension (q−1)(m−c)+1 for its binary semi-directed level-1 class, with root-suppressed, uncontracted edge count m and c cycles.

The [cyclic-inversion proof](inversion_proof.md) uses orbit-distance weights and the full characteristic-zero left kernel. It retains all q² consistent leaf-label coordinates while sharing edge parameters on the k=floor(q/2)+1 inversion orbits. Two orbit labels alone do not determine the orbit of their sum. The broader arbitrary-B dimension problem remains unresolved.

[Statement](statement.md), [sources](sources.md) and [review status](review-summary.md) record scope and credit. Finite matrices are diagnostics; the proofs establish the infinite quantifiers. Historical priority and external human confirmation remain pending.

## Reproduce

Use CPython3.13 and python-flint0.9.0 for all seven checkers. The scalar and elementary-two independent checkers also work with only the standard library.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -B verify_independent.py verify
python -B verify_elementary.py verify
python -B verify_scalar.py verify
python -B verify_elementary_independent.py verify
python -B verify_inversion.py verify
python -B verify_inversion_independent.py verify
```

Each command also works with python -O -B. Verification checks every archived input and compares recomputed certificates without writing files. The compute subcommand prints a certificate. Explicit regeneration proposes all seven certificates and the manifest in a new absolute directory outside the project:

```sh
python -B verify.py regenerate --output-dir /absolute/outside/path/new-sunlet-proposal
```

The output must not exist and its parent must exist. Mathematical changes need new review as well as new hashes. prior-r11-proof.md preserves the identical dependency copy referenced by the elementary-two proof. Original additions use the directory-scoped MIT license. No author paper, Julia source or private correspondence is redistributed.
