# A rational reversible matrix with a nonreversible embedding

**Exact counterexample with fresh AI mathematical and source review completed; external human review and priority confirmation remain pending.** This project audits only Theorem 3.16 (1)⇒(3) of Baake–Baake–Sumner, *Embedding of reversible Markov matrices*, arXiv:2511.22800v1. It does not claim to solve the general embedding problem.

The exact rational four-state transition matrix is strictly positive, symmetric and has positive eigenvalues 1, 2^-9, 2^-27, 2^-27. It is the exponential of a strict Markov generator, yet its principal logarithm has entry −3log(2)/4. Hence arbitrary continuous-time embeddability and reversible embeddability differ for this example. Reversibility of the sampled transition matrix does not impose detailed balance on every possible continuous-time generator.

The repeated positive eigenvalue is essential to this construction. It does not show a failure of a criterion whose input explicitly demands a reversible generator, a generic-spectrum failure, an open set of counterexamples, or a biological measurement error. General nonprincipal-generator examples are prior work. See [exact statement](statement.md), [complete proof](proof.md), [primary sources and limitations](sources.md), and [review status](review-summary.md).

## Reproduce

Python 3.12 or 3.13 with SymPy 1.13.3; CPU only. Keep environments and generated outputs outside this directory.

```sh
python -m pip install -r requirements.txt
python -B verify.py verify
python -O -B verify.py verify
python -B verify_independent.py verify
python -O -B verify_independent.py verify
```

The first implementation checks the orthogonal-projector construction. The second starts with the literal rational input matrix, reconstructs spectral projectors polynomially and checks the generator/exponential identities separately. They verify exact rational identities and symbolic prerequisites for the analytic exponential argument; no rounded eigensolver, matrix logarithm or matrix exponential is used. The proof supplies the elementary logarithm/pi bounds and exponential periodicity. These checks are not proof-assistant formalization or external human peer review.

Every verify command hashes all public project inputs and byte-compares its own recomputed certificate, without writing files. `compute` emits the certificate to stdout. Explicit regeneration proposes both certificates and a manifest in a new absolute directory outside the complete project:

```sh
python -B verify.py regenerate --output-dir /absolute/outside/path/new-embedding-proposal
```

The parent must exist and the output must not exist. No archive is overwritten; changed mathematics needs review, not only new hashes. Original additions have the directory-scoped [MIT license](LICENSE). [Third-party rights](THIRD_PARTY.md) remain separate; no cited paper or private correspondence is redistributed.
