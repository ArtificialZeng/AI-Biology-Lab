# A reversible cubic mass-action system with an algebraic limit cycle

This original code and proof packet gives one explicit two-species mass-action network with all ten complexes `(i,j)` satisfying `i+j <= 3` and all 90 directed edges between distinct complexes. Every rate is a fixed finite positive rational number. Its positive ellipse `(x-2)^2 + 4(y-5)^2 = 1` is one hyperbolic attracting limit cycle. All source and target complexes have molecularity at most three, so this is a witness in the source class W3.

The existential statement concerns this network and this cycle. It does not assert uniqueness of cycles in the whole phase plane or classify other W3 networks. Reversibility means that each directed reaction has its reverse; it does not mean detailed balance. There is no claim of a minimal reaction count or biological validation.

Read [proof.md](proof.md) for the analytic argument, [witness.json](witness.json) for all rational data, and [check.py](check.py) for the exact reconstruction and Sturm checks. The checker reconstructs the vector field from rates; stored polynomial identities are verified rather than accepted as inputs to the conclusion. The transverse multiplier argument is explained in the proof.

## Reproduce

Use Python with SymPy 1.13.3 (the CI template tests Python 3.12 and 3.13):

```sh
python -m pip install -r requirements.txt
python reproduce.py --work-dir /absolute/path/outside-this-project/replay
PYTHONOPTIMIZE=2 python -O reproduce.py --work-dir /absolute/path/outside-this-project/replay-optimized
```

Each run checks the manifest, starts an assertion-enabled child with `PYTHONOPTIMIZE` removed, verifies the pinned SymPy version, reconstructs the witness, and writes logs plus `result.json` only in the supplied external work directory. It works with a relocated read-only project. A passing hash check establishes byte consistency with the local manifest, not mathematical validity or trusted provenance; the scientific checker runs separately. The files themselves are not modified by replay.

## Source and status

Craciun and Erban, *Planar chemical reaction systems with algebraic and non-algebraic limit cycles*, Journal of Mathematical Biology 90, 64 (2025), [DOI](https://doi.org/10.1007/s00285-025-02221-0), [arXiv v2](https://arxiv.org/html/2406.05057v2), Definitions 2–3, Lemma 10 and the paragraph before Lemma 11, provide the model class and cubic weakly reversible algebraic-cycle question. See [sources.md](sources.md) for attribution and scope.

A fresh AI mathematical review reconstructed the field from all reaction rates and verified the oval and transverse-sign certificates using exact rational Bernstein bounds, separately from the supplied Sturm checks. It found no mathematical gap at the stated existence scope. A checker format strengthening was separately rechecked with valid and corrupted inputs. A separate proof-specific source/value review supported the contribution within the inspected primary literature. Literature priority and human expert verification are unconfirmed. AI assisted construction, coding, and checking; computational checks and AI reviews are not human peer review. No formal verification or manuscript release gate is claimed.

The MIT license applies only to the original files in this packet. No third-party paper, software source, or licensed dataset is redistributed. SymPy is an external dependency with its own license.
