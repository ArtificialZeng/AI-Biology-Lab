# A mass-action isola and the boundary-relay conjecture

An exact counterexample to the **given-parameter-family** reading of Conjecture 2 in Adenane, Avram and Halanay (2026), §6.3. The eight reactions in [proof.md](proof.md), with only the positive inflow rate lambda varying, give

\[
\dot x=x(-7+4x-x^2+4y-y^2),\qquad \dot y=\lambda-y,\quad\lambda>0.
\]

The entire positive equilibrium set is the closed circle

\[
y=\lambda,\qquad (x-2)^2+(\lambda-2)^2=1.
\]

It stays in the positive state and parameter domain. The only boundary equilibrium has eigenvalues \(-3-(\lambda-2)^2,-1\), so no boundary relay connects to the circle. The two folds are ordinary saddle-nodes. [Source and exact scope](source.md) · [Proof and reactions](proof.md) · [Exact evidence](evidence/exact-certificate.json).

**Scope matters.** This concerns a single smooth family in which the other rate constants are fixed. It does not give a counterexample in the full independently variable kinetic-rate space. Indeed, allowing the X-loss constant to vary gives an explicit positive-rate equilibrium path from this circle to a boundary transcritical point; the proof and verifier include that path. Isolas are classical. The proposed contribution is the source-specific counterexample, not a new bifurcation mechanism, a result about the Feng–Gavish model, or a failure of the local relay theorem.

Status: the mathematics in the stated single-family scope passed a fresh substantive AI review, including a separate integer-polynomial reconstruction. This is a manual AI review, not formal release-gate certification. The packaged checks were tested locally. No external human peer review, author confirmation, or priority guarantee is claimed. AI assisted source retrieval, construction, coding and review. See [review summary](review-summary.md).

## Reproduce

From this directory, use Python 3.12 or 3.13. No network or GPU is needed after installing the pinned dependency:

```sh
python -m pip install -r requirements.txt
python verify.py verify
python -O verify.py verify
python verify_independent.py
python -O verify_independent.py
```

The verifier checks the hashes of every public file in this project except the manifest itself, reconstructs all eight reactions, enumerates every siphon, recomputes the exact certificate and compares it byte-for-byte with the archive. It writes JSON only to stdout; ordinary verification never modifies the archive. All mathematical checks use explicit exceptions and remain active under Python optimization.

To review a deliberate regeneration, supply a **new, nonexistent directory outside this project**:

```sh
python verify.py regenerate --output-dir /tmp/relay-isola-regeneration
```

This exclusively creates a proposed certificate and manifest there. It never replaces the project archive. A mathematical change requires corresponding proof and evidence review; updating hashes alone is not validation.

Source papers are linked, not redistributed. [Third-party credit](THIRD_PARTY.md).
