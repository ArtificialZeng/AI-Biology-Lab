# Coexistence with arbitrary waning immunity

For the full eight-compartment, two-strain Feng–Gavish model, this project gives a **complete feasible-root and multiplicity classification** for all finite positive model parameters and all finite nonnegative waning rates, answering actual-model Question 2. A strictly convex population curve and an explicit cubic coefficient/discriminant table distinguish no equilibrium, one simple equilibrium, two simple equilibria, and one double scalar root, including coincident thresholds and polynomial degree cancellations. In the strict threshold order, ordered invasion gives a unique positive coexistence equilibrium; if both single-strain thresholds are met, it is also necessary. Fresh manual AI reviews accepted the stated mathematics; external human review and historical priority remain unconfirmed.

For the strict-order consequence, let `v_i=gamma_i+mu`, `p=v1/beta1<q=v2/beta2`, and `s0=Lambda/mu`. Assume `Lambda,mu,beta_i,gamma_i,eta_i,sigma_i>0` and `theta1,theta2,theta12>=0`. Define

```text
R2(E1) = [p + eta2*sigma2*gamma1*(s0-p)/(v1+theta1)] / q.
```

- If `R2(E1)>1`, exactly one strictly positive coexistence equilibrium exists on `N=s0`.
- If also `q<=s0`, positive coexistence exists **if and only if** `R2(E1)>1`.
- Exchange strain labels for the other strict ordering. Equal thresholds are handled separately by the full root-count table; zero values of the strictly positive parameters remain outside the theorem.

The proof writes the reduced susceptible equation as `F(s)=Lambda-mu*N(s)` on the positive reconstruction interval. Its left endpoint has `N(s)->infinity`; the invasion inequality gives `N(p)<s0`. Continuity then locates a feasible equilibrium. A population lower bound proves the conditional converse. An exact rational example explains why `q<=s0` cannot simply be dropped.

[The historical first-existence proof](proof.md) is preserved byte for byte from its audited research snapshot. The [complete root-count proof](root-count/evidence/proof.md) supplies the subsequent uniqueness and multiplicity classification. The proofs' opening candidate-status wording records their historical status; the subsequent fresh AI mathematical review accepted this exact scope. See [review summary](review-summary.md) for the audited SHA-256 and limitations. This is neither external human peer review nor journal acceptance or a priority determination.

The abstract A1–A5 **Conjecture 1 remains unresolved**. The subsequent [root-count extension](root-count/README.md) answers actual-model **Question 2** at the full positive-parameter/nonnegative-waning scope, including uniqueness under ordered invasion. No equilibrium stability, persistence, or periodic-orbit result is claimed. The permanent-immunity/no-ADE subcase and classical ADE-supported coexistence have prior literature; see [sources and nearest results](source.md).

## Reproduce from this directory

Python 3.12 or 3.13; CPU only; no model API or downloaded source paper is needed.

```sh
python -m pip install -r requirements.txt
python verify.py verify
python -O verify.py verify
python verify_independent.py verify
python -O verify_independent.py verify
```

Use an existing Python environment, or create and activate a virtual environment **outside this directory** before these commands. The project must contain only its public archive inputs; the manifest rejects extra files. Both checkers reconstruct 24 elementary reactions and the published eight equations, then check exact equilibrium, population, endpoint and converse identities. The independent checker derives the recovery amplitudes by matrix inversion without importing the original derivation. Each recomputes and byte-compares its own archived JSON certificate. Their shared helper implements archive I/O only. Analytic positivity, continuity, quantified scope and literature interpretation are explained in the proof and review summary; these scripts are not formal proof assistants.

Ordinary checks emit JSON on stdout, check the hashes of every public input, and do not write files. Python optimization must not change their outcome. The certificates are [original](evidence/exact-population-certificate.json) and [independent](evidence/independent-certificate.json).

Explicit regeneration computes all four certificates into a **new, absolute directory outside this project**; it never replaces the archive:

```sh
python verify.py regenerate --output-dir /absolute/outside/path/new-waning-proposal
```

The parent directory must exist and the requested output directory must not exist. Review any mathematical change and all four certificates before deliberately replacing an archive; updating hashes alone does not validate changed mathematics.

The directory is portable unchanged to another repository, including `notes/waning-coexistence`; commands stay relative to this directory. Original material here is covered by the directory-scoped [MIT license](LICENSE). [Third-party rights and credit](THIRD_PARTY.md) are separate. No source PDF, manuscript TeX, private prompts or raw audit conversations are included.


## Full root-count extension

The [root-count extension](root-count/README.md) gives the full actual-model Question 2 classification for all finite positive model parameters and nonnegative waning, including coincident thresholds and reduced polynomial degrees. Its strict-convex population argument gives uniqueness under ordered invasion and an explicit coefficient/discriminant table for zero, one, two, or double feasible roots. A fresh manual AI review accepted this new mathematical scope on 26 October 2026; this is not external human peer review, historical-priority certification or a formal release gate. The earlier accepted existence theorem and its original proof/certificates above retain their original bytes. Broad abstract Conjecture 1 remains unresolved; no stability, persistence, or eight-dimensional saddle-node claim is added.

From this project directory, the additional checks are:

```sh
python root-count/check.py verify
python -O root-count/check.py verify
python root-count/check_heldout.py verify
python -O root-count/check_heldout.py verify
```

All four entry points (the two earlier checkers and these two) now support `compute` to emit a recomputed certificate as JSON on stdout. `verify` first checks the complete parent manifest and then byte-compares the relevant certificate. Explicit `regenerate` now computes **all four certificates** into a new exclusive absolute directory outside the project. It does not replace any archive. Final JTB full-text access and the latest detailed SIAM HTML refresh remain limited; the literature comparison does not certify priority.
