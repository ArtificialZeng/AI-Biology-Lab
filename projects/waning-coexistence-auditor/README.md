# Coexistence with arbitrary waning immunity

For the full eight-compartment, two-strain Feng–Gavish model, ordered invasion implies a strictly positive coexistence equilibrium for **all finite positive model parameters and all finite nonnegative waning rates**. If both resident thresholds are met, the implication is an equivalence. This closes the additional feasible-root hypothesis in the source's full-model existence statement.

Let `v_i=gamma_i+mu`, `p=v1/beta1<q=v2/beta2`, and `s0=Lambda/mu`. Assume `Lambda,mu,beta_i,gamma_i,eta_i,sigma_i>0` and `theta1,theta2,theta12>=0`. Define

```text
R2(E1) = [p + eta2*sigma2*gamma1*(s0-p)/(v1+theta1)] / q.
```

- If `R2(E1)>1`, a strictly positive coexistence equilibrium exists on `N=s0`.
- If also `q<=s0`, positive coexistence exists **if and only if** `R2(E1)>1`.
- Exchange strain labels for the other strict ordering. Equal thresholds and zero values of the strictly positive parameters are outside this statement.

The proof writes the reduced susceptible equation as `F(s)=Lambda-mu*N(s)` on the positive reconstruction interval. Its left endpoint has `N(s)->infinity`; the invasion inequality gives `N(p)<s0`. Continuity then locates a feasible equilibrium. A population lower bound proves the conditional converse. An exact rational example explains why `q<=s0` cannot simply be dropped.

[The proof](proof.md) is preserved byte for byte from the audited research snapshot. Its opening “candidate awaiting review” sentence records its historical status; the subsequent fresh AI mathematical review accepted this exact scope. See [review summary](review-summary.md) for the audited SHA-256 and limitations. This is neither external human peer review nor journal acceptance or a priority determination.

The abstract A1–A5 **Conjecture 1 remains unresolved**. The source's formal **Question 2 root counts and multiple-root loci remain unresolved**. No uniqueness, equilibrium stability, persistence, or periodic-orbit result is claimed. The permanent-immunity/no-ADE subcase and classical ADE-supported coexistence have prior literature; see [sources and nearest results](source.md).

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

Explicit regeneration computes both certificates into a **new, absolute directory outside this project**; it never replaces the archive:

```sh
python verify.py regenerate --output-dir /absolute/outside/path/new-waning-proposal
```

The parent directory must exist and the requested output directory must not exist. Review any mathematical change and both certificates before deliberately replacing an archive; updating hashes alone does not validate changed mathematics.

The directory is portable unchanged to another repository, including `notes/waning-coexistence`; commands stay relative to this directory. Original material here is covered by the directory-scoped [MIT license](LICENSE). [Third-party rights and credit](THIRD_PARTY.md) are separate. No source PDF, manuscript TeX, private prompts or raw audit conversations are included.
