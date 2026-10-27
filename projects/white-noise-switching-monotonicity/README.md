# Fixed-target switching-rate example for clarification

**Shawn Zeng (Zijian Zeng), UCSI University**

This technical note supplies an exact fixed-target example for clarifying the switching-rate discussion in Browning and Hamis (2025). For one prescribed dormant proportion, doubling a positive switching rate strictly lowers stationary logarithmic growth. A companion slow-rate calculation explains this behavior on a natural parameter branch. The published prose conjecture remains ambiguous and unresolved outside this subsidiary fixed-target reading; a target reoptimized at every rate is not analyzed here.

The source mechanism is passive switching between proliferative and dormant/persister phenotypes under correlated environmental noise. The mathematical state is concentration and the objective is logarithmic growth; this is not a finite-population survival/extinction theorem or an experimentally calibrated bacterial prediction. Equivalent stochastic dispersal models, stationary laws and nonmonotone rate responses are prior work. The potential increment under review is a source-specific fixed-target rational clarification, not a new general discovery of nonmonotonicity. The source model, fast-switching optimum and prior theory are credited in [sources.md](sources.md).

For the exact parameters mu1=1, mu2=101/1000, sigma1=1, sigma2=1/10 and rho=0, the target remains q=1/10. The rational switching rates satisfy

    omega1=101/100000 < omega2=101/50000,
    F(omega2) < 499899/1000000 < F(omega1).

The [finite proof](proof.md) uses elementary integral inequalities and rational arithmetic, with no numerical integration in the decisive comparison. The [asymptotic theorem](asymptotics.md) shows F(omega)=phi1−q*omega+o(omega) for every source-admissible fixed optimum0<q<1/4, while the fast-rate limit exceeds phi1. Neither result says faster switching is always harmful.

**A fresh-context AI review found the frozen mathematics valid and required revisions to prior-work attribution and scientific framing.** Those attribution and framing revisions are incorporated here; this local review does not establish human peer review, priority or formal runner certification. AI assisted the derivation, coding and local checks. No external human peer review, public author endorsement, historical priority or formal runner-certified release pass is claimed. See [review-boundary.md](review-boundary.md).

## Reproduce

Python3.12 or3.13, CPU only, with pinned SymPy1.13.3. Keep the environment outside this project.

```sh
python -m pip install -r requirements.txt
python verify_model.py verify
python -O verify_model.py verify
python verify_bounds.py verify
python -O verify_bounds.py verify
```

The model checker independently applies the Itô generator to the original population equations, verifies the ratio/log-growth identities and the stationary density. The bounds checker uses standard-library exact fractions to check the explicit parameter, power and separating inequalities. These programs support reading the analytic proofs; they are not proof assistants. No Julia program or source-author code is executed or redistributed.

Ordinary verification checks the manifest against every project input and byte-compares a recomputed certificate. It writes no files. `compute` prints a certificate without changing the archive. Explicit regeneration recomputes both certificates and proposes a new manifest in a new absolute directory outside this complete project:

```sh
python verify_model.py regenerate --output-dir /absolute/outside/path/new-switching-proposal
```

The output's parent must exist, and the output directory must not exist. Regeneration never overwrites this archive; changed mathematics still needs review. The [MIT license](LICENSE) covers only this project's original additions; [third-party rights](THIRD_PARTY.md) remain separate.

This is a source-specific technical clarification, not a claimed standalone journal contribution or a new general model, stationary law, biological mechanism or discovery of nonmonotone dispersal.
