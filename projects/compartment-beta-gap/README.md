# Explosion in a polynomial-decay allocation mixture

**Status: partial theorem checked in a fresh substantive AI mathematical review.** [Review scope](review-summary.md). This package does not report external human peer review, journal acceptance, biological validation or a sharp transition.

For the published compartment model, fix `10/11 < beta < 1`, `0 < c < 1` and `0 < p < 1`. At a split of a compartment with `s` substrate molecules, independently allocating each molecule with bias `p` is chosen with probability `c/(1+s)^beta`; otherwise all molecules stay together in one daughter. From every finite initial state with positive total substrate, the theorem gives **strictly positive probability of finite-time explosion**. It does not assert almost-sure explosion. The interval `0 < beta <= 10/11` remains unresolved here.

The model receives enzymes at rate 1 per compartment, produces substrate at rate `es`, splits at rate `s`, and merges every unordered compartment pair at rate 2. It has no external compartment inflow or exit. [Full argument](proof.md), [precise statement](claim.json), [sources and overlap](source.md), [third-party credit](THIRD_PARTY.md).

Yoo's public September 2026 preprint already proves the empty-budding and independent-allocation endpoint outcomes, as well as the `beta >= 1` regime. Its stated intermediate gap motivates this narrower project. We do not claim a new resolution of the whole original Open Problem 3.8. A bounded primary-source search is not an absolute priority guarantee. No verified source-code GitHub upstream was located; the cited follow-up is distributed through Zenodo. This package does not invent an upstream PR destination.

## Reproduce the finite exact checks

Python 3.12 or 3.13, with the pinned SymPy version:

```sh
cd projects/compartment-beta-gap
python -m pip install -r requirements.txt
python -B verify.py verify
python -B -O verify.py verify
```

Normal and optimized Python must print identical deterministic JSON. Verification reads the archive, checks every project's input/evidence file against `archive_manifest.json`, recomputes the finite event counts and algebra, and compares them with `evidence/certificate.json`. It writes only to stdout and never updates project files. It uses explicit exceptions, so Python optimization does not remove checks.

The script confirms a product identity, exponent balances, rational/integer constants, and 70 finite enzyme-carrier configurations. **It does not prove explosion by simulation or finite enumeration.** The infinite-state theorem depends on the analytical drift, stopping, compensator, marked-process coupling and reachability arguments in `proof.md`.

Regeneration is an explicit separate action and requires a new absolute output filename outside this project:

```sh
python -B verify.py regenerate --output /absolute/external/directory/new-certificate.json
```

The parent directory must already exist. Existing destinations and destinations inside the project are rejected. Review mathematical changes before manually replacing archived evidence and rebuilding the manifest. A passing integrity hash is not mathematical certification. The manifest covers every regular project file except itself; it does not attempt a circular self-hash.

This directory includes its own [MIT license](LICENSE) so that it can be copied as a standalone package.

AI assisted derivation, code and separate model-context checking. Those checks are not external human peer review. Proof thresholds are deliberately large; no population or simulation at the displayed witness size was run.
