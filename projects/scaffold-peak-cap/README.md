# Peak bounds and one-vacancy peak construction for scaffold binding

This original code and proof packet studies **exclusive** occupancy in the independent-site equilibrium scaffold model. It combines two mathematically different quantifier statements:

- For **three sites**, exclusive occupancy of sites 1 and 2 has **at most two strict local maxima for every positive parameter vector**. The source's known two-peak example makes this cap sharp.
- For **each integer p>=2**, with **p+1 sites** and exclusive occupancy of sites 1 through p, an explicit **finite positive rational parameter vector** has **at least p strict local maxima**. The number of sites increases with p, and only one site is vacant. This is not an exact peak count or a claim of unbounded peaks at a fixed number of sites.

The three-site proof in `proof.md` reduces the stationary equation to a third-derivative sign. The checker reconstructs the resulting exact polynomial and verifies all 322 integer coefficients. Fraction-only checks certify the known two-peak parameter example. The analytic reductions and Rolle argument must also be read.

The all-p construction is proved in `one-vacancy-proof.md`. Rational cutpoints and an explicit finite error bound establish alternating derivative signs for a single finite parameter vector at each p. Real analyticity supplies a strict maximum in each of p disjoint intervals. Its proof is for every integer p>=2; the serialized p=2,...,6 radical certificates in `cases.json` are finite implementation checks, not an enumeration proof for all p. `generate_certificates.py` and `check_certificates.py` use only exact rational arithmetic for those cases.

The earlier three-site result received separate AI mathematical and source/value checks. A separate AI mathematical audit independently reconstructed the all-p proof and found no mathematical gap; it also certified held-out p=7 and p=11 cases by fresh Fraction bisection of the original quadratic. A separate proof-specific source/value audit recognized the finite one-vacancy construction as a structural contribution within the inspected sources. Formal manuscript-workflow certification remains unavailable, and no manuscript/PDF is released. These checks are not human peer review; no source-author endorsement or exhaustive literature-priority claim is made. This packet concerns equilibrium concentration curves, not time dynamics, clinical effects, or biological validation.

## Reproduce

Requires Python 3.12 or 3.13 and SymPy 1.13.3. Install with `python -m pip install -r requirements.txt`, then run:

```sh
python reproduce.py --work-dir /absolute/path/outside-this-project
PYTHONOPTIMIZE=2 python -O reproduce.py --work-dir /another/absolute/external/path
```

The wrapper verifies the flat file manifest and launches fresh scientific child interpreters with assertions enabled, stripping inherited `PYTHONOPTIMIZE`. It reconstructs the coefficient identity, checks both old rational sharpness certificates, and verifies the one-vacancy cases. It also regenerates the latter into the external work directory and compares their bytes with the frozen cases. All logs and generated data go outside the project; it can be read-only. The manifest is an integrity check, not mathematical evidence.

## Model and credit

The model and previously known numerical examples come from Eduardo D. Sontag, [Dynamics and dose response in scaffold ligand binding](https://arxiv.org/html/2508.06599v6), version 6 (29 June 2026), Theorem 2 and section 6.1. Remark 3/Figure 8 gives two exclusive m=3, I={1,2} peaks; Remark 5/Figure 10 gives three exclusive m=4, I={1,2} peaks. Those known examples are not claimed as new here. The proposed contributions are the full positive-parameter three-site cap and the explicit all-p one-vacancy construction. The latter uses I={1,...,p}, so its occupancy set differs from the known four-site example when p=3.

An earlier primary model is J. Yang and W. S. Hlavacek, [Scaffold-mediated nucleation of protein signaling complexes: elementary principles](https://arxiv.org/abs/1105.4688). Complete occupancy and marginal pair occupancy differ from the exclusive partial occupancy studied here. Bounded source searches do not establish exhaustive literature priority. Third-party papers are linked rather than redistributed.

The MIT license applies only to the original text, code, and certificates in this packet. It does not relicense the linked sources or SymPy.
